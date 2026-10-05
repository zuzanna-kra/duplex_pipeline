#!/usr/bin/env python3
"""
trace_alt_families.py - follow raw Horizon alt reads into their UMI families
and see what happened to each family at the SSCS and DCS steps.

Question: at the Horizon SNV positions, the raw BAM has alt reads for the exact
expected change, but in the 0.1% (and 1%) samples they are missing from SSCS/DCS.
Where do they go?

For every strand family (MI tag, e.g. 1234/A) covering each Horizon SNV, the
script records, from the UMI-grouped BAM:
  - family size (templates), alt / ref / other counts, UMI (RX) make-up
  - the opposite strand family (1234/B)
and from the consensus BAMs:
  - the base the SSCS / DCS consensus read has at that position, before
    filtering (…_mapped.bam) and after filter + clip (…_filtered_clipped_sorted.bam)

Each family carrying at least one alt template gets one category:
  retained_DCS                      alt survives to the final DCS
  SSCS_only:no_partner_strand       alt SSCS, but the other strand was never sequenced
  SSCS_only:partner_too_small       alt SSCS, other strand has < 3 reads
  SSCS_only:partner_strand_not_alt  alt SSCS, other strand says ref (single-strand event:
                                    DNA damage or first-cycle PCR error)
  SSCS_only:lost_in_DCS             both strands alt but not in final DCS  <- look at these
  small_family:all_alt              < 3 reads, all alt  (fails --min-reads 3; the typical
                                    signature of reads that came from ANOTHER library,
                                    e.g. index hopping / cross-contamination)
  small_family:mixed                < 3 reads, alt + ref
  alt_minority                      family >= 3 reads, alt is < 50% (PCR/sequencing error)
  alt_minority:umi_split            ...and the alt reads carry a different exact UMI from
                                    the ref reads (two molecules merged by --edits 1)
  alt_majority_lost_in_SSCS         family >= 3 reads, mostly alt, but no alt in final
                                    SSCS  <- look at these

Subcommands
  bed        write a 1-bp BED of the SNVs (for samtools view -L)
  trace      one sample -> <outdir>/<sample>_families.tsv
  summarise  all *_families.tsv in a directory -> summary tables, including a
             cross-sample check: do the alt molecules in one sample have the same
             fragment ends + UMI as a molecule in another sample?

Requires pysam. Run via run_trace.sh.
"""
import argparse
import collections
import csv
import glob
import os
import re
import sys

import pysam

MIN_SS_DEFAULT = 3          # legacy MUFS3: SSCS --min-reads 3; DCS 3 per strand (6 total)
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


# ---------------------------------------------------------------------------
# Variants
# ---------------------------------------------------------------------------
def read_variants(path):
    """Read a CSV/TSV truth table. Column names are matched case-insensitively."""
    with open(path, newline="") as fh:
        head = fh.read(8192)
        fh.seek(0)
        delim = "\t" if head.count("\t") >= head.count(",") else ","
        rdr = csv.DictReader(fh, delimiter=delim)
        cols = {c.lower().strip().lstrip("#"): c for c in rdr.fieldnames}

        def pick(*names):
            for n in names:
                if n in cols:
                    return cols[n]
            return None

        c_chr = pick("chrom", "chr", "chromosome", "contig")
        c_pos = pick("pos", "position", "start", "pos_hg38", "grch38_pos", "position_hg38")
        c_ref = pick("ref", "reference", "ref_allele")
        c_alt = pick("alt", "alternate", "alt_allele", "var")
        c_name = pick("name", "label", "variant", "id", "mutation", "protein", "hgvsp")
        c_gene = pick("gene", "gene_symbol", "symbol")
        missing = [n for n, c in [("chrom", c_chr), ("pos", c_pos), ("ref", c_ref), ("alt", c_alt)] if c is None]
        if missing:
            sys.exit(f"ERROR: {path}: could not find column(s) {missing}. Columns found: {rdr.fieldnames}")
        out, skipped = [], []
        for row in rdr:
            ref, alt = row[c_ref].strip().upper(), row[c_alt].strip().upper()
            chrom, pos = row[c_chr].strip(), int(float(row[c_pos]))
            name = " ".join(x for x in [row.get(c_gene, "") if c_gene else "",
                                        row.get(c_name, "") if c_name else ""] if x).strip()
            name = name or f"{chrom}:{pos}{ref}>{alt}"
            if len(ref) == 1 and len(alt) == 1 and ref in "ACGT" and alt in "ACGT":
                out.append(dict(name=name, chrom=chrom, pos=pos, ref=ref, alt=alt))
            else:
                skipped.append(name)
    if skipped:
        print(f"[variants] skipping {len(skipped)} non-SNV entries: {', '.join(skipped)}", file=sys.stderr)
    return out


def match_contig(chrom, references):
    if chrom in references:
        return chrom
    alt = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
    return alt if alt in references else None


# ---------------------------------------------------------------------------
# Read helpers
# ---------------------------------------------------------------------------
def base_at(read, pos0):
    """(base, qual, qpos) at 0-based ref pos; ('-', None, None) if deleted; None if not covered."""
    if read.cigartuples is None:
        return None
    rpos, q = read.reference_start, 0
    for op, ln in read.cigartuples:
        if op in (0, 7, 8):                      # M = X
            if rpos <= pos0 < rpos + ln:
                qp = q + pos0 - rpos
                quals = read.query_qualities
                return read.query_sequence[qp], (quals[qp] if quals is not None else None), qp
            rpos += ln
            q += ln
        elif op in (1, 4):                       # I S
            q += ln
        elif op in (2, 3):                       # D N
            if rpos <= pos0 < rpos + ln:
                return ("-", None, None) if op == 2 else None
            rpos += ln
        # H, P consume nothing
    return None


def _clips(ops):
    """leading and trailing clip (S+H) lengths from [(op_char, len)]"""
    lead = trail = 0
    for op, ln in ops:
        if op in "SH":
            lead += ln
        else:
            break
    for op, ln in reversed(ops):
        if op in "SH":
            trail += ln
        else:
            break
    return lead, trail


def unclipped_5p(read):
    ops = [("MIDNSHP=X"[op], ln) for op, ln in read.cigartuples]
    lead, trail = _clips(ops)
    if read.is_reverse:
        return read.reference_end + trail
    return read.reference_start - lead


def mate_unclipped_5p(read):
    if read.has_tag("MC"):
        ops = [(op, int(ln)) for ln, op in CIGAR_RE.findall(read.get_tag("MC"))]
        lead, trail = _clips(ops)
        if read.mate_is_reverse:
            rlen = sum(ln for op, ln in ops if op in "MDN=X")
            return read.next_reference_start + rlen + trail
        return read.next_reference_start - lead
    return read.next_reference_start          # approximation without MC


def molecule_key(read, chrom):
    """Library-independent molecule signature: chrom, both unclipped 5' ends, UMI halves (order-free)."""
    rx = read.get_tag("RX") if read.has_tag("RX") else ""
    umi = "-".join(sorted(rx.split("-"))) if "-" in rx else rx
    ends = sorted([unclipped_5p(read), mate_unclipped_5p(read)])
    return f"{chrom}:{ends[0]}-{ends[1]}:{umi}"


def split_mi(mi):
    mi = str(mi)
    if "/" in mi:
        b, s = mi.rsplit("/", 1)
        return b, s
    return mi, ""


def combine(calls):
    """Combine per-read calls of one template / one consensus molecule into one call."""
    good = {c for c in calls if c not in ("lowq", "N")}
    if len(good) == 1:
        return good.pop()
    if len(good) > 1:
        return "discord"
    if "lowq" in calls:
        return "lowq"
    return "N" if calls else "none"


# ---------------------------------------------------------------------------
# Loading layers at one position
# ---------------------------------------------------------------------------
def raw_families(bam, chrom, pos0, min_mapq, min_bq, max_depth):
    """strand-family MI -> {templates: {qname: [calls]}, rx: {qname: RX}, key}"""
    fams = collections.defaultdict(lambda: {"templates": collections.defaultdict(list), "rx": {}, "key": None})
    n = 0
    for read in bam.fetch(chrom, pos0, pos0 + 1):
        if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail:
            continue
        if read.mapping_quality < min_mapq or not read.has_tag("MI"):
            continue
        b = base_at(read, pos0)
        if b is None:
            continue
        base, qual, _ = b
        call = base if base == "-" else ("lowq" if (qual is not None and qual < min_bq) else base)
        f = fams[str(read.get_tag("MI"))]
        f["templates"][read.query_name].append(call)
        f["rx"][read.query_name] = read.get_tag("RX") if read.has_tag("RX") else ""
        if f["key"] is None and read.is_paired and not read.mate_is_unmapped:
            f["key"] = molecule_key(read, chrom)
        n += 1
        if max_depth and n >= max_depth:
            print(f"[warn] {chrom}:{pos0+1} hit --max-depth {max_depth}", file=sys.stderr)
            break
    return fams


def consensus_calls(bam, chrom, pos0, min_mapq, with_strand_tags=False):
    """consensus MI -> dict(call, detail). Low-MAPQ reads are reported as 'lowMAPQ:<base>'."""
    if bam is None:
        return None
    per = collections.defaultdict(lambda: {"calls": [], "lowmq": [], "cd": [], "ce": [], "ac": [], "bc": []})
    for read in bam.fetch(chrom, pos0, pos0 + 1):
        if read.is_unmapped or read.is_secondary or read.is_supplementary or not read.has_tag("MI"):
            continue
        b = base_at(read, pos0)
        if b is None:
            continue
        base, qual, qp = b
        d = per[str(read.get_tag("MI"))]
        (d["calls"] if read.mapping_quality >= min_mapq else d["lowmq"]).append(base)
        if qp is not None:
            for tag, key in (("cd", "cd"), ("ce", "ce")):
                if read.has_tag(tag):
                    try:
                        d[key].append(int(read.get_tag(tag)[qp]))
                    except Exception:
                        pass
            if with_strand_tags:
                for tag in ("ac", "bc"):
                    if read.has_tag(tag):
                        try:
                            d[tag].append(read.get_tag(tag)[qp])
                        except Exception:
                            pass
    out = {}
    for mi, d in per.items():
        call = combine(d["calls"]) if d["calls"] else ("lowMAPQ:" + combine(d["lowmq"]) if d["lowmq"] else "none")
        out[mi] = {
            "call": call,
            "depth": max(d["cd"]) if d["cd"] else "",
            "err": max(d["ce"]) if d["ce"] else "",
            "ac": combine(d["ac"]) if d["ac"] else "",
            "bc": combine(d["bc"]) if d["bc"] else "",
        }
    return out


# ---------------------------------------------------------------------------
# Family stats and categories
# ---------------------------------------------------------------------------
def family_stats(f, ref, alt):
    tcalls = {q: combine(c) for q, c in f["templates"].items()}
    n_alt = sum(1 for c in tcalls.values() if c == alt)
    n_ref = sum(1 for c in tcalls.values() if c == ref)
    n_lowq = sum(1 for c in tcalls.values() if c == "lowq")
    n_tot = len(tcalls)
    n_other = n_tot - n_alt - n_ref - n_lowq
    informative = n_alt + n_ref + n_other
    alt_frac = n_alt / informative if informative else 0.0
    # UMI make-up
    rx_all = collections.Counter(f["rx"].values())
    rx_alt = collections.Counter(f["rx"][q] for q, c in tcalls.items() if c == alt)
    rx_ref = collections.Counter(f["rx"][q] for q, c in tcalls.items() if c == ref)
    top_alt_rx, top_alt_share, ref_share_same = "", "", ""
    umi_split = False
    if rx_alt:
        top_alt_rx, k = rx_alt.most_common(1)[0]
        top_alt_share = k / n_alt
        ref_share_same = (rx_ref[top_alt_rx] / n_ref) if n_ref else ""
        umi_split = (n_alt >= 2 and n_ref >= 2 and top_alt_share >= 0.8
                     and ref_share_same != "" and ref_share_same <= 0.2)
    return dict(n_templates=n_tot, n_alt=n_alt, n_ref=n_ref, n_other=n_other, n_lowq=n_lowq,
                alt_frac=alt_frac, n_rx=len(rx_all), top_alt_rx=top_alt_rx,
                top_alt_rx_share=top_alt_share, ref_share_with_top_alt_rx=ref_share_same,
                umi_split=umi_split)


def categorise(st, partner, sscs_final, dcs_final, alt, min_ss):
    if st["n_alt"] == 0:
        return ""
    if sscs_final == alt:
        if dcs_final == alt:
            return "retained_DCS"
        if partner is None:
            return "SSCS_only:no_partner_strand"
        if partner["n_templates"] < min_ss:
            return "SSCS_only:partner_too_small"
        if partner["alt_frac"] < 0.5:
            return "SSCS_only:partner_strand_not_alt"
        return "SSCS_only:lost_in_DCS"
    if st["n_templates"] < min_ss:
        return "small_family:all_alt" if st["n_alt"] == st["n_templates"] else "small_family:mixed"
    if st["alt_frac"] < 0.5:
        return "alt_minority:umi_split" if st["umi_split"] else "alt_minority"
    return "alt_majority_lost_in_SSCS"


FAMILY_COLS = ["sample", "variant", "chrom", "pos", "ref", "alt", "mi", "strand", "mol_key",
               "n_templates", "n_alt", "n_ref", "n_other", "n_lowq", "alt_frac",
               "n_rx", "top_alt_rx", "top_alt_rx_share", "ref_share_with_top_alt_rx", "umi_split",
               "partner_n", "partner_alt", "partner_ref",
               "sscs_called", "sscs_called_depth", "sscs_called_err", "sscs_final",
               "dcs_called", "dcs_called_a", "dcs_called_b", "dcs_final", "category"]


def open_bam(path, required=True):
    if not path:
        if required:
            sys.exit("ERROR: missing BAM path")
        return None
    if not os.path.exists(path):
        if required:
            sys.exit(f"ERROR: not found: {path}")
        print(f"[warn] not found, skipping layer: {path}", file=sys.stderr)
        return None
    bam = pysam.AlignmentFile(path, "rb")
    if not bam.has_index():
        sys.exit(f"ERROR: {path} has no index (samtools index it)")
    return bam


def fmt(x):
    if isinstance(x, float):
        return f"{x:.4g}"
    if isinstance(x, bool):
        return "yes" if x else "no"
    return str(x)


def cmd_trace(a):
    variants = read_variants(a.variants)
    raw = open_bam(a.grouped)
    s_called = open_bam(a.sscs_called, required=False)
    s_final = open_bam(a.sscs_final)
    d_called = open_bam(a.dcs_called, required=False)
    d_final = open_bam(a.dcs_final)
    os.makedirs(a.outdir, exist_ok=True)
    out_path = os.path.join(a.outdir, f"{a.sample}_families.tsv")
    n_alt_fams = 0
    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(FAMILY_COLS)
        for v in variants:
            chrom = match_contig(v["chrom"], raw.references)
            if chrom is None:
                print(f"[warn] contig {v['chrom']} not in BAM, skipping {v['name']}", file=sys.stderr)
                continue
            pos0 = v["pos"] - 1
            fams = raw_families(raw, chrom, pos0, a.min_mapq, a.min_bq, a.max_depth)
            stats = {mi: family_stats(f, v["ref"], v["alt"]) for mi, f in fams.items()}
            sc = consensus_calls(s_called, chrom, pos0, a.min_mapq) if s_called else {}
            sf = consensus_calls(s_final, chrom, pos0, a.min_mapq)
            dc = consensus_calls(d_called, chrom, pos0, a.min_mapq, with_strand_tags=True) if d_called else {}
            df = consensus_calls(d_final, chrom, pos0, a.min_mapq)
            for mi, st in stats.items():
                base, strand = split_mi(mi)
                pmi = f"{base}/{'B' if strand == 'A' else 'A'}" if strand in ("A", "B") else None
                partner = stats.get(pmi) if pmi else None
                sfin = sf.get(mi, {}).get("call", "absent")
                dfin = df.get(base, {}).get("call", "absent")
                cat = categorise(st, partner, sfin, dfin, v["alt"], a.min_ss)
                n_alt_fams += bool(cat)
                scd = sc.get(mi, {}) if s_called else {}
                dcd = dc.get(base, {}) if d_called else {}
                row = dict(sample=a.sample, variant=v["name"], chrom=chrom, pos=v["pos"],
                           ref=v["ref"], alt=v["alt"], mi=mi, strand=strand,
                           mol_key=fams[mi]["key"] or "", **st,
                           partner_n=partner["n_templates"] if partner else 0,
                           partner_alt=partner["n_alt"] if partner else 0,
                           partner_ref=partner["n_ref"] if partner else 0,
                           sscs_called=scd.get("call", "absent" if s_called else "NA"),
                           sscs_called_depth=scd.get("depth", ""), sscs_called_err=scd.get("err", ""),
                           sscs_final=sfin,
                           dcs_called=dcd.get("call", "absent" if d_called else "NA"),
                           dcs_called_a=dcd.get("ac", ""), dcs_called_b=dcd.get("bc", ""),
                           dcs_final=dfin, category=cat)
                w.writerow([fmt(row[c]) for c in FAMILY_COLS])
            print(f"[{a.sample}] {v['name']}: {len(stats)} strand families, "
                  f"{sum(1 for s in stats.values() if s['n_alt'])} with alt", file=sys.stderr)
    print(f"[{a.sample}] wrote {out_path} ({n_alt_fams} alt-carrying families)", file=sys.stderr)


# ---------------------------------------------------------------------------
# Summarise
# ---------------------------------------------------------------------------
CATEGORY_ORDER = ["retained_DCS", "SSCS_only:no_partner_strand", "SSCS_only:partner_too_small",
                  "SSCS_only:partner_strand_not_alt", "SSCS_only:lost_in_DCS",
                  "small_family:all_alt", "small_family:mixed",
                  "alt_minority", "alt_minority:umi_split", "alt_majority_lost_in_SSCS"]


def cmd_summarise(a):
    files = sorted(glob.glob(os.path.join(a.indir, "*_families.tsv")))
    if not files:
        sys.exit(f"ERROR: no *_families.tsv in {a.indir}")
    rows = []
    for fp in files:
        with open(fp) as fh:
            rows.extend(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        for k in ("n_templates", "n_alt", "n_ref", "n_other", "n_lowq"):
            r[k] = int(r[k])
    samples = sorted({r["sample"] for r in rows})

    # --- 1. where do the alt templates go? (per sample, all SNVs pooled)
    by_sample = collections.defaultdict(collections.Counter)      # templates
    by_sample_f = collections.defaultdict(collections.Counter)    # families
    for r in rows:
        if r["category"]:
            by_sample[r["sample"]][r["category"]] += r["n_alt"]
            by_sample_f[r["sample"]][r["category"]] += 1
    cats = [c for c in CATEGORY_ORDER if any(by_sample[s][c] for s in samples)]
    p1 = os.path.join(a.indir, "summary_alt_templates_by_category.tsv")
    with open(p1, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "total_alt_templates"] + [f"{c} (templates|families)" for c in cats])
        for s in samples:
            tot = sum(by_sample[s].values())
            w.writerow([s, tot] + [f"{by_sample[s][c]} ({100*by_sample[s][c]/tot:.1f}%) | {by_sample_f[s][c]}"
                                   if tot else "0" for c in cats])

    # --- 2. per sample x variant, incl. allele specificity of raw non-ref bases
    p2 = os.path.join(a.indir, "summary_per_variant.tsv")
    agg = collections.defaultdict(lambda: collections.Counter())
    for r in rows:
        k = (r["sample"], r["variant"])
        g = agg[k]
        g["families"] += 1
        g["templates"] += r["n_templates"]
        g["alt"] += r["n_alt"]
        g["other_nonref"] += r["n_other"]
        if r["category"]:
            g["alt_families"] += 1
            g["cat:" + r["category"]] += r["n_alt"]
    with open(p2, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "variant", "strand_families", "raw_templates", "raw_alt_templates",
                    "raw_alt_VAF", "raw_other_nonref_templates", "alt_to_other_ratio",
                    "alt_families"] + [f"alt_templates:{c}" for c in cats])
        for (s, v), g in sorted(agg.items()):
            vaf = g["alt"] / g["templates"] if g["templates"] else 0
            ratio = (g["alt"] / g["other_nonref"]) if g["other_nonref"] else ("inf" if g["alt"] else "")
            w.writerow([s, v, g["families"], g["templates"], g["alt"], f"{vaf:.3g}", g["other_nonref"],
                        fmt(ratio) if ratio != "" else "", g["alt_families"]] + [g["cat:" + c] for c in cats])

    # --- 3. cross-sample molecule matching
    # molecule = (variant, mol_key); a molecule can be present in several samples
    present = collections.defaultdict(dict)   # (variant, key) -> {sample: (n_templates, n_alt)}
    for r in rows:
        if not r["mol_key"]:
            continue
        k = (r["variant"], r["mol_key"])
        n, al = present[k].get(r["sample"], (0, 0))
        present[k][r["sample"]] = (n + r["n_templates"], al + r["n_alt"])
    p3 = os.path.join(a.indir, "summary_cross_sample.tsv")
    detail = os.path.join(a.indir, "cross_sample_alt_matches.tsv")
    with open(p3, "w", newline="") as fh, open(detail, "w", newline="") as fd:
        w = csv.writer(fh, delimiter="\t")
        wd = csv.writer(fd, delimiter="\t")
        w.writerow(["sample", "group", "molecules", "matched_in_other_sample", "pct_matched",
                    "matched_in_other_sample_as_alt"])
        wd.writerow(["sample", "variant", "mol_key", "category", "n_templates_here", "n_alt_here",
                     "other_samples(n_templates,n_alt)"])
        for s in samples:
            groups = collections.defaultdict(lambda: [0, 0, 0])
            seen = set()
            for r in rows:
                if r["sample"] != s or not r["mol_key"]:
                    continue
                k = (r["variant"], r["mol_key"])
                if (k, r["category"]) in seen:
                    continue
                seen.add((k, r["category"]))
                others = {t: v for t, v in present[k].items() if t != s}
                grp = r["category"] if r["category"] else "no_alt (background)"
                g = groups[grp]
                g[0] += 1
                if others:
                    g[1] += 1
                    if any(al > 0 for _, al in others.values()):
                        g[2] += 1
                if r["category"] and others:
                    wd.writerow([s, r["variant"], r["mol_key"], r["category"], r["n_templates"], r["n_alt"],
                                 "; ".join(f"{t}({n},{al})" for t, (n, al) in sorted(others.items()))])
            for grp in ["no_alt (background)"] + [c for c in CATEGORY_ORDER if c in groups]:
                if grp not in groups:
                    continue
                m, mt, ma = groups[grp]
                w.writerow([s, grp, m, mt, f"{100*mt/m:.2f}" if m else "", ma])

    # --- console summary
    print("\nWhere the raw alt templates end up (all Horizon SNVs pooled):\n")
    colw = max(len(c) for c in cats) if cats else 10
    print("category".ljust(colw) + "".join(s[:22].rjust(24) for s in samples))
    for c in cats:
        line = c.ljust(colw)
        for s in samples:
            tot = sum(by_sample[s].values())
            line += (f"{by_sample[s][c]} ({100*by_sample[s][c]/tot:.0f}%)" if tot else "-").rjust(24)
        print(line)
    print("\nWrote:\n  " + "\n  ".join([p1, p2, p3, detail]))


def cmd_bed(a):
    vs = read_variants(a.variants)
    refs = pysam.AlignmentFile(a.bam, "rb", check_sq=False).references if a.bam else None
    with open(a.out, "w") as fh:
        for v in vs:
            chrom = match_contig(v["chrom"], refs) if refs else v["chrom"]
            if chrom is None:
                sys.exit(f"ERROR: contig {v['chrom']} not in {a.bam}")
            fh.write(f"{chrom}\t{v['pos']-1}\t{v['pos']}\t{v['name'].replace(' ', '_')}\n")
    print(f"wrote {len(vs)} SNVs to {a.out}", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("bed")
    b.add_argument("--variants", required=True)
    b.add_argument("--out", required=True)
    b.add_argument("--bam", help="BAM whose contig names (chr1 vs 1) the BED should use")
    t = sub.add_parser("trace")
    t.add_argument("--variants", required=True, help="Horizon truth table (CSV/TSV with chrom,pos,ref,alt[,gene,name])")
    t.add_argument("--sample", required=True)
    t.add_argument("--grouped", required=True, help="UMI-grouped BAM (subset), sorted + indexed")
    t.add_argument("--sscs-called", help="SSCS consensus before FilterConsensusReads (…_SSCS_mapped.bam subset), indexed")
    t.add_argument("--sscs-final", required=True, help="…_SSCS_mapped_filtered_clipped_sorted.bam")
    t.add_argument("--dcs-called", help="DCS consensus before filtering (…_DCS_mapped.bam subset), indexed")
    t.add_argument("--dcs-final", required=True, help="…_DCS_mapped_filtered_clipped_sorted.bam")
    t.add_argument("--outdir", required=True)
    t.add_argument("--min-mapq", type=int, default=20)
    t.add_argument("--min-bq", type=int, default=20)
    t.add_argument("--min-ss", type=int, default=MIN_SS_DEFAULT, help="min reads per strand family (default 3)")
    t.add_argument("--max-depth", type=int, default=0, help="stop after N raw reads per position (0 = no limit)")
    s = sub.add_parser("summarise")
    s.add_argument("--indir", required=True)
    a = ap.parse_args()
    {"trace": cmd_trace, "summarise": cmd_summarise, "bed": cmd_bed}[a.cmd](a)


if __name__ == "__main__":
    main()
