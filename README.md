# Duplex Sequencing Pipeline

Snakemake pipeline that turns paired-end duplex sequencing reads with UMIs into
single-strand (SSCS) and duplex (DCS) consensus BAMs at several consensus
stringencies, and calls variants on each with two callers:

- **VarDict-Java** (original path) → ANNOVAR → annotated TSV
- **duplexcaller** (custom, transparent caller) → ANNOVAR → `final.tsv` with every
  allele, its evidence, and the protein change

```
FASTQ → unmapped BAM (UMI in RX) → align + ZipperBams → GroupReadsByUmi
      → consensus called once, permissively (SSCS + DCS)
      → for each filter setting (e.g. 6-3-3, 4-2-2, 3-2-1, 2-1-1) and layer (SSCS, DCS):
            filter → clip → BAM
              ├─ VarDict-Java → ANNOVAR → annotated VCF + TSV
              ├─ duplexcaller  → ANNOVAR → annotated VCF + final.tsv   (optional, see below)
              └─ samtools mpileup at panel positions
```

**Read structure:** `3M2S+T 3M2S+T` (3 nt UMI + 2 nt skip + template, both reads)

**Tools:** fgbio, BWA MEM, samtools, Picard, VarDict-Java, duplexcaller, ANNOVAR,
bcftools, FastQC

---

## Where things live

| What | Repository | Branch |
|---|---|---|
| This pipeline (Snakefile, configs, envs) | `github.com/zuzanna-kra/duplex_pipeline` | `duplexcaller` (custom caller rules; to be merged into `main`) |
| The custom caller **duplexcaller** (source, tests, docs) | `github.com/zuzanna-kra/variant_caller` | `restructure` (to be merged into `main`) |

The caller is developed in its own repository. The pipeline does not install
it: it runs it straight from a local clone (`PYTHONPATH=<clone>/src python -m
duplexcaller …`), whose path you give in the config (`duplexcaller.src`). To
update the caller, `git pull` in that clone; the pipeline picks up the new code
on the next run.

---

## Quick start

### 1. Prerequisites

- Conda (or Mamba)
- Snakemake ≥ 7.8 (tested with 8.30): `conda create -n snakemake -c conda-forge -c bioconda snakemake`
- ANNOVAR, installed manually (see [ANNOVAR setup](#annovar-setup))
- hg38 reference with BWA index, `.fai` and `.dict`
- For the custom caller: a clone of the caller repository
  ```bash
  git clone https://github.com/zuzanna-kra/variant_caller.git
  cd variant_caller && git switch restructure
  # optional check: conda env create -f environment.yml && conda activate duplexcaller && PYTHONPATH=src pytest -q
  ```

### 2. Clone this repository

```bash
git clone https://github.com/zuzanna-kra/duplex_pipeline.git
cd duplex_pipeline
git switch duplexcaller
```

### 3. Sample sheet

One row per sample in a TSV (e.g. `config/samples.tsv`; lines starting with `#` are ignored):

| Column | Description |
|---|---|
| `sample_udi` | Unique ID (sample + UDI); used in every output path |
| `udi` | UDI barcode name (e.g. `xGenUDI1`) |
| `sample` | Sample name |
| `library` | Library / sequencing run ID |
| `fq1`, `fq2` | Absolute paths to R1/R2 FASTQs (may be empty when starting from existing BAMs) |

### 4. Config files

Settings are split over two (or more) YAML files given together with `--configfile`:

- **`config/config.local.yaml`** — server-specific paths (reference, panel BED,
  tool jars, ANNOVAR, the caller clone). Not committed (git-ignored). Start from
  `config/config.yaml`, which documents every key.
- **an experiment config**, e.g. `config/config_horizon_umifix.yaml` — sample
  sheet, `results_dir`, consensus settings, filter settings, and optional start
  points (`existing_unmapped_bam`, `existing_grouped_bam`).

Later files override earlier ones key by key.

Key settings:

| Key | Meaning |
|---|---|
| `results_dir` | Output folder (keep outside the repo) |
| `reference.fasta` | hg38 FASTA |
| `snv_bed` | Panel BED (regions for calling and mpileup) |
| `tools.fgbio_jar`, `tools.picard_jar` | Jar paths (find them with `find .snakemake/conda -name "fgbio.jar"` after a first `--use-conda` run) |
| `tools.annovar_annotate`, `tools.annovar_humandb` | ANNOVAR `table_annovar.pl` and `humandb/` |
| `consensus.*` | Permissive consensus calling and quality settings (see `config/config.yaml`) |
| `filter_settings` | Consensus filter settings to sweep, e.g. `"6-3-3": {duplex: "6 3 3", sscs: 3}` |
| `existing_unmapped_bam` / `existing_grouped_bam` | Start from existing BAMs instead of FASTQs (pattern with `{sample_udi}`, `{sample}`) |
| `resources` | Override per-rule threads / RAM (defaults in `DEFAULT_RESOURCES`, `workflow/Snakefile`) |
| `duplexcaller.*` | Custom caller, see below |

### 5. Run

Always dry-run first (`-n`) and check the job list:

```bash
snakemake -s workflow/Snakefile \
    --configfile config/config.local.yaml config/config_horizon_umifix.yaml \
    --use-conda --conda-frontend conda --rerun-triggers mtime -n
```

then run with resource limits:

```bash
snakemake -s workflow/Snakefile \
    --configfile config/config.local.yaml config/config_horizon_umifix.yaml \
    --use-conda --conda-frontend conda --rerun-triggers mtime \
    --cores 80 --resources mem_mb=190000 annovar=4
```

Notes:
- **Server limits (RADT-OUILLE):** ≤ 80 threads and < 200 GB RAM in total at
  once. `--cores` / `--resources mem_mb` are those totals; each rule declares its
  own threads and `mem_mb`. `annovar=4` caps concurrent ANNOVAR jobs. Split the
  totals if several runs share the server.
- **`--rerun-triggers mtime`:** only rerun steps whose inputs are newer than their
  outputs. Without it, Snakemake also reruns steps whose rule code or parameters
  changed, which after a pipeline update can mean redoing alignment and
  consensus for every sample.
- **Targets go first.** `--configfile` and `--resources` accept several values
  and swallow anything after them, so put file targets right after
  `snakemake`, e.g. `snakemake <file1> <file2> -s workflow/Snakefile …`.
- **Running only some rules:** `--until <rule> [<rule> …]` runs only what is
  needed for those rules (e.g. only the custom caller, without VarDict):
  `--until merge_duplexcaller_annovar compress_index_duplexcaller_annotated_vcf`.
- The first run builds the conda environments (5–10 min).
- Use `tmux`/`screen` for long runs.

---

## Custom caller (duplexcaller)

A transparent SNV caller for duplex consensus BAMs. It reports **every**
non-reference allele with its evidence and applies no filters on evidence;
the only exclusions are MAPQ, base quality and SAM flags, and every excluded
read is counted by reason. Thresholds such as AS−XS only add alert labels.
Full documentation (design, every output column, parameters) is in the caller
repository's README; `duplexcaller schema alleles` prints every column.

**Status:** SNVs only (indels and MNVs to come). Runs alongside VarDict.

### Enable it

Add to `config/config.local.yaml`:

```yaml
duplexcaller:
  enabled: true
  src: "/path/to/variant_caller/src"   # the clone's src/ folder
  params: ""                           # optional YAML of caller parameters ("" = defaults)
```

Default caller parameters (print with `PYTHONPATH=<src> python -m duplexcaller params`):
MAPQ ≥ 20, base quality ≥ 20, SAM flags 0xF04 excluded (unmapped, secondary,
QC-fail, duplicate, supplementary); AS−XS alert below 50.

### Rules

| Rule | Does |
|---|---|
| `duplexcaller_pileup` | Counts every allele at every panel position of the clipped BAM |
| `annotate_duplexcaller` | ANNOVAR (`table_annovar.pl -vcfinput`, same databases as the VarDict rule) on the caller's VCF |
| `compress_index_duplexcaller_annotated_vcf` | bgzip + tabix of the annotated VCF |
| `merge_duplexcaller_annovar` | Joins ANNOVAR's annotation onto the allele table → `final.tsv` |

Single-threaded, ~4 GB per job; ~3 min per DCS file on the Horizon panel, longer for SSCS.

### Outputs (`{setting}/DuplexCaller/`)

| File | Content |
|---|---|
| `{sample_udi}_{layer}.final.tsv` | **Main output.** Every allele + evidence + alert labels + gene, protein/cDNA change, consequence, gnomAD, ClinVar, COSMIC |
| `{sample_udi}_{layer}.alleles.tsv` | Every allele with its evidence (before annotation) |
| `{sample_udi}_{layer}.depth.tsv.gz` | Every panel position: fragment counts per base and every exclusion by reason |
| `{sample_udi}_{layer}.reads.tsv.gz` | One row per read supporting a non-reference base |
| `{sample_udi}_{layer}.vcf.gz` (+ `.tbi`) | The allele table as VCF (input to ANNOVAR); FILTER is `.` (nothing filtered) |
| `{sample_udi}_{layer}.hg38_multianno.vcf.gz` (+ `.tbi`) | ANNOVAR-annotated VCF |
| `{sample_udi}_{layer}.hg38_multianno.txt` | ANNOVAR table |
| `{sample_udi}_{layer}.run.yaml` | Caller version, git commit, command, parameters, run summary |

### Run only the custom caller

```bash
snakemake -s workflow/Snakefile \
    --configfile config/config.local.yaml config/config_horizon_umifix.yaml \
    --use-conda --conda-frontend conda --rerun-triggers mtime \
    --until merge_duplexcaller_annovar compress_index_duplexcaller_annotated_vcf \
    --cores 40 --resources mem_mb=150000 annovar=4
```

---

## Output structure

```
{results_dir}/
└── {sample_udi}/
    ├── fastqc/
    ├── Metrics_and_images/            insert size, family size, UMI metrics
    ├── {sample_udi}_unmapped_bam.bam
    ├── {sample_udi}_mapped_merged_sorted.bam
    ├── {sample_udi}_UMI_grouped.bam
    ├── consensus_permissive/          SSCS/DCS consensus, called once
    └── {setting}/                     one folder per filter setting, e.g. 6-3-3
        ├── SSCS/  {sample_udi}_SSCS_mapped_filtered_clipped_sorted.bam (+ .bai)
        ├── DCS/   {sample_udi}_DCS_mapped_filtered_clipped_sorted.bam (+ .bai)
        ├── Metrics/       flagstat, overlap metrics
        ├── Variants/      VarDict path
        │   ├── {sample_udi}_{layer}_variants.vcf.gz
        │   ├── {sample_udi}_{layer}_variants_annotated.hg38_multianno.vcf.gz
        │   ├── {sample_udi}_{layer}_variants_annotated.tsv
        │   └── {sample_udi}_{layer}_pileup.txt
        └── DuplexCaller/  custom caller path (if enabled; see above)
```

Settings that share the same SSCS min-reads (e.g. 3-2-1 and 2-1-1) give identical
SSCS BAMs; the later one is a symlink to the first (`alias_duplicate_sscs`).

### VarDict path: what to know

- The VarDict TSV (`*_variants_annotated.tsv`, no header) has: CHROM, POS, REF,
  ALT, DP, VD, AF, Func.refGene, Gene.refGene, ExonicFunc.refGene,
  gnomAD_genome_ALL, CLNSIG. It has **no FILTER column**: VarDict's own filter
  labels (PASS, v2, p8, pSTD, NM5.25, …) are only in the VCF.
- `var2vcf_valid.pl` is run without `-A`, so only the highest-frequency allele
  per position is reported.
- VarDict runs with local realignment off (`-k 0`).

---

## Shell script alternative

`shell/data_processing_v13.sh` (current; v12 kept for reference) runs the
VarDict path for one sample without Snakemake, at one consensus setting. Edit
its `USER CONFIGURATION` block, activate an environment with all tools, then:

```bash
bash shell/data_processing_v13.sh -s sample_indexes.csv \
    -f /path/to/R1.fq.gz -g /path/to/R2.fq.gz -r my_run \
    [-t threads] [-m min_reads_SSCS] [-a min_reads_DCS] [-c clip_bases]
```

Run `bash shell/data_processing_v13.sh` without arguments for all options.

The shell script does not include the custom caller.

---

## ANNOVAR setup

ANNOVAR is not on conda (licensing). Register and download it from
https://annovar.openbioinformatics.org, then build the databases:

```bash
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad_genome humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar clinvar_20240917 humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar exac03 humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar cosmic102_coding humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar cosmic101_noncoding humandb/
```

The protocol list is defined once in the Snakefile (`ANNOVAR_PROTOCOL`) and
used by both callers' annotation rules.

---

## Environments

Per-rule conda environments in `workflow/envs/` (created automatically with `--use-conda`):

| File | Tools |
|---|---|
| `fgbio.yaml` | fgbio, bwa, samtools |
| `picard.yaml` | picard, samtools |
| `vardict.yaml` | vardict-java, R |
| `bcftools.yaml` | bcftools, bgzip, tabix |
| `fastqc.yaml` | fastqc |
| `duplexcaller.yaml` | python, pysam, numpy, pandas, pyyaml (for the custom caller) |

For the shell script, one environment with everything:

```bash
conda create -n duplex-seq -c conda-forge -c bioconda \
    fgbio bwa samtools picard vardict-java bcftools fastqc r-base
```

---

## Citation

Please cite the underlying tools:

- **fgbio**: https://github.com/fulcrumgenomics/fgbio
- **BWA**: Li H. & Durbin R. (2009) *Bioinformatics* 25:1754-60
- **VarDict**: Lai Z. et al. (2016) *Nucleic Acids Research* 44:e108
- **ANNOVAR**: Wang K. et al. (2010) *Nucleic Acids Research* 38:e164
- **pysam / htslib** (used by duplexcaller): https://github.com/pysam-developers/pysam
