#!/usr/bin/env python3
"""Build samples.tsv for the duplex sequencing Snakemake pipeline."""

import argparse
import csv
import sys
from pathlib import Path
import os

# Import sample_reader from the same directory
sys.path.insert(0, str(Path(__file__).parent))
import sample_reader


def collect_fastq_r1_files(inputs: list[str]) -> list[Path]:
    """
    Given a list of paths (files or directories), return a list of R1 FASTQs.
    For directories, glob *1.fq.gz inside.
    For files, include them directly if they match the R1 pattern.
    """
    r1_files = []
    for input_path in inputs:
        if os.path.isdir(input_path):
            # Glob *1.fq.gz inside the directory
            r1_files.extend(Path(fq) for fq in sorted(Path(input_path).glob("*1.fq.gz")))
        elif os.path.isfile(input_path) and input_path.endswith("1.fq.gz"):
            # Include the file directly
            r1_files.append(Path(input_path))
        else:
            # Skip unknown files
            pass
    return r1_files


def build_rows(r1_files: list[Path], csv_manifest: str) -> list[dict]:
    """
    For each R1 FASTQ, find the mate, extract metadata, return a list of dicts.
    Skip and warn on broken pairs or unknown UDIs.
    """
    rows = []
    for fastq_r1 in r1_files:
        # Use sample_reader to find the mate and extract metadata
        fastq_r2 = Path(str(fastq_r1).replace("1.fq.gz", "2.fq.gz"))
        if not fastq_r2:
            print(f"Warning: no mate found for {fastq_r1}, skipping", file=sys.stderr)
            continue
        sample_name = sample_reader.extract_sample_name(csv_manifest,fastq_r1)
        udi = sample_reader.extract_udi(fastq_r1)
        library_name = sample_reader.extract_library_name(fastq_r1)
        sample_name_with_udi = sample_reader.extract_sample_name_with_UDI(csv_manifest, fastq_r1)
        if not sample_name or not udi or not sample_name_with_udi:
            print(f"Warning: could not extract sample name, UDI or sample name with UDI from {fastq_r1}, skipping", file=sys.stderr)
            continue
        rows.append({
            "sample_udi": sample_name_with_udi,
            "udi": udi,
            "sample": sample_name,
            "library": library_name,
            "fq1": str(fastq_r1),
            "fq2": str(fastq_r2)
        })

    return rows


def write_tsv(rows: list[dict], output_path: str):
    """Write rows to a TSV with the standard column order."""
    fieldnames = ["sample_udi", "udi", "sample", "library", "fq1", "fq2"]
    with open(output_path, "w", newline="") as tsvfile:
        writer = csv.DictWriter(tsvfile, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+",
                        help="One or more FASTQ files or directories containing them")
    parser.add_argument("--csv", required=True,
                        help="CSV manifest mapping sample names to UDI")
    parser.add_argument("--out", required=True,
                        help="Output TSV path")
    args = parser.parse_args()

    r1_files = collect_fastq_r1_files(args.inputs)
    if not r1_files:
        sys.exit("Error: no R1 FASTQ files found in input paths")

    rows = build_rows(r1_files, args.csv)
    if not rows:
        sys.exit("Error: no valid sample pairs could be built")

    write_tsv(rows, args.out)
    print(f"Wrote {len(rows)} samples to {args.out}")


if __name__ == "__main__":
    main()