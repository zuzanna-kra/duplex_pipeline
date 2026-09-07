# Duplex Sequencing Pipeline

A pipeline for processing duplex sequencing data — paired-end reads with UMI barcodes — into high-confidence SSCS (Single-Strand Consensus Sequence) and DCS (Duplex Consensus Sequence) variant calls with ANNOVAR annotation.

## Overview

```
FASTQs → unmapped BAM → UMI extraction → alignment → UMI grouping
       → SSCS consensus → filter → clip → variant calling (VarDict) → ANNOVAR annotation → TSV
       → DCS consensus  → filter → clip → variant calling (VarDict) → ANNOVAR annotation → TSV
```

**Read structure:** `3M2S+T 3M2S+T` (3 nt UMI + 2 nt skip + template, both reads)

**Tools used:** fgbio, BWA MEM, samtools, Picard, VarDict-Java, ANNOVAR, bcftools, FastQC

---

## Quick start (Snakemake — recommended)

### 1. Prerequisites

- [Conda](https://docs.conda.io/en/latest/miniconda.html) or [Mamba](https://github.com/mamba-org/mamba)
- Snakemake ≥ 7.0: `conda install -n base -c conda-forge -c bioconda snakemake`
- ANNOVAR (see [ANNOVAR setup](#annovar-setup) below — manual install required)
- hg38 reference genome with BWA index and samtools `.fai`

### 2. Clone the repo

```bash
git clone https://github.com/YOUR_USERNAME/duplex-seq-pipeline.git
cd duplex-seq-pipeline
```

### 3. Fill in the sample sheet

Edit `config/samples.tsv` — one row per sample:

| Column | Description |
|---|---|
| `sample_udi` | Unique ID combining sample name and UDI (used as wildcard throughout) |
| `udi` | UDI barcode name (e.g. `xGenUDI1`) |
| `sample` | Sample name |
| `library` | Library / sequencing run ID |
| `fq1` | Absolute path to R1 FASTQ |
| `fq2` | Absolute path to R2 FASTQ |

### 4. Edit the config

Copy and edit `config/config.yaml` — all paths marked `# <-- CHANGE THIS` must be updated:

```bash
cp config/config.yaml config/config.local.yaml   # optional: keep a local copy
```

Key settings:
- `results_dir` — where outputs go (keep outside the repo)
- `reference.fasta` — path to hg38.fa
- `snv_bed` — BED file of target regions
- `tools.fgbio_jar` / `tools.picard_jar` — jar file paths (see note below)
- `tools.annovar_annotate` / `tools.annovar_humandb` — ANNOVAR paths

**Finding jar paths after `--use-conda`:** Run once with `--use-conda`, then:
```bash
find .snakemake/conda -name "fgbio.jar"   # copy this path to fgbio_jar
find .snakemake/conda -name "picard.jar"  # copy this path to picard_jar
```

### 5. Run

```bash
snakemake --use-conda --conda-frontend conda --cores 16
```

On first run, Snakemake builds the conda environments (5–10 min). Subsequent runs reuse them.

**Dry run** (check the DAG without executing):
```bash
snakemake --use-conda --cores 16 -n
```

---

## Shell script alternative

A self-contained bash script is provided for users who prefer not to use Snakemake.

```bash
shell/data_processing_v12.sh
```

### Setup

1. Edit the `USER CONFIGURATION` block at the top of the script (clearly marked)
2. Activate your conda environment containing all tools
3. Run:

```bash
conda activate YOUR_ENV
bash shell/data_processing_v12.sh \
    --sample MySample_xGenUDI1 \
    --fq1 /path/to/R1.fq.gz \
    --fq2 /path/to/R2.fq.gz \
    --outdir /path/to/results
```

The shell script requires a single conda environment with all tools installed. See [Environment setup](#environment-setup) below.

---

## ANNOVAR setup

ANNOVAR is not available via conda due to licensing restrictions. You must register and download it manually from https://annovar.openbioinformatics.org.

After downloading, build the required databases:

```bash
# Required databases
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad_genome humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar clinvar_20240917 humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar exac03 humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar cosmic102_coding humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar cosmic101_noncoding humandb/
```

---

## Environment setup

### Snakemake (automatic with `--use-conda`)

Per-rule conda environments are defined in `workflow/envs/`:

| File | Tools |
|---|---|
| `fgbio.yaml` | fgbio, bwa, samtools |
| `picard.yaml` | picard, samtools |
| `vardict.yaml` | vardict-java, R |
| `bcftools.yaml` | bcftools |
| `fastqc.yaml` | fastqc |

Snakemake creates and manages these automatically.

### Shell script (manual)

Create a conda environment with all tools:

```bash
conda create -n duplex-seq \
    -c conda-forge -c bioconda \
    fgbio bwa samtools picard vardict-java bcftools fastqc r-base
conda activate duplex-seq
```

---

## Output structure

```
results/
└── {sample_udi}/
    ├── fastqc/                          # FastQC reports
    ├── Metrics_and_images/
    │   ├── *_insert_size_metrics.txt
    │   ├── *_insert_size_histogram.pdf
    │   ├── *_family_size_histogram.txt
    │   ├── *_SSCS_UMI_metrics.txt
    │   ├── *_SSCS_overlap_metrics.txt
    │   └── *_DCS_overlap_metrics.txt
    ├── SSCS_files_MUFS{call_min}/
    │   └── *_SSCS_mapped_filtered_clipped_sorted.bam
    ├── DCS_files_MUFS{call_min}/
    │   └── *_DCS_mapped_filtered_clipped_sorted.bam
    └── Variants/
        ├── *_SSCS_variants.vcf.gz
        ├── *_SSCS_variants_annotated.hg38_multianno.vcf.gz
        ├── *_SSCS_variants_annotated.tsv       # ← main output
        ├── *_SSCS_pileup.txt
        ├── *_DCS_variants.vcf.gz
        ├── *_DCS_variants_annotated.hg38_multianno.vcf.gz
        ├── *_DCS_variants_annotated.tsv        # ← main output
        └── *_DCS_pileup.txt
```

The annotated TSV files are the primary outputs, with columns: CHROM, POS, REF, ALT, DP, VD, AF, Func.refGene, Gene.refGene, ExonicFunc.refGene, gnomAD_genome_ALL, CLNSIG.

---

## Key parameters

| Parameter | Default | Description |
|---|---|---|
| `consensus.call_min` | 3 | Min reads per strand for SSCS |
| `consensus.call_min_duplex` | 3 | Min SSCS family size for DCS |
| `consensus.base_qual_min` | 20 | Min input base quality |
| `consensus.cons_qual_min` | 45 | Min consensus base quality |
| `consensus.max_base_error_rate` | 0.2 | Max per-base error rate |
| `consensus.clip_bases` | 3 | Bases hard-clipped from each read end |
| `variants.all_freq` | 0.000001 | Min allele frequency for VarDict |

---

## Citation

If you use this pipeline, please cite the underlying tools:

- **fgbio**: https://github.com/fulcrumgenomics/fgbio
- **BWA**: Li H. & Durbin R. (2009) *Bioinformatics* 25:1754-60
- **VarDict**: Lai Z. et al. (2016) *Nucleic Acids Research* 44:e108
- **ANNOVAR**: Wang K. et al. (2010) *Nucleic Acids Research* 38:e164
