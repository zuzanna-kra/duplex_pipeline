#!/bin/bash

# Exit on error
set -euo pipefail

################################################################################
# === USER CONFIGURATION — EDIT THESE LINES BEFORE RUNNING ==================
################################################################################
# All paths that need to be set for your environment are listed here.
# Nothing else in this script should need to be changed.

# --- Conda environment ---
# Name of the conda environment that contains all pipeline tools
# (fgbio, bwa, samtools, fastqc, vardict-java, bcftools, R with teststrandbias)
CONDA_ENV_NAME="nanoseq"                          # <-- CHANGE IF NEEDED

# --- Tool paths (jar files) ---
# These are typically inside your conda environment's share/ directory
FGBIO_JAR="$HOME/miniconda3/envs/${CONDA_ENV_NAME}/share/fgbio/fgbio.jar"          # <-- CHANGE IF NEEDED
PICARD_JAR="$HOME/miniconda3/envs/${CONDA_ENV_NAME}/share/picard-3.4.0-0/picard.jar" # <-- CHANGE IF NEEDED

# --- ANNOVAR ---
ANNOVAR_ANNOTATE="$HOME/tools/annovar/table_annovar.pl"   # <-- CHANGE IF NEEDED
ANNOVAR_HUMANDB="$HOME/tools/annovar/humandb"             # <-- CHANGE IF NEEDED

# --- Reference genome ---
# Must have .fai index and BWA index files alongside it
REF="$HOME/references/hg38/hg38.fa"                      # <-- CHANGE IF NEEDED

# --- Panel BED file ---
# BED file defining the target regions for variant calling
SNV_BED="$HOME/panels/MDS_screening.bed"                  # <-- CHANGE IF NEEDED

# --- Helper script ---
# sample_reader.py maps sample names to UDI indices from the CSV manifest
SAMPLE_READER="$HOME/sequencing/scripts/sequencing_pipeline/sample_reader.py" # <-- CHANGE IF NEEDED

################################################################################
# === PIPELINE PARAMETERS — DEFAULTS (override with flags) ===================
################################################################################

CSV=""
FQ1=""
FQ2=""
THREADS=16
CALLMIN=3
CALLMINDUPLEX=3
ALLFREQ=0.000001
BASEQUALMIN=20
MAPQUAL=23
MAXBASEERRORRATE=0.2
CONSQUALMIN=45
CLIPBASES=3
RUN_NAME=""

usage() {
  cat << EOF
Usage: $0 [options] -s <sample_indexes.csv> -f <read1.fq> -g <read2.fq> -r <run_name>

Duplex sequencing pipeline: FASTQ → SSCS/DCS BAMs → VarDict variant calls → ANNOVAR annotation.

OPTIONS:
    -s <sample_indexes.csv>  CSV mapping sample names to UDI indices (required)
    -f <read1.fq>            R1 FASTQ file (required)
    -g <read2.fq>            R2 FASTQ file (required)
    -r <run_name>            Name used for log/error files (required)
    -t <cpus>                Threads (default: $THREADS)
    -m <min-reads-SSCS>      Min reads to form an SSCS read (default: $CALLMIN)
    -a <min-reads-DCS>       Min SSCS family size for DCS calling (default: $CALLMINDUPLEX)
    -b <min-base-quality>    Min base quality for consensus calling (default: $BASEQUALMIN)
    -e <max-base-error-rate> Max base error rate for consensus calling (default: $MAXBASEERRORRATE)
    -q <min-cons-quality>    Min consensus quality (default: $CONSQUALMIN)
    -c <clip-bases>          Bases to clip from 5' and 3' ends (default: $CLIPBASES)
EOF
}

################################################################################
# === PARSE ARGUMENTS =========================================================
################################################################################

while getopts "s:m:a:b:c:t:f:g:r:" opt; do
  case $opt in
    s) CSV="$OPTARG" ;;
    m) CALLMIN="$OPTARG" ;;
    a) CALLMINDUPLEX="$OPTARG" ;;
    b) BASEQUALMIN="$OPTARG" ;;
    c) CLIPBASES="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    f) FQ1="$OPTARG" ;;
    g) FQ2="$OPTARG" ;;
    r) RUN_NAME="$OPTARG" ;;
    *) usage; exit 1 ;;
  esac
done

# Mandatory argument checks
if [[ -z "$CSV" ]];      then usage; echo; echo "Error: -s (CSV manifest) is required."; exit 1; fi
if [[ -z "$FQ1" ]];      then usage; echo; echo "Error: -f (R1 FASTQ) is required."; exit 1; fi
if [[ -z "$FQ2" ]];      then usage; echo; echo "Error: -g (R2 FASTQ) is required."; exit 1; fi
if [[ -z "$RUN_NAME" ]]; then usage; echo; echo "Error: -r (run name) is required."; exit 1; fi

# Redirect stdout and stderr to log files
LOG_FILE="${RUN_NAME}_analysis.log"
ERROR_FILE="${RUN_NAME}_analysis.err"
exec 1> >(tee -a "$LOG_FILE") 2> >(tee -a "$ERROR_FILE" >&2)

################################################################################
# === ACTIVATE CONDA ENVIRONMENT ==============================================
################################################################################

set +u  # Conda hooks need this temporarily
if ! command -v conda &>/dev/null; then
    if [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
        source "$HOME/miniconda3/etc/profile.d/conda.sh"
    elif [[ -f "/opt/conda/etc/profile.d/conda.sh" ]]; then
        source "/opt/conda/etc/profile.d/conda.sh"
    else
        echo "Error: conda not found." >&2; exit 1
    fi
fi
eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV_NAME"
set -u

################################################################################
# === RETRIEVE SAMPLE METADATA ================================================
################################################################################

echo "Extracting sample name from UDI index"
SAMPLE=$(python "$SAMPLE_READER" --mode sample --samples_csv_file "$CSV" --fastq "$FQ1")
sample_name_UDI=$(python "$SAMPLE_READER" --mode sample_UDI --samples_csv_file "$CSV" --fastq "$FQ1")
library_name=$(python "$SAMPLE_READER" --mode library --fastq "$FQ1")
echo "Sample: $SAMPLE  |  UDI label: $sample_name_UDI  |  Library: $library_name"

################################################################################
# === OUTPUT FILE PATHS =======================================================
################################################################################

unmapped_bam="$sample_name_UDI/${SAMPLE}_unmapped_bam.bam"
unmapped_bam_UMI="$sample_name_UDI/${SAMPLE}_unmapped_bam_UMI.bam"
mapped_merged_BAM="$sample_name_UDI/${SAMPLE}_mapped_merged.bam"
mapped_merged_BAM_sorted="$sample_name_UDI/${SAMPLE}_mapped_merged_sorted.bam"
grouped_BAM="$sample_name_UDI/${SAMPLE}_UMI_grouped.bam"

# SSCS
SSCS_bam_unmapped="$sample_name_UDI/SSCS_files_MUFS${CALLMIN}/${SAMPLE}_SSCS_unmapped.bam"
SSCS_bam_mapped="$sample_name_UDI/SSCS_files_MUFS${CALLMIN}/${SAMPLE}_SSCS_mapped.bam"
SSCS_bam_filtered="$sample_name_UDI/SSCS_files_MUFS${CALLMIN}/${SAMPLE}_SSCS_mapped_filtered.bam"
SSCS_clipped="$sample_name_UDI/SSCS_files_MUFS${CALLMIN}/${SAMPLE}_SSCS_mapped_filtered_clipped.bam"
SSCS_clipped_sorted="$sample_name_UDI/SSCS_files_MUFS${CALLMIN}/${SAMPLE}_SSCS_mapped_filtered_clipped_sorted.bam"
pileup_SSCS="$sample_name_UDI/Variant_calling_MUFS${CALLMIN}/${SAMPLE}_SSCS.pileup"
VarDict_SSCS_vcf="$sample_name_UDI/Variant_calling_MUFS${CALLMIN}/${SAMPLE}_SSCS_VarDictJava.vcf.gz"
annovar_annotated_VarDict_SSCS="$sample_name_UDI/Variant_calling_MUFS${CALLMIN}/${SAMPLE}_SSCS_variants_annovar"
annotated_VarDict_SSCS_tsv="$sample_name_UDI/Variant_calling_MUFS${CALLMIN}/${SAMPLE}_SSCS_annotated_variants.tsv"

# DCS
DCS_bam_unmapped="$sample_name_UDI/DCS_files_MUFS${CALLMINDUPLEX}/${SAMPLE}_DCS_unmapped.bam"
DCS_bam_mapped="$sample_name_UDI/DCS_files_MUFS${CALLMINDUPLEX}/${SAMPLE}_DCS_mapped.bam"
DCS_bam_filtered="$sample_name_UDI/DCS_files_MUFS${CALLMINDUPLEX}/${SAMPLE}_DCS_mapped_filtered.bam"
DCS_clipped="$sample_name_UDI/DCS_files_MUFS${CALLMINDUPLEX}/${SAMPLE}_DCS_mapped_filtered_clipped.bam"
DCS_clipped_sorted="$sample_name_UDI/DCS_files_MUFS${CALLMINDUPLEX}/${SAMPLE}_DCS_mapped_filtered_clipped_sorted.bam"
pileup_DCS="$sample_name_UDI/Variant_calling_MUFS${CALLMINDUPLEX}/${SAMPLE}_DCS.pileup"
VarDict_DCS_vcf="$sample_name_UDI/Variant_calling_MUFS${CALLMINDUPLEX}/${SAMPLE}_DCS_VarDictJava.vcf.gz"
annovar_annotated_VarDict_DCS="$sample_name_UDI/Variant_calling_MUFS${CALLMINDUPLEX}/${SAMPLE}_DCS_variants_annovar"
annotated_VarDict_DCS_tsv="$sample_name_UDI/Variant_calling_MUFS${CALLMINDUPLEX}/${SAMPLE}_DCS_annotated_variants.tsv"

# Metrics
insert_size_metrics="$sample_name_UDI/Metrics_and_images/${SAMPLE}_insert_size_metrics.txt"
insert_size_histogram="$sample_name_UDI/Metrics_and_images/${SAMPLE}_insert_size_histogram.pdf"
SSCS_overlap_metrics="$sample_name_UDI/Metrics_and_images/${SAMPLE}_SSCS_overlap_metrics.txt"
DCS_overlap_metrics="$sample_name_UDI/Metrics_and_images/${SAMPLE}_DCS_overlap_metrics.txt"
family_size_histogram="$sample_name_UDI/Metrics_and_images/${SAMPLE}_family_size_histogram.txt"
SSCS_UMI_metrics="$sample_name_UDI/Metrics_and_images/${SAMPLE}_SSCS_UMI_metrics.txt"
TEMP="$sample_name_UDI/TEMP"

################################################################################
# === CREATE OUTPUT DIRECTORIES ===============================================
################################################################################

mkdir -p \
    "$sample_name_UDI" \
    "$TEMP" \
    "$sample_name_UDI/SSCS_files_MUFS${CALLMIN}" \
    "$sample_name_UDI/DCS_files_MUFS${CALLMINDUPLEX}" \
    "$sample_name_UDI/Variant_calling_MUFS${CALLMIN}" \
    "$sample_name_UDI/Variant_calling_MUFS${CALLMINDUPLEX}" \
    "$sample_name_UDI/Metrics_and_images"

################################################################################
# === PIPELINE ================================================================
################################################################################

echo "=== [1/12] FastQC on raw FASTQs ==="
fastqc "$FQ1" "$FQ2" --outdir="$sample_name_UDI" --threads "$THREADS" --quiet

echo "=== [2/12] FASTQ → unmapped BAM (fgbio FastqToBam) ==="
java -Xmx32G -jar "$FGBIO_JAR" --compression 1 --async-io FastqToBam \
  --input "$FQ1" "$FQ2" \
  --read-structures 3M2S+T 3M2S+T \
  --sample "$SAMPLE" \
  --library "$library_name" \
  --platform-unit FLOWCELL.LANE \
  --output "$unmapped_bam"

echo "=== [3/12] Extract UMIs (fgbio ExtractUmisFromBam) ==="
java -Xmx32G -jar "$FGBIO_JAR" --tmp-dir "$TEMP" \
  ExtractUmisFromBam \
  --input "$unmapped_bam" \
  --output "$unmapped_bam_UMI" \
  --read-structure 3M2S+T \
  --read-structure 3M2S+T \
  --molecular-index-tags ZA ZB \
  --single-tag RX

echo "=== [4/12] Align (BWA MEM) + merge (fgbio ZipperBams) ==="
samtools fastq -n -1 "$TEMP/r1.fq.gz" -2 "$TEMP/r2.fq.gz" "$unmapped_bam_UMI"
bwa mem -t "$THREADS" "$REF" "$TEMP/r1.fq.gz" "$TEMP/r2.fq.gz" \
  | samtools view -b -o "$TEMP/aligned.bam"
java -Xmx32G -jar "$FGBIO_JAR" --tmp-dir "$TEMP" ZipperBams \
  --unmapped "$unmapped_bam_UMI" \
  --input "$TEMP/aligned.bam" \
  --output "$mapped_merged_BAM" \
  --ref "$REF"

echo "=== [5/12] Sort by template coordinate ==="
samtools sort --template-coordinate -@ "$THREADS" -o "$mapped_merged_BAM_sorted" "$mapped_merged_BAM"

echo "=== [5b] Insert size metrics (Picard) ==="
java -jar "$PICARD_JAR" CollectInsertSizeMetrics \
  I="$mapped_merged_BAM_sorted" O="$insert_size_metrics" H="$insert_size_histogram"

echo "=== [6/12] Group reads by UMI (fgbio GroupReadsByUmi) ==="
java -Xmx32G -jar "$FGBIO_JAR" --compression 1 --async-io GroupReadsByUmi \
  --input "$mapped_merged_BAM_sorted" \
  --strategy Paired \
  --edits 1 \
  --output "$grouped_BAM" \
  --family-size-histogram "$family_size_histogram" \
  --grouping-metrics "$SSCS_UMI_metrics" \
  --threads "$THREADS"

##############################################################################
# SSCS
##############################################################################

echo "=== [7/12] Call SSCS (min reads = $CALLMIN) ==="
java -Xmx32G -jar "$FGBIO_JAR" --compression 1 --async-io CallMolecularConsensusReads \
  --input "$grouped_BAM" \
  --output "$SSCS_bam_unmapped" \
  --min-reads "$CALLMIN" \
  --min-input-base-quality "$BASEQUALMIN" \
  --threads "$THREADS"

echo "=== [8/12] Remap + merge SSCS ==="
samtools fastq "$SSCS_bam_unmapped" \
  | bwa mem -t "$THREADS" -p -K 150000000 -Y "$REF" - \
  | java -Xmx32G -jar "$FGBIO_JAR" --compression 1 --async-io ZipperBams \
      --unmapped "$SSCS_bam_unmapped" \
      --ref "$REF" \
      --tags-to-reverse Consensus \
      --tags-to-revcomp Consensus \
      --output "$SSCS_bam_mapped"

echo "=== [9/12] Filter + sort SSCS ==="
java -Xmx32G -jar "$FGBIO_JAR" --compression 0 FilterConsensusReads \
  --input "$SSCS_bam_mapped" \
  --output /dev/stdout \
  --ref "$REF" \
  --min-reads "$CALLMIN" \
  --min-base-quality "$CONSQUALMIN" \
  --max-base-error-rate "$MAXBASEERRORRATE" \
  | samtools sort -n --threads "$THREADS" -o "$SSCS_bam_filtered" -

echo "=== [10/12] Clip overlapping reads + ends (SSCS) ==="
java -Xmx32G -jar "$FGBIO_JAR" --compression 0 ClipBam \
  --input "$SSCS_bam_filtered" \
  --output "$SSCS_clipped" \
  --ref "$REF" \
  --clipping-mode Hard \
  --clip-overlapping-reads true \
  --read-one-five-prime "$CLIPBASES" \
  --read-one-three-prime "$CLIPBASES" \
  --read-two-five-prime "$CLIPBASES" \
  --read-two-three-prime "$CLIPBASES" \
  --metrics "$SSCS_overlap_metrics"

echo "=== Sort + index SSCS ==="
samtools sort -@ "$THREADS" -o "$SSCS_clipped_sorted" "$SSCS_clipped"
samtools index "$SSCS_clipped_sorted"

echo "=== SSCS variant calling (VarDictJava) ==="
vardict-java \
  -G "$REF" -N "$SAMPLE" -b "$SSCS_clipped_sorted" \
  -z 1 -c 1 -S 2 -E 3 -g 4 \
  -f "$ALLFREQ" -r 1 -Q "$BASEQUALMIN" -q 0 -o 0 -k 0 -y \
  -th "$THREADS" "$SNV_BED" \
  | teststrandbias.R \
  | var2vcf_valid.pl -N "$SAMPLE" -E -f "$ALLFREQ" \
  | bgzip -c > "$VarDict_SSCS_vcf"

echo "=== SSCS annotation (ANNOVAR) ==="
perl "$ANNOVAR_ANNOTATE" "$VarDict_SSCS_vcf" "$ANNOVAR_HUMANDB/" \
  -buildver hg38 -vcfinput \
  -out "$annovar_annotated_VarDict_SSCS" \
  -remove \
  -protocol refGene,cosmic102_coding,cosmic101_noncoding,exac03,gnomad_genome,clinvar_20240917 \
  -operation g,f,f,f,f,f \
  -nastring . -thread "$THREADS"

bgzip -c "${annovar_annotated_VarDict_SSCS}.hg38_multianno.vcf" \
  > "${annovar_annotated_VarDict_SSCS}.hg38_multianno.vcf.gz"
tabix -p vcf "${annovar_annotated_VarDict_SSCS}.hg38_multianno.vcf.gz"

bcftools query \
  -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP\t%INFO/VD\t%INFO/AF\t%INFO/Func.refGene\t%INFO/Gene.refGene\t%INFO/ExonicFunc.refGene\t%INFO/gnomAD_genome_ALL\t%INFO/CLNSIG\n' \
  "${annovar_annotated_VarDict_SSCS}.hg38_multianno.vcf.gz" \
  > "$annotated_VarDict_SSCS_tsv"

echo "=== SSCS pileup ==="
samtools mpileup -A -q 0 -Q 0 -B -a -f "$REF" -d 1000000 \
  -l "$SNV_BED" "$SSCS_clipped_sorted" > "$pileup_SSCS"

##############################################################################
# DCS
##############################################################################

echo "=== [7/12] Call DCS (min reads = $CALLMINDUPLEX) ==="
java -Xmx32G -jar "$FGBIO_JAR" --compression 1 --async-io CallDuplexConsensusReads \
  --input "$grouped_BAM" \
  --output "$DCS_bam_unmapped" \
  --min-reads "$CALLMIN" \
  --min-input-base-quality "$BASEQUALMIN" \
  --threads "$THREADS"

echo "=== [8/12] Remap + merge DCS ==="
samtools fastq "$DCS_bam_unmapped" \
  | bwa mem -t "$THREADS" -p -K 150000000 -Y "$REF" - \
  | java -Xmx32G -jar "$FGBIO_JAR" --compression 1 --async-io ZipperBams \
      --unmapped "$DCS_bam_unmapped" \
      --ref "$REF" \
      --tags-to-reverse Consensus \
      --tags-to-revcomp Consensus \
      --output "$DCS_bam_mapped"

echo "=== [9/12] Filter + sort DCS ==="
java -Xmx32G -jar "$FGBIO_JAR" --compression 0 FilterConsensusReads \
  --input "$DCS_bam_mapped" \
  --output /dev/stdout \
  --ref "$REF" \
  --min-reads "$CALLMIN" \
  --min-base-quality "$CONSQUALMIN" \
  --max-base-error-rate "$MAXBASEERRORRATE" \
  | samtools sort -n --threads "$THREADS" -o "$DCS_bam_filtered" -

echo "=== [10/12] Clip overlapping reads + ends (DCS) ==="
java -Xmx32G -jar "$FGBIO_JAR" --compression 0 ClipBam \
  --input "$DCS_bam_filtered" \
  --output "$DCS_clipped" \
  --ref "$REF" \
  --clipping-mode Hard \
  --clip-overlapping-reads true \
  --read-one-five-prime "$CLIPBASES" \
  --read-one-three-prime "$CLIPBASES" \
  --read-two-five-prime "$CLIPBASES" \
  --read-two-three-prime "$CLIPBASES" \
  --metrics "$DCS_overlap_metrics"

echo "=== Sort + index DCS ==="
samtools sort -@ "$THREADS" -o "$DCS_clipped_sorted" "$DCS_clipped"
samtools index "$DCS_clipped_sorted"

echo "=== DCS variant calling (VarDictJava) ==="
vardict-java \
  -G "$REF" -N "$SAMPLE" -b "$DCS_clipped_sorted" \
  -z 1 -c 1 -S 2 -E 3 -g 4 \
  -f "$ALLFREQ" -r 1 -Q 0 -q 0 -o 0 -k 0 -y \
  -th "$THREADS" "$SNV_BED" \
  | teststrandbias.R \
  | var2vcf_valid.pl -N "$SAMPLE" -E -f "$ALLFREQ" -r 1 \
  | bgzip -c > "$VarDict_DCS_vcf"

echo "=== DCS annotation (ANNOVAR) ==="
perl "$ANNOVAR_ANNOTATE" "$VarDict_DCS_vcf" "$ANNOVAR_HUMANDB/" \
  -buildver hg38 -vcfinput \
  -out "$annovar_annotated_VarDict_DCS" \
  -remove \
  -protocol refGene,cosmic102_coding,cosmic101_noncoding,exac03,gnomad_genome,clinvar_20240917 \
  -operation g,f,f,f,f,f \
  -nastring . -thread "$THREADS"

bgzip -c "${annovar_annotated_VarDict_DCS}.hg38_multianno.vcf" \
  > "${annovar_annotated_VarDict_DCS}.hg38_multianno.vcf.gz"
tabix -p vcf "${annovar_annotated_VarDict_DCS}.hg38_multianno.vcf.gz"

bcftools query \
  -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP\t%INFO/VD\t%INFO/AF\t%INFO/Func.refGene\t%INFO/Gene.refGene\t%INFO/ExonicFunc.refGene\t%INFO/gnomAD_genome_ALL\t%INFO/CLNSIG\n' \
  "${annovar_annotated_VarDict_DCS}.hg38_multianno.vcf.gz" \
  > "$annotated_VarDict_DCS_tsv"

echo "=== DCS pileup ==="
samtools mpileup -A -a -q 0 -Q 0 -d 1000000 -f "$REF" \
  -l "$SNV_BED" "$DCS_clipped_sorted" > "$pileup_DCS"

##############################################################################
# CLEANUP
##############################################################################

rm -rf "$TEMP"
echo "Pipeline finished for sample $SAMPLE"
