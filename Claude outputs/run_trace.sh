#!/bin/bash
# =====================================================================
# run_trace.sh - trace raw Horizon alt reads through SSCS / DCS (legacy MUFS3 = 6-3-3)
#
# Usage (on the server, in the env that has pysam + samtools):
#   bash run_trace.sh                   # all 8 samples, then summary
#   bash run_trace.sh 999_01_rep1 ...   # only samples whose name contains these strings
#   bash run_trace.sh --summary-only
#
# Needs the neat + 10% samples too: they are the positive control and the
# reference for the cross-sample (index hopping / contamination) check.
# =====================================================================
set -euo pipefail

# ---- EDIT THESE -------------------------------------------------------
BASE=/powervault/zuzanna/horizon_legacy_duplex/standard_hybridisation
VARIANTS=$HOME/sequencing/horizon_variants_lookup/horizon_check/horizon_variants.tsv   # <-- the truth table horizon_pileup.py uses
OUT=$HOME/sequencing/horizon_variants_lookup/horizon_check/alt_family_trace
SCRIPT="$(dirname "$0")/trace_alt_families.py"
MUFS=3
THREADS=8
# ----------------------------------------------------------------------

# sample_udi directory  :  sample name (file prefix)
SAMPLES=(
  "Horizon_rep1_xGenUDI1:Horizon_rep1"
  "Horizon_rep2_xGenUDI2:Horizon_rep2"
  "LEG016_Horizon_90_10_rep1_xGenUDI11:LEG016_Horizon_90_10_rep1"
  "LEG016_Horizon_90_10_rep2_xGenUDI4:LEG016_Horizon_90_10_rep2"
  "LEG016_Horizon_99_1_rep1_xGenUDI13:LEG016_Horizon_99_1_rep1"
  "LEG016_Horizon_99_1_rep2_xGenUDI6:LEG016_Horizon_99_1_rep2"
  "LEG016_Horizon_999_01_rep1_xGenUDI7:LEG016_Horizon_999_01_rep1"
  "LEG016_Horizon_999_01_rep2_xGenUDI16:LEG016_Horizon_999_01_rep2"
)

mkdir -p "$OUT/subsets"
SUMMARY_ONLY=0
FILTERS=()
for a in "$@"; do
  if [[ "$a" == "--summary-only" ]]; then SUMMARY_ONLY=1; else FILTERS+=("$a"); fi
done

# Subset a BAM to reads overlapping the SNVs. samtools view -L works on
# unindexed BAMs (it streams the whole file), so this also works for the
# template-coordinate-sorted grouped BAM and the unsorted *_mapped.bam.
subset() {   # in out
  local in=$1 out=$2
  if [[ -s "$out" && -s "$out.bai" && "$out" -nt "$in" ]]; then echo "  reuse $(basename "$out")"; return; fi
  echo "  subsetting $(basename "$in")"
  samtools view -@ "$THREADS" -u -L "$BED" "$in" | samtools sort -@ "$THREADS" -o "$out" -
  samtools index "$out"
}

if [[ $SUMMARY_ONLY -eq 0 ]]; then
  first_grouped="$BASE/${SAMPLES[0]%%:*}/${SAMPLES[0]#*:}_UMI_grouped.bam"
  BED="$OUT/horizon_snvs.bed"
  python "$SCRIPT" bed --variants "$VARIANTS" --bam "$first_grouped" --out "$BED"

  for entry in "${SAMPLES[@]}"; do
    udi=${entry%%:*}; s=${entry#*:}
    if [[ ${#FILTERS[@]} -gt 0 ]]; then
      keep=0; for f in "${FILTERS[@]}"; do [[ "$udi" == *"$f"* ]] && keep=1; done
      [[ $keep -eq 1 ]] || continue
    fi
    echo "== $udi"
    d="$BASE/$udi"
    sub="$OUT/subsets/$udi"
    subset "$d/${s}_UMI_grouped.bam"                          "${sub}_grouped.bam"
    subset "$d/SSCS_files_MUFS${MUFS}/${s}_SSCS_mapped.bam"   "${sub}_SSCS_called.bam"
    subset "$d/DCS_files_MUFS${MUFS}/${s}_DCS_mapped.bam"     "${sub}_DCS_called.bam"
    python "$SCRIPT" trace \
      --variants "$VARIANTS" --sample "$udi" \
      --grouped     "${sub}_grouped.bam" \
      --sscs-called "${sub}_SSCS_called.bam" \
      --sscs-final  "$d/SSCS_files_MUFS${MUFS}/${s}_SSCS_mapped_filtered_clipped_sorted.bam" \
      --dcs-called  "${sub}_DCS_called.bam" \
      --dcs-final   "$d/DCS_files_MUFS${MUFS}/${s}_DCS_mapped_filtered_clipped_sorted.bam" \
      --outdir "$OUT" --min-ss "$MUFS"
  done
fi

python "$SCRIPT" summarise --indir "$OUT"
