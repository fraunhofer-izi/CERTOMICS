#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

#############################################
# Usage:
# ./compare_all_CARs_against_reference.sh \
#   <BAM> <SAMPLE> <REF_TG> <KMER_DIR>
#############################################

# ---------------------------
# Input arguments
# ---------------------------
BAM="$1"
SAMPLE="$2"
REF_TG="$3"
KMER_DIR="$4"

#############################################
# Output files
#############################################

TG_BAM="${SAMPLE}.${REF_TG}.bam"
TG_FASTQ="${SAMPLE}.${REF_TG}.fastq"
TG_SEQS="${SAMPLE}.${REF_TG}.seqs.txt"
SUMMARY="${SAMPLE}.CAR_kmer_summary.txt"

#############################################
# 0. Check reference transgene
#############################################

echo "Checking if reference transgene is present in BAM..."

REF_READS=$(samtools view -c "$BAM" "$REF_TG" || true)

if [[ "$REF_READS" -eq 0 ]]; then
    echo "❌ Reference transgene $REF_TG NOT found in BAM"
    FOUND_REF="NO"
else
    echo "✔ Reference transgene $REF_TG found ($REF_READS reads)"
    FOUND_REF="YES"
fi

#############################################
# 1. Write header
#############################################

echo -e "Sample\tReference_transgene\tRef_found\tRef_reads\tTested_CAR\tKmer_matches\tReads_tested" > "$SUMMARY"

#############################################
# 2. If REF_TG NOT found → write zeros and exit
#############################################

if [[ "$FOUND_REF" == "NO" ]]; then
    echo "▶ No reference reads found – setting all k-mer counts to 0"

    for kmer in "$KMER_DIR"/*unique_31mers.txt; do
        CAR=$(basename "$kmer" _unique_31mers.txt)
        echo -e "${SAMPLE}\t${REF_TG}\tNO\t0\t${CAR}\t0\t0" >> "$SUMMARY"
    done

    echo "✅ Done (no reference transgene reads)"
    echo "Summary written to: $SUMMARY"
    exit 0
fi

#############################################
# 3. Extract reads mapped to reference transgene
#############################################

echo "Extracting reads mapped to $REF_TG..."
samtools view -b "$BAM" "$REF_TG" > "$TG_BAM"
samtools fastq "$TG_BAM" > "$TG_FASTQ"

#############################################
# 4. Extract sequences
#############################################

awk 'NR % 4 == 2' "$TG_FASTQ" > "$TG_SEQS"
TOTAL_READS=$(wc -l < "$TG_SEQS")

#############################################
# 5. Test ALL CAR kmers
#############################################

for kmer in "$KMER_DIR"/*unique_31mers.txt; do
    CAR=$(basename "$kmer" _unique_31mers.txt)
    MATCHES=$(grep -F -c -f "$kmer" "$TG_SEQS" || true)

    echo -e "${SAMPLE}\t${REF_TG}\tYES\t${REF_READS}\t${CAR}\t${MATCHES}\t${TOTAL_READS}" >> "$SUMMARY"
done

#############################################
# Done
#############################################

echo "✅ Done"
echo "Summary written to: $SUMMARY"
