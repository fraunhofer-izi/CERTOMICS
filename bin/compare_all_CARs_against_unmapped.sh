#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

#############################################
# Usage:
# ./compare_all_CARs_against_unmapped.sh \
#   <BAM> <SAMPLE> <KMER_DIR> <K> <THREADS>
#############################################

BAM="$1"
SAMPLE="$2"
KMER_DIR="$3"
K="$4"
THREADS="$5"

#############################################
# Paths
#############################################

UNMAPPED_FASTA="${SAMPLE}.unmapped.fasta"
UNMAPPED_DB="${SAMPLE}.unmapped_k${K}"

TMPDIR="${SAMPLE}.kmc_tmp"
CAR_DB_DIR="${SAMPLE}.car_kmc_dbs"
CAR_FASTA_DIR="${SAMPLE}.car_kmer_fastas"

SUMMARY="${SAMPLE}.unmapped_kmc_intersect_summary.txt"

mkdir -p "$TMPDIR" "$CAR_DB_DIR" "$CAR_FASTA_DIR"

#############################################
# 1. Extract unmapped reads (FASTA)
#############################################

echo "▶ Extracting unmapped reads from BAM..."
samtools view -@ "$THREADS" -f 4 "$BAM" \
| awk '{print ">"$1"\n"$10}' \
> "$UNMAPPED_FASTA"

READS_TESTED=$(grep -c "^>" "$UNMAPPED_FASTA")
echo "▶ Unmapped reads: $READS_TESTED"

#############################################
# 2. Build KMC DB from unmapped reads (if needed)
#############################################


echo "▶ Building KMC DB from unmapped reads..."

kmc \
    -k"$K" \
    -t"$THREADS" \
    -m64 \
    -ci1 \
    -fm \
    "$UNMAPPED_FASTA" \
    "$UNMAPPED_DB" \
    "$TMPDIR"

#############################################
# 3. Convert CAR kmer TXT → FASTA (if needed)
#############################################

echo "▶ Preparing CAR kmer FASTA files..."

for kmer_txt in "$KMER_DIR"/*unique_${K}mers.txt; do
    CAR=$(basename "$kmer_txt" _unique_${K}mers.txt)
    CAR_FASTA="${CAR_FASTA_DIR}/${CAR}.fasta"

    if [[ ! -f "$CAR_FASTA" ]]; then
        echo "  → Converting $CAR kmers to FASTA"
        awk '{print ">kmer_"NR"\n"$0}' "$kmer_txt" > "$CAR_FASTA"
    else
        echo "  → FASTA for $CAR already exists – skipping"
    fi
done

#############################################
# 4. Build KMC DBs for CAR kmer sets (if needed)
#############################################

echo "▶ Building CAR KMC DBs..."

for CAR_FASTA in "$CAR_FASTA_DIR"/*.fasta; do
    CAR=$(basename "$CAR_FASTA" .fasta)
    CAR_DB="${CAR_DB_DIR}/${CAR}.k${K}"

    if [[ ! -f "${CAR_DB}.kmc_pre" ]]; then
        echo "  → Building DB for $CAR"

        kmc \
          -k"$K" \
          -t"$THREADS" \
          -m4 \
          -ci1 \
          -fm \
          "$CAR_FASTA" \
          "$CAR_DB" \
          "$TMPDIR"
    else
        echo "  → DB for $CAR already exists – skipping"
    fi
done

#############################################
# 5. Intersect unmapped DB with each CAR DB
#############################################

echo -e "Sample\tReads_tested\tTested_CAR\tKmer_hits" > "$SUMMARY"

for PRE in "$CAR_DB_DIR"/*.k${K}.kmc_pre; do
    CAR_DB="${PRE%.kmc_pre}"
    CAR=$(basename "$CAR_DB" .k${K})
    HIT_DB="${CAR}.hits.tmp"

    echo "▶ Intersecting with $CAR"

    kmc_tools simple \
        "$UNMAPPED_DB" \
        "$CAR_DB" \
        intersect \
        "$HIT_DB"

    HITS=$(kmc_dump "$HIT_DB" /dev/stdout \
           | awk '{sum += $2} END {print sum+0}')

    echo -e "${SAMPLE}\t${READS_TESTED}\t${CAR}\t${HITS}" >> "$SUMMARY"

    rm -f "$HIT_DB".*
done

#############################################
# Done
#############################################

echo "✅ Done"
echo "Summary written to: $SUMMARY"
