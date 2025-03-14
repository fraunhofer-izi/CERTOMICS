# Genome metadata
genome="GRCh38"
version="2024-A"

# Set up source and build directories
build="GRCh38-GENCODEv44_build"
mkdir -p "$build"

# Download source files if they do not exist in reference_sources/ folder
source="reference_sources"
mkdir -p "$source"
mv "!{src_fa}" $source
mv "!{src_gtf}" $source
fasta_in="${source}/!{src_fa}"
gtf_in="${source}/!{src_gtf}"

# Modify sequence headers in the Ensembl FASTA to match the file
# "GRCh38.primary_assembly.genome.fa" from GENCODE. Unplaced and unlocalized
# sequences such as "KI270728.1" have the same names in both versions.
#
# Input FASTA:
#   >1 dna:chromosome chromosome:GRCh38:1:1:248956422:1 REF
#
# Output FASTA:
#   >chr1 1
fasta_modified="$build/$(basename "$fasta_in").modified"
# sed commands:
# 1. Replace metadata after space with original contig name, as in GENCODE
# 2. Add "chr" to names of autosomes and sex chromosomes
# 3. Handle the mitochrondrial chromosome
cat "$fasta_in" \
    | sed -E 's/^>(\S+).*/>\1 \1/' \
    | sed -E 's/^>([0-9]+|[XY]) />chr\1 /' \
    | sed -E 's/^>MT />chrM /' \
    > "$fasta_modified"


# Remove version suffix from transcript, gene, and exon IDs in order to match
# previous Cell Ranger reference packages
#
# Input GTF:
#     ... gene_id "ENSG00000223972.5"; ...
# Output GTF:
#     ... gene_id "ENSG00000223972"; gene_version "5"; ...
gtf_modified="$build/$(basename "$gtf_in").modified"
# Pattern matches Ensembl gene, transcript, and exon IDs for human or mouse:
ID="(ENS(MUS)?[GTE][0-9]+)\.([0-9]+)"
cat "$gtf_in" \
    | sed -E 's/gene_id "'"$ID"'";/gene_id "\1"; gene_version "\3";/' \
    | sed -E 's/transcript_id "'"$ID"'";/transcript_id "\1"; transcript_version "\3";/' \
    | sed -E 's/exon_id "'"$ID"'";/exon_id "\1"; exon_version "\3";/' \
    > "$gtf_modified"


# Define string patterns for GTF tags
# NOTES:
# Since Ensembl 110, polymorphic pseudogenes are now just protein_coding.
# Readthrough genes are annotated with the readthrough_transcript tag.
BIOTYPE_PATTERN=\
"(protein_coding|protein_coding_LoF|lncRNA|\
IG_C_gene|IG_D_gene|IG_J_gene|IG_LV_gene|IG_V_gene|\
IG_V_pseudogene|IG_J_pseudogene|IG_C_pseudogene|\
TR_C_gene|TR_D_gene|TR_J_gene|TR_V_gene|\
TR_V_pseudogene|TR_J_pseudogene)"
GENE_PATTERN="gene_type \"${BIOTYPE_PATTERN}\""
TX_PATTERN="transcript_type \"${BIOTYPE_PATTERN}\""
READTHROUGH_PATTERN="tag \"readthrough_transcript\""


# Construct the gene ID allowlist. We filter the list of all transcripts
# based on these criteria:
#   - allowable gene_type (biotype)
#   - allowable transcript_type (biotype)
#   - no "readthrough_transcript" tag
# We then collect the list of gene IDs that have at least one associated
# transcript passing the filters.
cat "$gtf_modified" \
    | awk '$3 == "transcript"' \
    | grep -E "$GENE_PATTERN" \
    | grep -E "$TX_PATTERN" \
    | grep -Ev "$READTHROUGH_PATTERN" \
    | sed -E 's/.*(gene_id "[^"]+").*/\1/' \
    | sort \
    | uniq \
    > "${build}/gene_allowlist"


# NOTES:
# Since Ensembl 110, the PAR locus genes are included on chrY as copies of chrX
# Using the GRCh38.p13 assembly hard masks these regions on chrY, but removing the
# chrY PAR genes is still desirable so they do not end up as extra entries in the output.
# The awk command below excludes all PAR_Y genes, including XGY2.
# The non-coding gene XGY2 straddles the PAR1 boundary on chrY, and is homologous to XG on chrX.
# GRCh38-2024-A excludes XGY2, but includes SRY and ENSG00000286130, which are in an intron of XGY2,
# and RPS4Y1, which overlaps XGY2.

# Filter the GTF file based on the gene allowlist
gtf_filtered="${build}/$(basename "$gtf_in").filtered"
# Copy header lines beginning with "#"
grep -E "^#" "$gtf_modified" > "$gtf_filtered"
# Filter to the gene allowlist, and then remove PAR_Y genes
grep -Ff "${build}/gene_allowlist" "$gtf_modified" \
    | awk -F "\t" '$1 != "chrY" || $1 == "chrY" && $4 >= 2752083 && $4 < 56887903 && !/ENSG00000290840/' \
    >> "$gtf_filtered"

# Concatenate CAR files
fasta_final="$fasta_modified"
gtf_final="$gtf_filtered"
if [ -f "!{car_fa}" ] && [ -f "!{car_gtf}" ]; then
  gtf_car="${build}/$(basename "$gtf_final")_car.gtf"
  cat "$gtf_final" "!{car_gtf}" > "$gtf_car"
  gtf_final="$gtf_car"

  fasta_car="${build}/$(basename "$fasta_final")_car.fa"
  cat "$fasta_final" "!{car_fa}" > "$fasta_car"
  fasta_final="$fasta_car"
else
  [[ ! -f "!{car_fa}" ]] && echo "Missing: !{car_fa}"
  [[ ! -f "!{car_gtf}" ]] && echo "Missing: !{car_gtf}"
  echo "Missing at least one car file. Not concatenating."
fi

# Create reference package
cellranger mkref \
  --ref-version="$version" \
  --genome="$genome" \
  --fasta="$fasta_final" \
  --genes="$gtf_final" \
  --output-dir="gex_reference" \
  !{task.cpus ? '--nthreads=' + task.cpus : ''} \
  !{task.memory ? '--memgb=' + task.memory.toGiga() : ''}
