#!/bin/bash

# List of biotypes

# Genome metadata
genome="GRCh38"
version="p13-Gencode_v33"


# Set up source and build directories
build="GRCh38-GENCODE-V33-resources"
mkdir -p "$build"


# Download source files if they do not exist in reference_sources/ folder
sourcee="reference_sources"
mkdir -p "$sourcee"

# https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
# https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz
#fasta_url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.p14.genome.fa.gz" #1
#fasta_url="https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" #2
fasta_url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/GRCh38.p13.genome.fa.gz"  #3
fasta_in="${sourcee}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

## Mathcing GTF files to FASTA files (#number indicates which FASTA-GTF pairs are matching)
#gtf_url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz" #1
#gtf_url="https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz" #2
gtf_url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.primary_assembly.annotation.gtf.gz" #3
gtf_in="${sourcee}/gencode.v33.primary_assembly.annotation.gtf"

if [ ! -f "$fasta_in" ]; then
    curl -sS "$fasta_url" | zcat > "$fasta_in"
fi

if [ ! -f "$gtf_in" ]; then
    curl -sS "$gtf_url" | zcat > "$gtf_in"
fi


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

## Define the neame of the filtered (or modified) GTF file, that will be filtered for the selcted biotypes
gtf_modified="$build/$(basename "$gtf_in").modified"

biotypes=("protein_coding" "lncRNA" "antisense" "IG_LV_gene" "IG_V_gene" "IG_V_pseudogene" "IG_D_gene" "IG_J_gene" "IG_J_pseudogene" "IG_C_gene" "IG_C_pseudogene" "TR_V_gene" "TR_V_pseudogene" "TR_D_gene" "TR_J_gene" "TR_J_pseudogene" "TR_C_gene")

# Initialize the command
command="cellranger mkgtf $gtf_in $gtf_modified"

# Iterate through the biotypes and add --attribute options to the command
for biotype in "${biotypes[@]}"; do
    command+=" --attribute=gene_biotype:$biotype"
done



# Execute the command
echo "Executing the following command:"
echo "$command"

eval $command

echo "start runnung STAR reference genom genaration" 
#STAR --runThreadN 16 \
#--runMode genomeGenerate \
#--genomeDir "$genome" \
#--genomeFastaFiles $fasta_in \
#--sjdbGTFfile $gtf_modified \
#--sjdbOverhang 99


echo "Starting to run mkref"
# Create reference package
cellranger mkref --ref-version="$version" \
    --genome="$genome-$version" --fasta="$fasta_modified" --genes="$gtf_modified"
