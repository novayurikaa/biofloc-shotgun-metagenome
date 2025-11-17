
# ==========================================
# Step 3: Taxonomic and Functional Profiling
# ==========================================

#checi if trimmed reads are available
ls 1_quality_control/trimmed_reads/ 
#make subfolder
mkdir -p 2_taxonomy_functional/{kraken2,bracken,humann}

#install kraken2 and bracken
mamba install -c conda-forge -c bioconda kraken2 bracken


cd 2_taxonomy_functional
wget https://genome-idx.s3.amazonaws.com/kraken/minikraken2_v2_8GB_201904.tgz
tar -xvzf minikraken2_v2_8GB_201904.tgz
mv minikraken2_v2_8GB_201904 kraken2_db #download database
tar -xvzf minikraken2_v2_8GB_201904.tgz  #etxract database

mkdir -p kraken2_db
mv minikraken2_v2_8GB_201904_UPDATE kraken2_db
ls kraken2_db
mv kraken2_db/minikraken2_v2_8GB_201904_UPDATE/* kraken2_db/
rmdir kraken2_db/minikraken2_v2_8GB_201904_UPDATE


=====
#install kraken2 and bracken in phyton 3.8 enviornment
=======
mkdir -p kraken2_output bracken_output kraken2_db
conda install mamba -n base -c conda-forge -y
mamba install bracken #not working because phyton 3.13

conda create -n bioinfo38 python=3.8 -c defaults -c conda-forge -c bioconda -y #this is for phyton 3.8, safely to install bracken
mamba install kraken2 sra-tools #download other tools too
#it works. now we have two environments; bioinfo for the latest phyton and bioinfo 38 for phyton 3.8

#next time always specify the phytonversion when creating new environment
#conda create -n bioinfo38 python=3.8 -c defaults -c conda-forge -c bioconda -y
#bioinfo38 is the environment name.
#Channels (conda-forge and bioconda) ensure you can install bioinformatics tools.

#-y automatically confirms installation.
#After that, you can activate it:
#conda activate bioinfo38
#Then install Bracken:
#mamba install bracken=2.6.1


=============================================================
# Taxonomic classification and functional profiling pipeline
# Nova Yurika | Biofloc metagenomics project
=========================================================

# Activate environment (adjust if needed)
conda activate bioinfo38

# Paths
PROJECT_DIR="/Users/mm/Documents/project/metagenome/2_taxonomy_functional"
DB_DIR="$PROJECT_DIR/kraken2_db"
KRAKEN_DIR="$PROJECT_DIR/kraken2"
BRACKEN_DIR="$PROJECT_DIR/bracken"

# Input reads (trimmed FASTQ files)
READS_DIR="/Users/mm/Documents/project/metagenome/1_quality_control/trimmed_reads"
READS_DIR="/Users/mm/Documents/project/metagenome/1_quality_control/trimmed_reads"


# Create output folders if missing
mkdir -p "$KRAKEN_DIR" "$BRACKEN_DIR"

# List of samples (adjust if different)
SAMPLES=("SRR24442555" "SRR24442556" "SRR24442557")

# Run Kraken2 + Bracken

for SAMPLE in "${SAMPLES[@]}"; do
    READ1="$READS_DIR/${SAMPLE}_1.trimmed.fastq"
    READ2="$READS_DIR/${SAMPLE}_2.trimmed.fastq"

    kraken2 \
      --db "$DB_DIR" \
      --paired "$READS_DIR/${SAMPLE}_1.trimmed.fastq" "$READS_DIR/${SAMPLE}_2.trimmed.fastq" \
      --threads 16 \
      --report "$KRAKEN_DIR/${SAMPLE}_kraken2.report" \
      --output "$KRAKEN_DIR/${SAMPLE}_kraken2.output"

    bracken \
      -d "$DB_DIR" \
      -i "$KRAKEN_DIR/${SAMPLE}_kraken2.report" \
      -o "$BRACKEN_DIR/${SAMPLE}_bracken.txt" \
      -r 150 -l S
done

#look at kraken2 output
# Example for one sample
less /Users/mm/Documents/project/metagenome/2_taxonomy_functional/kraken2/SRR24442557_kraken2.report

#% reads – proportion of total reads assigned to that taxon
# reads – number of reads assigned to this taxon
#Rank – taxonomic rank (D = Domain, P = Phylum, C = Class, O = Order, F = Family, G = Genus, S = Species)
#Name – name of the taxon

## Top 10 species
grep "S" SRR24442557_kraken2.report | sort -nrk2 | head
#2. Look at the full output file
#This shows the classification for each read:
head /Users/mm/Documents/project/metagenome/2_taxonomy_functional/kraken2/SRR24442557_kraken2.output



#check all bracken output for each sample
ls -lh /Users/mm/Documents/project/metagenome/2_taxonomy_functional/bracken/

#Optional: merge all species tables into a single matrix for downstream analysis (e.g., R, Python):
#bracken_combine.py -i /path/to/bracken/ -o bracken_combined.tsv
#Proceed to downstream analysis like diversity metrics, differential abundance, or visualization.




=============================================================
summarize kraken2 output summary
======================================
PROJECT_DIR="/Users/mm/Documents/project/metagenome/2_taxonomy_functional"
KRAKEN_DIR="$PROJECT_DIR/kraken2"
OUTPUT_SUMMARY="$PROJECT_DIR/kraken2_summary_species.tsv"

# Create header
echo -e "Sample\tReads\tTaxID\tRank\tTaxon" > "$OUTPUT_SUMMARY"

# Loop through Kraken2 reports
for REPORT in "$KRAKEN_DIR"/*_kraken2.report; do
    SAMPLE=$(basename "$REPORT" _kraken2.report)
    
    awk -F'\t' -v sample="$SAMPLE" '$4=="S"{print sample "\t" $2 "\t" $5 "\t" $4 "\t" $6}' "$REPORT" >> "$OUTPUT_SUMMARY"
done
#check
head "$OUTPUT_SUMMARY" #this is only the first sample

#lets summarize for all samples

PROJECT_DIR="/Users/mm/Documents/project/metagenome/2_taxonomy_functional"
KRAKEN_DIR="$PROJECT_DIR/kraken2"
OUTPUT_SUMMARY="$PROJECT_DIR/kraken2_summary_species.tsv"

# Create header
echo -e "Sample\tReads\tTaxID\tRank\tTaxon" > "$OUTPUT_SUMMARY"


# Loop through all Kraken2 report files
for REPORT in "$KRAKEN_DIR"/*_kraken2.report; do
    SAMPLE=$(basename "$REPORT" _kraken2.report)
    awk -F'\t' -v sample="$SAMPLE" '$4=="S"{print sample "\t" $2 "\t" $5 "\t" $4 "\t" $6}' "$REPORT" >> "$OUTPUT_SUMMARY"
done


==========================
#LETS CREATE BRACKEN SUMMARY
=========================

PROJECT_DIR="/Users/mm/Documents/project/metagenome/2_taxonomy_functional"
BRACKEN_DIR="$PROJECT_DIR/bracken"
OUTPUT_SUMMARY="$PROJECT_DIR/bracken_summary_species.tsv"

# Create header
echo -e "Sample\tReads\tTaxID\tRank\tTaxon" > "$OUTPUT_SUMMARY"

# Loop over Bracken outputs
for FILE in "$BRACKEN_DIR"/*_bracken.txt; do
    SAMPLE=$(basename "$FILE" _bracken.txt)

    awk -F'\t' -v sample="$SAMPLE" '$3=="S"{print sample "\t"$5"\t"$2"\t"$3"\t"$1}' "$FILE" >> "$OUTPUT_SUMMARY"
done




#OPEN
open -a "Microsoft Excel" "$OUTPUT_SUMMARY"
open -a "Microsoft Excel" "$OUTPUT_SUMMARY"

=================================================================
#lets visualize kracken2 output using pavian
==========================================================
#go to https://fbreitwieser.shinyapps.io/pavian/ 
#uplaod
#SRR24442555_kraken2.report
#SRR24442556_kraken2.report
#SRR24442557_kraken2.report

=================================================================
#lets visualize kracken2 and bracken output using krona
==========================================================

cd /Users/mm/Documents/project/metagenome/2_taxonomy_functional
brew install krona
./updateTaxonomy.sh
cd $CONDA_PREFIX/share/krona
./updateTaxonomy.sh


cd /Users/mm/miniconda/envs/bioinfo38/opt/krona
./updateTaxonomy.sh

cd /Users/mm/Documents/project/metagenome/2_taxonomy_functional
ktImportTaxonomy kraken2_summary_species.tsv -o krona_kraken2.html


=====
#not finished yet, use galaxy instead
=================





conda activate bioinfo
touch MAG.sh  #create new script for mag and gene clustering
