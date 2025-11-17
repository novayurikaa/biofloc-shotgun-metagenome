============================================================================================
######################functional annotation using eggnog mapper and prokka#######################
===========================================================================================
pwd
cd /Users/mm/Documents/project/metagenome/
mkdir 4_functional_annotation
cd /Users/mm/Documents/project/metagenome/4_functional_annotation

#install prokka and eggnog mapper using conda
conda activate bioinfo
mamba install -y -c conda-forge -c bioconda prokka eggnog-mapper diamond hmmer prodigal
#verify installation
prokka --version  #there is problem.
===========================
#lets run PROKKA first
=========================
cd /Users/mm/Documents/project/metagenome/4_functional_annotation
mkdir -p bins prokka eggnog

BINS_DIR=$(pwd)/bins
PROKKA_OUT=$(pwd)/prokka

for bin in $BINS_DIR/*.fa; do
    bin_name=$(basename $bin .fa)
    echo "Annotating $bin_name with Prokka..."
    prokka --outdir $PROKKA_OUT/$bin_name \
           --prefix $bin_name \
           --cpus 8 \
           $bin
done
#Creates one folder per bin in prokka/
#Outputs: .faa (protein), .ffn (CDS nucleotide), .gff, .gbk


=================================================
#NEXT STEP: FUNCTIONAL ANNTOTATION
#Goal: assign functions (e.g., COG, GO, KEGG) to your proteins using EggNOG-mapper
==================================================
mkdir -p /Users/mm/Documents/project/metagenome/4_functional_annotation/eggnog_db
#create directory eggnog database
cd eggnog_db

# 1. Annotation database
wget http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.db.gz
gunzip eggnog.db.gz 

# 2. DIAMOND database
wget http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz
gunzip eggnog_proteins.dmnd.gz 

# 3. Taxonomic database
wget http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz
tar -zxf eggnog.taxa.tar.gz 
=========================================================================
#step 3:run functional annotation for each bin using eggNOG-mapper
========================================================================
cd /users/mm/Documents/project/metagenome/4_functional_annotation

#1. set path
# Folder with Prokka bins
PROKKA_DIR=/Users/mm/Documents/project/metagenome/4_functional_annotation/prokka

# Folder with EggNOG databases (eggnog.db, eggnog_proteins.dmnd, eggnog.taxa)
EGGNOG_DB=/Users/mm/Documents/project/metagenome/4_functional_annotation/eggnog_db

# Output folder for EggNOG annotations
EGGNOG_OUT=/Users/mm/Documents/project/metagenome/4_functional_annotation/eggnog_annotations
mkdir -p $EGGNOG_OUT

#2.run eggnog mapper via docker, for all bin in prokka folder
#biocontainers
for bin_dir in /Users/mm/Documents/project/metagenome/4_functional_annotation/prokka/*; do
    bin_name=$(basename $bin_dir)
    docker run --rm -v /Users/mm/Documents/project/metagenome/4_functional_annotation:/data \
               quay.io/biocontainers/eggnog-mapper:2.1.9--pyhdfd78af_0 \
               emapper.py -i /data/prokka/$bin_name/$bin_name.faa \
                          --output /data/eggnog_annotations/$bin_name \
                          --cpu 8 \
                          --itype proteins \
                          --data_dir /data/eggnog_db
done


###############
#steo 4: check output files


#A. For each bin, you should have:
#bin.X.emapper.hits → raw DIAMOND/sequence search hits
#bin.X.emapper.seed_orthologs → selected orthologs
#bin.X.emapper.annotations → final functional annotation
#The .annotations file is the main output you’ll use for downstream analyses.

#B. Inspect the .annotations file
#The .annotations file contains:
#Query gene IDs (from your Prokka .faa)
#Predicted orthologs in eggNOG
#Functional categories: COG, KEGG, GO terms, EC numbers
#Gene description and annotations
#You can quickly inspect it with:
head /Users/mm/Documents/project/metagenome/4_functional_annotation/eggnog_annotations/bin.1.emapper.annotations


#C. If you want a single summary table for all bins:
cat /Users/mm/Documents/project/metagenome/4_functional_annotation/eggnog_annotations/*.emapper.annotations > all_bins_annotations.tsv
#You can later filter by functional category, GO term, or KEGG pathway.

#D. Optional: Mapping to functional categories
#EggNOG-mapper outputs have these fields for functional analysis:
#COG_category → general functional class
#KEGG_ko → metabolic pathways
#GO_terms → gene ontology (molecular function, biological process, cellular component)
#EC → enzyme classification

#E.You can use R, Python (pandas), or Excel to summarize:
#Count genes per functional category
#Compare metabolic potential between bins
#Plot KEGG pathways or GO terms

#F.Downstream ideas
#Now that you have annotations:
#Functional comparison of MAGs/bins
#Metabolic reconstruction (based on KEGG modules/EC numbers)
#Identify genes of interest (e.g., carbon cycling, antibiotic resistance, etc.)
#Visualize with heatmaps, barplots, or pathway maps



#EggNOG-mapper (you will ran)
#Purpose:
#EggNOG-mapper is for global functional annotation based on orthology.
#Assigns general gene function:
#COG functional categories (general functions like metabolism, transcription, replication)
#KEGG pathways (metabolism, energy, transport, etc.)
#GO terms (molecular function, biological process, cellular component)
#EC numbers (enzymes)
#Input:
#Protein sequences (your Prokka .faa files).
#Output:
#Annotations for every gene (or most genes) in the genome/bin.
#Gives a comprehensive functional landscape, not just BGCs.
#Key point:
#Broad functional overview, including housekeeping genes and metabolic genes.
#Not specialized for secondary metabolites like AntiSMASH.



# Key point:
#The .annotations files are the “finished product” from EggNOG-mapper. 
#Your next steps are mainly data parsing, functional summarization, and visualization.


#in bioflock context:

#1.AntiSMASH (secondary metabolite potential)
#Why it matters for biofloc:
#Biofloc systems rely on microbial communities to improve water quality and suppress pathogens.
#AntiSMASH tells you which microbes in your MAGs can produce bioactive compounds, such as:
#-Antimicrobials → natural control against pathogenic bacteria in the tank.
#-Siderophores → help in nutrient cycling (iron scavenging).
#-Pigments or signaling molecules → may influence microbial interactions.
#-How to use it:
#-Identify microbes with potential probiotic traits or bioactive metabolites.
#-Focus on MAGs that have clusters producing compounds beneficial for water quality or shrimp health.
#2. EggNOG-mapper (general functional profile)
#Why it matters for biofloc:
#Gives a genome-wide view of microbial metabolic capabilities, which is essential in a biofloc system where microbes drive nutrient recycling.
#Examples of what you can see:
#Nitrogen metabolism genes → ammonia oxidation, nitrate reduction → key for water quality.
#Carbon metabolism genes → help understand how microbes break down organic matter in biofloc.
#Stress response genes → indicate which microbes are resilient under biofloc conditions (high density, low oxygen, variable pH).
#Transporters and secretion systems → show how microbes interact with each other and the host.
#How to use it:
#Quantify which functional categories dominate in your MAGs.
#Identify microbes that can contribute to nutrient cycling or biofloc stability.
#Compare bins to see which MAGs have complementary functions — e.g., one fixes nitrogen, another degrades carbon-rich waste.
#3. Combined interpretation
#Tool	What it shows	Biofloc relevance
#AntiSMASH	BGCs (secondary metabolites)	Microbe-microbe interactions, natural pathogen suppression
#EggNOG-mapper	Functional genes (all metabolism)	Nutrient cycling, resilience, overall biofloc ecosystem function