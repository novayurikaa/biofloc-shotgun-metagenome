pwd
cd /Users/mm/Documents/project/metagenome/3_MAG_gene_cluster
==================================================================================================
##################.        MAG analysis using sumsampled reads###############################
=========================================================================================
#the fastq file is too big so we need to cut into small chunks and do the analysis on them

#create subsample folder and fastq for trimmed reads for SRR24442557  
TRIMMED_DIR="/Users/mm/Documents/project/metagenome/1_quality_control/trimmed_reads"
OUTDIR="$TRIMMED_DIR/SRR24442557_subsample"
R1="$TRIMMED_DIR/SRR24442557_1.trimmed.fastq"
R2="$TRIMMED_DIR/SRR24442557_2.trimmed.fastq"

READS=1000000 

head -n $((READS * 4)) $R1 > $OUTDIR/SRR24442557_1.subsample.fastq
head -n $((READS * 4)) $R2 > $OUTDIR/SRR24442557_2.subsample.fastq

ls -lh $OUTDIR  #check subsample files
=================================
#step 1: assembly
============================

#run megahit using docker
docker run -v /Users/mm/Documents/project/metagenome:/data quay.io/biocontainers/megahit:1.2.9--h8b12597_0 \
  megahit \
  -1 /data/1_quality_control/trimmed_reads/SRR24442557_subsample/SRR24442557_1.subsample.fastq \
  -2 /data/1_quality_control/trimmed_reads/SRR24442557_subsample/SRR24442557_2.subsample.fastq \
  -o /data/3_MAG_gene_cluster/mag_SRR24442557/assembly \
  --min-contig-len 1000 \
  -t 8 \
  --k-list 21,29,39


#now i have  MEGAHIT assembly (final.contigs.fa)

=============================
# Step 2: Mapping reads
============================

# Base folder for mapping
MAPPING_DIR=/Users/mm/Documents/project/metagenome/3_MAG_gene_cluster/mag_SRR24442557/mapping
mkdir -p $MAPPING_DIR


ASSEMBLY=/Users/mm/Documents/project/metagenome/3_MAG_gene_cluster/mag_SRR24442557/assembly/final.contigs.fa

# Index assembly
bwa index $ASSEMBLY


R1=/Users/mm/Documents/project/metagenome/1_quality_control/trimmed_reads/SRR24442557_subsample/SRR24442557_1.subsample.fastq
R2=/Users/mm/Documents/project/metagenome/1_quality_control/trimmed_reads/SRR24442557_subsample/SRR24442557_2.subsample.fastq

# Map reads
bwa mem -t 8 $ASSEMBLY $R1 $R2 | samtools view -b - > $MAPPING_DIR/aln.bam


# Sort BAM
samtools sort -@ 8 -o $MAPPING_DIR/aln.sorted.bam $MAPPING_DIR/aln.bam

# Index BAM
samtools index $MAPPING_DIR/aln.sorted.bam


#check mapping summary
samtools flagstat $MAPPING_DIR/aln.sorted.bam
#~14% of reads mapped to the assembly.
#This is expected for a subsampled dataset (2 million reads) and a de novo assembly with possibly fragmented contigs.
#This analysis (subsample) only used 2 million reads out of the full dataset (which is 20–25 GB per read file). That’s ~10% of the total reads, so only a fraction of  assembly is “seen” by this subset.
# In short: The result is valid for testing the workflow and getting an idea of how MAGs might be recovered, but it is not representative of the full metagenome. For high-confidence MAGs, we need:
#Full dataset (all reads), sufficient coverage per genome, possibly multiple k-mer assemblies or larger contig minimums

=========================
# Step 3: Calculate depth
=========================
#set path
ASSEMBLY=/Users/mm/Documents/project/metagenome/3_MAG_gene_cluster/mag_SRR24442557/assembly/final.contigs.fa
MAPPING_DIR=/Users/mm/Documents/project/metagenome/3_MAG_gene_cluster/mag_SRR24442557/mapping
OUTDIR=/Users/mm/Documents/project/metagenome/3_MAG_gene_cluster/mag_SRR24442557

# Output depth file
DEPTH_FILE=$OUTDIR/depth.txt


# Check if BAM file exists
if [ ! -f $MAPPING_DIR/aln.sorted.bam ]; then
    echo "Error: BAM file not found at $MAPPING_DIR/aln.sorted.bam"
    exit 1
fi


# Run depth calculation using MetaBAT2 script
echo "Calculating contig depth..."
jgi_summarize_bam_contig_depths --outputDepth $DEPTH_FILE $MAPPING_DIR/aln.sorted.bam


# Make sure OUTDIR is defined
OUTDIR=/Users/mm/Documents/project/metagenome/3_MAG_gene_cluster/mag_SRR24442557

# Path to depth file
DEPTH_FILE="$OUTDIR/depth.txt"

# Check if depth file exists
if [-f "$DEPTH_FILE"]; then
    echo "Depth calculation completed successfully!"
    echo "Preview of depth file:"
    head -n 10 "$DEPTH_FILE"
else
    echo "Error: Depth calculation failed."
    exit 1
fi

# if it not working, check chek_depth.sh script, run using nano

#successfully created binning contigs. File name: depth.txt
#Folder: /Users/mm/Documents/project/metagenome/3_MAG_gene_cluster/mag_SRR24442557

#to view it completely

cat /Users/mm/Documents/project/metagenome/3_MAG_gene_cluster/mag_SRR24442557/depth.txt


======================================================
# Step 4: Binning contigs
=====================================================

cd /Users/mm/Documents/project/metagenome/3_MAG_gene_cluster
conda activate bioinfo

#1. prepare directory
OUTDIR=/Users/mm/Documents/project/metagenome/3_MAG_gene_cluster/mag_SRR24442557
ASSEMBLY="$OUTDIR/assembly/final.contigs.fa"
DEPTH_FILE="$OUTDIR/depth.txt"
BINS_DIR="$OUTDIR/bins"

mkdir -p $BINS_DIR

#2. run metabat2 for binning

metabat2 -i $ASSEMBLY -a $DEPTH_FILE -o $BINS_DIR/bin
#-i → input contigs
#-a → depth/coverage file
#-o → output prefix for bins
#Output: bin.1.fa, bin.2.fa, etc. in $BINS_DIR
#Each file represents a MAG (draft genome).
#It grouped them into 3 distinct bins (bin.1.fa, bin.2.fa, bin.3.fa).
#Each bin likely represents a different microbial genome (or partial genome).

#Contigs	Short assembled sequences from reads
#Binning	Grouping contigs belonging to the same genome
#MAG (bin)	One reconstructed microbial genome from your sample



#3. check quality, but we skip this step for now because we cannot install checkm
checkm lineage_wf -x fa $BINS_DIR $OUTDIR/checkm_out
#-x fa → extension of your bins (FASTA)
#$OUTDIR/checkm_out → folder for CheckM results
#Ideally, Completeness > 90%, Contamination < 5% for high-quality MAGs.

################################4. gene cluster analysis#######################
#Once we have MAGs: Predict genes in each MAG (e.g., using Prokka); Identify biosynthetic gene clusters (BGCs) with antiSMASH:

#lets install via docker
docker pull antismash/standalone:7.1.0

#######################.  antismash for bin.1.fa ##########################
docker run --rm -u $(id -u):$(id -g) \
  -v "/Users/mm/Documents/project/metagenome/3_MAG_gene_cluster/mag_SRR24442557:/input" \
  -v "/Users/mm/Documents/project/metagenome/3_MAG_gene_cluster/mag_SRR24442557:/output" \
  antismash/standalone:7.1.0 \
  bins/bin.1.fa \
  --output-dir /output/antismash_bin1 \
  --genefinding-tool prodigal


#check
ls -lh /Users/mm/Documents/project/metagenome/3_MAG_gene_cluster/mag_SRR24442557/antismash_bin1
#File / Folder	Description
#index.html	Main interactive report — open this in your browser to view antiSMASH results (cluster maps, domains, similarities, etc.).
#bin.1.gbk	GenBank file with all annotated clusters and their predicted functions — can be opened in tools like Geneious or Artemis.
#bin.1.json	JSON format of the results (useful for scripting or downstream bioinformatics).
#k39_XXXXX.region001.gbk	Each of these represents one predicted biosynthetic gene cluster (BGC) — for example, a NRPS, PKS, terpene, etc.
#regions.js	Used internally by index.html for visualization.
#bin.1.zip	Compressed version of all the above (for easy sharing or archiving).

#open interactive report
open /Users/mm/Documents/project/metagenome/3_MAG_gene_cluster/mag_SRR24442557/antismash_bin1/index.html

#summarize BCG
#1. SUMMARY BY FILENAME
ls /Users/mm/Documents/project/metagenome/3_MAG_gene_cluster/mag_SRR24442557/antismash_bin1/*.region*.gbk
#that means antiSMASH successfully detected 8 biosynthetic gene clusters (BGCs) in your first MAG (bin.1.fa).
#summarize cluster rype
grep "product" /Users/mm/Documents/project/metagenome/3_MAG_gene_cluster/mag_SRR24442557/antismash_bin1/*.region*.gbk | sed 's/.*\/\(.*region[0-9]*\.gbk\):.*product="\([^"]*\)"/\1\t\2/' 
#create tsv file
grep "product" /Users/mm/Documents/project/metagenome/3_MAG_gene_cluster/mag_SRR24442557/antismash_bin1/*.region*.gbk \
| sed 's/.*\/\(.*region[0-9]*\.gbk\):.*product="\([^"]*\)"/\1\t\2/' \
> /Users/mm/Documents/project/metagenome/3_MAG_gene_cluster/mag_SRR24442557/antismash_bin1/antismash_cluster_bin1_summary.tsv
#open summary tsv
cat /Users/mm/Documents/project/metagenome/3_MAG_gene_cluster/mag_SRR24442557/antismash_bin1/antismash_cluster_bin1_summary.tsv




========================== Antismash for bin.2.fa ==========================

docker run --rm -u $(id -u):$(id -g) \
  -v "/Users/mm/Documents/project/metagenome/3_MAG_gene_cluster/mag_SRR24442557:/input" \
  -v "/Users/mm/Documents/project/metagenome/3_MAG_gene_cluster/mag_SRR24442557:/output" \
  antismash/standalone:7.1.0 \
  bins/bin.2.fa \
  --output-dir /output/antismash_bin2 \
  --genefinding-tool prodigal


#open interactive report
open /Users/mm/Documents/project/metagenome/3_MAG_gene_cluster/mag_SRR24442557/antismash_bin2/index.html
#no bcg found in bin.fa2
#If it still returns nothing, that confirms bin.2 has no detected biosynthetic gene clusters (BGCs) — which can happen if:
#The genome bin is too small (your bin.2.fa is only ~212 KB vs. bin.1.fa being 2.9 MB).
#It doesn’t contain any genes resembling known secondary metabolite pathways.

========================== Antismash for bin.3.fa ==========================

docker run --rm -u $(id -u):$(id -g) \
  -v "/Users/mm/Documents/project/metagenome/3_MAG_gene_cluster/mag_SRR24442557:/input" \
  -v "/Users/mm/Documents/project/metagenome/3_MAG_gene_cluster/mag_SRR24442557:/output" \
  antismash/standalone:7.1.0 \
  bins/bin.3.fa \
  --output-dir /output/antismash_bin3 \
  --genefinding-tool prodigal

#check
ls -lh /Users/mm/Documents/project/metagenome/3_MAG_gene_cluster/mag_SRR24442557/antismash_bin3

#check if clusters found
grep "cluster type" /Users/mm/Documents/project/metagenome/3_MAG_gene_cluster/mag_SRR24442557/antismash_bin3/*.gbk
#show no result. why?
#Short answer: your Docker command is fine. The empty grep output is because antiSMASH either found no BGCs in bin.3 or the JSON/GBK fields you searched for differ. Inspect the JSON structure and extract cluster info with these checks.

# show top-level keys
jq 'keys' /Users/mm/Documents/project/metagenome/3_MAG_gene_cluster/mag_SRR24442557/antismash_bin3/bin.3.json

# show .regions type and content (first ~200 chars)
jq -r '.regions | type, . | tostring[0:200]' /Users/mm/Documents/project/metagenome/3_MAG_gene_cluster/mag_SRR24442557/antismash_bin3/bin.3.json

# search JSON for cluster-related keys
grep -i -n '"cluster_types"\|"regions"\|"product"\|"region_id"' /Users/mm/Documents/project/metagenome/3_MAG_gene_cluster/mag_SRR24442557/antismash_bin3/bin.3.json || true

#Reason: your bin.3 JSON has no regions — antiSMASH produced no BGC entries for that input. The top-level keys show "records" and .regions is null, so grep for "cluster type" returns nothing.

#quick check

#show record count and keys
jq '.records | length, .records[] | keys' /Users/mm/Documents/project/metagenome/3_MAG_gene_cluster/mag_SRR24442557/antismash_bin3/bin.3.json

# show whether each record has regions and how many
jq -r '.records[] | {input_file: .input_file, regions_count: (if .regions then (.regions|length) else 0 end)}' /Users/mm/Documents/project/metagenome/3_MAG_gene_cluster/mag_SRR24442557/antismash_bin3/bin.3.json

# if regions_count > 0, list their ids, types and product
jq -r '.records[] | if (.regions and (.regions|length>0)) then (.regions[] | "\(.region_id)\t\((.cluster_types//[])|join(","))\t\(.product//"")") else empty end' /Users/mm/Documents/project/metagenome/3_MAG_gene_cluster/mag_SRR24442557/antismash_bin3/bin.3.json


#antiSMASH found no BGCs in bin.3 — evidence: bin.3.json .regions is null and your grep/jq returned nothing. Prodigal produced proteins (237), so gene calling worked; the absence of regions means antiSMASH simply didn’t detect any recognizable biosynthetic gene clusters in that bin (common: clusters can be missing, partial, or split across contigs).

#If you want to be certain or to look deeper, run these checks:
docker run --rm -u $(id -u):$(id -g) \
  -v "/Users/mm/Documents/project/metagenome/3_MAG_gene_cluster/mag_SRR24442557:/data" \
  antismash/standalone:7.1.0 \
  /data/bins/bin.3.fa \
  --output-dir /data/antismash_bin3 \
  --genefinding-tool prodigal 2>&1 | tee /tmp/antismash_bin3.log
tail -n 120 /tmp/antismash_bin3.log

#open
open /Users/mm/Documents/project/metagenome/3_MAG_gene_cluster/mag_SRR24442557/antismash_bin3/index.html
#result: no bcg found in bin 3.fa

#summary
#Why only bin.1 shows BGCs

#Bin size / content:  script notes bin.1 is much larger than bin.2 (bin.2 ≈ 212 KB vs bin.1 ≈ 2.9 MB). Smaller bins often lack full BGCs.
#Fragmentation / completeness: partial BGCs split across contigs or missing genes → no region called.
#Low coverage in subsample: subsampling reduces read depth and assembly contiguity, hurting BGC detection.
#Gene calling differences: Prodigal may call fewer/short proteins on short/fragmented contigs.
#antiSMASH sensitivity and thresholds: some clusters are borderline and need more context to be reported.


# 1) Check bin sizes
ls -lh mag_SRR24442557/bins

# 2) Inspect JSON for regions (already in MAG.sh, but run manually)
jq '.[].regions | type, .[0:5]' 3_MAG_gene_cluster/mag_SRR24442557/antismash_bin3/bin.3.json || true

# 3) Count predicted proteins (run Prodigal locally to confirm gene calling)
prodigal -i 3_MAG_gene_cluster/mag_SRR24442557/bins/bin.2.fa -a /tmp/bin2.proteins.faa -p single

# 4) Re-run antiSMASH on full dataset / larger bin set (avoid subsample) or relax detection params
# and re-run binning with different parameters (Metabat2, CONCOCT, VAMB) to get less fragmented bins.

#Recommendation: verify bin sizes and gene counts first (commands above). If bins are small/fragmented, run full reads (no subsample) or re-bin/merge contigs before re-running antiSMASH.

#conclusion
#Conclusion — why only bin.1 has BGCs

#bin.1 is ~2.9 MB and contains long, contiguous contigs so antiSMASH can detect multi-gene BGCs.
#bin.2 and bin.3 are ~200 KB each — too small/fragmented for most BGCs (BGCs often span ≥10–20 kb and many genes).
#we used a subsampled read set (1M reads per file). Subsampling lowers coverage and contig lengths → fewer/fragmented BGCs recovered.




#conclusion:

#AntiSMASH Purpose:
#AntiSMASH is specialized in secondary metabolite biosynthesis.
#It looks for biosynthetic gene clusters (BGCs), e.g., for antibiotics, pigments, toxins, siderophores, etc.
#Input:
#Annotated genomes or bins (your Prokka outputs).
#Output:
#Identified gene clusters (BGCs) with predicted chemical products.
#Reports the type of cluster (NRPS, PKS, RiPPs, etc.) and sometimes similarity to known compounds.
#Key point:
#Very specific, focused on specialized metabolism (secondary metabolites).
#Doesn’t give a global functional profile of all genes.


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


#AntiSMASH tells you: “These bins can potentially make antibiotics or other secondary metabolites.”
#EggNOG-mapper tells you: “These bins have genes for metabolism, transport, replication, stress response, etc











==================================================================================
###############END OF ANALYSIS FOR GENE CLUSTER FROM SUBSAMPLED READS########################
===================================================================================

#Lets do functional annotation
pwd
cd /Users/mm/Documents/project/metagenome/

touch functional_annotation.sh