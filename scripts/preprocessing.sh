#!/bin/bash
# =======================================
# Metagenomic analysis of biofloc data, using VS code
# Step 1: Create folder structure
# =======================================

echo "Creating main folders..."

mkdir -p 0_raw_data/raw
mkdir -p 1_quality_control/trimmed_reads
mkdir -p 1_quality_control/qc_reports
mkdir -p 1_quality_control/qc_reports_trimmed

mkdir -p 2_taxonomy_functional
mkdir -p 3_MAG_gene_cluster

chmod +x preprocessing. sh
./preprocessing.sh
===================================
# step 2 install  toolkit
===================================
cd ~/Documents/project/metagenome/0_rawdata
===================================
#1. download tools
==================================
conda install mamba -n base -c conda-forge
mamba --version
mamba install -c bioconda sra-tools
conda activate bioinfo #we are inside bioinfo
conda install -c bioconda sra-tools
====================================================
#step 3 lets download raw data from SRA accession number
==================================================

for sra in SRR24442555 SRR24442556 SRR24442557; do
  prefetch $sra
  fasterq-dump $sra -O .
done
#I forgot to cd,  install it in raw data folder, so lets move it manually to raw data folder using cursor on left page of vs code

#finish, lets move to next step, quality control
touch quality_control.sh #create new script for quality control
