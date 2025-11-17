

# ==============================================
# Quality Control of Biofloc Metagenomic Data
# Step 2: Run FastQC and MultiQC
# ==============================================

# Define directories
RAW_DIR="0_rawdata"
QC_DIR="1_quality_control/qc_reports"
TRIMMED_DIR="1_quality_control/trimmed_reads"

# Create output folders if not exist
mkdir -p $QC_DIR
mkdir -p $TRIMMED_DIR


# Run FastQC on all FASTQ files
fastqc $RAW_DIR/*.fastq -o $QC_DIR

# Combine all FastQC results using MultiQC

multiqc $QC_DIR -o $QC_DIR

#lets do trimming 
touch trimming.sh
