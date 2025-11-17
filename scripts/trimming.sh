# =========================================================
# Step 2:  trimming and quality control of trimming reads using fastp
# =========================================================

# Define input and output directories
RAW_DIR="0_rawdata"
TRIMMED_DIR="1_quality_control/trimmed_reads"
REPORT_DIR="1_quality_control/qc_reports_trimmed"

# Create output folders if not exist
mkdir -p ${TRIMMED_DIR}
mkdir -p ${REPORT_DIR}

# Loop through all FASTQ pairs in the raw data folder
for SAMPLE in $(ls ${RAW_DIR}/*_1.fastq | sed 's/_1.fastq//' | xargs -n 1 basename); do
    echo "Processing sample: ${SAMPLE}"
    
    fastp \
        -i ${RAW_DIR}/${SAMPLE}_1.fastq \
        -I ${RAW_DIR}/${SAMPLE}_2.fastq \
        -o ${TRIMMED_DIR}/${SAMPLE}_1.trimmed.fastq \
        -O ${TRIMMED_DIR}/${SAMPLE}_2.trimmed.fastq \
        --html ${REPORT_DIR}/${SAMPLE}_fastp.html \
        --json ${REPORT_DIR}/${SAMPLE}_fastp.json \
        --detect_adapter_for_pe \
        --thread 4
done

#lets summarize
multiqc 1_quality_control/qc_reports_trimmed -o 1_quality_control/qc_reports_trimmed





############################
#lets upload step 1 and 2 (raw data dan quality control) to github
############################

#create a new repository on github named biofloc-metagenome-project

cd ~/Documents/project/metagenome
nano .gitignore

# Ignore macOS system files
.DS_Store

# Ignore raw data files
0_rawdata/*.fastq
0_rawdata/*.fastq.gz
0_rawdata/*.sra
0_rawdata/*.zip
0_rawdata/*.tar.gz

# Ignore trimmed read files
1_quality_control/trimmed_reads/*.fastq
1_quality_control/trimmed_reads/*.fastq.gz

# Ignore QC report zips
1_quality_control/qc_reports/*.zip

# Ignore temporary and backup files
*.tmp
*.bak

git init
git add. 
git commit -m "Initial commit: full project structure with preprocessing, QC, and trimming scripts"

git remote add origin https://github.com/novayurikaa/biofloc-metagenome-project.git
git branch -M main
git push -u origin main
git push origin main --force #lets cancel, too many big files


echo "0_rawdata/" >> .gitignore
echo "*.zip" >> .gitignore
git add .gitignore
git commit -m "Ignore large raw data and zip files"
git push origin main --force #next 


touch taxonomy_functional.sh

