
**Biofloc Metagenomics Shotgun Analysis Project**

**A. Background Statement**

Biofloc technology (BFT) enhances aquaculture sustainability by leveraging microbial food webs to convert nitrogenous waste into beneficial biomass while reducing water exchange. Understanding the successional shifts in microbial taxa and their encoded functions is key to designing interventions (e.g., probiotics) that promote system stability and pathogen resistance.

This project investigates the microbial community structure and functional potential of a biofloc aquaculture system using three samples representing early, mid, and late biofloc stage, without biological replication, using shotgun metagenomic sequencing data from Rajeev et al. (2023).

**B. Unix-Based Shotgun Metagenomics Workflow**

**1. Quality Control:** FastQC and fastp for read trimming and filtering.

**2. Taxonomic Profiling:** Kraken2 + Bracken for species-level quantification across successional stages.

**3. Genome Assembly: ** De novo assembly of metagenomic reads using MEGAHIT.

**4. MAG Recovery:** Draft Metagenome-Assembled Genome (MAG) recovery using MetaBAT2.

**5. Functional Annotation ** (MAG-Centric):

**6. EggNOG-mapper:** Evolutionary genealogy of genes: Non-supervised Orthologous Groups for COG/KEGG categorization.

**7. Prokka: ** Rapid prokaryotic genome annotation.

****8. Secondary Metabolism: ** Biosynthetic Gene Cluster (BGC) discovery using antiSMASH.

**Note:** Due to local computational limits, MAG recovery and BGC discovery were focused on a representative subsample of the Late-stage biofloc (SRR24442557). Full-scale implementation for all samples is planned for future High-Performance Computing (HPC) environments.

**C. Research Questions & Findings**
**1. Microbial Succession and Composition**

**Question: How does the microbial community evolve from setup to maturity?**

Finding: The system follows a clear ecological succession. Early stages are dominated by pioneer heterotrophs (Rhodobacter sp.), transitioning into a highly diverse climax community in the late stage (Sample 557).

Diversity: The late stage achieved a peak Shannon Diversity Index > 6.0 and richness of >3,000 species.

Dominant Taxa: Phaeobacter gallaeciensis was identified as the primary biological driver, consistently maintaining high abundance.

**2. Functional Potential of the Climax Community**


**Question: What metabolic strategies allow dominant microbes to maintain the biofloc?
**
Finding: Functional profiling of the dominant MAG (Bin 1: ECFFIKGN) revealed a high density of TonB-dependent receptors. This indicates a specialized ecological niche focused on high-efficiency nutrient scavenging and iron acquisition from the water column.

**3. Biosynthetic Gene Clusters (BGCs)**

Question: Does the biofloc microbiome provide chemical defense for aquaculture species?

Finding: Bin 1 acts as a "biosynthetic powerhouse," harboring over 80 BGCs. The dominance of NRPS-like (Non-Ribosomal Peptide Synthetase) clusters suggests a high capacity for producing antimicrobial compounds and siderophores, contributing to the system's biosecurity.

**D. Downstream Analysis in R**

The following analyses are available in the repository scripts, translating bioinformatics outputs into biological insights:

Taxonomic & Diversity Analysis

Alpha Diversity: Visualization of Richness and Shannon Index across stages.

Beta Diversity: Hierarchical Clustering (UPGMA) and Heatmaps based on Bray-Curtis dissimilarity to identify significant community shifts.

Dominant Taxa Profiling: Targeted heatmaps of the Top 20 Species tracking the bloom of beneficial taxa.

Functional & MAG Analysis (Bin 1 Focus)

Secondary Metabolite Profiling: Donut charts and bar plots of BGCs, highlighting NRPS-like cluster dominance.

Core Functional Gene Analysis: Frequency bar plots of cleaned functional annotations identifying TonB-dependent receptors.

Ecological Niche Mapping: Stacked bar plots identifying genes related to Nitrogen cycling, Phosphate remediation, and Oxidative stress defense.

**E. Data Source**

The project utilizes three shotgun metagenomic samples from BioProject PRJNA967453.

Sample	Accession	Succession Stage:

Biofloc S1	SRR24442555	Early stage microbiome, 

Biofloc S2	SRR24442556	Mid development microbiome, 

Biofloc S3	SRR24442557	Late stage (Climax community), 

Reference: Rajeev et al. (2023) :  Metagenome sequencing and recovery of 444 MAGs from a biofloc aquaculture system. DOI: 10.1038/s41597-023-02622-0


**Author
Nova Yurika
Marine Science & Bioinformatics
GitHub Profile**

