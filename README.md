# Biofloc metagenomics shotgun analysis project


## Background Statement
Biofloc technology enhances aquaculture sustainability by leveraging microbial
food webs to convert nitrogenous waste into beneficial biomass while reducing
water exchange. Understanding microbial taxa and their encoded functions is key
to designing interventions (e.g., probiotics) that promote desirable traits.

This project investigates the microbial community structure and functional potential of biofloc aquaculture using shotgun metagenomic sequencing from shotgun metagenomic analysis of biofloc aquaculture microbiome (Rajeev et al., 2023).

### This project applies a full unix-based shotgun metagenomics workflow, including: 

1. Quality control and trimming (FastQC, fastp)
2. Taxonomic profiling (Kraken2 + Bracken)
3. Functional annotation (HUMAnN3)
4. Genome assembly+draft MAG recovery, small scale (MEGAHIT, MetaBAT2, CheckM)
   1 Sample for demonstration due to local compute limits
6. Biosynthetic gene cluster (BGC) discovery (antiSMASH/DeepBCG)

   NOTE: MAG and BCG run at reduced scale on local macOS. Full scale implementation could be expanded using HPC during future PhD work. 

### Research Questions
### 1. (community composition)
 #### How do microbial community composition differ between multiple biofloc?
 #### which microbial taxa dominate the biofloc microbiome?

  
### 2. (Functional potential)
#### Which functional pathway relevant to water quality and host health are enriched in different biofloc samples?

### 3. Can we recover draft genomes of key ecological players?
### 4. What are key biosynthetic gene clusters that may benefit aquaculture performance?

### Data source


We use three real shotgun metagenomic samples from a biofloc aquaculture system
Rajeev et al. 2023 — Metagenome sequencing and recovery of 444 MAGs from a biofloc aquaculture system. DOI: https://doi.org/10.1038/s41597-023-02622-0

### BioProject PRJNA967453:

| Sample | Accession | Description |
|--------|-----------|-------------|
| Biofloc S1 | SRR24442555 | Early stage microbiome |
| Biofloc S2 | SRR24442556 | Mid development microbiome |
| Biofloc S3 | SRR24442557 | Later stage / treatment microbiome |


### Progress Checklist
Setup
- [x] GitHub repository
- [x] Project directory structure

Analysis
- [x] Download SRA metadata & raw reads
- [x] QC + MultiQC report
- [ ] Taxonomic profiles + plots
- [ ] Functional profiles + pathway visualization

 MAG + BGC 
- [ ] Single-sample assembly + binning
- [ ] CheckM quality
- [ ] antiSMASH / DeepBGC on recovered MAGs
- [ ] Report BGC traits

Interpretation
- [ ] Ecological insights summary in `docs/interpretation.md`

### Nova Yurika 
Marine Science + Bioinformactics
GitHub: https://github.com/novayurikaa
