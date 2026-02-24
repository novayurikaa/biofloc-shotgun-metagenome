# === 1. Define project paths ===
project_dir <- "Biofloc_Metagenome_Analysis"
scripts_dir <- file.path(project_dir, "scripts")
results_dir <- file.path(project_dir, "results")
data_dir <- file.path(project_dir, "data")
metadata_dir <- file.path(project_dir, "metadata")




# Create folders if they don't exist
dirs <- c(project_dir, scripts_dir, results_dir, data_dir, metadata_dir)
for(d in dirs) if(!dir.exists(d)) dir.create(d, recursive = TRUE)

# === 2. Create metadata file ===
metadata_file <- file.path(metadata_dir, "samples.csv")
metadata_content <- "Sample,Phase
SRR24442555,Early
SRR24442556,Mid
SRR24442557,Late"
writeLines(metadata_content, metadata_file)

# === 3. Write the analysis script (bracken.R) ===
script_file <- file.path(scripts_dir, "bracken.R")
script_content <- "
# Load required packages
library(tidyverse)
library(pheatmap)
library(vegan)

# Load Bracken data
bracken <- read_tsv('bracken_species.tsv')


read.delim("~/Documents/project/rbiofloc/rbiofloc/bracken_species.tsv", sep = "\t", header = TRUE)

# Load metadata
metadata <- read_csv('metadata/samples.csv')

# Keep only species-level data
bracken <- bracken %>%
  filter(Rank == 'S') %>%
  select(Sample, Taxon, Reads)

# Compute relative abundance
bracken <- bracken %>%
  group_by(Sample) %>%
  mutate(rel_abundance = Reads / sum(Reads)) %>%
  ungroup()

# Merge with metadata
bracken <- bracken %>%
  left_join(metadata, by='Sample')

# Top 10 species
top_species <- bracken %>%
  group_by(Taxon) %>%
  summarise(mean_abundance = mean(rel_abundance)) %>%
  top_n(10, mean_abundance) %>%
  pull(Taxon)

bracken_top <- bracken %>% filter(Taxon %in% top_species)

# 1. Stacked bar plot
ggplot(bracken_top, aes(x=Sample, y=rel_abundance, fill=Taxon)) +
  geom_bar(stat='identity') +
  facet_wrap(~Phase, scales='free_x') +
  theme_bw() +
  ylab('Relative abundance') +
  xlab('Sample') +
  ggtitle('Top species composition in biofloc ecosystem') +
  ggsave('results/species_composition_barplot.png', width=8, height=6)

# 2. Heatmap

library(tidyverse)

species_matrix <- bracken %>%
  # Sum abundance for duplicate Taxon/Sample pairs
  group_by(Sample, Taxon) %>%
  summarize(rel_abundance = sum(rel_abundance), .groups = "drop") %>%
  # Now pivot
  pivot_wider(names_from = Sample, values_from = rel_abundance, values_fill = 0) %>%
  column_to_rownames('Taxon')


species_matrix <- bracken %>%
  pivot_wider(names_from = Sample, values_from = rel_abundance, values_fill = 0) %>%
  column_to_rownames('Taxon')

annotation_col <- metadata %>% column_to_rownames('Sample')
pheatmap(species_matrix,
         cluster_rows = TRUE, cluster_cols = TRUE,
         annotation_col = annotation_col,
         filename = 'results/species_heatmap.png',
         width = 8, height = 6)


#gagal, once again




library(tidyverse)

# 1. Summarize: Sum abundances for duplicate Taxon names within each Sample
cleaned_bracken <- bracken %>%
  group_by(Sample, Taxon) %>%
  summarise(rel_abundance = sum(rel_abundance), .groups = "drop")

# 2. Pivot: Wider format
species_matrix <- cleaned_bracken %>%
  pivot_wider(names_from = Sample, 
              values_from = rel_abundance, 
              values_fill = 0)

# 3. Finalize: Make Taxon row names
species_matrix <- species_matrix %>%
  column_to_rownames('Taxon')

# 4. Optional: Only show top 30-50 species for a clearer heatmap
# top_species <- rowSums(species_matrix) %>% sort(decreasing = TRUE) %>% head(50) %>% names()
# species_matrix <- species_matrix[top_species, ]

# 5. Plot heatmap (e.g., using pheatmap)
# pheatmap::pheatmap(species_matrix, scale = "row", show_rownames = TRUE)
#gagal
===
# 1. Remove rows where the standard deviation is 0
species_matrix_filtered <- species_matrix[apply(species_matrix, 1, sd) > 0, ]

# 2. Now run pheatmap
pheatmap::pheatmap(species_matrix_filtered, 
                   scale = "row", 
                   show_rownames = TRUE)
==
# Calculate total abundance per taxon and pick top 50
top_50_taxa <- rowSums(species_matrix_filtered) %>% 
  sort(decreasing = TRUE) %>% 
  head(50) %>% 
  names()

species_top <- species_matrix_filtered[top_50_taxa, ]

pheatmap::pheatmap(species_top, 
                   scale = "row", 
                   clustering_distance_rows = "correlation", # Often better for bio data
                   main = "Top 50 Most Abundant Taxa")
                   
                  #we did it
                  

========================
# 3. Shannon diversity

# Calculate Shannon Diversity using your 'Reads' column
diversity_data <- bracken %>%
    group_by(Sample) %>%
    summarise(
      shannon = vegan::diversity(pmax(Reads, 0), index = 'shannon'),
      richness = sum(Reads > 0) # Count how many taxa have > 0 reads
    ) %>%
    left_join(metadata, by = 'Sample')

# View the result to make sure it worked
print(diversity_data)
#we did it


#visualise
library(ggplot2)

# Ensure Phase is in the correct chronological order
diversity_data$Phase <- factor(diversity_data$Phase, levels = c("Early", "Mid", "Late"))

ggplot(diversity_data, aes(x = Phase, y = shannon, group = 1)) +
  geom_line(color = "steelblue", size = 1.2) + # Connects the phases
  geom_point(aes(color = Phase), size = 5) +    # Highlights the data points
  geom_text(aes(label = round(shannon, 2)), vjust = -1.5) + # Shows exact value
  theme_classic() +
  ylim(min(diversity_data$shannon)-0.5, max(diversity_data$shannon)+0.5) +
  labs(title = "Microbial Diversity Succession",
       subtitle = "Shannon Index across Biofloc Development Phases",
       y = "Shannon Index (H')",
       x = "Sampling Phase") +
  theme(legend.position = "none")


#shanon vs richness


library(patchwork) # For combining plots
install.packages("patchwork")
library(patchwork)

# Plot 1: Shannon
p1 <- ggplot(diversity_data, aes(x = Phase, y = shannon, fill = Phase)) +
  geom_col(alpha = 0.8) +
  labs(y = "Shannon Index") +
  theme_minimal() + theme(legend.position = "none")

# Plot 2: Richness
p2 <- ggplot(diversity_data, aes(x = Phase, y = richness, fill = Phase)) +
  geom_col(alpha = 0.8) +
  labs(y = "Species Richness (Count)") +
  theme_minimal() + theme(legend.position = "none")

# Display side-by-side
p1 + p2

###see if these biofloc have eukariote
library(dplyr)
library(ggplot2)

# 1. Classify Taxa into broad groups
kingdom_summary <- bracken %>%
  mutate(Group = case_when(
    # Search for common Algal identifiers
    grepl("Chlorella|Chlamydomonas|Bacillariophyta|Phaeophyceae|Chlorophyta", Taxon, ignore.case = TRUE) ~ "Algae",
    # Search for Fungal identifiers
    grepl("Saccharomyces|Aspergillus|Penicillium|Candida|Fungi", Taxon, ignore.case = TRUE) ~ "Fungi",
    # Most Bracken files are predominantly Bacteria
    TRUE ~ "Bacteria" 
  )) %>%
  group_by(Sample, Phase, Group) %>%
  summarize(TotalReads = sum(Reads), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(Percentage = (TotalReads / sum(TotalReads)) * 100)

# 2. View the result
print(kingdom_summary)
# Ensure Phase order
kingdom_summary$Phase <- factor(kingdom_summary$Phase, levels = c("Early", "Mid", "Late"))

ggplot(kingdom_summary, aes(x = Phase, y = Percentage, fill = Group)) +
  geom_bar(stat = "identity", width = 0.7) +
  theme_minimal() +
  scale_fill_manual(values = c("Bacteria" = "#1f77b4", "Algae" = "#2ca02c", "Fungi" = "#d62728")) +
  labs(title = "Biofloc Composition: Bacteria vs. Algae vs. Fungi",
       subtitle = "Relative Abundance based on Estimated Reads",
       y = "Percentage (%)",
       x = "Sampling Phase")



####result

SRR24442555 (Early): Likely dominated by fast-growing Bacteria.
SRR24442556 (Mid): The slight dip in Shannon diversity often indicates a "Bloom" where one specific group (could be a specific Algae or Bacteria) has become very dominant.
SRR24442557 (Late): The high diversity suggests a complex "Mature" floc with a balance of many different organisms.


