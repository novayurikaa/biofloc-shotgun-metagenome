=====================================================
####### visualisation and interpretation of results ########
===================================================
mkdir -p /Users/mm/Documents/project/metagenome/5_visualization_and_interpretation

#copy antismash results to 4_functional_annotation folder too fist, along with eggnog_annotations

cd /Users/mm/Documents/project/metagenome/5_visualization_and_interpretation

#activate R environment
#OPEN NEW TERMINAL, type R

==========================
# Load required packages
if (!require("tidyverse")) install.packages("tidyverse", repos="https://cran.rstudio.com/")
library(tidyverse)

# ----------------------------
# Set paths
# ----------------------------
antismash_dir <- "/Users/mm/Documents/project/metagenome/4_functional_annotation/antismash"
eggnog_dir <- "/Users/mm/Documents/project/metagenome/4_functional_annotation/eggnog_annotations"
output_dir <- "/Users/mm/Documents/project/metagenome/visualization_and_interpretation"

# Create output folder if it does not exist
dir.create(output_dir, showWarnings = FALSE)

# ----------------------------
# 1. Antismash summary
# ----------------------------
# List all .gbk files in subfolders
antismash_files <- list.files(antismash_dir, pattern = "\\.gbk$", recursive = TRUE, full.names = TRUE)

antismash_summary <- map_dfr(antismash_files, function(f) {
  bin_name <- basename(dirname(f)) # get folder name as bin
  tibble(
    bin = bin_name,
    gbk_file = basename(f),
    n_clusters = length(antismash_files) # or extract from file if needed
  )
})

# Write Antismash summary CSV
write_csv(antismash_summary, file.path(output_dir, "antismash_summary.csv"))

##Here’s what each column means:
#bin → which MAG/bin the file came from (antismash_bin1, antismash_bin2, …).
#gbk_file → the individual GenBank output files produced by Antismash for that bin.
#n_clusters → the number of biosynthetic gene clusters (BGCs) predicted in that file.
#So for example, in your file:
#antismash_bin1  bin.1.gbk  9
#means bin 1 has a main GenBank file (bin.1.gbk) containing 9 BGCs, and the other rows are region-specific GenBank files (subclusters) also reporting 9 clusters — sometimes Antismash splits clusters into regions.


# ----------------------------
# 2. EggNOG-mapper summary
# ----------------------------
# List all annotation files
annotation_files <- list.files(eggnog_dir, pattern = "*.emapper.annotations$", full.names = TRUE)

eggnog_summary <- map_dfr(annotation_files, function(f) {
  bin_name <- basename(f) %>% str_remove(".emapper.annotations")
  df <- read_tsv(f, comment = "#") # skip comment lines
  tibble(
    bin = bin_name,
    n_genes = nrow(df)
  )
})

# Write EggNOG summary CSV
write_csv(eggnog_summary, file.path(output_dir, "eggnog_summary.csv"))

# ----------------------------
# Finished
# ----------------------------
cat("Summary CSV files created in:", output_dir, "\n")


====================================================================
library(tidyverse)

# EggNOG functional annotation
egg_df <- read_tsv("/path/to/all_bins_annotation.tsv")


# AntiSMASH cluster summary
antismash_df <- read_tsv("/path/to/antismash_cluster_bin1_summary.tsv")
# Repeat for other bins if needed or combine into one dataframe

egg_summary <- egg_df %>%
  group_by(bin, COG_category) %>%
  summarise(n_genes = n(), .groups = "drop")

antismash_summary <- antismash_df %>%
  group_by(bin, bgc_type) %>%
  summarise(n_clusters = n(), .groups = "drop")

write_csv(egg_summary, "/path/to/visualization_and_interpretation/eggnog_summary.csv")
write_csv(antismash_summary, "/path/to/visualization_and_interpretation/antismash_summary.csv")


ggplot(egg_summary, aes(x = bin, y = n_genes, fill = COG_category)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Functional categories per MAG", y = "Number of genes")


ggplot(antismash_summary, aes(x = bin, y = n_clusters, fill = bgc_type)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "BGC types per MAG", y = "Number of clusters")


