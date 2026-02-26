setwd("~/Documents/project/rbiofloc/rbiofloc")

#we focus on
#1.⁠ ⁠Functional capabilities of dominant microbes (EggNOG functional categories).
#2.	Biosynthetic gene clusters and their ecological/biotech relevance (antiSMASH).

==========
  #1. load libraries
=================

library (vegan)
library (readxl)
library (readr)
library (dplyr)
library (tidyr)
library (tibble)
library (ggplot2)
library (pheatmap)
library (reshape2)

#I only analysed sub sample of late biofloc (....57) due to computational limit
#we cannot analyse three samples and make community scale interpertation of its functional profiling
#however we can still analyse:
#1. Potential secondaty metabolite (analysis):AntiSMASH
#we have files input: antismash_summary_csv , antismash_cluster_bin1_summary.tsv
#we aim to get the  visualisation of numbers of BCG per bin, types of clusers in tables or plot
#2. gene fucntional: EggNOG
#we have files eggnog_summary.tsv , all_bins_annotations.tsv
#we expect to create  gene funtional profile per bin (KGG, GO, COG) total per category. by heatmap and stacked bar
#3. MAG completeness/bin size 
#using file bins/bin.S.fa


# ============================================
# 2. AntiSMASH: BGC summary per bin
# ============================================
# Input: antismash_summary.csv
# Columns example: Bin, Num_BGC, Dominant_Type
antismash <- read_csv("~/Documents/project/rbiofloc/rbiofloc/antismash_summary.csv")

# Quick preview
head(antismash)

# Make a summary per bin
# Sum all clusters (if n_clusters already reflects total per file)
antismash_summary_per_bin <- antismash %>%
  group_by(bin) %>%
  summarise(Total_BGCs = sum(n_clusters), .groups = "drop")


# Check
antismash_summary_per_bin

# Plot per bin, this shows which bins have the most BGC. That is BIN 1
ggplot(antismash_summary_per_bin, aes(x = bin, y = Total_BGCs, fill = bin)) +
  geom_bar(stat = "identity") +
  labs(title = "Total BGCs per bin", y = "Number of BGCs", x = "MAG Bin") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")



antismash <- tibble(
  Bin = c("antismash_bin1","antismash_bin2","antismash_bin3"),
  n_clusters = c(9,1,1)
)

# Quick look
antismash
ggplot(antismash, aes(x = Bin, y = n_clusters, fill = Bin)) +
  geom_bar(stat = "identity") +
  labs(title = "Number of BGCs per MAG bin",
       x = "MAG Bin",
       y = "Number of biosynthetic gene clusters (BGCs)") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = "none")

#pie chart

antismash %>%
  mutate(perc = n_clusters / sum(n_clusters) * 100) %>%
  ggplot(aes(x = "", y = perc, fill = Bin)) +
  geom_col() +
  coord_polar(theta = "y") +
  labs(title = "Proportion of BGCs per MAG bin") +
  theme_void() +
  scale_fill_brewer(palette = "Set2")


# AntiSMASH detailed clusters for bin.1

antismash_bin1 <- read_tsv("~/Documents/project/rbiofloc/rbiofloc/antismash_cluster_bin1_summary.tsv") # Read in the bin1 detailed clusters


# Rename columns for easier handling
colnames(antismash_bin1) <- c("Region_ID", "Cluster_type")

# Quick check
head(antismash_bin1)

# Summarize number of clusters per type
cluster_summary <- antismash_bin1 %>%
  group_by(Cluster_type) %>%
  summarise(Num_clusters = n(), .groups = "drop")

cluster_summary

# Plot BGC cluster types in bin 1
ggplot(cluster_summary, aes(x = Cluster_type, y = Num_clusters, fill = Cluster_type)) +
  geom_bar(stat = "identity") +
  labs(title = "BGC Cluster Types in bin.1", y = "Number of Clusters", x = "Cluster Type") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# ============================================
# 3. EggNOG: functional profile per bin
# ============================================
# Input: eggnog_summary.tsv
#  Read EggNOG summary as CSV (comma-separated)
eggnog <- read_csv("~/Documents/project/rbiofloc/rbiofloc/eggnog_summary.csv", show_col_types = FALSE)

# Quick check
eggnog


# rename columns for clarity
colnames(eggnog) <- c("Bin", "Num_genes")

#  Plot total genes per bin
ggplot(eggnog, aes(x = Bin, y = Num_genes, fill = Bin)) +
  geom_bar(stat = "identity") +
  labs(title = "Total genes per bin (EggNOG)", y = "Number of genes", x = "MAG Bin") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")  #i think i dont need this



#EggNOG annotations

# Load all bins annotations
egg_annotations <- read_tsv("~/Documents/project/rbiofloc/rbiofloc/all_bins_annotations.tsv") #FAIL


egg_annotations <- read_tsv(
  "~/Documents/project/rbiofloc/rbiofloc/all_bins_annotations.tsv",
  comment = "##",
  show_col_types = FALSE
)

# Check first few rows
head(egg_annotations)
egg_annotation

cog_summary <- egg_annotations %>%
  group_by(Bin, COG_category) %>%
  summarise(Num_genes = n(), .groups = "drop")

# Check summary
head(cog_summary)







####summarize EggNOG gene functions

# Load EggNOG annotations
egg_annotations <- read_tsv(
  "~/Documents/project/rbiofloc/rbiofloc/all_bins_annotations.tsv",
  comment = "##",
  show_col_types = FALSE
)

# Select only the columns we need
egg_fun <- egg_annotations %>%
  select(COG_category, KEGG_ko, Description)

# Example: summarize number of genes per COG category
cog_summary <- egg_fun %>%
  filter(!is.na(COG_category)) %>%
  group_by(COG_category) %>%
  summarise(Num_genes = n(), .groups = "drop")

# Make a simple matrix for heatmap
# Here, since we only have one sample, it will be 1 column
heatmat <- cog_summary %>%
  column_to_rownames("COG_category") %>%
  as.matrix()

# Plot heatmap
pheatmap(
  heatmat,
  cluster_rows = TRUE,
  cluster_cols = FALSE, # only one sample
  display_numbers = TRUE,
  main = "Functional gene abundance (COG categories) in subsample late biofloc"
)
#not informative


######heatmap by kegg pathway

kegg_summary <- egg_fun %>%
  filter(!is.na(KEGG_ko)) %>%
  group_by(KEGG_ko) %>%
  summarise(Num_genes = n(), .groups = "drop")

heatmat_kegg <- kegg_summary %>%
  column_to_rownames("KEGG_ko") %>%
  as.matrix()

pheatmap(
  heatmat_kegg,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  main = "Functional gene abundance (KEGG pathways) in late biofloc"
)
#not informatvie




egg_annotations <- read_tsv(
  "~/Documents/project/rbiofloc/rbiofloc/all_bins_annotations.tsv",
  comment = "##",
  show_col_types = FALSE
)

# Check columns
colnames(egg_annotations)
# Key: COG_category, KEGG_ko, KEGG_Pathway, GOs, EC, Description

# ============================================
#  Summarize by functional categories
# ============================================

# A. COG categories
cog_summary <- egg_annotations %>%
  group_by(COG_category) %>%
  summarise(Num_genes = n(), .groups = "drop") %>%
  arrange(desc(Num_genes))

write_csv(cog_summary, "~/Documents/project/rbiofloc/rbiofloc/COG_summary_genelevel.csv")

# Barplot for COG categories
ggplot(cog_summary, aes(x = reorder(COG_category, Num_genes), y = Num_genes, fill = COG_category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "COG Functional Categories (Late Biofloc)",
       x = "COG Category", y = "Number of genes") +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() #not informative

# B. KEGG pathways
kegg_summary <- egg_annotations %>%
  separate_rows(KEGG_Pathway, sep = ";") %>%
  group_by(KEGG_Pathway) %>%
  summarise(Num_genes = n(), .groups = "drop") %>%
  arrange(desc(Num_genes))

write_csv(kegg_summary, "~/Documents/project/rbiofloc/rbiofloc/KEGG_summary_genelevel.csv")

# Top 20 KEGG pathways barplot
top_kegg <- kegg_summary %>% slice_max(order_by = Num_genes, n = 20)

ggplot(top_kegg, aes(x = reorder(KEGG_Pathway, Num_genes), y = Num_genes, fill = Num_genes)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top 20 KEGG Pathways (Late Biofloc)",
       x = "KEGG Pathway", y = "Number of genes") +
  scale_fill_gradient(low = "skyblue", high = "darkblue") +
  theme_minimal() #not informative


         




##############dont use heatmap.it is for comparingsamples.
#we need to focus on BIN 1
#WE NEED TO improve the eggng visualisation


#A. Fix the COG categories

# Create a lookup table for COG descriptions
cog_lookup <- data.frame(
  COG_category = c("J", "A", "K", "L", "B", "D", "Y", "V", "T", "M", "N", "Z", "W", "U", "O", "C", "G", "E", "F", "H", "I", "P", "Q", "S"),
  Description = c("Translation", "RNA processing", "Transcription", "Replication/Repair", "Chromatin structure", "Cell cycle", "Nuclear structure", "Defense mechanisms", "Signal transduction", "Cell wall/membrane", "Cell motility", "Cytoskeleton", "Extracellular structures", "Intracellular trafficking", "Posttranslational mod", "Energy production", "Carbohydrate metabolism", "Amino acid metabolism", "Nucleotide metabolism", "Coenzyme metabolism", "Lipid metabolism", "Inorganic ion metabolism", "Secondary metabolites", "Function unknown")
)

# Join and Plot
cog_summary_named <- cog_summary %>%
  left_join(cog_lookup, by = "COG_category") %>%
  filter(!is.na(Description), Description != "Function unknown") # Remove noise

ggplot(cog_summary_named, aes(x = reorder(Description, Num_genes), y = Num_genes, fill = Description)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  guides(fill = "none") +
  labs(title = "Functional Potential of Late Biofloc MAGs", x = "Functional Category", y = "Gene Count")


B. Visualizing the biotech relevance
#In biofloc, you are likely looking for Siderophores (iron scavenging) or Bacteriocins (antimicrobial activity).



# Data for Bin 1 specifically
ggplot(cluster_summary, aes(x = 2, y = Num_clusters, fill = Cluster_type)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) + # Creates the "hole" in the middle
  theme_void() +
  labs(title = "BGC Diversity in Bin 1", subtitle = "Dominant Secondary Metabolites")



# Using your cluster_summary from antiSMASH_bin1
ggplot(cluster_summary, aes(x = reorder(Cluster_type, Num_clusters), y = Num_clusters, fill = Cluster_type)) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  labs(
    title = "Biosynthetic Potential of Bin 1",
    subtitle = "Identification of Secondary Metabolite Clusters",
    x = "BGC Type",
    y = "Number of Clusters"
  ) +
  theme_minimal() +
  scale_fill_brewer(palette = "Dark2") +
  theme(legend.position = "none")




#ecological subset of  bin 1, gene nitrogen

# Extract Bin names correctly from the #query column

# 1. Fix the Bin extraction (Note the closing bracket after `#query`)
egg_annotations_fixed <- egg_annotations %>%
  mutate(Bin = sub("_.*", "", `#query`)) 

# 2. Now create the summary table
mag_summary <- egg_annotations_fixed %>%
  group_by(Bin) %>%
  summarise(
    Total_Genes = n(),
    .groups = "drop"
  )

# View the result
print(mag_summary)




####Ecological Niche (EggNOG) nitrogen cycle 

# Filter for Nitrogen and Phosphorus specific roles
biofloc_functions <- egg_annotations_fixed %>%
  filter(grepl("ammonia|nitrite|nitrate|nitrogen|phosphate", Description, ignore.case = TRUE)) %>%
  mutate(Role = ifelse(grepl("phosph", Description, ignore.case = TRUE), "P-Remediation", "N-Cycling"))

# Plot the ecological potential
ggplot(biofloc_functions, aes(x = Bin, fill = Role)) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = c("P-Remediation" = "#E69F00", "N-Cycling" = "#56B4E9")) +
  labs(title = "Key Ecological Services per MAG",
       subtitle = "Based on EggNOG functional annotation",
       x = "MAG Bin", y = "Gene Count") +
  theme_minimal()

#### oxygen

# Filter for Oxygen-related metabolism and defense
oxygen_story <- egg_annotations_fixed %>%
  filter(grepl("oxygen|aerobic|oxidase|catalase|superoxide|peroxidase", Description, ignore.case = TRUE)) %>%
  mutate(Oxygen_Role = case_when(
    grepl("catalase|superoxide|peroxidase|reductase", Description, ignore.case = TRUE) ~ "Oxidative Stress Defense",
    grepl("oxidase|cytochrome", Description, ignore.case = TRUE) ~ "Aerobic Respiration",
    TRUE ~ "Other Oxygen-related"
  ))

# Plot the distribution across Bins
ggplot(oxygen_story, aes(x = Bin, fill = Oxygen_Role)) +
  geom_bar(position = "stack") +
  scale_fill_brewer(palette = "Paired") +
  labs(title = "Oxygen Metabolism & Defense Profiles",
       subtitle = "Adaptation of MAGs to aerated Biofloc conditions",
       x = "MAG Bin", y = "Number of Genes") +
  theme_minimal()

# Optional: Rename Bins for better labels
oxygen_story_clean <- oxygen_story %>%
  mutate(Bin_Label = case_when(
    Bin == "ECFFIKGN" ~ "Bin 1 (Dominant)", #
    Bin == "GNLOJFBN" ~ "Bin 2",
    Bin == "DMFPEKNM" ~ "Bin 3",
    TRUE ~ Bin
  ))


# This will show you the unique names of all your bins
unique(egg_annotations_fixed$Bin)
#if you look at your first screenshot, your #query names look like ECFFIKGN_00001. The code sub("_.*", "", "#query") tells R: "Take everything before the underscore and call it the Bin". Therefore, for that specific bin, the ID is indeed ECFFIKGN.



###identify top gene in bin 1

# Filter for your dominant bin and count specific gene descriptions
top_genes_bin1 <- egg_annotations_fixed %>%
  filter(Bin == "ECFFIKGN") %>%
  group_by(Description) %>%
  summarise(Count = n()) %>%
  arrange(desc(Count))

# View the top 20 most frequent genes
head(top_genes_bin1, 20)




#visualize


# 1. Prepare data: Get top 15 genes and shorten long descriptions
top_plot_data <- top_genes_bin1 %>%
  slice_max(order_by = Count, n = 15) %>%
  mutate(Description = ifelse(nchar(Description) > 50, 
                              paste0(substr(Description, 1, 47), "..."), 
                              Description))

# 2. Create the plot
ggplot(top_plot_data, aes(x = reorder(Description, Count), y = Count, fill = Count)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # Makes it horizontal for readability
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  theme_minimal() +
  labs(
    title = "Dominant Functional Genes in Bin: ECFFIKGN",
    subtitle = "Top 15 functions based on gene count",
    x = "Gene Description",
    y = "Number of Genes"
  ) +
  theme(
    axis.text.y = element_text(size = 9),
    legend.position = "none"
  )
#the blue is annotted genes




# Filter out the "-" (unannotated) genes before plotting
top_plot_data_clean <- top_genes_bin1 %>%
  filter(Description != "-") %>%  # This removes the weird empty bar
  slice_max(order_by = Count, n = 15) %>%
  mutate(Description = ifelse(nchar(Description) > 50, 
                              paste0(substr(Description, 1, 47), "..."), 
                              Description))

# Re-run your ggplot code with 'top_plot_data_clean'
ggplot(top_plot_data_clean, aes(x = reorder(Description, Count), y = Count, fill = Count)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  theme_minimal() +
  labs(title = "Top Functional Genes in Bin 1: ECFFIKGN (Cleaned)",
       x = "Gene Description", y = "Number of Genes")





===
#identify the most dominant taxa in gene function in bin 1
  
  # Filter for Bin 1 and remove unannotated results
  top_genes_bin1 <- egg_annotations_fixed %>%
  filter(Bin == "ECFFIKGN") %>%
  filter(Description != "-") %>% 
  group_by(Description) %>%
  summarise(Gene_Count = n(), .groups = "drop") %>%
  arrange(desc(Gene_Count))

# View the Top 10 Dominant Genes
head(top_genes_bin1, 10)

# Visualizing the functional dominance
ggplot(head(top_genes_bin1, 15), aes(x = reorder(Description, Gene_Count), y = Gene_Count, fill = Gene_Count)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  theme_minimal() +
  labs(title = "Dominant Functional Genes: Bin 1 (ECFFIKGN)",
       subtitle = "Cleaned functional profile for late-stage biofloc",
       x = "Gene Function", y = "Total Count")


