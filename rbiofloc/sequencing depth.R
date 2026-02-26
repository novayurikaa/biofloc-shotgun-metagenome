##before doing the analysis we should count the sequencing depth
# if depth differs you should rarefy to the minimum read depth

library (vegan)
library (readxl)
library (dplyr)
library (tidyr)
library (tibble)

#1. calculate sequencing depth with kraken2 data
abundance_table=read.csv('kraken2_summary_species.csv')
dim (abundance_table) #check structure
head(abundance_table)

#all of our data is in one column, we need to make multiple columns
library(tidyr)
library(dplyr)

# Split the single column into multiple columns
abundance_table_clean <- abundance_table %>%
  separate(Sample.Reads.TaxID.Rank.Taxon, 
           into = c("Sample","Reads","TaxID","Rank","Taxon"), 
           sep = "\t", convert = TRUE)

head (abundance_table_clean) #chek

depth <- tapply(abundance_table_clean$Reads, abundance_table_clean$Sample, sum)
depth



#unfortunately it differs quite big
#SRR24442555 SRR24442556 SRR24442557 
9386416    16513700    10602368 




####sequencing depth on bracken file

abundance_table_bracken=read.csv('bracken_species.csv')
dim (abundance_table_bracken) #check structure
head(abundance_table)
# Split the single column into multiple columns
abundance_table_bracken_clean <- abundance_table_bracken %>%
  separate(Sample.Reads.TaxID.Rank.Taxon, 
           into = c("Sample","Reads","TaxID","Rank","Taxon"), 
           sep = "\t", convert = TRUE)


head (abundance_table_clean) #check

depth_bracken <- tapply(abundance_table_clean$Reads, abundance_table_clean$Sample, sum)
depth_bracken #this is sequencing depth per sample

# depth_bracken
#SRR24442555 SRR24442556 SRR24442557 
#9386416    16513700    10602368



#####step 2: rarefraction

min_depth_original<- min(depth_bracken)
min_depth_original
# 9,386,416


#create otu

# Pivot the cleaned table to OTU format: rows = Taxon, columns = Sample
otu_mat <- abundance_table_bracken_clean %>%
  select(Sample, Taxon, Reads) %>%
  pivot_wider(
    names_from = Sample,
    values_from = Reads,
    values_fill = 0  # fill missing taxa with 0 reads
  ) %>%
  column_to_rownames("Taxon")  # taxa as row names



# Check OTU matrix
dim(otu_mat)
head(otu_mat) #value read count numeric,

#prepare otu matrix for rarefraction
otu_mat_rarefy <- t(otu_mat)  # now rows = samples, columns = taxa
dim(otu_mat_rarefy)           # should be 3 × 3436
otu_mat_rarefy #total reads in your otu table per sample



 #rarefraction to minimum sequencing depth


# Find minimum sequencing depth across your samples
min_depth <- min(rowSums(otu_mat_rarefy))
min_depth # 1242141


# Rarefy to minimum depth present in OTU table
min_depth <- min(rowSums(otu_mat_rarefy))  # 1,242,141
otu_rarefied <- rrarefy(otu_mat_rarefy, sample = min_depth)

# Quick check
rowSums(otu_rarefied)  # all rows should equal min_depth
head(otu_rarefied[,1:5])  # first few taxa



#############################################################################
########################### 1.  calculate alpha diversity######################
#we want to see if there are succession.  by ciuntry alpha diversity
#never filter your OTU for alpha diversity.

# Shannon diversity
shannon <- diversity(otu_rarefied, index = "shannon")

# Species richness
richness <- specnumber(otu_rarefied)

# Combine into a table
alpha_div <- data.frame(
  Sample = rownames(otu_rarefied),
  Shannon = shannon,
  Richness = richness
)
alpha_div #gives normalised shannon diversity and richness per sample

####create single figure to show both shanon index and species richness

#prepare the data
# alpha_div contains Shannon and Richness
# Sample   Shannon   Richness
# SRR24442555  x       y

alpha_div #check

alpha_long <- alpha_div %>%
  pivot_longer(
    cols = c("Shannon", "Richness"),   # must match column names exactly
    names_to = "Metric",
    values_to = "Value"
  )

alpha_long

ggplot(alpha_long, aes(x = Sample, y = Value, fill = Sample)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Metric, scales = "free_y") +  # separate y-axis for Shannon vs Richness
  labs(title = "Alpha Diversity Across Samples",
       x = "Sample", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")



########################### 2. count bray curtis#####################
#purpose: we want to see overall community differences, using rarefied counts.
#first analysysi s is using all otu
# Bray-Curtis distance between samples
dist_bc <- vegdist(otu_rel, method = "bray")
dist_bc

# Perform hierarchical clustering bray curtis
hc <- hclust(dist_bc, method = "average")  # UPGMA
plot(hc, main = "Hierarchical Clustering of Samples",
     xlab = "", sub = "", cex = 1.2)
#the result show that sample 57 cluster away from 55 and 56, because 57 has many rare species or different composition
#bray curtis considers all species and their abundance
###############################

library(pheatmap)
# Clear any existing graphics devices
dev.off()  # run once; if it says "null device", that’s fine

############################## Heatmap of relative abundance
pheatmap(otu_rel, 
         cluster_rows = TRUE,   # cluster samples
         cluster_cols = TRUE,   # cluster taxa
         scale = "row",         # optional: scale per sample
         main = "OTU Clustering Heatmap")


#the heatmap cannot be read, we neef to filter otu FOR HEATMAP ONLY

# Filter taxa for visualization, based on relative abundance treeshold, LESS or bigger than 1% relative abundance in at least one sample
#means only taxa that aresomewaht abundant in at least 1 sample are shown. rare taxa ( more than !%) are removed for readabilty.
otu_heatmap <- otu_rel[, apply(otu_rel, 2, max) >= 0.01]
dim(otu_heatmap)  # should be smaller


# Reorder samples according to hc dendrogram
otu_heatmap_ordered <- otu_heatmap[hc$labels[hc$order], ]

#create heatmap of most dominant taxa only from  filtered otu
pheatmap(otu_heatmap_ordered,
         cluster_rows = FALSE,  # rows already ordered by Bray-Curtis
         cluster_cols = TRUE,   # cluster taxa for readability
         scale = "row",         # optional: scale each row for visual contrast
         main = "Bray-Curtis Heatmap (Dominant Taxa)",
         angle_col = 45)




####################### 3. relative abundance##########################
#convert count to relative abundance, for clustering
#1. calculate relative abundance per sample:
otu_rel<- sweep(otu_rarefied, 1, rowSums(otu_rarefied), FUN = "/")
head(otu_rel[,1:5])


# 1. Calculate total relative abundance of each species across all samples
taxa_total <- colSums(otu_rel)

# 2. Select top 20 species globally
top20_species <- names(sort(taxa_total, decreasing = TRUE))[1:20]

# 3. Subset OTU table for only these top 20 species
otu_top20 <- otu_rel[, top20_species]

# 4. Convert to long format for ggplot
otu_long20 <- otu_top20 %>%
  as.data.frame() %>%
  mutate(Sample = rownames(.)) %>%
  pivot_longer(
    cols = -Sample,
    names_to = "Species",
    values_to = "RelAbundance"
  )

# 5. Order species by total abundance for better visualization
otu_long20$Species <- factor(otu_long20$Species, levels = top20_species)

# 6. Plot heatmap with ggplot
ggplot(otu_long20, aes(x = Sample, y = Species, fill = RelAbundance)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "red", labels = percent_format(accuracy = 1)) +
  labs(
    title = "Heatmap of Top 20 Species Across Samples",
    x = "Sample",
    y = "Species",
    fill = "Relative Abundance"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8)
  )





#2. find the most dominant species per sample

dominant_species_sample <- apply(otu_rel, 1, function(x) {
  sorted <- sort(x, decreasing = TRUE)
  data.frame(
    Top1 = names(sorted)[1],
    Abundance1 = sorted[1],
    Top2 = names(sorted)[2],
    Abundance2 = sorted[2],
    Top3 = names(sorted)[3],
    Abundance3 = sorted[3]
  )
})

# Convert list to data frame
dominant_species_sample_df <- do.call(rbind, dominant_species_sample)
dominant_species_sample_df <- cbind(Sample = rownames(otu_rel), dominant_species_sample_df)
dominant_species_sample_df



####2. the most dominant species globally

# Sum across all samples
taxa_global_sum <- colSums(otu_rel)

# Sort descending
taxa_global_sorted <- sort(taxa_global_sum, decreasing = TRUE)

# Top 10 globally
top_global_species <- data.frame(
  Species = names(taxa_global_sorted)[1:10],
  TotalAbundance = taxa_global_sorted[1:10]
)
top_global_species



# Prepare top 10 taxa globally for plotting
top10_taxa <- names(taxa_global_sorted)[1:10]

# Subset relative abundance for top 10 species
otu_top10 <- otu_rel[, top10_taxa]

# Convert to long format
otu_top10_long <- otu_top10 %>%
  as.data.frame() %>%
  mutate(Sample = rownames(.)) %>%
  pivot_longer(
    cols = -Sample,
    names_to = "Species",
    values_to = "RelAbundance"
  )

# Stacked bar plot
library(ggplot2)

ggplot(otu_top10_long, aes(x = Sample, y = RelAbundance, fill = Species)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Top 10 Most Abundant Species Across Samples",
    x = "Sample",
    y = "Relative Abundance (%)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )





########create stacked barplot of relative abundance of  20 top sepcies
#1. prepare the data

# Convert otu_top to long format for ggplot
otu_top_long <- otu_top %>%
  as.data.frame() %>%
  mutate(Sample = rownames(.)) %>%   # keep sample names
  pivot_longer(
    cols = -Sample,
    names_to = "Species",
    values_to = "RelAbundance"
  )

#make ggplot

ggplot(otu_top_long, aes(x = Sample, y = RelAbundance, fill = Species)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Relative Abundance of Top 20 Species Across Samples",
    x = "Sample",
    y = "Relative Abundance (%)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )



#######create top 50 abundance taxa

# Convert rarefied OTU counts to relative abundance
otu_rel <- sweep(otu_rarefied, 1, rowSums(otu_rarefied), FUN = "/")


# Sum abundance across all samples per taxon
taxa_sums <- colSums(otu_rel)

# Get names of top 50 taxa
top50_taxa <- names(sort(taxa_sums, decreasing = TRUE))[1:50]

# Subset OTU table for top 50 taxa
otu_top50 <- otu_rel[, top50_taxa]



# Clear any existing graphics device
if(!is.null(dev.list())) dev.off()

# Heatmap
pheatmap(otu_top50, 
         cluster_rows = TRUE,   # cluster samples
         cluster_cols = TRUE,   # cluster taxa
         scale = "row",         # scale per sample for better visualization
         main = "Top 50 Most Abundant Taxa Across Samples")



