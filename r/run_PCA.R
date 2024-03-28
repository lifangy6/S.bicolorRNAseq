# Install necessary packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
install.packages("ggplot2")
install.packages("plotly")  
install.packages("gplots")

library(DESeq2)
library(ggplot2)
library(plotly)
library(gplots)

# Load DESeq2 dds result
dds <- readRDS("E:/S.bicolorRNAseq/data/rds/dds.rds")

# Extract count data
countData <- counts(dds)

# Normalize count data
normalizedCounts <- counts(dds, normalized=TRUE)

# Perform PCA
pcaData <- prcomp(t(normalizedCounts))

# Make metadata
colData <- data.frame(condition = factor(rep(c("seed_dry",
                                               "seed_plasma",
                                               "seed_wt",
                                               "shoot_plasma", 
                                               "shoot_wt"), each = 3)))

# Combine PCA data with metadata
pca_df <- data.frame(PC1 = pcaData$x[,1], PC2 = pcaData$x[,2], Group = colData)

# Define colors for sample groups
group_colors <- c("seed_dry" = "blue", "seed_plasma" = "red", "seed_wt" = "green", "shoot_plasma" = "orange", "shoot_wt" = "purple")

# Create PCA plot with ggplot2
ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  labs(x = "PC1", y = "PC2", title = "PCA of DESeq2 Normalized Counts") +
  scale_color_manual(values = group_colors) +
  theme_minimal()

