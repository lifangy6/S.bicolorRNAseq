# Install necessary packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
install.packages("ggplot2")
install.packages("plotly")  

library(DESeq2)
library(ggplot2)
library(plotly)


# Load DESeq2 dds result
dds <- readRDS("E:/S.bicolorRNAseq/data/rds/dds.rds")

# Extract count data
countData <- counts(dds)

# Normalize count data
normalizedCounts <- counts(dds, normalized=TRUE)

# Perform PCA
pcaData <- prcomp(t(normalizedCounts))


### 5 Condition PCA ############################################################
# Make metadata for 5 conditions
colData_5 <- data.frame(condition = factor(rep(c("seed_dry", "seed_plasma", "seed_wt", "shoot_plasma", "shoot_wt"), each = 3)))

# Combine PCA data with 5 conditions metadata
pca_df_5 <- data.frame(PC1 = pcaData$x[,1], PC2 = pcaData$x[,2], Group = colData_5)

# Define colors for sample groups (5 conditions)
group_colors_5 <- c("seed_dry" = "blue", "seed_plasma" = "red", "seed_wt" = "green", "shoot_plasma" = "orange", "shoot_wt" = "purple")

# Create PCA plot with ggplot2 for 5 conditions
pca_plot_5 <- ggplot(pca_df_5, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  labs(x = "PC1", y = "PC2", title = "PCA of DESeq2 Normalized Counts (5 Conditions)") +
  scale_color_manual(values = group_colors_5) +
  theme_minimal()

# Plot the PCA plot for 5 conditions
print(pca_plot_5)


### 4 Condition PCA ############################################################
# Make metadata for 4 conditions
colData_4 <- data.frame(condition = factor(rep(c("seed_plasma", "seed_wt", "shoot_plasma", "shoot_wt"), each = 3)))

# Subset PCA data to match the number of rows in metadata for 4 conditions
pca_df_4 <- data.frame(PC1 = pcaData$x[4:15,1], PC2 = pcaData$x[4:15,2], Group = colData_4)

# Define colors for sample groups (4 conditions)
group_colors_4 <- c("seed_plasma" = "red", "seed_wt" = "green", "shoot_plasma" = "orange", "shoot_wt" = "purple")

# Create PCA plot with ggplot2 for 4 conditions
pca_plot_4 <- ggplot(pca_df_4, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  labs(x = "PC1", y = "PC2", title = "PCA of DESeq2 Normalized Counts (4 Conditions)") +
  scale_color_manual(values = group_colors_4) +
  theme_minimal()

# Plot the PCA plot for 4 conditions
print(pca_plot_4)


### 3D PCA #####################################################################
# Make metadata for 5 conditions and 3 principal components
colData_pc3 <- data.frame(condition = factor(rep(c("seed_dry", "seed_plasma", "seed_wt", "shoot_plasma", "shoot_wt"), each = 3)))

# Combine PCA data with 5 conditions and 3 principal components metadata
pca_df_pc3 <- data.frame(PC1 = pcaData$x[,1], PC2 = pcaData$x[,2], PC3 = pcaData$x[,3], Group = colData_pc3)

# Define colors for sample groups
group_colors_pc3 <- c("seed_dry" = "blue", "seed_plasma" = "red", "seed_wt" = "green", "shoot_plasma" = "orange", "shoot_wt" = "purple")

# Create 3D PCA plot with plotly
pca_3d <- plot_ly(data = pca_df_pc3, 
             x = ~PC1, y = ~PC2, z = ~PC3, 
             color = ~condition, colors = group_colors_pc3,
             type = "scatter3d", mode = "markers", marker = list(size = 5)) %>%
  layout(scene = list(xaxis = list(title = "PC1"),
                      yaxis = list(title = "PC2"),
                      zaxis = list(title = "PC3"),
                      aspectmode = "cube"),
         title = "3D PCA of DESeq2 Normalized Counts")

# Plot the 3D PCA plot
pca_3d
