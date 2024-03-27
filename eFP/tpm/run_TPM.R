# Install and load necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")
install.packages("reshape2")
library(edgeR)
library(stats)
library(reshape2)
library(corrplot)
library(ggplot2)
library(hrbrthemes)


#' Calculate TPM values from raw counts and gene length
#' 
#' https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/
#' 
#' @param counts The matrix of raw count data, where rows represent genes and 
#' columns represent samples.
#' @param len A numeric vector of gene lengths corresponding to each gene.
#' @return A matrix of TPM values, with genes in rows and samples in columns.
run_TPM <- function(counts, len) {
  # Divide counts by gene length
  x <- counts / (len / 1000)
  
  # Calculate scale factor
  scale_factor <- colSums(counts) / 1e6
  
  # Calculate TPM values
  tpm <- x / scale_factor
  
  return(tpm)
}


# Load rds and get counts & len
fc_result <- readRDS("E:/S.bicolorRNAseq/data/rds/fc_result.rds")
raw_counts <- readRDS("E:/S.bicolorRNAseq/data/rds/raw_counts.rds")
gene_length <- fc_result$annotation$Length

# Calculate TPM counts
tpm_counts <- run_TPM(raw_counts, gene_length)

# Save TPM values
saveRDS(tpm_counts, file = "E:/S.bicolorRNAseq/data/rds/tpm_counts.rds")


### Correlation Analysis #######################################################

#' Process and save databased TPM
#' 
#' @param mean_tpm_file Path to the file containing mean TPM values.
#' @param all_tpm_file Path to the file containing all TPM values.
#' @param output_file Path to save the databased TPM.
#' @return NULL
process_databased_tpm <- function(mean_tpm_file, all_tpm_file, output_file) {
  # Process Mean TPM
  mean_tpm <- read.table(mean_tpm_file, header = FALSE, sep = "|", strip.white = TRUE)[, c(2, 3)]
  colnames(mean_tpm) <- c("gene_id", "mean_tpm")
  mean_tpm$gene_id <- trimws(mean_tpm$gene_id)
  
  # Process ALL TPM VALUES
  all_tpm <- read.table(all_tpm_file, header = FALSE, sep = "|", strip.white = TRUE)[, c(4, 5, 6)]
  colnames(all_tpm) <- c("gene_id", "tpm", "sample")
  
  # Reshape and rename columns
  all_tpm_reshape <- reshape2::dcast(all_tpm, gene_id ~ sample, value.var = "tpm")
  colnames(all_tpm_reshape)[2:4] <- c("seed_dry_1", "seed_dry_2", "seed_dry_3")
  
  # Combine data
  databased_tpm <- all_tpm_reshape[, -1]
  rownames(databased_tpm) <- all_tpm_reshape[, 1]
  databased_tpm$mean_tpm <- mean_tpm$mean_tpm
  
  # Save data
  saveRDS(databased_tpm, file = output_file)
}


# Get databased TPM
process_databased_tpm("E:/S.bicolorRNAseq/eFP/tpm/databased_dry_seed_mean_tpm.txt",
                      "E:/S.bicolorRNAseq/eFP/tpm/databased_dry_seed_tpm.txt",
                      "E:/S.bicolorRNAseq/eFP/tpm/databased_tpm.rds")

# Load data
fc_result <- readRDS("E:/S.bicolorRNAseq/data/rds/fc_result.rds")
raw_counts <- readRDS("E:/S.bicolorRNAseq/data/rds/raw_counts.rds")
databased_tpm <- readRDS("E:/S.bicolorRNAseq/eFP/tpm/databased_tpm.rds")

# Preprocess raw_counts
raw_counts_edit <- raw_counts
rownames(raw_counts_edit) <- sub("\\.v3\\.2$", "", rownames(raw_counts_edit))


# APPROACH 1: Calculate then Subset
set1_tpm <- run_TPM(raw_counts_edit, fc_result$annotation$Length)
set1_tpm <- as.data.frame(set1_tpm)
set1_tpm <- set1_tpm[rownames(databased_tpm), , drop = FALSE][, 1:3]
set1_tpm$raw_mean_tpm <- rowMeans(set1_tpm)

cm1 <- cor(databased_tpm, set1_tpm)
print(cm1)

# Approach 1's correlation coefficients are generally weak (<0.5). 
# The approach must be corrected since the number of genes in the gene list affects the scale factor. 
# The subsection should be done first, followed by the TPM calculation.


# APPROACH 2: Subset then calculate
set2_tpm <- raw_counts_edit[, 1:3]
set2_tpm <- as.data.frame(set2_tpm)
set2_gene_length <- fc_result$annotation[, c(1, 6)]
rownames(set2_gene_length) <- sub("\\.v3\\.2$", "", set2_gene_length[, 1])
set2_gene_length <- set2_gene_length[rownames(databased_tpm), , drop = FALSE]
set2_tpm <- run_TPM(set2_tpm[rownames(databased_tpm), , drop = FALSE], set2_gene_length$Length)
set2_tpm <- cbind(set2_tpm, raw_mean_tpm = rowMeans(set2_tpm))

cm2 <- cor(databased_tpm, set2_tpm)
print(cm2)

# Approach 2's correlation coefficients slightly increased but are still weak (~0.5). 
# There's an obvious outlier that should be removed.


# APPROACH 3: Remove outlier then subset then calculation
set3 <- raw_counts_edit[, 1:3]
set3 <- as.data.frame(set3)
set3_gene_length <- fc_result$annotation[, c(1, 6)]
rownames(set3_gene_length) <- sub("\\.v3\\.2$", "", set3_gene_length[, 1])
set3_gene_length <- set3_gene_length[rownames(databased_tpm), , drop = FALSE]

# Remove an obvious outlier
outlier_code <- "Sobic.K044407"
set3_databased_tpm <- databased_tpm[rownames(databased_tpm) != outlier_code, ]
set3_gene_length <- set3_gene_length[rownames(set3_gene_length) != outlier_code, ]
set3 <- set3[rownames(set3) != outlier_code, ]

# Run TPM calculation
set3_tpm <- run_TPM(set3[rownames(set3_databased_tpm), , drop = FALSE], set3_gene_length$Length)
set3_tpm <- cbind(set3_tpm, raw_mean_tpm = rowMeans(set3_tpm))

cm3 <- cor(set3_databased_tpm, set3_tpm)
print(cm3)

# Approach 3's correlation coefficients are all above 0.8.

corrplot(cm3, method = "color", addCoef.col = "black", tl.col = "black", tl.srt = 45)



### Scatter Plot ###############################################################

df1 <- set2_tpm
df2 <- databased_tpm

# Convert row names to a column for plotting
df1$GeneID <- rownames(df1)
df2$GeneID <- rownames(df2)

# Merge data frames to compare TPM values
merged_df <- merge(df1, df2, by = "GeneID", suffixes = c("_raw", "_databased"))
colnames(merged_df)[which(names(merged_df) == "mean_tpm")] <- "databased_mean_tpm"

p1 <- ggplot(merged_df, aes(x=seed_dry_1_raw, y=seed_dry_1_databased)) +
  ggtitle("Raw vs. Databased TPM - seed_dry_1") +
  geom_point() +
  theme_ipsum()

p2 <- ggplot(merged_df, aes(x=seed_dry_2_raw, y=seed_dry_2_databased)) +
  ggtitle("Raw vs. Databased TPM - seed_dry_2") +
  geom_point() +
  theme_ipsum()

p3 <- ggplot(merged_df, aes(x=seed_dry_3_raw, y=seed_dry_3_databased)) +
  ggtitle("Raw vs. Databased TPM - seed_dry_3") +
  geom_point() +
  theme_ipsum()

p4 <- ggplot(merged_df, aes(x=raw_mean_tpm, y=databased_mean_tpm)) +
  ggtitle("Raw vs. Databased TPM - Mean TPM") +
  geom_point() +
  theme_ipsum()

# View the plots
p1
p2
p3
p4
