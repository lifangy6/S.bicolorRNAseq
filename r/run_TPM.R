# Install and load necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")
install.packages("reshape2")
library(edgeR)
library(stats)
library(reshape2)

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
process_databased_tpm("E:/S.bicolorRNAseq/data/tpm/dry_seed_all_genes_mean_tpm.txt",
                      "E:/S.bicolorRNAseq/data/tpm/dry_seed_all_genes.txt",
                      "E:/S.bicolorRNAseq/data/tpm/databased_tpm.rds")

# Load data
fc_result <- readRDS("E:/S.bicolorRNAseq/data/rds/fc_result.rds")
raw_counts <- readRDS("E:/S.bicolorRNAseq/data/rds/raw_counts.rds")
databased_tpm <- readRDS("E:/S.bicolorRNAseq/data/tpm/databased_tpm.rds")

# Preprocess raw_counts
raw_counts_edit <- raw_counts
rownames(raw_counts_edit) <- sub("\\.v3\\.2$", "", rownames(raw_counts_edit))

# Calculate then Subset
cal_then_sub <- run_TPM(raw_counts_edit, fc_result$annotation$Length)
cal_then_sub <- as.data.frame(cal_then_sub)
cal_then_sub <- cal_then_sub[rownames(databased_tpm), , drop = FALSE][, 1:3]
cal_then_sub$raw_mean_tpm <- rowMeans(cal_then_sub)

cm1 <- cor(databased_tpm, cal_then_sub)

# Subset then calculate
sub_then_cal <- raw_counts_edit[, 1:3]
sub_then_cal <- as.data.frame(sub_then_cal)
gene_length_subset <- fc_result$annotation[, c(1, 6)]
rownames(gene_length_subset) <- sub("\\.v3\\.2$", "", gene_length_subset[, 1])
gene_length_subset <- gene_length_subset[rownames(databased_tpm), , drop = FALSE]
sub_then_cal <- run_TPM(sub_then_cal[rownames(databased_tpm), , drop = FALSE], gene_length_subset$Length)
sub_then_cal <- cbind(sub_then_cal, raw_mean_tpm = rowMeans(sub_then_cal))

cm2 <- cor(databased_tpm, sub_then_cal)



