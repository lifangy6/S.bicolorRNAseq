# Install and load necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")
library(edgeR)
library(stats)

#' Calculate TPM values from raw counts and gene length
#' 
#' https://support.bioconductor.org/p/91218/
#' 
#' @param counts The matrix of raw count data, where rows represent genes and 
#' columns represent samples.
#' @param len A numeric vector of gene lengths corresponding to each gene.
#' @return A matrix of TPM values, with genes in rows and samples in columns.
run_TPM <- function(counts, len) {
  # Divide counts by gene length
  x <- counts / len
  
  # Calculate TPM values
  tpm <- t(t(x) * 1e6 / colSums(x))
  
  return(tpm)
}

# Usage
# Load featureCounts data
fc_result <- readRDS("E:/S.bicolorRNAseq/data/rds/fc_result.rds")

# Calculate TPM counts
tpm_counts <- run_TPM(fc_result$counts, fc_result$annotation$Length)

# Save TPM values
saveRDS(tpm_counts, file = "E:/S.bicolorRNAseq/data/rds/tpm_counts.rds")
