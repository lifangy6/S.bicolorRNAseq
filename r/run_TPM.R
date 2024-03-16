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
  x <- counts/ (len / 1000)
  
  # Calculate scale factor
  scale_factor <- colSums(counts) / 1e6
  
  # Calculate TPM values
  tpm <- x / scale_factor
  
  return(tpm)
}

# Usage
# Load rds and get counts & len
fc_result <- readRDS("E:/S.bicolorRNAseq/data/rds/fc_result.rds")
raw_counts <- readRDS("E:/S.bicolorRNAseq/data/rds/raw_counts.rds")
dds <- readRDS("E:/S.bicolorRNAseq/data/rds/dds.rds")
norm_counts <- DESeq2::counts(dds, normalized=TRUE) 
gene_length <- fc_result$annotation$Length

# Calculate TPM counts
tpm_counts <- run_TPM(raw_counts, gene_length)

# Save TPM values
saveRDS(tpm_counts, file = "E:/S.bicolorRNAseq/data/rds/tpm_counts.rds")

################################################################################

x <- raw_counts / (fc_result$annotation$Length / 1000)

y <- colSums(raw_counts)

y <- y / 1e6

x <- x / y 


################################################################################

norm_counts <- DESeq2::counts(dds, normalized=TRUE)

z <- norm_counts / (fc_result$annotation$Length / 1000)

sf <- colSums(norm_counts)

sf <- sf / 1e6

z <- z / sf
