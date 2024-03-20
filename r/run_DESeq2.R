# Install necessary packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)

#' Perform differential expression analysis using DESeq2 package.
#'
#' @param raw_counts A matrix of raw counts where rows represent genes and columns represent samples.
#' @param colData A data frame containing sample metadata. Each row corresponds to a sample and columns represent different experimental conditions.
#' @return A DESeqDataSet object containing normalized counts and differential expression results.
run_DESeq2 <- function(raw_counts, colData) {
  # Create a DeSeqDataSet object
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = raw_counts,
                                        colData = colData,
                                        design = ~ condition)
  
  # Pre-filtering: removing rows with low gene counts and keeping rows that have at least 10 reads total
  keep <- rowSums(DESeq2::counts(dds)) >= 10
  dds <- dds[keep,]
  
  # Set the factor level
  dds$condition <- relevel(dds$condition, ref = "seed_wt")
  
  # Run DESeq
  dds <- DESeq2::DESeq(dds)
  return(dds)
}

# Load featureCounts raw counts
raw_counts <- readRDS("E:/S.bicolorRNAseq/data/rds/raw_counts.rds")

# Prepare the metadata
colData <- data.frame(condition = factor(rep(c("seed_dry",
                                               "seed_plasma",
                                               "seed_wt",
                                               "shoot_plasma", 
                                               "shoot_wt"), each = 3)))

rownames(colData) <- c(paste("seed_dry", 1:3, sep = "_"),
                       paste("seed_plasma", 1:3, sep = "_"),
                       paste("seed_wt", 1:3, sep = "_"),
                       paste("shoot_plasma", 1:3, sep = "_"),
                       paste("shoot_wt", 1:3, sep = "_"))

# Check formatting
if (all(colnames(raw_counts) == rownames(colData))) {
  # Process DESeq
  dds <- run_DESeq2(raw_counts, colData)
  
  # Save result
  saveRDS(dds, file = "E:/S.bicolorRNAseq/data/rds/dds.rds")
} else {
  print("Error: column names of raw_counts do not match row names of colData.")
}


### Save DESeq2 Separate Results ###############################################

# Load dds
dds <- readRDS("E:/S.bicolorRNAseq/data/rds/dds.rds")

# Extract results under different comparisons
res_seed_wt_dry <- DESeq2::results(dds, contrast = c("condition", "seed_wt", "seed_dry"))
res_seed_wt_plasma <- DESeq2::results(dds, contrast = c("condition", "seed_wt", "seed_plasma"))
res_shoot_wt_plasma <- DESeq2::results(dds, contrast = c("condition", "shoot_wt", "shoot_plasma"))

# Save RDS
saveRDS(res_seed_wt_dry, file = "E:/S.bicolorRNAseq/data/rds/res_seed_wt_dry.rds")
saveRDS(res_seed_wt_plasma, file = "E:/S.bicolorRNAseq/data/rds/res_seed_wt_plasma.rds")
saveRDS(res_shoot_wt_plasma, file = "E:/S.bicolorRNAseq/data/rds/res_shoot_wt_plasma.rds")
