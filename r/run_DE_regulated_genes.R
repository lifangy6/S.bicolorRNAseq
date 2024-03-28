# Load DESeq2 results
dds_seed_wt_dry <- readRDS("E:/S.bicolorRNAseq/data/rds/dds_seed_wt_dry.rds")
dds_seed_wt_plasma <- readRDS("E:/S.bicolorRNAseq/data/rds/dds_seed_wt_plasma.rds")
dds_shoot_wt_plasma <- readRDS("E:/S.bicolorRNAseq/data/rds/dds_shoot_wt_plasma.rds")

#' Write Differential Expression (DE) Analysis Results to Files
#'
#' @param res A data frame containing DE analysis results. It should contain columns for 'padj' (adjusted p-values) and 'log2FoldChange' (log2 fold changes).
#' @param padj_threshold The threshold for adjusted p-value (default is 0.05).
#' @param log2FC_threshold The threshold for absolute log2 fold change (default is 1).
#' @param prefix A prefix to be added to the output file names (default is "DESeq2").
#' @return This function does not explicitly return any value, but writes output files.
write_DE_results <- function(res, padj_threshold = 0.05, log2FC_threshold = 1, prefix = "DESeq2") {
  # Subset up-regulated genes
  up_regulated <- subset(res, padj < padj_threshold & log2FoldChange > log2FC_threshold)
  up_regulated_genes <- rownames(up_regulated)
  write.csv(up_regulated, file = paste0(prefix, "_up.csv"))
  write.table(up_regulated_genes, file = paste0(prefix, "_up_genes.txt"), quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)
  
  # Subset down-regulated genes
  down_regulated <- subset(res, padj < padj_threshold & log2FoldChange < -log2FC_threshold)
  down_regulated_genes <- rownames(down_regulated)
  write.csv(down_regulated, file = paste0(prefix, "_down.csv"))
  write.table(down_regulated_genes, file = paste0(prefix, "_down_genes.txt"), quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)
  
  # Edit gene names (remove .v3.2)
  up_regulated_genes_edited <- gsub("\\.v3\\.2", "", up_regulated_genes)
  write.table(up_regulated_genes_edited, file = paste0(prefix, "_up_genes_edited.txt"), 
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  down_regulated_genes_edited <- gsub("\\.v3\\.2", "", down_regulated_genes)
  write.table(down_regulated_genes_edited, file = paste0(prefix, "_down_genes_edited.txt"), 
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

# Usage
write_DE_results(dds_seed_wt_dry, padj_threshold = 0.05, log2FC_threshold = 1, prefix = "DESeq2_seed_wt_dry")
write_DE_results(dds_seed_wt_plasma, padj_threshold = 0.05, log2FC_threshold = 1, prefix = "DESeq2_seed_wt_plasma")
write_DE_results(dds_shoot_wt_plasma, padj_threshold = 0.05, log2FC_threshold = 1, prefix = "DESeq2_shoot_wt_plasma")
