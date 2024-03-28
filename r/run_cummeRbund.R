if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("cummeRbund")
library(cummeRbund)

# Read Cufflinks data
cuff_4c <- cummeRbund::readCufflinks(dir = "E:/S.bicolorRNAseq/data/cufflinks/cufflinks_4c")
cuff_5c <- cummeRbund::readCufflinks(dir = "E:/S.bicolorRNAseq/data/cufflinks/cufflinks_5c")
cuff_tak <- cummeRbund::readCufflinks(dir = "E:/S.bicolorRNAseq/data/cufflinks/cufflinks_tak")
cuff_tak_replicates <- cummeRbund::readCufflinks(dir = "E:/S.bicolorRNAseq/data/cufflinks/cufflinks_tak_replicates")

# Draw dendrograms
par(mar=c(10.2, 4.1, 4.1, 2.1))

dendro_4c <- cummeRbund::csDendro(cummeRbund::genes(cuff_4c))
dendro_4c <- cummeRbund::csDendro(cummeRbund::genes(cuff_4c), replicates = T)

dendro_5c <- cummeRbund::csDendro(cummeRbund::genes(cuff_5c))
dendro_5c <- cummeRbund::csDendro(cummeRbund::genes(cuff_5c), replicates = T)

dendro_tak <- cummeRbund::csDendro(cummeRbund::genes(cuff_tak))
dendro_tak <- cummeRbund::csDendro(cummeRbund::genes(cuff_tak), replicates = T)

dendro_tak_replicates <- cummeRbund::csDendro(cummeRbund::genes(cuff_tak_replicates))

# Extracting gene expression data
diff_4c <- cummeRbund::diffData(cummeRbund::genes(cuff_4c))
diff_5c <- cummeRbund::diffData(cummeRbund::genes(cuff_5c))
diff_tak <- cummeRbund::diffData(cummeRbund::genes(cuff_tak))

# Subset
cuff_seed_wt_dry <- subset(diff_5c, sample_1 == "seed_wt" & sample_2 == "seed_dry")
cuff_seed_wt_plasma <- subset(diff_5c, sample_1 == "seed_wt" & sample_2 == "seed_plasma")
cuff_shoot_wt_plasma <- subset(diff_5c, sample_1 == "shoot_wt" & sample_2 == "shoot_plasma")


#' Write Cufflinks Results to Files
#' 
#' @param data cufflinks dataframe subset
#' @param middle_part the name for the data
#' @param up_threshold threshold for up log2fold value
#' @param down_threshold threshold for down log2fold value
#' @param p_value_threshold threshold for p-value
write_cuff_data <- function(data, middle_part, up_threshold = 1, down_threshold = -1, p_value_threshold = 0.05) {
  cuff_up <- subset(data, p_value < p_value_threshold & log2_fold_change > up_threshold)
  cuff_down <- subset(data, p_value < p_value_threshold & log2_fold_change < down_threshold)
  
  write.csv(cuff_up, file = paste0("cuff_", middle_part, "_up.csv"))
  write.csv(cuff_down, file = paste0("cuff_", middle_part, "_down.csv"))
  
  write.table(cuff_up$gene_id, file = paste0("cuff_", middle_part, "_up_genes.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  write.table(cuff_down$gene_id, file = paste0("cuff_", middle_part, "_down_genes.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  cuff_up_genes_edited <- gsub("\\.v3\\.2", "", cuff_up$gene_id)
  write.table(cuff_up_genes_edited, file = paste0("cuff_", middle_part, "_up_genes_edited.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  cuff_down_genes_edited <- gsub("\\.v3\\.2", "", cuff_down$gene_id)
  write.table(cuff_down_genes_edited, file = paste0("cuff_", middle_part, "_down_genes_edited.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

# Write files
write_cuff_data(cuff_seed_wt_dry, "seed_wt_dry")
write_cuff_data(cuff_seed_wt_plasma, "cuff_seed_wt_plasma")
write_cuff_data(cuff_shoot_wt_plasma, "cuff_shoot_wt_plasma")
