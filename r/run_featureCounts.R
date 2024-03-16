# Install and load necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Rsubread")
library(Rsubread)

#' Run featureCounts on BAM files and save the results
#' 
#' @param bam_path Character string specifying the path to BAM files.
#' @param gtf Character string specifying the path to the GTF file.
#' @return A list containing featureCounts result and raw counts.
run_featureCounts <- function(bam_path, gtf) {
  # Perform featureCounts
  fc_result <- Rsubread::featureCounts(files = list.files(path = bam_path, pattern = "\\.bam$", full.names = TRUE),
                                       annot.ext = gtf,
                                       isGTFAnnotationFile = TRUE,
                                       GTF.featureType = "exon",
                                       GTF.attrType = "gene_id",
                                       useMetaFeatures = TRUE,
                                       isPairedEnd = TRUE)
  
  # Extract and rename raw counts
  raw_counts <- fc_result$counts
  colnames(raw_counts) <- gsub("\\.bam", "", colnames(raw_counts))
  
  return(list(fc_result = fc_result, raw_counts = raw_counts))
}

# Usage
bam_path <- "E:/S.bicolorRNAseq/data/misc/bam_files"
gtf <- "E:/S.bicolorRNAseq/data/misc/Sbicolor_annotations.gtf"

result <- run_featureCounts(bam_path, gtf)

# Save featureCounts result and raw counts
saveRDS(result$fc_result, file = "E:/S.bicolorRNAseq/data/rds/fc_result.rds")
saveRDS(result$raw_counts, file = "E:/S.bicolorRNAseq/data/rds/raw_counts.rds")
