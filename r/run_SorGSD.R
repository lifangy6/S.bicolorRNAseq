#' Use SorGSD Conversion CSV to Convert Gene IDs
#' 
#' @param conversion_sheet_path Path to the conversion sheet in CSV format
#' @param gene_path Path to the regulated gene lists in txt format
#' @param filename Filename to output
#' @return NULL
write_converted_genes <- function(conversion_sheet_path, gene_path, filename) {
  # Read the conversion sheet CSV file
  conversion_sheet <- read.csv(conversion_sheet_path)
  # Extract relevant columns from the conversion sheet
  conversion_sheet <- conversion_sheet[ , 3:6]
  
  # Read the gene IDs text file
  gene_ids <- readLines(gene_path)
  
  # Merge gene IDs with conversion sheet based on current_version and Version2.1
  merged_df <- merge(data.frame(current_version = gene_ids), conversion_sheet, by.x = "current_version", by.y = "Version2.1", all.x = TRUE)
  
  # Filter out rows with NA or "null" UniProt values
  filtered_df <- merged_df[!(is.na(merged_df$UniProt) | merged_df$UniProt == "null"), ]
  
  # Extract converted gene IDs
  converted_gene_ids <- filtered_df$UniProt
  
  # Write converted gene IDs to output file
  write.table(converted_gene_ids, file = filename, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
}

# Usage
write_converted_genes(conversion_sheet_path = "E:/S.bicolorRNAseq/data/enrichment/gProfiler/sorgsd_seed_wt_plasma_up.csv",
                      gene_path = "E:/S.bicolorRNAseq/data/regulated_genes/genes_edited/DESeq2_seed_wt_plasma_up_genes_edited.txt",
                      filename = "E:/S.bicolorRNAseq/data/enrichment/gProfiler/DESeq2_seed_wt_plasma_up_genes_converted.txt")

write_converted_genes(conversion_sheet_path = "E:/S.bicolorRNAseq/data/enrichment/gProfiler/sorgsd_seed_wt_plasma_down.csv",
                      gene_path = "E:/S.bicolorRNAseq/data/regulated_genes/genes_edited/DESeq2_seed_wt_plasma_down_genes_edited.txt",
                      filename = "E:/S.bicolorRNAseq/data/enrichment/gProfiler/DESeq2_seed_wt_plasma_down_genes_converted.txt")
