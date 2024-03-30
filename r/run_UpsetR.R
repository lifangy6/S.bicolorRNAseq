# Load packages
install.packages("UpSetR")
library(UpSetR)

### All 14 #####################################################################
# Set the directory where gene lists are stored
gene_list_directory <- "E:/S.bicolorRNAseq/data/regulated_genes/genes_edited"

# List the files in the directory
gene_list_files <- list.files(gene_list_directory, full.names = TRUE)

# Initialize an empty list to store the gene sets
gene_sets <- list()

# Read gene sets from files and store them in the list
for (file in gene_list_files) {
  # Extract the condition name from the file name 
  condition_name <- gsub("_genes_edited.txt$", "", basename(file))
  
  # Read the gene IDs from the file
  gene_ids <- readLines(file)
  
  # Store the gene IDs in the list with the condition name as the key
  gene_sets[[condition_name]] <- gene_ids
}

# Create an UpSetR plot
upSetData <- UpSetR::fromList(gene_sets)
UpSetR::upset(upSetData, nsets = 14, order.by = "freq")


### DESeq2 Only ################################################################
DESeq2_gene_sets <- gene_sets[9:14]
UpSetData_DESeq2 <- UpSetR::fromList(DESeq2_gene_sets)
UpSetR::upset(UpSetData_DESeq2, nsets = 6, order.by = "freq")
