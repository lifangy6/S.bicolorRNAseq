# Install necessary packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)

# Load dds
dds <- readRDS("E:/S.bicolorRNAseq/data/rds/dds.rds")

# Extract results under different comparisons
res_seed_wt_dry <- DESeq2::results(dds, contrast = c("condition", "seed_wt", "seed_dry"))
res_seed_wt_plasma <- DESeq2::results(dds, contrast = c("condition", "seed_wt", "seed_plasma"))
res_shoot_wt_plasma <- DESeq2::results(dds, contrast = c("condition", "shoot_wt", "shoot_plasma"))

# Save
saveRDS(res_seed_wt_dry, file = "E:/S.bicolorRNAseq/data/rds/res_seed_wt_dry.rds")
saveRDS(res_seed_wt_plasma, file = "E:/S.bicolorRNAseq/data/rds/res_seed_wt_plasma.rds")
saveRDS(res_shoot_wt_plasma, file = "E:/S.bicolorRNAseq/data/rds/res_shoot_wt_plasma.rds")