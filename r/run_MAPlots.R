# Install necessary packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)

# Load DESeq2 results
res_seed_wt_dry <- readRDS("E:/S.bicolorRNAseq/data/rds/res_seed_wt_dry.rds")
res_seed_wt_plasma <- readRDS("E:/S.bicolorRNAseq/data/rds/res_seed_wt_plasma.rds")
res_shoot_wt_plasma <- readRDS("E:/S.bicolorRNAseq/data/rds/res_shoot_wt_plasma.rds")

# MA plots
DESeq2::plotMA(res_seed_wt_plasma, alpha = 0.1, main = "MA Plot: seed_wt vs. seed_dry (alpha = 0.1)")
DESeq2::plotMA(res_seed_wt_dry, alpha = 0.1, main = "MA Plot: seed_wt vs. seed_plasma (alpha = 0.1)")
DESeq2::plotMA(res_shoot_wt_plasma, alpha = 0.1, main = "MA Plot: shoot_wt vs. shoot_plasma (alpha = 0.1)")
