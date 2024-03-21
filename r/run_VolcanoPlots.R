# Install necessary packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)

# Load DESeq2 results
dds_seed_wt_dry <- readRDS("E:/S.bicolorRNAseq/data/rds/dds_seed_wt_dry.rds")
dds_seed_wt_plasma <- readRDS("E:/S.bicolorRNAseq/data/rds/dds_seed_wt_plasma.rds")
dds_shoot_wt_plasma <- readRDS("E:/S.bicolorRNAseq/data/rds/dds_shoot_wt_plasma.rds")

# Cutoff: +1/-1 log2FC & <0.05 adjusted p-value

# seed_wt vs. seed_dry
par(mfrow=c(1,1))

# Draw volcano plot
with(dds_seed_wt_dry, plot(log2FoldChange, -log10(pvalue), pch=20, 
                           main="Volcano Plot - Seed Wt vs. Dry (DESeq2)", 
                           sub="+1/-1 log2FC & <0.05 adjusted p-value", 
                           xlim=c(-10,10), ylim=c(0,30)))

# Add colored points: blue if padj<0.05, red if log2FC>1 and padj<0.05)
with(subset(dds_seed_wt_dry, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(dds_seed_wt_dry, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))


# seed_wt vs. seed_plasma
par(mfrow=c(1,1))

with(dds_seed_wt_plasma, plot(log2FoldChange, -log10(pvalue), pch=20, 
                              main="Volcano Plot - Seed Wt vs. Plasma (DESeq2)", 
                              sub="+1/-1 log2FC & <0.05 adjusted p-value", 
                              xlim=c(-10,10), ylim=c(0,30)))

with(subset(dds_seed_wt_plasma, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(dds_seed_wt_plasma, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))


# shoot_wt vs. shoot_plasma
par(mfrow=c(1,1))

with(dds_shoot_wt_plasma, plot(log2FoldChange, -log10(pvalue), pch=20, 
                               main="Volcano Plot - Shoot Wt vs. Plasma (DESeq2)", 
                               sub="+1/-1 log2FC & <0.05 adjusted p-value", 
                               xlim=c(-10,10), ylim=c(0,30)))

with(subset(dds_shoot_wt_plasma, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(dds_shoot_wt_plasma, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
