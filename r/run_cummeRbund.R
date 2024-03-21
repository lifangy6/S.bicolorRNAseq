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
