# *Sorghum bicolor* RNA-Seq Analysis

Check repository [Wiki](https://github.com/lifangy6/S.bicolorRNAseq/wiki) for details

## Tree Structure

```         
.
├── bash                            # Bash scripts
|   ├── cufflinks_4c.sh
|   ├── cufflinks_5c.sh
|   ├── hisat2.sh          
├── data                            # Data (BAM, GTF, rds...)
|   ├── cufflinks
|   ├── misc                
|   ├── rds
|   ├── regulated_genes
├── eFP                             # eFP Browser
|   ├── configuration
|   ├── image
|   ├── tpm
├── figure                          # Figures produced
|   ├── CufflinksDendro
|   ├── DeepVenn
|   ├── MAPlot
|   ├── PCA
|   ├── QualityControl
|   ├── VolcanoPlot
├── r                               # R scripts
|   ├── run_DESeq2.R
|   ├── run_DE_regulated_genes.R
|   ├── run_MAPlots.R
|   ├── run_PCA.R
|   ├── run_VolcanoPlots.R
|   ├── run_cummeRbund.R
|   ├── run_featureCounts.R
├── .gitignore
├── LICENSE
├── README.md
└── S.bicolorRNAseq.Rproj
```
