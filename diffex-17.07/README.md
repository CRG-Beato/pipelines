# README

Pipeline to perform differential expression analysis.



## Requirements

The following R packages need to be installed:
```R
source("http://bioconductor.org/biocLite.R")
biocLite("rhdf5")
install.packages("devtools")
devtools::install_github("pachterlab/sleuth")
install.packages("RSQLite")
biocLite("biomaRt")
install.packages("DESeq2")
install.packages("pheatmap")
```