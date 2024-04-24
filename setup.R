install.packages("R.utils")
install.packages("RCurl")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(version = '3.16', ask = FALSE)

BiocManager::install("limma")
BiocManager::install("DeSeq2")
BiocManager::install("edgeR")
install.packages("statmod")

install.packages("tidyverse")
install.packages("ComplexHeatmap")
install.packages("circlize")

# verify
require(limma)

