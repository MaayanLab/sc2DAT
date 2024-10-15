install.packages("R.utils")
install.packages("RCurl")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(version = '3.16', ask = FALSE)

BiocManager::install("limma")
BiocManager::install("DESeq2")
BiocManager::install("edgeR")
BiocManager::install("Biobase")
install.packages("statmod")
install.packages("devtools")
devtools::install_github("humengying0907/InstaPrism")
# verify
require(limma)
require(InstaPrism)



