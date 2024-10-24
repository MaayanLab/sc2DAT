install.packages("R.utils")
install.packages("RCurl")
install.packages("statmod")
install.packages("remotes")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(version = '3.16', ask = FALSE)
BiocManager::install("Biobase")
remotes::install_github("humengying0907/InstaPrism")