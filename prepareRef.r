# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# if (!require("Biobase", quietly = TRUE))
#   BiocManager::install("Biobase")
# 
# if (!require("devtools", quietly = TRUE))
#   install.packages("devtools")
# 
# devtools::install_github("humengying0907/InstaPrism")



library(InstaPrism)

# load reference data (cellxgene rds format)
library(Seurat)
ref = readRDS('example_cellxgene.rds')
meta_data = ref@meta.data

# Save the reference object which can be added to the sc2DAT referece scRNA-seq options
ref_obj = refPrepare(sc_Expr = ref[["RNA"]]@counts, cell.type.labels = meta_data[["cell_type"]], cell.state.labels = meta_data[["cell_type"]])
saveRDS(ref_obj, "example_ref.rdf")