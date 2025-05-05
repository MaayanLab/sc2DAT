import anndata as ad
import pandas as pd
from rpy2.robjects import pandas2ri, numpy2ri
import rpy2.robjects as ro
import numpy as np

pandas2ri.activate()
numpy2ri.activate()

adata = ad.read_h5ad("data/76c942bd-45c3-47ec-b290-1e695ec9c177.h5ad")

# Filter for any specific cell types or tissues if needed
#adata = adata[adata.obs['tissue'] == 'liver'].copy()

# Ensure counts are in dense matrix format
if not isinstance(adata.raw.X, np.ndarray):
    counts = adata.raw.X.toarray()
else:
    counts = adata.raw.X

# Push counts and cell_type to R
ro.globalenv['counts'] = numpy2ri.py2rpy(counts.T)  # Transpose to genes x cells
ro.globalenv['cell_types'] = pandas2ri.py2rpy(adata.obs['cell_type'].reset_index(drop=True))

# Send gene names to R
gene_names = adata.var_names.to_list()
ro.globalenv['gene_names'] = ro.StrVector(gene_names)


ro.r('''
library(Matrix)
library(InstaPrism)

# Create a sparse matrix
counts_mat <- Matrix(counts, sparse = TRUE)

# Add gene names as rownames
rownames(counts_mat) <- gene_names

# Prepare the reference object
ref_obj <- refPrepare(
    sc_Expr = counts_mat,
    cell.type.labels = cell_types,
    cell.state.labels = cell_types
)

saveRDS(ref_obj, "guimaraes-et-al.rds")
''')

