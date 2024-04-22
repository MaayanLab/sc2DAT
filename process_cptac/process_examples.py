#%%
import json
import pandas as pd
from maayanlab_bioinformatics.harmonization import ncbi_genes_lookup

lookup = ncbi_genes_lookup(organism='Mammalia/Homo_sapiens')
#%%
metadata = pd.read_csv('sample_descriptions.tsv', sep='\t', index_col=0)
metadata_HNSCC = metadata[(metadata['cohort'] == 'HNSCC') & (metadata['tissue_type'] == 'Tumor')]
HNSCC_tumors = {i: row['GDC_id'] for i, row in metadata_HNSCC.iterrows()}

# %%
rna_counts_df = pd.read_csv('ALL.rsem_genes_expected_count.txt', sep='\t', index_col=0)
HNSCC_tumor_counts = rna_counts_df[list(HNSCC_tumors.keys())]

HNSCC_tumor_counts.columns = HNSCC_tumor_counts.columns.map(lambda id: HNSCC_tumors[id])
HNSCC_tumor_counts.index = HNSCC_tumor_counts.index.map(lambda g: lookup(g.split('.')[0]) if lookup(g.split('.')[0]) else g.split('.')[0])

# %%
HNSCC_tumor_counts.to_csv('CPTAC3_HNSCC_tumor_counts.tsv', sep='\t')
# %%

phospho_df = pd.read_csv('phospho_isoform_adjusted_v1.1_HNSCC_CB_obs.tsv', sep='\t')
available_samps = list(set(phospho_df.columns.values).intersection(list(map(lambda s: s + '-T', HNSCC_tumors.values()))))
phospho_df_HSNCC = phospho_df[available_samps]
phospho_df_HSNCC.columns = phospho_df_HSNCC.columns.map(lambda s: s.replace('-T', ''))
phospho_df_HSNCC.to_csv('HSNCC_phospho_observed.tsv', sep='\t')

# %%
prot_df = pd.read_csv('HNSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt', sep='\t', index_col=0)
prot_df.index = prot_df.index.map(lambda g: lookup(g.split('.')[0]) if lookup(g.split('.')[0]) else g.split('.')[0])
prot_df.to_csv('HNSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.tsv', sep='\t')
# %%
