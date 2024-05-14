#%%
import json
import pandas as pd
from maayanlab_bioinformatics.harmonization import ncbi_genes_lookup

lookup = ncbi_genes_lookup(organism='Mammalia/Homo_sapiens')

#%%
for ct in ['BRCA', 'CCRCC', 'Colon', 'GBM', 'LSCC', 'LUAD', 'Ovary', 'PDAC', 'UCEC']:
    metadata = pd.read_csv('data/sample_descriptions.tsv', sep='\t', index_col=0)
    metadata_ct = metadata[(metadata['cohort'] == ct) & (metadata['tissue_type'] == 'Tumor')]
    ct_tumors = {i: row['GDC_id'] for i, row in metadata_ct.iterrows()}

    rna_counts_df = pd.read_csv('data/ALL.rsem_genes_expected_count.txt', sep='\t', index_col=0)
    ct_tumor_counts = rna_counts_df[list(ct_tumors.keys())]

    ct_tumor_counts.columns = ct_tumor_counts.columns.map(lambda id: ct_tumors[id])
    ct_tumor_counts.index = ct_tumor_counts.index.map(lambda g: lookup(g.split('.')[0]) if lookup(g.split('.')[0]) else g.split('.')[0])

    ct_tumor_counts.to_csv(f'out/CPTAC3_{ct}_tumor_counts.tsv', sep='\t')

#%%
#Phosphoproteome_single_site_isoform_adjusted_SinaiPreprocessed
import os
for file in os.listdir('data/observed_data_phospho'):
    phospho_df = pd.read_csv(f'data/observed_data_phospho/{file}', sep='\t')
    phospho_df.columns = phospho_df.columns.map(lambda s: s.replace('-T', '').replace('.T', ''))
    cols_to_drop = [c for c in phospho_df.columns if '-N' in c]
    phospho_df.drop(columns=cols_to_drop).to_csv(f'out/{file.split("_")[4]}_phospho_observed.tsv', sep='\t')
# %%
for file in os.listdir('data/observed_data_proteomics'):
    if file.split('_')[-1].split('.')[0] == 'Tumor':
        ct = file.split('_')[0]
        prot_df = pd.read_csv(f'data/observed_data_proteomics/{file}', sep='\t', index_col=0)
        prot_df.index = prot_df.index.map(lambda g: lookup(g.split('.')[0]) if lookup(g.split('.')[0]) else g.split('.')[0])
        prot_df.to_csv(f'out/{ct}_proteomics_tumor.tsv', sep='\t')
# %%
