import os
from IPython.display import HTML
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
from scipy.sparse import csr_matrix
import numpy as np
import seaborn as sns
import requests
import math
import anndata
from scipy.stats import norm
import scipy.stats
from maayanlab_bioinformatics.normalization import log2_normalize
import numpy as np
import json
import pandas as pd
from rpy2 import robjects as ro
from dotenv import load_dotenv
import statsmodels
import statsmodels.stats
import statsmodels.stats.multitest
from scipy.io import mmread
import gzip
from os import fspath
load_dotenv()

def load_seurat_files(mtx_filename, gene_filename, barcodes_filename):
    
    adata = anndata.read_mtx(mtx_filename).T
    with open(barcodes_filename, "r") as f:
        cells = f.readlines()
        cells = [x.strip() for x in cells]
    genes = pd.read_csv(
        gene_filename,
        header=None,
        sep='\t',
    )
    
    adata.var['gene_ids'] = genes.iloc[:, 0].values    
    adata.var['gene_symbols'] = genes.iloc[:, 1].values
    adata.var_names = adata.var['gene_symbols']
    adata.var_names_make_unique(join="-")
    
    
    adata.obs['barcode'] = cells
    adata.obs_names = cells
    adata.obs_names_make_unique(join="-")
    return adata


def load_mtx(mtx_filename, barcodes_filename, gene_filename):
    
    if mtx_filename.endswith(".gz"):
        with gzip.open(mtx_filename, 'rb') as f:
            matrix = mmread(f)
        X = csr_matrix(X)
        adata = anndata.AnnData(X).T
    else:
        adata = anndata.read_mtx(mtx_filename).T
        
    with open(barcodes_filename, "r") as f:
        cells = f.readlines()
        cells = [x.strip() for x in cells]
        
    genes = pd.read_csv(
        gene_filename,
        header=None,
        sep='\t',
    )
    
    adata.var['gene_ids'] = genes.iloc[:, 0].values    
    adata.var['gene_symbols'] = genes.iloc[:, 1].values
    adata.var_names = adata.var['gene_symbols']
    adata.var_names_make_unique(join="-")
    
    
    adata.obs['barcode'] = cells
    adata.obs_names = cells
    adata.obs_names_make_unique(join="-")
    return adata

def read_sc_data(sc_data_file: str, sc_metadata_file: str, type: str):
    if type == 'plain':
        if sc_data_file.endswith('.csv') or sc_data_file.endswith('.csz.gz'):
            sc_data = pd.read_csv(sc_data_file, index_col=0, compression='gzip' if sc_data_file.endswith('.gz') else None)
        elif sc_data_file.endswith('.txt') or sc_data_file.endswith('.txt.gz') or sc_data_file.endswith('.tsv') or sc_data_file.endswith('.tsv.gz'):
            sc_data = pd.read_csv(sc_data_file, index_col=0, sep='\t', compression='gzip' if sc_data_file.endswith('.gz') else None)
        else:
            raise ValueError('File type for scRNA-seq control profile not supported (.csv, .tsv, .txt)')

        if sc_metadata_file.endswith('.csv') or sc_metadata_file.endswith('.csz.gz'):
            sc_metadata = pd.read_csv(sc_metadata_file, index_col=0, compression='gzip' if sc_metadata_file.endswith('.gz') else None)
        elif sc_data_file.endswith('.txt') or sc_metadata_file.endswith('.txt.gz') or sc_data_file.endswith('.tsv') or sc_data_file.endswith('.tsv.gz'):
            sc_metadata= pd.read_csv(sc_metadata_file, index_col=0, sep='\t', compression='gzip' if sc_metadata_file.endswith('.gz') else None)
        else:
            raise ValueError('File type for scRNA-seq control profile not supported (.csv, .tsv, .txt)')

        adata = sc.AnnData(sc_data.T.values, obs=sc_metadata)
        adata.var['gene_names'] = sc_data.index.values
        adata.var.set_index('gene_names', drop=False, inplace=True)
        adata.obs['samples'] = sc_metadata.index.values
        return adata
    
def read_bulk_data(filename: str):
    if filename.endswith('.csv') or filename.endswith('.csz.gz'):
        return pd.read_csv(filename, index_col=0, compression='gzip' if filename.endswith('.gz') else None)
    elif filename.endswith('.txt') or filename.endswith('.txt.gz') or filename.endswith('.tsv') or filename.endswith('.tsv.gz'):
        return pd.read_csv(filename, index_col=0, sep='\t', compression='gzip' if filename.endswith('.gz') else None)
    else:
        raise ValueError('File type for bulk data not supported (.csv, .tsv, .txt)')

def deconvolution_insta_prism(ref_url: str, bulk_expr: pd.DataFrame, output_dir: str, convert_dict: dict):
    """
    Perform deconvolution using InstaPrism and save results.

    Parameters:
    - ref_obj_path: str, path to the reference RDS file.
    - bulk_expr_path: str, path to the bulk expression file.
    - output_dir: str, directory to save output files.
    - convert_dict: dict, mapping of gene symbols to Ensembl IDs.

    Returns:
    - estimated_frac: pd.DataFrame, cell type fractions.
    """
    # Example conversion dictionary (you should populate this with actual mappings)
    os.makedirs(name='temp', exist_ok=True)
    # Step 1: Read the bulk expression data and convert gene symbols to Ensembl IDs
    bulk_expr.index = bulk_expr.index.map(lambda x: convert_dict.get(x, x))  # Convert gene symbols
    print(ref_url)
    response = requests.get(ref_url)
    print(response.status_code)
    if response.status_code == 200:
        # Write the content to a file
        with open('temp/ref_obj.rds', 'wb') as file:
            file.write(response.content)
    else:
        print('Failed to retrieve reference object.', response.status_code)
    
    # Step 2: Write out the modified bulk expression DataFrame
    bulk_expr_ensembl_path = f'temp/bulk_expr_ensembl.tsv'
    bulk_expr = bulk_expr.loc[~bulk_expr.index.duplicated(keep='first')]
    bulk_expr.to_csv(bulk_expr_ensembl_path, sep='\t')

    # Step 3: Create the R code string
    r_code = f"""
    library(InstaPrism)

    # Read the reference object
    print('reading reference object')
    ref_obj <- readRDS('temp/ref_obj.rds')

    # Read the bulk expression data
    print('reading bulk expression')
    bulk_expr <- read.csv('{bulk_expr_ensembl_path}', sep='\\t', row.names=1)

    # Print diagnostic information
    cat("Bulk expression dimensions:", dim(bulk_expr), "\\n")
    cat("Reference object class:", class(ref_obj), "\\n")

    # Perform deconvolution
    deconv_res <- InstaPrism::InstaPrism(bulk_Expr=bulk_expr, refPhi_cs=ref_obj)

    # Get estimated fractions
    estimated_frac <- t(deconv_res@Post.ini.ct@theta)

    # Save estimated fractions
    write.csv(estimated_frac, file='{output_dir}/estimated_frac.csv', row.names=TRUE)

    # Get Z array and save for each cell type
    Z <- get_Z_array(deconv_res)
    cell_types <- dimnames(Z)[[3]]  # Get cell type names

    for (cell_type in cell_types) {{
        DCT_Z <- Z[, , cell_type]
        write.csv(DCT_Z, file=paste0('{output_dir}/', gsub("/", "", cell_type), '_Z.csv'), row.names=TRUE)
    }}

    # Return estimated fractions
    """

    # Step 4: Execute the R code
    try:
        ro.r(r_code)  # Execute the R code and get estimated fractions
    except Exception as e:
        print(f"Error during execution: {str(e)}")
        raise
    
    estimated_frac_df = pd.read_csv(f'{output_dir}/estimated_frac.csv', index_col=0)
    cell_type_dfs = {}
    for f in os.listdir(output_dir):
        if f.endswith('_Z.csv'):
            cell_type_dfs[f[:-4]] = pd.read_csv(f'{output_dir}/{f}', index_col=0)
    return estimated_frac_df, cell_type_dfs

import scanpy as sc
import decoupler as dc
def normalize_sc_data(adata, n_neighbors, min_dist, resolution):
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    # Annotate the group of mitochondrial genes as 'mt'
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # Filter cells following standard QC criteria.
    adata = adata[adata.obs.n_genes_by_counts < 2500, :]
    adata = adata[adata.obs.pct_counts_mt < 5, :]

    # Normalize the data
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers['log_norm'] = adata.X.copy()
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')
    dc.swap_layer(adata, 'log_norm', X_layer_key=None, inplace=True)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=40)
    sc.tl.umap(adata, min_dist=min_dist)
    sc.tl.leiden(adata, resolution=resolution)
    return adata


def mannwhitney3(n1, n2, ranksum):
	meanRankExpected = (n1*n2)/2
	sigma = math.sqrt(n1)*math.sqrt(n2/12)*math.sqrt(n1+n2+1)
	U = ranksum - n1*(n1+1)/2
	z = (U - meanRankExpected)/sigma
	pvals = np.minimum(norm.cdf(z), norm.sf(z))*2
	return z, pvals



metadata_api = "https://maayanlab.cloud/sigcom-lincs/metadata-api"
data_api = "https://maayanlab.cloud/sigcom-lincs/data-api/api/v1"


input_gene_set = {
    "up_entities": ["TARBP1", "APP", "RAP1GAP", "UFM1", "DNAJA3", "PCBD1", "CSRP1"],
    "down_entities": ["CEBPA", "STAT5B", "DSE", "EIF4EBP1", "CARD8", "HLA-DMA", "SERPINE1"]
}

meta_columns = ["pert_name", "pert_dose", "pert_time"]
set_columns = ['zscore', 'p-value', 'type']
up_down_columns = ['z-up', 'p-up', 'z-down', 'p-down', 'z-sum', 'type']

def get_genes_metadata():
	print("Getting genes metadata")
	if not os.path.isfile("data/genes_meta.json"):
		if not os.path.exists('data/'):
			os.makedirs('data/')
		print("Fetching signature metadata...")
		url = "https://s3.dev.maayanlab.cloud/sigcom-lincs/ranker/genes_meta.json"
		res = requests.get(url)
		with open("data/genes_meta.json", "w") as o:
			o.write(json.dumps(res.json()))
		return res.json()
	else:
		with open("data/genes_meta.json") as o:
			return json.loads(o.read())


def get_metadata():
	if not os.path.isfile("data/signatures_meta.json"):
		if not os.path.exists('data/'):
			os.makedirs('data/')
		print("Fetching signature metadata...")
		url = "https://s3.dev.maayanlab.cloud/sigcom-lincs/ranker/signatures_meta.json"
		res = requests.get(url)
		with open("data/signatures_meta.json", "w") as o:
			o.write(json.dumps(res.json()))
		return res.json()
	else:
		with open("data/signatures_meta.json") as o:
			return json.loads(o.read())


metadata_api = "https://maayanlab.cloud/sigcom-lincs/metadata-api"
data_api = "https://maayanlab.cloud/sigcom-lincs/data-api/api/v1"

datasets = {
	'l1000_shRNA': '8f1ff550-ece8-591d-a213-2763f854c008',
	'l1000_oe': 'ef9389a8-53d3-50db-90cc-57e7d150b76c',
	'l1000_cp': '54198d6e-fe17-5ef8-91ac-02b425761653',
	'l1000_xpr': '96c7b8c5-1eca-5764-88e4-e4ccaee6603f',
	'l1000_mean_xpr': '98c6a1ef-e3b7-5a96-9913-480acf78a577',
	'l1000_lig': 'e24989ba-5258-511c-9c6a-00634f74857b',
	'l1000_mean_cp': '42cd56da-0ad8-5dad-b27c-fe1d135401b2',
	'l1000_aby': 'b759a653-0877-549a-a76a-3a8914b73106',
	'l1000_siRNA': 'b953025a-4356-5cc8-b6e3-dcf2f4f85420'
}


def convert_entities(genes, genes_meta=None):
	metadata_api = "https://maayanlab.cloud/sigcom-lincs/metadata-api"
	if not genes_meta:
		payload = {
			"filter": {
				"where": {
					"meta.symbol": {"inq": [g.upper() for g in genes]}
				},
				"fields": ["id"]
			}
		}
		res = requests.post(metadata_api + "/entities/find", json=payload)
		if res.ok:
			return [i["id"] for i in res.json()]
		else:
			print(res.text)
	else:
		return [genes_meta[i] for i in genes if i in genes_meta]

def get_sigcom_link(input, convert=False):
	if convert:
		if 'entities' in input:
			input['entities'] = convert_entities(input['entities'])
		if 'up_entities' in input:
			input['up_entities'] = convert_entities(input['up_entities'])
		if 'down_entities' in input and 'up_entities' in input:
			input['down_entities'] = convert_entities(input['down_entities'])
	payload = {
		"meta": {
			**input,
			"$validator": "/dcic/signature-commons-schema/v6/meta/user_input/user_input.json"
		},
		"type": "signatures"
	}
	res = requests.post(metadata_api + "/user_input", json=payload)
	if res.ok:
		uid = res.json()["id"]
		endpoint = "Set" if "entities" in input else "UpDown"
		link = "https://maayanlab.cloud/sigcom-lincs/#/SignatureSearch/%s/%s"%(endpoint, uid)
		return input.get('description', ''), link
	else: 
		print(res.text)
	# print("%s: %s"%(input.get('description', ''), link))


def get_sigcom_data(input: dict, database: str, enrichment_id: str, convert=False):
	DATA_API = "https://maayanlab.cloud/sigcom-lincs/data-api/api/v1/"
	METADATA_API = "https://maayanlab.cloud/sigcom-lincs/metadata-api/"

	if convert:
		if 'entities' in input:
			input['entities'] = convert_entities(input['entities'])
		if 'up_entities' in input:
			input['up_entities'] = convert_entities(input['up_entities'])
		if 'down_entities' in input and 'up_entities' in input:
			input['down_entities'] = convert_entities(input['down_entities'])

	query = {
		**input,
		"limit": 10,
		"database": database,
		"enrichment_id": enrichment_id
	}

	res = requests.post(DATA_API + "enrich/ranktwosided", json=query)
	results = res.json()

	# Optional, multiply z-down and direction-down with -1
	for i in results["results"]:
		i["z-down"] = -i["z-down"]
		i["direction-down"] = -i["direction-down"]

	sigids = {i["uuid"]: i for i in results["results"]}

	payload = {
		"filter": {
			"where": {
				"id": {
					"inq": list(sigids.keys())
				}
			}
		}
	}

	res = requests.post(METADATA_API + "signatures/find", json=payload)
	signatures = res.json()

	## Merge the scores and the metadata
	for sig in signatures:
		uid = sig["id"]
		scores = sigids[uid]
		scores.pop("uuid")
		sig["scores"] = scores
	bar_chat_data = []
	for sig in signatures:
		if sig["scores"]["type"] == "reversers":
			name = f"{sig['meta']['pert_name']} {sig['meta']['pert_type']}"
			bar_chat_data.append([name, abs(sig["scores"]["z-sum"])])

	sigs_df = pd.DataFrame(bar_chat_data, columns=["name", "abs(Z-sum)"]).sort_values("abs(Z-sum)", ascending=False)[:10]
	return sigs_df


def ttest_differential_expression(controls_mat: pd.DataFrame, cases_mat: pd.DataFrame, equal_var=False, alternative='two-sided', log2norm=True):
  ''' Given two separate dataframes (controls, cases) with a shared index (genes),
  we compute the ttest differential expression for all genes. Benjamini-Hochberg Adjusted p-value.

  :param controls_mat: (pd.DataFrame) the control samples (samples as columns and genes as rows)
  :param cases_mat: (pd.DataFrame) the case samples (samples as columns and genes as rows)
  :param equal_var: (bool) Should t-test assume equal variance (default: False)
  :param alternative: (str) Alternative hypothesis (see scipy.stats.ttest_ind) (default: two-sided)
  :param log2norm: (bool) Apply log2norm, typically keep with raw counts but disable if you have normalized data (default: True)
  :return: A data frame with the results
  '''
  assert (controls_mat.index == cases_mat.index).all(), 'Index between controls and cases must be the same'
  if log2norm:
    cases_mat = log2_normalize(cases_mat)
    controls_mat = log2_normalize(controls_mat)
  results = scipy.stats.ttest_ind(cases_mat.T, controls_mat.T, equal_var=equal_var, alternative=alternative)
  df_results = pd.DataFrame({
    'Statistic': results.statistic,
    'Pval': results.pvalue,
  }, index=controls_mat.index)
  adj_pval = statsmodels.stats.multitest.multipletests(df_results['Pval'].fillna(1.), method='fdr_bh')
  df_results['AdjPval'] = adj_pval[1]
  df_results.sort_values('AdjPval', inplace=True)
  return df_results
		







