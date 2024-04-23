import pandas as pd
from IPython.display import HTML

import requests
import json
import rpy2
import rpy2.rinterface

def read_table(filename):
    if filename.endswith('.tsv') or filename.endswith('.tsv.gz') or filename.endswith('.txt') or filename.endswith('.txt.gz'):
        return pd.read_csv(filename, sep='\t', index_col=0)
    elif filename.endswith('.csv') or filename.endswith('.csv.gz'):
        return pd.read_csv(filename, sep=',', index_col=0)
    elif filename.endswith('.gct') or filename.endswith('.gct.gz'):
        return pd.read_csv(filename, sep='\t', index_col=0, skiprows=2)
    else:
        return pd.read_table(filename, sep=None, engine='python', index_col=0)

def create_download_link(df, title = "Download CSV file: {}", filename = "data.csv"):  
    df.to_csv(filename)
    html = "<a href=\"./{}\" target='_blank'>{}</a>".format(filename, title.format(filename))
    return HTML(html)

def get_chea3_results(gene_set, query_name='', db='Integrated--meanRank'):
    """ Inputs:
    gene_set: list[]
    query_name: str
    db: ['Integrated--meanRank', 'Integrated--topRank'], ranking method for ChEA3 
    Return: list[] of TFs by decreasing score (higher match) """ 

    ADDLIST_URL = 'https://maayanlab.cloud/chea3/api/enrich/'
    payload = {
        'gene_set': gene_set,
        'query_name': query_name
    }
    
    # Retrieve response
    response = requests.post(ADDLIST_URL, data=json.dumps(payload))
    if not response.ok:
        raise Exception('Error analyzing gene list')
    chea = json.loads(response.text)

    # Save nodes
    ## Lower scores indicate more relevancy to the transcription factor
    chea_results = pd.DataFrame(chea[db])
    chea_results = chea_results[['TF', 'Score']].set_index('TF')
    return chea_results.astype(float)

def get_kea3_results(gene_set, query_name='', db='Integrated--meanRank'):
    """ Inputs:
    gene_set: list[]
    query_name: str
    db: ['Integrated--meanRank', 'Integrated--topRank'], ranking method for KEA3 
    Return: list[] of kinases by decreasing score (higher match) """ 
    
    ADDLIST_URL = 'https://amp.pharm.mssm.edu/kea3/api/enrich/'
    payload = {
        'gene_set': gene_set,
        'query_name': query_name
    }
    
    # Retrieve response
    response = requests.post(ADDLIST_URL, data=json.dumps(payload))
    if not response.ok:
        raise Exception('Error analyzing gene list')
    data = json.loads(response.text)

    df = pd.DataFrame([[k['TF'], k['Score']] for k in data[db]], columns=['Kinase','Score'])
    df = df.set_index('Kinase')
    return df.astype(float)

def geneshot_set_augment(gene_list, similarity_type='coexpression', n_genes=100):
    """ Inputs:
    gene_list: list[]
    similarity_type: ['generif', 'tagger', 'coexpression', 'enrichr'], resource to do the augmentation
    n_genes: int, number of the most similar genes to return
    Return: list[] of augmented genes """

    GENESHOT_URL = 'https://maayanlab.cloud/geneshot/api/associate'
    payload = {
        "gene_list": gene_list,
        "similarity": similarity_type 
        }
    
    # Retrieve response
    response = requests.post(GENESHOT_URL, json=payload)
    data = json.loads(response.text)
  
    # Return augmented genes
    geneshot_dataframe = pd.DataFrame.from_dict(data['association'], orient="index").drop(['topGenes', 'topScores'], axis=1).sort_values(by='simScore', ascending=False)
    geneshot_dataframe = geneshot_dataframe.iloc[:n_genes]
    return geneshot_dataframe.index.to_list()

def get_g2n_results(geneset, subgraph_size=30, path_length=2):
    # Inputs: 
    # 'geneset' : ['gene name'] - needs a list of strings
    # 'databases' : ['database name'] - needs a list of strings,  
    # options: iid, bioGRID, STRING, bioPlex 3.0, default - all of them 
    # 'subgraph_size' : int, - max size of the subgraph you want to produce
    # 'path_length' : int, number of edges between target proteins 

    payload = {'geneset' : geneset, 'subgraph_size': str(subgraph_size), 'path_length': path_length} 

    # Retrieve response
    res = requests.post("https://g2nkg.dev.maayanlab.cloud/api/knowledge_graph/ppi_kg", json = payload)
    results = res.json()

    # Format results
    nodes_data = {}
    for i in results:
        if i["data"]["kind"] == 'Relation':
            continue
        else:
            nodes_data[i["data"]["id"]] = i["data"]
    
    # Return nodes
    nodes = pd.DataFrame.from_dict(nodes_data, orient="index")
    return nodes['label'].to_list()


def heatmap(dt: str, title='', width=800, height=800):
    """ Inputs:
    dt: string x2k, chea, kea"""

    rpy2.robjects.r("""
    function create_heatmap(data) {
        # User parameters
        # data = dt # x2k, chea, kea
        tracks <- c('baseline.tumor_stage_pathological', 'binned_age', 'Sex') # sample (column) annotation tracks
        N = 50 # number of genes (rows) to show in heatmap
        column_clust_method = 'average'
        row_clust_method = 'average'
        low_color = "blue"
        high_color = "red"

        # Data to plot
        top <- read.csv(sprintf("results/up_%s_zscores.csv", data), header=T, row.names=1, check.names = F)
        bottom <- read.csv(sprintf("results/dn_%s_zscores.csv", data), header=T, row.names=1, check.names = F)

        metadata <- read.csv("data/metadata.csv", header=T, row.names=1, check.names = F)

        # Align sample order
        col_order <- rownames(metadata)[rownames(metadata) %in% names(top)]
        top <- top %>% select(all_of(col_order))
        bottom <- bottom %>% select(all_of(col_order))
        metadata <- metadata[col_order, ]

        # Heatmap annotation dataframe
        anno_df <- metadata[tracks]
        ha = HeatmapAnnotation(df = anno_df,
                            simple_anno_size = unit(0.3, "cm"),
                            na_col = "white",
                            annotation_name_gp = gpar(fontsize = 8))

        # Rows to plot
        row_vars <- apply(top, 1, sd)
        sorted_indices <- order(row_vars, decreasing = TRUE)
        top_rows <- row.names(top)[sorted_indices[1:N]]
        row_vars <- apply(bottom, 1, sd)
        sorted_indices <- order(row_vars, decreasing = TRUE)
        bottom_rows <- row.names(bottom)[sorted_indices[1:N]]

        mat1 = data.matrix(top[rownames(top) %in% top_rows, ])
        mat2 = data.matrix(bottom[rownames(bottom) %in% bottom_rows, ])

        if (data == "kea") {
        split_df = c(rep('Enriched in up-phosphosites', length(top_rows)), rep('Enriched in down-phosphosites', length(bottom_rows)))
        } else {
        split_df = c(rep('Enriched in up-genes', length(top_rows)), rep('Enriched in down-genes', length(bottom_rows)))
        }

        colors = colorRamp2(seq(min(min(mat1), min(mat2)), max(max(mat1), max(mat2)), length = 3), c(low_color, "#EEEEEE", high_color))

        hm = Heatmap(rbind(mat1, mat2), 
                        name = "Enrichment score",
                        clustering_method_columns = column_clust_method,
                        clustering_method_rows = row_clust_method,
                        top_annotation = ha, 
                        row_split = split_df,
                        show_column_names = F,
                        column_title = NULL,
                        col = colors,
                        row_names_gp = grid::gpar(fontsize = 8),
                        heatmap_legend_param = list(legend_direction = "horizontal", title_position = "lefttop", legend_width = unit(5, "cm")),
                        row_gap = unit(3.2, "mm"))
        png(file = sprintf("results/%s_heatmap.png", data), width = 225, height = 205, units='mm', res = 300)
        draw(hm, heatmap_legend_side = "bottom")
        dev.off()
    }""")
