import pandas as pd
from IPython.display import HTML

import requests
import json

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
    return chea_results

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
    return df

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
