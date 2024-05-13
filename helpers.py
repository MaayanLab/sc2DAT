import os
from IPython.display import HTML
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import requests
import json
import time
from tqdm import tqdm
import openai
from dotenv import load_dotenv
load_dotenv()

client = openai.Client(api_key=os.getenv("OPENAI_API_KEY"))

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



def enrichr_figure(res_list: list): 
    all_terms,all_pvalues, all_adjusted_pvalues, all_libraries = res_list
    # Bar colors
    bar_color_not_sig = 'lightgrey'
    bar_color = 'lightblue'
    edgecolor=None
    linewidth=0
    fig, axes = plt.subplots(nrows=len(all_libraries), ncols=1)
        

    for i, library_name in enumerate(all_libraries):
        bar_colors = [bar_color if (x < 0.05) else bar_color_not_sig for x in all_pvalues[i]]
        sns.barplot(x=np.log10(all_pvalues[i])*-1, y=all_terms[i],ax=axes[i], palette=bar_colors, edgecolor=edgecolor, linewidth=1)
        axes[i].axes.get_yaxis().set_visible(False)
        axes[i].set_title(library_name.replace('_',' '),fontsize=30)
        if i == len(all_libraries)-1:
            axes[i].set_xlabel('-Log10(p-value)',fontsize=30)
        axes[i].xaxis.set_major_locator(MaxNLocator(integer=True))
        axes[i].tick_params(axis='x', which='major', labelsize=20)
        if max(np.log10(all_pvalues[i])*-1)<1:
            axes[i].xaxis.set_ticks(np.arange(0, max(np.log10(all_pvalues[i])*-1), 0.1))
        for ii,annot in enumerate(all_terms[i]):
            if all_adjusted_pvalues[i][ii] < 0.05:
                annot = '  *'.join([annot, str(str(np.format_float_scientific(all_pvalues[i][ii],precision=2)))]) 
            else:
                annot = '  '.join([annot, str(str(np.format_float_scientific(all_pvalues[i][ii],precision=2)))])

            title_start= max(axes[i].axes.get_xlim())/200
            axes[i].text(title_start,ii,annot,ha='left',wrap = True, fontsize = 30)
    plt.subplots_adjust(top=4.8, right = 4.7,wspace = 0.03,hspace = 0.2)
    plt.show()
    return fig


def enrich_libraries(user_list_id: str, all_libraries: list = ['WikiPathway_2023_Human', 'GO_Biological_Process_2023', 'MGI_Mammalian_Phenotype_Level_4_2021']):

    all_terms = []
    all_pvalues =[] 
    all_adjusted_pvalues = []
    library_success = []

    for library_name in all_libraries: 
        
        ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/enrich'
        query_string = '?userListId=%s&backgroundType=%s'
        gene_set_library = library_name
        response = requests.get(
            ENRICHR_URL + query_string % (user_list_id, gene_set_library)
         )
        if not response.ok:
            raise Exception('Error fetching enrichment results')
        try:
            data = json.loads(response.text)
            results_df  = pd.DataFrame(data[library_name][0:5])
            all_terms.append(list(results_df[1]))
            all_pvalues.append(list(results_df[2]))
            all_adjusted_pvalues.append(list(results_df[6]))
            library_success.append(library_name)
        except:
            print('Error for ' + library_name + ' library')
        time.sleep(1)

    return [all_terms, all_pvalues, all_adjusted_pvalues, library_success]


def label_clusters(cluster_enrichments: dict):
    """ Inputs:
    cluster_enrichments: dict, dictionary of cluster enrichments
    Return: dict, dictionary of labeled clusters (consensus with GPT-4) """
    
    cluster_labels = {}

    for c in tqdm(cluster_enrichments):
        # Get top 5 enrichments
        up_enrichments = []
        down_enrichments = []
        for i in range(3): # libraries
            terms = cluster_enrichments[c]['up'][0][i]
            p_values = cluster_enrichments[c]['up'][1][i]
            up_enrichments.extend([terms[j] for j in range(len(terms)) if p_values[j] < 0.01])

            terms = cluster_enrichments[c]['down'][0][i]
            p_values = cluster_enrichments[c]['down'][1][i]
            down_enrichments.extend([terms[j] for j in range(len(terms)) if p_values[j] < 0.01])
        # Get GPT-4 labels
        
        response = client.chat.completions.create(
            model='gpt-4-turbo',
            messages=[{'role': 'system', 'content': 'You are an assistant to a biologist who has performed enrichment analysis on a set of up-regulated and down-regulated genes for a certain cluster of patient samples. You must determine a consensus label for a given cluster based on the signfigantly enriched terms which relate to biological processes and phenotypes.'},
                      {'role': 'user', 'content': f"The most significantly enriched terms for the upregulated genes of cluster {str(c)} are: {', '.join(up_enrichments)}. For the down genes the significantly enriched terms are: {', '.join(down_enrichments)}. Please provide a consensus label for this cluster with no other reasoning. The label should be at maximum 5 words in length."}],
            max_tokens=50,
            temperature=0.3
        )
        cluster_labels[c] = response.choices[0].message.content
    return cluster_labels

def convert_to_string(info_dict: dict):
    res_str = ""
    for key in info_dict:
        if isinstance(info_dict[key], list):
            res_str += f"{key}: {','.join(info_dict[key])}\n"
        elif isinstance(info_dict[key], dict): res_str += convert_to_string(info_dict[key])
        else: res_str += f"{key}: {info_dict[key]}\n" 
    return res_str

def create_results_text(results: dict):
    """ Inputs:
    results: dict, dictionary of results -- keys represent 'section' and values are relevant data to be included in the discussion
    Return: str, discussion section in markdown
    """
    prompt = ""
    for k in results:
        prompt += f"## {k}\n"
        prompt += results[k] + "\n"
    system_prompt = 'You are an assistant to a biologist who has performed various analysis on gene, protein and phosphoprotein data related to tumor expression. A certain result or set of results will be provided with a section header that describes the analysis. Do not include the headers in your response. Write a discussion of the results that mainly describes the results without interpretation. You may specifically denote differences between clusters or specific samples/patients.'
    response = client.chat.completions.create(
            model='gpt-4-turbo',
            messages=[{'role': 'system', 'content': system_prompt},
                      {'role': 'user', 'content': prompt}],
            temperature=0.5
        )
    text = response.choices[0].message.content
    return text

def create_results_text_prompt(results, desc):
    """ Inputs:
    results: dict, dictionary of results -- keys represent 'section' and values are relevant data to be included in the discussion
    Return: str, discussion section in markdown
    """
    prompt = desc + '\n' + str(results)
    system_prompt = 'You are an assistant to a biologist who has performed various analysis on gene, protein and phosphoprotein data related to tumor expression. A certain result or set of results will be provided with a section header that describes the analysis. Do not include the headers in your response. Write a discussion of the results that mainly describes the results without interpretation. You may specifically denote differences between clusters or specific samples/patients.'
    response = client.chat.completions.create(
            model='gpt-4-turbo',
            messages=[{'role': 'system', 'content': system_prompt},
                      {'role': 'user', 'content': prompt}],
            temperature=0.5
        )
    text = response.choices[0].message.content
    return text

def clean_enrichr_lables(enrichr_labels: dict):
    enrichr_labels_clean = {}    
    for clus_num, enrichments in enrichr_labels.items():
            cluster_name = f'Cluster {clus_num}'
            enrichr_labels_clean[cluster_name] = {}
            
            for dir, values in enrichments.items():
                terms = values[0]
                pvals = values[1]
                adjpvals = values[2]
                libraries = values[3]
                for i, lib in enumerate(libraries):
                    if lib not in enrichr_labels_clean[cluster_name]:
                        enrichr_labels_clean[cluster_name][lib] = {}
                    if dir not in enrichr_labels_clean[cluster_name][lib]:
                        enrichr_labels_clean[cluster_name][lib][dir] = []
                    lib_terms = terms[i]
                    lib_adjpvals = adjpvals[i]
                    for j, term in enumerate(lib_terms):
                        if lib_adjpvals[j] < 0.05:
                            enrichr_labels_clean[cluster_name][lib][dir].append((term, "{:.4g}".format(lib_adjpvals[j])))
    return enrichr_labels_clean

def describe_clusters(results: dict):
    """ Inputs:
    results: dict, dictionary of results -- keys represent clusters and values are libraries and direction of signficant labels to be included in the discussion
    Return: str, discussion section in markdown
    """
    prompt = str(results)
    system_prompt = 'You are an assistant to a biologist who has performed enrichment analysis on significant terms for a set of clusters. Each cluster has been analyzed for significant terms in different libraries and directions. A certain result or set of results will be provided with a section header that describes the analysis. Do not include the headers in your response. Write a discussion of the results that mainly describes the common themes of the enriched terms from the up and down genes. You may specifically denote differences between clusters or specific samples/patients. When referencing a term please use the full term name (with NO quotes) followed by the adj. pvalue as follows: (term, adj. pvalue = 0.00123)'
    response = client.chat.completions.create(
            model='gpt-4-turbo',
            messages=[{'role': 'system', 'content': system_prompt},
                      {'role': 'user', 'content': prompt}],
            temperature=0.5
        )
    text = response.choices[0].message.content
    return text
    
