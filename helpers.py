import pandas as pd
from IPython.display import HTML


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
