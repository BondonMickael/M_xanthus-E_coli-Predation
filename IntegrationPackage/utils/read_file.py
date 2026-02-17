import pandas as pd

def read_expression_file(path):
    if path.endswith('.tsv'):
        df = pd.read_csv(path, sep='\t')
    else:
        df = pd.read_csv(path)
    return df