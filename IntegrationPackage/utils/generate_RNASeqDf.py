
import pandas as pd

def generate_RNASeqDf(metabolicModel, expression_df,geneColName, expressionColName):
    '''
    Generates a processed RNA-seq Dataframe aligned with the metabolic Model. Additionally,
    genes present in the metabolicModel but not in the expression_df are added as new rows with NaN 
    entries.
    
    :param metabolicModel: cobra.core.Model
    :param expression_df: pandas.df
        dataframe created from input file
    :param geneColName: str
        Column name containing gene IDs
    :param expressionColName: str
        Column name containing expression values
    :param ignore_human: bool, optional
    
    Returns
    -------
    RNASeqData : pandas.DataFrame
    DataFrame indexed by gene IDs, containing a single column of expression values,
    with NaN entries for genes present in the model but missing in the data file.
    '''
    expression_df.set_index(geneColName, inplace = True)
    if expression_df.index[0].startswith('MXAN'):
        species_prefix = 'MXAN'
        
    elif expression_df.index[0].startswith('MXAN'):
        species_prefix = 'ENSG'
    else: 
        raise SyntaxError('Verify gene ids in gene expressions. Only human or mice gene IDs allowed.')
    
    genes = []
    missingGenes = []
    for g in metabolicModel.genes:
        if not g.id.startswith(species_prefix):
            continue
        (genes if g.id in expression_df.index else missingGenes).append(g.id)
    expression_df = expression_df.loc[genes, [expressionColName]]
    new_entries = pd.DataFrame(index = missingGenes, columns = [expressionColName])
    expression_df = pd.concat([expression_df, new_entries])
    return expression_df