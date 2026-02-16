"""
Utility functions for (weighted) iMAT integration method.
"""

import pandas as pd
import cobra
import numpy as np
from cobra.io import read_sbml_model


def read_expression_file(path):
    if path.endswith(".tsv"):
        df = pd.read_csv(path, sep="\t")
    else:
        df = pd.read_csv(path)
    return df


def generate_RNASeqDataFrame(
    metabolicModel, expressionFile, geneColName, expressionColName, ignore_human
):
    """
    Generate a processed RNA-seq DataFrame aligned with a metabolic model.

    This function reads a gene expression file (CSV or TSV) and filters it to include
    only genes present in the provided COBRApy metabolic model. It also ensures that
    all model genes (of the relevant species) are represented, adding missing genes
    as new rows with NaN values for their expression.

    Parameters
    ----------
    metabolicModel : cobra.core.Model
        Cobra model
    expressionFile : str
        Path to a .csv or .tsv file containing gene expression measurements.
    geneColName : str
        Column name containing gene IDs.
    expressionColName : str
        Column name containing expression values.
    ignore_human : bool, optional
        If True (default), ignore human genes (ENSG...) and focus on mouse (ENSMUS...).
        If False, use human genes instead.

    Returns
    -------
    RNASeqData : pandas.DataFrame
        DataFrame indexed by gene IDs, containing a single column of expression values,
        with NaN entries for genes present in the model but missing in the data file.
    """
    genes_in_model = []
    for gene in metabolicModel.genes:
        genes_in_model.append(gene.id)

    RNASeqData = read_expression_file(expressionFile)
    RNASeqData.set_index(geneColName, inplace=True)

    genes = []
    missingGenes = []
    species_prefix = "MXAN" if ignore_human else "ENSG"

    for i in metabolicModel.genes:
        if not i.id.startswith(species_prefix):
            continue
        (genes if i.id in RNASeqData.index else missingGenes).append(i.id)
    print(f"There are {len(missingGenes)} missing genes: {';'.join(missingGenes)}")
    # df with intersection of genes in model and in RNASeqData and only gene expression values of interest
    RNASeqData = RNASeqData.loc[genes, [expressionColName]]
    # Add missingGenes with NaN entries
    new_entries = pd.DataFrame(index=missingGenes, columns=[expressionColName])
    RNASeqData = pd.concat([RNASeqData, new_entries])

    return RNASeqData


def discretization(RNASeq_Df, expressionColName, discretization_method, quantiles=None):
    """
    Depending on the method chosen, this function adds a column 'scaled_expression' and 'discretization' to
    the dataframe.
    Args
    -----
    RNASeq_Df: pandas.DataFrame
        dataframe with gene ids and expression values
    expressionColName : str
        Column name containing expression values.
    discretization_method: str
        string of discretization method chosen. Options are 'mean' or 'quantile'
    quantiles: list of int
        two integer for lower and upper quantile

    Returns
    -------
    RNASeq_Df: pandas.DataFrame
        with added scaled_expression and discretization column
    [lower_threshold, upper_threshold]: list of float
        threshold used for discretization (scaled and unscaled thresholds)

    """
    RNASeq_Df["scaled_expression"] = (
        RNASeq_Df[expressionColName] - RNASeq_Df[expressionColName].min()
    ) / (RNASeq_Df[expressionColName].max() - RNASeq_Df[expressionColName].min())

    if discretization_method == "mean":
        mean = RNASeq_Df[expressionColName].mean()
        sd = RNASeq_Df[expressionColName].std()
        upper_threshold = mean + 1 / 2 * sd
        lower_threshold = mean - 1 / 2 * sd
        upper_threshold_scaled = (
            RNASeq_Df["scaled_expression"].mean()
            + 1 / 2 * RNASeq_Df["scaled_expression"].std()
        )
        lower_threshold_scaled = (
            RNASeq_Df["scaled_expression"].mean()
            - 1 / 2 * RNASeq_Df["scaled_expression"].std()
        )

        print(f"Thresholds for mean discretization:")
        classification_rules = [
            RNASeq_Df[expressionColName] > upper_threshold,
            RNASeq_Df[expressionColName] < lower_threshold,
        ]
        classes = [1, -1]
        RNASeq_Df["discretization"] = np.select(
            classification_rules, classes, default=0
        )

    elif discretization_method == "quantile":
        # Quantile discretization
        q1, q2 = quantiles
        lower_threshold = RNASeq_Df[expressionColName].quantile(q1 / 100)
        upper_threshold = RNASeq_Df[expressionColName].quantile(q2 / 100)
        lower_threshold_scaled = RNASeq_Df["scaled_expression"].quantile(q1 / 100)
        upper_threshold_scaled = RNASeq_Df["scaled_expression"].quantile(q2 / 100)
        print(f"Thresholds for quantile discretization")
        print(f"quantiles: {quantiles}")
        RNASeq_Df["discretization"] = 0
        RNASeq_Df.loc[
            RNASeq_Df[expressionColName] < lower_threshold, "discretization"
        ] = -1
        RNASeq_Df.loc[
            RNASeq_Df[expressionColName] > upper_threshold, "discretization"
        ] = 1

    print(
        f"Mean: {RNASeq_Df[expressionColName].mean()}\n"
        f"Median: {np.nanquantile(RNASeq_Df[expressionColName],0.5)}\n"
        f"Lower threshold: {lower_threshold}\n"
        f"Upper threshold: {upper_threshold}\n"
        f"Scaled thresholds: ({lower_threshold_scaled,upper_threshold_scaled})"
    )

    return (
        [lower_threshold, upper_threshold],
        [lower_threshold_scaled, upper_threshold_scaled],
        RNASeq_Df,
    )


def filter_gpr(gpr_rule, ignore_human=False):
    """'
    Filters gpr rules from human/mouse genes/

    Args:
    ------
    gpr_rule: cobra.core.gene.GPR
        gpr rule from cobra model
    ignore_hum: Boolean
        if True drops ENSMUS genes.

    Returns:
    ------
    filtered_gpr: str
        Filtered gpr expression
    """
    gpr_rule = str(gpr_rule)
    tokens = gpr_rule.split(" ")
    if ignore_human:
        prefix = "ENSG"
    else:
        prefix = "ENSMUS"
    indices_to_remove = set()
    for i in range(len(tokens) - 1, -1, -1):
        if prefix in tokens[i]:
            indices_to_remove.update({i - 1, i, i + 1})

    tokens = [tok for idx, tok in enumerate(tokens) if idx not in indices_to_remove]
    if tokens and tokens[-1].lower() == "or":
        tokens = tokens[:-1]
    filtered_gpr = " ".join(tokens)

    return filtered_gpr


def recursive_evaluation_iMAT(op_list, depth=0, debbug=False):
    """
    Recursively evaluates a parsed GPR expression represented as a list.

    Args:
        op_list (list): Expression list with values and operations ('min', 'max', '(', ')').
        depth (int): Internal use for recursion depth tracking.
        debbug (bool): If True, prints debbug info.

    Returns:
        float: Final evaluated expression value.
    """
    if debbug:
        print("  " * depth + f"Evaluating: {op_list}")
    if len(op_list) == 1:
        return op_list[0]
    # Loops from outside to the innermost parenthesis
    i = 0
    while i < len(op_list):
        if op_list[i] == "(":
            level = 1
            j = i + 1
            # Find matching closing parenthesis
            while j < len(op_list) and level > 0:
                # next token
                if op_list[j] == "(":
                    level += 1
                elif op_list[j] == ")":
                    level -= 1
                j += 1

            inner = op_list[i + 1 : j - 1]
            value = recursive_evaluation_iMAT(inner, depth=depth + 1, debbug=debbug)
            op_list = op_list[:i] + [value] + op_list[j:]
            if debbug:
                print("  " * depth + f"After eval: {op_list}")
            return recursive_evaluation_iMAT(op_list, depth=depth, debbug=debbug)
        i += 1
    result = None
    i = 0
    while i < len(op_list):
        token = op_list[i]
        if token == "min":
            result = (
                min(result, op_list[i + 1])
                if result is not None
                else min(op_list[i - 1], op_list[i + 1])
            )
            i += 2
        elif token == "max":
            result = (
                max(result, op_list[i + 1])
                if result is not None
                else max(op_list[i - 1], op_list[i + 1])
            )
            i += 2
        else:
            if result is None and not isinstance(token, str):
                result = token
            i += 1
    if debbug:
        print("  " * depth + f"Returning: {result}")
    return result


def parse_token(token, df_expression, default_value=0.5):
    """
    Replaces gene id token with expression values from the dataframe and replaces logical
    operators with min/max

    Args
    -----
    token: str
        token from the gpr rule
    df_expression: pandas.DataFrame
        DataFrame with gene ids as index and a column 'scaled_expression' with scaled expression value
    default_value: float
        default value to return if gene id is not found in the dataframe (e.g., gene was not measured)

    Returns
    -------
    exchanged value: float/str
    """
    if token == "AND":
        return "min"
    elif token == "OR":
        return "max"
    elif token in ("(", ")"):
        return token
    else:
        expr_series = df_expression["scaled_expression"]
        value = expr_series.get(token, default_value)
        if pd.isna(value):
            # Remove the possibility of NaN values when scaled_expression is generated
            return default_value
        return value
