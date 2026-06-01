import numpy as np
import re
from os.path import join
from cobra.util.solver import linear_reaction_coefficients


def read_csv_data(
    path, id_gene, id_val, head, list_of_genes, sep=",", comment="_", quantile=np.nan
):
    """This function is used to read a csv file and return a dictionnary with gene id as key and expression value. The expression value is also divided by it's maximum to scale it between 0 and 1.

    Args:
        path (str): Path to the txt file.
        id_gene (int): Column number of the gene id in the data
        id_val (int): Column number of the expression value
        head (bool): Bool to remove the first line or not (header)
        list_of_genes (list): a list of genes used to only select the genes that are inside the model for the scaling.
        sep (char): Separator used in the file (default is ",")
        comment (char): character used to comment in the data file, generally _ or # (the line that start with this char will be ignored)

        Returns:
            dict : a dict that can be used for the eflux method
    """
    dict_exp = {}
    with open(path) as file:
        for line in file:
            if head:
                head = False
                continue
            if line.startswith(comment):  # If the line is a comment ignore it
                continue
            line = line.replace('"', "").split(sep)
            if line[id_gene] in list_of_genes:
                dict_exp[line[id_gene]] = float(
                    line[id_val]
                )  # create the entry in the dictionnary with the gene id and its value
    if quantile == quantile:
        quantile = np.quantile(list(dict_exp.values()), quantile)
        for i in dict_exp:
            dict_exp[i] = dict_exp[i] / quantile
            if dict_exp[i] > 1:
                dict_exp[i] = 1
    else:
        maximum = max(dict_exp.values())
        for i in dict_exp:  # Divide each entry by the maximum value
            dict_exp[i] = dict_exp[i] / maximum

    return dict_exp  # return the dictionnary


def recursive_evaluation(op_list, depth=0, debbug=False):
    """
    Recursively evaluates a parsed GPR expression represented as a list.

    Args:
        op_list (list): Expression list with values and operations ('min', 'sum', '(', ')').
        depth (int): Internal use for recursion depth tracking.
        debbug (bool): If True, prints debbug info.

    Returns:
        float: Final evaluated expression value.
    """
    if debbug:
        print("  " * depth + f"Evaluating: {op_list}")
    if len(op_list) == 1:
        return op_list[0]
    i = 0
    while i < len(op_list):
        if op_list[i] == "(":
            # Find the matching closing parenthesis
            level = 1
            j = i + 1
            while j < len(op_list) and level > 0:
                if op_list[j] == "(":
                    level += 1
                elif op_list[j] == ")":
                    level -= 1
                j += 1
            # Evaluate the sub-expression recursively
            inner = op_list[i + 1 : j - 1]
            value = recursive_evaluation(inner, depth=depth + 1, debbug=debbug)

            # Replace the whole subexpression with its result
            op_list = op_list[:i] + [value] + op_list[j:]
            if debbug:
                print("  " * depth + f"After eval: {op_list}")
            return recursive_evaluation(op_list, depth=depth, debbug=debbug)
        i += 1
    # No parentheses left, apply operators
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
        elif token == "sum":
            result = (
                result + op_list[i + 1]
                if result is not None
                else op_list[i - 1] + op_list[i + 1]
            )
            i += 2
        else:
            if result is None and not isinstance(token, str):
                result = token
            i += 1
    if debbug:
        print("  " * depth + f"Returning: {result}")
    return result


def read_gpr_from_reaction(
    model,
    identifier,
    expr_dict,
    const,
    default_expr_val,
    ignore_mus=False,
    ignore_human=False,
    debbug=False,
):
    """
    Applies the eFlux method to a specific reaction in a COBRA model by adjusting flux bounds
    based on gene expression data and the reaction's GPR (Gene-Protein-Reaction) rule.

    Args:
        model: COBRApy metabolic model. This function modifies bounds directly in this model.
        identifier (str): Reaction ID used to access model.reactions[identifier].
        expr_dict (dict): Gene expression dictionary (gene_id -> expression value scaled between 0 and 1).
        const (float): Max value used to scale reaction bounds.
        default_expr_val (float): Expression value used when no GPR or expression data is available.
        ignore_mus (bool): If True, ignore Mus musculus (ENSMUS*) genes.
        ignore_human (bool): If True, ignore human (ENSG*) genes.
        debbug (bool): Print debbug messages if True.

    Returns:
        int: 1 if the reaction was constrained, 0 otherwise.
    """
    reaction = model.reactions[identifier]
    if reaction in linear_reaction_coefficients(model):
        return 0  # Skip objective function reactions

    gpr = str(reaction._gpr)
    # Part only used for Mitomammal
    if ignore_mus:
        gpr = remove_genes(gpr, "ENSMUS")
    if ignore_human:
        gpr = remove_genes(gpr, "ENSG")

    if gpr.strip() == "":
        final_value = default_expr_val
    else:
        final_value = evaluate_gpr_expression(gpr, expr_dict, default_expr_val, debbug)

    if final_value == -1:
        final_value = 1
        constrained = 0
        if debbug:
            print(f"Reaction '{identifier}' has no specific constraint.")
    else:
        constrained = 1
    upper = float(final_value * const)
    lower = 0 if reaction.lower_bound == 0 else -upper
    reaction.bounds = (lower, upper)
    if debbug:
        print(
            f"Final expression value: {final_value}, bounds set to: ({lower}, {upper})"
        )

    return constrained


def remove_genes(gpr_str, gene_prefix):
    """Remove genes with a specific prefix (e.g., 'ENSMUS' or 'ENSG') from the GPR string."""
    tokens = gpr_str.split(" ")
    indices_to_remove = set()

    for i in range(len(tokens) - 1, -1, -1):
        if gene_prefix in tokens[i]:
            indices_to_remove.update({i - 1, i, i + 1})

    tokens = [tok for idx, tok in enumerate(tokens) if idx not in indices_to_remove]
    return " ".join(tokens)


def evaluate_gpr_expression(gpr, expr_dict, default_val, debbug=False):
    """
    Evaluate the expression in a GPR string using gene expression values.

    Args:
        gpr (str): The GPR string.
        expr_dict (dict): Expression values (gene_id -> expression).
        default_val (float): Default expression value for unknown genes.
        debbug (bool): If True, prints debbug info.

    Returns:
        float: Evaluated expression value.
    """
    # Replace logical operators with Pythonic versions
    gpr = gpr.replace(" or ", " OR ").replace(" and ", " AND ")

    # Tokenize while keeping parentheses and operators
    tokens = re.findall(r"\(|\)|AND|OR|[^\s()]+", gpr)
    parsed_tokens = [parse_token(tok, expr_dict, default_val) for tok in tokens]

    if debbug:
        print("Parsed GPR expression:")
        print(parsed_tokens)

    return recursive_evaluation(parsed_tokens, debbug=debbug)


def parse_token(token, expr_dict, default_val):
    if token == "AND":
        return "min"
    elif token == "OR":
        return "sum"
    elif token in ("(", ")"):
        return token
    else:
        return expr_dict.get(token, default_val)


def constrain_exchanges(model, quantile):
    lb = []
    for rxn in model.reactions:
        if rxn in model.exchanges or rxn in model.sinks or rxn in model.demands:
            continue
        lb.append(rxn._upper_bound)
    constraint = np.quantile(lb, quantile)
    for ex in model.exchanges:
        ex.bounds = {-constraint, constraint}
    for dm in model.demands:
        dm.bounds = {0, constraint}
    for sk in model.sinks:
        sk.bounds = {-constraint, 0}


def Eflux(
    model,
    dict_val,
    const,
    default_exp_val,
    ignore_mus=False,
    ignore_human=False,
    debbug=False,
):
    """This function is used to apply the eflux method to a cobra model.

    Args:
        model (): Metabolic model opened with cobrapy. The function will change the upper and lower bound in this model.
        dict_val (dict): Dictionnary with the gene id as key and the expression as value (scaled between 0 and 1)
        const (int): This is the value selected to set the flux, They will have a range between 0 and this value.
        default_exp_val (float): Default value used to set the reactions without GPR or for the genes without expression data
        ignore_mus (bool): value set to ignore the mus musculus genes (default value: False)
        ignore_human (bool): value set to ignore the human genes (default value: False)
    """
    constrained = 0
    for i in range(len(model.reactions)):
        if (
            model.reactions[i] in model.exchanges
            or model.reactions[i] in model.demands
            or model.reactions[i] in model.sinks
        ):  # If the reaction is a sink, an exchange or a demand it continue without changing it's boundaries.
            continue
        if debbug:
            print(f"\nReaction: {model.reactions[i].id}\tGPR:{model.reactions[i].gpr}")
        constrained += read_gpr_from_reaction(
            model,
            i,
            dict_val,
            const,
            default_exp_val,
            ignore_mus=ignore_mus,
            ignore_human=ignore_human,
            debbug=debbug,
        )
