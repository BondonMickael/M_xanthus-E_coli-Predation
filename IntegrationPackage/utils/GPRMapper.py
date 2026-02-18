from dataclasses import dataclass, field
import pandas as pd
from cobra import Model
from cobra.core.gene import GPR
from typing import List
import re


@dataclass
class GPRMapperOutput:
    RL: List[str]
    RM: List[str]
    RH: List[str]


@dataclass
class GPRMapper:
    metabolicModel: Model
    expression_df: pd.DataFrame
    ignore_human: bool

    RL: list = field(init=False, default_factory=list)
    RM: list = field(init=False, default_factory=list)
    RH: list = field(init=False, default_factory=list)

    def get_reaction_expression(self, gpr: str):
        """
        For a given reaction, finds the leading gene according to the GPR rule and returns it's scaled gene expression.
        And evaluates the GPR rule.

        Args
        -----
        rct: cobra.core.reaction.Reaction
            reaction from the cobra model
        gene_expression: pandas.DataFrame
            DataFrame with gene ids as index and a column 'scaled_expression' with scaled expression value (between 0 and 1)
            and a column 'discretization' with discretized values (-1, 0 ,1)

        Returns
        -------
        leading_gene_expression: float
            Evaluated expression value (if OR, max value, if AND min value of leading gene)
        """
        gpr_rule = self._filter_gpr(gpr)
        gpr_rule = gpr_rule.replace(" or ", " OR ").replace(" and ", " AND ")
        # tokenize while keeping parentheses and operators, all elements of rules correspond to one list entry
        tokens = re.findall(r"\(|\)|AND|OR|[^\s()]+", gpr_rule)
        # Exchange gene ids and logical expressions (and, or) with expression values and (min, max)
        parsed_tokens = [self._parse_token(tok) for tok in tokens]
        val = self._recursive_evaluation(parsed_tokens)
        return val

    def create_reaction_classes(self) -> GPRMapperOutput:
        """
        Maps model GPR rules to gene expression dataframe and returns list of reactions IDs
        within the classes lowly-, moderate-, and highly active reactions, RL, RM, and RH respectively.
        """
        active_genes = self.expression_df.index[
            self.expression_df["discretization"] == 1
        ].tolist()
        nonActive_genes = self.expression_df.index[
            self.expression_df["discretization"] == -1
        ].tolist()
        moderateActive_genes = self.expression_df.index[
            self.expression_df["discretization"] == 0
        ].tolist()

        for rct in self.metabolicModel.reactions:
            gpr_rule = self._filter_gpr(str(rct.gpr))
            gpr_rule = GPR().from_string(gpr_rule)

            if not GPR(gpr_rule).eval(nonActive_genes):
                self.RL.append(rct.id)
            elif not GPR(gpr_rule).eval(nonActive_genes + moderateActive_genes):
                self.RM.append(rct.id)
            elif not GPR(gpr_rule).eval(
                active_genes + moderateActive_genes + nonActive_genes
            ):
                self.RH.append(rct.id)

        return GPRMapperOutput(RL=self.RL, RM=self.RM, RH=self.RH)

    def _filter_gpr(self, gpr_rule: str):
        """'
        Filters gpr rules from human/mouse genes/
        """
        tokens = gpr_rule.split(" ")
        if self.ignore_human:
            prefix = "ENSG"
        else:
            prefix = "MXAN"
        indices_to_remove = set()
        for i in range(len(tokens) - 1, -1, -1):
            if prefix in tokens[i]:
                indices_to_remove.update({i - 1, i, i + 1})

        tokens = [tok for idx, tok in enumerate(tokens) if idx not in indices_to_remove]
        if tokens and tokens[-1].lower() == "or":
            tokens = tokens[:-1]
        filtered_gpr = " ".join(tokens)
        return filtered_gpr

    def _parse_token(self, token: str, default_value=0.5):
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
            value = self.expression_df.loc[token, "scaled_expression"]
            if pd.isna(value):
                # Remove the possibility of NaN values when scaled_expression is generated
                return default_value
        return value

    def _recursive_evaluation(self, op_list: List, depth=0, debug=False):
        """
        Recursively evaluates a parsed GPR expression represented as a list.

        Args:
            op_list (list): gpr rule list with gene ids and operations ('min', 'max', '(', ')').
            depth (int): Internal use for recursion depth tracking.
            debug (bool): If True, prints debug info.
        Returns:
            float: Final evaluated expression value.
        """
        if debug:
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
                # i and j will be the parentesises of the inner most rule
                inner = op_list[i + 1 : j - 1]
                value = self._recursive_evaluation(inner, depth=depth + 1, debug=debug)

                op_list = op_list[:i] + [value] + op_list[j:]
                if debug:
                    print("  " * depth + f"After eval: {op_list}")
                return self._recursive_evaluation(op_list, depth=depth, debug=debug)
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
        if debug:
            print("  " * depth + f"Returning: {result}")
        return result
