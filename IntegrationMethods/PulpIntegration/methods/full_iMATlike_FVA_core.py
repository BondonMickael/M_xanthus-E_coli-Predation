"""
Module for the core functions of full_iMATlike_FVA. Goal is to run
regular FVA, but constrain all reactions in set of highly active reactions
(RH) to carry at least a flux of +/-epsilon.

"""

import numpy as np
import pandas as pd
import pulp
import cobra
from cobra.io import read_sbml_model, write_sbml_model
import methods.utils.iMAT_utils as iMAT_utils
from cobra.core.gene import GPR
import methods.utils.iMAT_FVA_utils as iMAT_FVA_utils


def run_full_iMAT_FVA(args):
    """
    rct_of_interest, model, classification_df, fraction_of_optimum
    """
    metabolicModel = read_sbml_model(args.model)
    # Check if given reaction exist
    if args.fva_reactionOfInterest not in metabolicModel.reactions:
        raise ValueError(
            f"Reaction {args.fva_reactionOfInterest} does not exist in given model."
        )
    RNASeqData = iMAT_utils.generate_RNASeqDataFrame(
        args.model,
        args.expressionFile,
        args.geneColName,
        args.expressionColName,
        args.ignore_human,
    )
    # Discetization on scaled expression data

    (
        [lower_threshold, upper_threshold],
        [lower_threshold_scaled, upper_threshold_scaled],
        RNASeqData,
    ) = iMAT_utils.discretization(
        RNASeqData, args.expressionColName, args.discretization, args.quantiles
    )

    # map GPR to reaction
    penalizedReactions = []
    rewardedReactions = []
    active_genes = RNASeqData.index[RNASeqData["discretization"] == 1].tolist()
    nonActive_genes = RNASeqData.index[RNASeqData["discretization"] == -1].tolist()
    moderate_genes = RNASeqData.index[RNASeqData["discretization"] == 0].tolist()

    for rct in metabolicModel.reactions:
        # Remove mouse/human genes from gpr rules
        gpr_rule = iMAT_utils.filter_gpr(rct.gpr, ignore_human=args.ignore_human)
        gpr_rule = GPR().from_string(gpr_rule)
        if not GPR(gpr_rule).eval(nonActive_genes):
            penalizedReactions.append(rct.id)
        elif not GPR(gpr_rule).eval(nonActive_genes + moderate_genes):
            continue
        elif not GPR(gpr_rule).eval(active_genes + moderate_genes + nonActive_genes):
            rewardedReactions.append(rct.id)

    print(f"There are {len(rewardedReactions)} highly active reactions.")

    # Check feasibility
    infeasible_reactions = iMAT_FVA_utils.get_infeasible_reactions(
        rewardedReactions, metabolicModel, args.epsilon
    )

    # Create optimization Model
    objectiveValue = metabolicModel.optimize().objective_value
    obj_reactions = [
        rct for rct in metabolicModel.reactions if rct.objective_coefficient != 0.0
    ]
    objective_ids = [rct.id for rct in obj_reactions]
    if len(objective_ids) == 1:
        objective_id = obj_reactions[0].id
        print(f"Objective reaction: {objective_id}")
    else:
        print(f"Multiple objective reactions: {objective_ids}")
        raise ValueError(
            f"Multiple objective reactions cannot be handled fo now, sorry!"
        )
    print(
        f"The optimal solution is {objectiveValue}, the model will ask a minimum of {objectiveValue *args.optimum/100} for the objective reaction."
    )

    prob_max = pulp.LpProblem("full_iMAT_FVA_maximize", pulp.LpMaximize)
    prob_min = pulp.LpProblem("full_iMAT_FVA_min", pulp.LpMinimize)

    v_vars = {}
    y_var = {}
    for rct in metabolicModel.reactions:
        v_vars[rct.id] = pulp.LpVariable(
            f"v_{rct.id}",
            lowBound=rct.bounds[0],
            upBound=rct.bounds[1],
            cat="Continuous",
        )
    # Binary variable for highly active reversible reactions
    for rct in [rxt for rxt in rewardedReactions if rxt not in infeasible_reactions]:
        rct = metabolicModel.reactions.get_by_id(rct)
        if rct.lower_bound != 0 and rct.upper_bound != 0:
            # reversible reactions
            bigM = 1000
            y = pulp.LpVariable(f"y_{rct}", cat="Binary")
            y_var[rct.id] = y
            prob_max += v_vars[rct.id] >= args.epsilon - bigM * (1 - y)
            prob_max += v_vars[rct.id] <= -args.epsilon + bigM * y
            prob_min += v_vars[rct.id] >= args.epsilon - bigM * (1 - y)
            prob_min += v_vars[rct.id] <= -args.epsilon + bigM * y
        elif rct.lower_bound >= 0:
            # forward reaction
            prob_max += v_vars[rct.id] >= args.epsilon
            prob_min += v_vars[rct.id] >= args.epsilon
        else:
            # backward reactions
            prob_max += v_vars[rct.id] <= args.epsilon
            prob_min += v_vars[rct.id] <= args.epsilon

    # Mass balance constraint
    for met in metabolicModel.metabolites:
        constraint = (
            pulp.lpSum(
                rct.metabolites.get(met, 0) * v_vars[rct.id] for rct in met.reactions
            )
            == 0
        )
        prob_max += constraint, f"mass_balance_{met.id}"
        prob_min += constraint, f"mass_balance_{met.id}"
    # objective reaction
    prob_max += v_vars[objective_id] >= objectiveValue * args.optimum / 100
    prob_min += v_vars[objective_id] >= objectiveValue * args.optimum / 100

    prob_max += v_vars[args.fva_reactionOfInterest]
    prob_min += v_vars[args.fva_reactionOfInterest]

    prob_min.solve()
    if pulp.LpStatus[prob_min.status] == "Infeasible":
        raise ValueError(f"The minimization problem is infeasible")

    fva_min = v_vars[args.fva_reactionOfInterest].varValue
    prob_max.solve()
    if pulp.LpStatus[prob_max.status] == "Infeasible":
        raise ValueError(f"The maximization problem is infeasible")
    fva_max = v_vars[args.fva_reactionOfInterest].varValue

    print(f"{args.fva_reactionOfInterest}\nMin: {fva_min}\tMax: {fva_max}")
    # for rct in rewardedReactions:
    #     print(f'{rct}: {v_vars[rct].varValue}')
