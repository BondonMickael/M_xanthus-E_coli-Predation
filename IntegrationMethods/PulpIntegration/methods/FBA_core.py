'''
Module for the core functions of the FBA integration method with PulP.
'''

import numpy as np
import pandas as pd
import pulp
from cobra.io import read_sbml_model, write_sbml_model, load_json_model
import cobra

from IntegrationMethods.PulpIntegration.methods.utils.iMAT_utils import *

def run_classical_FBA(args,only_FBA=False):
    '''
    Objective reaction is hardcoded!!
    Args:
    ------
    only_FBA: bool
        decides whether result should be printed
    '''
    metabolicModel = read_sbml_model(args.model)
    metabolites = [met.id for met in metabolicModel.metabolites]
    reactions = [rxn.id for rxn in metabolicModel.reactions]

    prob = pulp.LpProblem("FBA", pulp.LpMaximize)
    
    flux={}
    for rxn in metabolicModel.reactions:
        flux[rxn.id] = pulp.LpVariable(f"flux_{rxn.id}", lowBound=rxn.bounds[0], upBound=rxn.bounds[1], cat='Continuous')
    
    S = cobra.util.array.create_stoichiometric_matrix(metabolicModel)
    num_mets, num_rxns = S.shape
    for i, met in enumerate(metabolites):
        dp = pulp.lpSum([S[i, j] * flux[reactions[j]] for j in range(num_rxns)])
        prob += (dp == 0), f"MassBalance_{met}"
    
    #Objecive function
    obj_reactions =[rct for rct in metabolicModel.reactions if rct.objective_coefficient !=0.0]
    objective_ids = [rct.id for rct in obj_reactions]
    if len(objective_ids) ==1:
        objective_id = obj_reactions[0].id
        print(f'Objective reaction: {objective_id}')
    else:
        print(f'Multiple objective reactions: {objective_ids}')
        raise ValueError (f'Not implemented yet. Sorry!')

    prob += flux[objective_id], f'Objective Function'

    prob.solve()
    cobraSol = metabolicModel.optimize()

    if only_FBA:
        print(f"\nStatus:{pulp.LpStatus[prob.status]}\nObjective (Total Penalty):, {pulp.value(prob.objective)}\nFlux values:")
        for rxn in flux:
            if flux[rxn].varValue != 0.0:
                print(f"\t{rxn}: {flux[rxn].varValue}\t{cobraSol.fluxes[rxn]}")
                
    objectiveValue = flux["OF_ATP_mitoMap"].varValue
    prob = pulp.LpProblem("remove free cycle FBA", pulp.LpMaximize)
    flux={}
    for rxn in metabolicModel.reactions:
        flux[rxn.id] = pulp.LpVariable(f"flux_{rxn.id}", lowBound=rxn.bounds[0], upBound=rxn.bounds[1], cat='Continuous')
    flux[objective_id] = pulp.LpVariable(f"flux_{objective_id}", lowBound=objectiveValue, upBound=rxn.bounds[1], cat='Continuous')
    
    S = cobra.util.array.create_stoichiometric_matrix(metabolicModel)
    num_mets, num_rxns = S.shape
    for i, met in enumerate(metabolites):
        dp = pulp.lpSum([S[i, j] * flux[reactions[j]] for j in range(num_rxns)])
        prob += (dp == 0), f"MassBalance_{met}"
    prob += flux[objective_id]
    prob.solve()
    if only_FBA:
        print(f"Cycle free FBA\nStatus:{pulp.LpStatus[prob.status]}\nObjective (Total Penalty):, {pulp.value(prob.objective)}\nFlux values:")
        for rxn in flux:
            if flux[rxn].varValue != 0.0:
                print(f"\t{rxn}: {flux[rxn].varValue}\t{cobraSol.fluxes[rxn]}")
    
#Could return objective reaction for iMAT