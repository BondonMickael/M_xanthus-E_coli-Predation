import pulp
from pathlib import Path
import os
import numpy as np
import pandas as pd
from pulp import LpStatus

def print_output(lpProblem, v_vars, rewardedReactions, penalizedReactions):
    '''
    Solves a linear Programming problem, collects flux distribution and prints result to stdout.

    Args
    -------
    lpProblem: pulp.LpProblem

    rewardedReactions/penalizedReactions: list
        list of reaction ids of penalized/rewarded reactions
    v_vars, y_vars: dict
        Dictionary with entries per reaction id for iMAT problems. 
    '''
    status = lpProblem.solve(pulp.PULP_CBC_CMD(msg=True))
    
    if pulp.LpStatus[lpProblem.status] == 'Infeasible':
        raise ValueError(f'Problem is infeasible')
    fluxes = {rid.replace('v_', ''): v for rid, v in v_vars.items()}
    print(f'\nStatus: {pulp.LpStatus[lpProblem.status]}\nObjective: {lpProblem.objective}\nFluxes:')
    
    for rct in fluxes:
        if fluxes[rct].varValue !=0.0:
            print(f'\t{rct}: {fluxes[rct].varValue}')
        

    count1=0
    for rid in rewardedReactions:
        if fluxes[rid].varValue!=0.0:
            count1+=1

    print(f'From {len(rewardedReactions)} highly active reactions, {count1} are active.')
    count2=0
    for rid in penalizedReactions:
        if fluxes[rid].varValue==0.0:
            count2+=1
    print(f'From {len(penalizedReactions)} lowly active reactions, {count2} reactions are inactive.') 


def create_outputFiles(lpProblem, path, dataframe, disc_method, expressionColName, penalizedRcts, rewardedRcts, moderateRcts, uT, lT, epsilon, v_vars, y_vars, c_vars =None, triple=None, o2Constraint= None ):
    '''
    Creates 3 different output files from Linear Programming Problem:
    - expressionData_classification.tsv with gene id, gene expression, scaled expression, and discretization value
    - reactionData_classification.tsv with reaction ids of reaction in rH or rL, y_vars, flux value, and c_vars
    - parameter_summary.txt with Parameter summary

    Args
    -------
    lpProblem: pulp.LpProblem
    path: str
        path to output directory
    dataframe: pandas.DataFrame
        with gene id, gene expression data, discretization values
    disc_method: str
        discretization method chosen
    expressionColName: str
        Name of the column with expression values
    penalizedRct, rewardedRct: list
        List with reaction id of penalized/rewarded reactions
    uT, lT: float
        upper and lower thresholds from discretization
    epsilon: float
        epsilon value chosen for iMAT 
    v_vars, y_vars, c_vars: dict
        variables from Linear Programming Problem
    triple: bool
        if true, triple-weighted-iMAT was used
    '''
    status = lpProblem.solve(pulp.PULP_CBC_CMD(msg=True))
    if pulp.LpStatus[lpProblem.status] == 'Infeasible':
        print(f"⚠️ Problem {lpProblem.name} for {expressionColName} is infeasible — skipping.")
        if o2Constraint:
            print(f'o2 Constraint given" {o2Constraint}')
        return None  # exit early


    cwd = Path.cwd()
    output_path = cwd / path  
    
    if not os.path.exists(output_path):
        print(f' Output path {output_path} does not exist. Please verify.')
    else:
        # expressionData_clasification
        output_file = Path(output_path) / 'expressionData_classification.tsv'
        output_file.parent.mkdir(parents=True, exist_ok=True)
        dataframe.to_csv(output_file, index=True, sep='\t')
        print(f'Expression data file has been saved to {output_file}')

        # reactionData_classification
        output_file = Path(output_path) / 'reactionData_classification.tsv'
        output_file.parent.mkdir(parents=True, exist_ok=True)
        rows = []  # collect rows here
        if c_vars is None:
            #regular iMAT
            for rid in v_vars.keys():
                rid.replace('v_','')
                if rid in y_vars.keys():
                    ytot, y_f, y_r = y_vars[rid]
                    classification = 'high' if rid in (rewardedRcts) else 'low'
                    rows.append({'reaction_id': rid, 
                                 'flux':v_vars[rid].value(), 
                                 'y_f':y_f.value(), 
                                 'y_r': y_r.value(),
                                 'classification':classification})
                elif rid in moderateRcts:
                    rows.append({'reaction_id': rid, 
                                 'flux':v_vars[rid].value(), 
                                 'classification':'moderate'})
                else:
                    rows.append({'reaction_id': rid, 
                                 'flux':v_vars[rid].value()})

        else:
            #weighted iMATs
            for rid in v_vars.keys():
                rid.replace('v_','')
                if rid in y_vars.keys():
                    ytot, y_f, y_r = y_vars[rid]
                    classification = 'high' if rid in rewardedRcts else 'moderate' if rid in moderateRcts else 'low'
                    if (rid in rewardedRcts) or (rid in moderateRcts and triple) or (rid in penalizedRcts and not triple):
                        rows.append({'reaction_id': rid, 
                                     'flux':v_vars[rid].value(), 
                                     'y_f':y_f.value(), 'y_r': y_r.value(), 
                                     'classification':classification, 
                                     'c_value':c_vars[rid]})
                    elif rid in penalizedRcts and triple:
                        rows.append({'reaction_id': rid, 
                                     'flux':v_vars[rid].value(), 
                                     'y_f':y_f.value(), 
                                     'y_r': y_r.value(), 
                                     'classification':classification})
                    elif rid in moderateRcts and not triple:
                        rows.append({'reaction_id': rid, 
                                     'flux':v_vars[rid].value(),
                                     'classification':classification})             
                #Reaction is classified but no y
                elif rid in (rewardedRcts+moderateRcts+penalizedRcts):
                    classification = 'high' if rid in rewardedRcts else 'moderate' if rid in moderateRcts else 'low'
                    rows.append({'reaction_id': rid, 
                                 'flux':v_vars[rid].value(),
                                 'classification':classification})             
                #Reactions is not classified 
                else:
                    rows.append({'reaction_id': rid, 'flux':v_vars[rid].value()})
        reactionData_df = pd.DataFrame(rows)      
        reactionData_df.to_csv(output_file, index=False, sep='\t')
        print(f'Reaction data file has been saved to {output_file}')

        # Save .txt file with all parameters
        output_file = Path(output_path) / 'Parameters_classification.txt'
        with open(output_file,'w') as f:
            if c_vars is not None:
                f.write('==== weighted iMAT Parameters ====\n\n')
            else:
                f.write('==== iMAT Parameters ====\n\n')

            if disc_method =='mean':
                f.write('Mean Discretization\n')
                f.write(f'Mean: {dataframe[expressionColName].mean()}\n' 
                        f'Standard Deviation: {dataframe[expressionColName].std()}\n'
                        f'Upper Threshold: {uT}\n'
                        f'Lower Threshold: {lT}\n\n')
                f.write('Parameters:\n'
                        f'Epsilon: {epsilon}\n\n')
            else:
                f.write('Quantile Discretization\n')
                f.write(f'Mean: {dataframe[expressionColName].mean()}\n' 
                        f'Median: {np.nanquantile(dataframe[expressionColName],0.5)}\n'
                        f'Upper Threshold: {uT}\n'
                        f'Lower Threshold: {lT}\n\n')
                f.write('Parameters:\n'
                        f'Epsilon: {epsilon}\n\n')
            f.write('Linear Programming Problem information:\n')
            f.write(f'Solver status:\t{LpStatus[status]}\n'
                    f'Objective Value:\t{pulp.value(lpProblem.objective)}')
        print(f"Text file saved at: {output_file.resolve()}")