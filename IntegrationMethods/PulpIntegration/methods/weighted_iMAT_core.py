import numpy as np
import pandas as pd
import pulp
from cobra.io import read_sbml_model
import cobra
from cobra.core.gene import GPR
import methods.utils.iMAT_utils as iMAT_utils
import re
import methods.utils.output_utils as output_utils


def find_leading_expression(rct, gene_expression, ignore_human = False, debbug=False):
    '''
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
    '''
    
    gpr_rule = iMAT_utils.filter_gpr(rct.gpr, ignore_human=ignore_human)
    if gpr_rule is None:
        #This will not happen because reactions without gpr are not included in rH or rL
        raise ValueError(f'Reactions {rct.id} has no GPR rule and should therefore not be classified into rH or rL')
    else:
        gpr = gpr_rule.replace(" or ", " OR ").replace(" and ", " AND ")
        # Tokenize while keeping parentheses and operators
        tokens = re.findall(r'\(|\)|AND|OR|[^\s()]+', gpr)
        if debbug:
            print(f'tokens: {tokens}')
        #Exchange gene ids with expression values 
        parsed_tokens = [iMAT_utils.parse_token(tok,gene_expression) for tok in tokens]
    
        val = iMAT_utils.recursive_evaluation_iMAT(parsed_tokens, debbug=False)
    return val

def create_expression_variable_full(model, rH, rL, gene_expression,fraction_of_optimum, epsilon, threshold_low, threshold_high, o2_constraint=None, ignore_human = False):
    '''
    Using the gene expression values, GPR rules and reaction weights, this function creates individual 
    weight variables per reaction so that reactions with associated gene expression value. Reactions in the 
    lowly expressed class, will receive higher penalty, whereas reactions of 
    the highly expressed set will 

    Args
    --------
    model: cobra.Model
        a cobrapy model
    rH: list
        list of highly active reactions
    rL: list
        list of lowly active reactions
    gene_expression: pandas.DataFrame
        DataFrame with gene ids as index and a column 'scaled_expression' with scaled expression values between 0 and 1
        and a column 'discretization' with discretized values (-1,0,1)
    fraction_of_optimum: float
        Fraction of optimum to ask for objective reaction.
    epsilon: float
        Activation threshold for highly active reactions
    threshold_low, threshold_high: float
        thresholds for lowly and highly expressed genes to calculate penalty/reward
    o2_constraint: float
        if given, oxygen exchange reaction (EX_o2_e) will be set to that value

    Returns 
    -------
    pulp variables for weighted_iMAT formulation
    '''    
    objectiveValue = model.optimize().objective_value
    obj_reactions =[rct for rct in model.reactions if rct.objective_coefficient !=0.0]
    objective_ids = [rct.id for rct in obj_reactions]
    if len(objective_ids) ==1:
        objective_id = obj_reactions[0].id
        print(f'Objective reaction: {objective_id}')
    else:
        print(f'Multiple objective reactions: {objective_ids}')
    print(f'The optimal solution is {objectiveValue}, the model will ask a minimum of {objectiveValue *fraction_of_optimum/100} for the objective reaction.')

    prob = pulp.LpProblem('weighted_iMAT_full', pulp.LpMaximize)
    #v_min <= v <= v_max
    v_vars ={}
    if o2_constraint is not None:
        oxygen_constrain = True
    else:
        oxygen_constrain = False
    
    for rct in model.reactions:
        rid = rct.id
        lb = rct.lower_bound
        ub = rct.upper_bound
        #Flux constraints (all reactions)
        if oxygen_constrain and rid == 'EX_o2_e':
            v = pulp.LpVariable(f'v_{rid}', o2_constraint, o2_constraint)
        # CII only runs forward 
        elif rid == 'CII':
            v = pulp.LpVariable(f'v_{rid}',0.0,1000.0, cat='Continuous')    
        # Glucose transporter constraint
        elif rid == 'GLCt1r':
            v = pulp.LpVariable(f'v_{rid}',lb,5.0, cat='Continuous')
        # Creatine/ Phosphocreatine exchange
        elif rid =='r0942':
            v = pulp.LpVariable(f'v_{rid}',-0.05,1000.0, cat='Continuous')
        elif rid == 'r0942b_mitoMap':
            v = pulp.LpVariable(f'v_{rid}',-0.05,1000.0, cat='Continuous')

        elif lb == np.inf or ub == np.inf:  
            v = pulp.LpVariable(f'v_{rid}', cat='Continuous') 
        else:
            v = pulp.LpVariable(f'v_{rid}', lb, ub, cat='Continuous')
        v_vars[rid] = v
    # Constrain CKc and CK to run in the same direction to avoid loop
    y_creatine = pulp.LpVariable('y_CK_sign', cat = 'Binary')
    v_CK = v_vars['CK']
    v_CKc = v_vars['CKc']
    M_ck = 1000.0
    prob += v_CK   <=  M_ck * y_creatine, f"CK_pos_upper_if_y1"
    prob += v_CK   >= -M_ck * (1 - y_creatine), f"CK_pos_lower_if_y1"
    prob += v_CKc  <=  M_ck * y_creatine, f"CKc_pos_upper_if_y1"
    prob += v_CKc  >= -M_ck * (1 - y_creatine), f"CKc_pos_lower_if_y1"
        
    # constrain objective reaction 
    prob += v_vars[objective_id]>= objectiveValue * fraction_of_optimum/100
    #Mass balance constraint
    for met in model.metabolites:
        constraint = pulp.lpSum(rct.metabolites.get(met,0) * v_vars[rct.id] for rct in met.reactions) == 0
        prob += constraint, f'mass_balance_{met.id}'

    #flux constraint on classified reactions and their weights for objective function
    y_vars = {}
    c_vars = {}
    # Binary activity variables
    for rid in (rH + rL):
        y_tot =pulp.LpVariable(f'ytot_{rid}', cat='Binary')
        y_f = pulp.LpVariable(f'yf_{rid}', cat = 'Binary')
        y_r = pulp.LpVariable(f'yr_{rid}', cat='Binary')
        prob += y_tot == y_f + y_r, f'ytot_def_{rid}'
        y_vars[rid] = (y_tot, y_f, y_r)
        
        v  = v_vars[rid]
        lb = model.reactions.get_by_id(rid).lower_bound
        ub = model.reactions.get_by_id(rid).upper_bound
        if rid in rH:
            #Highly expressed reactions
            prob += v + y_f*(lb-epsilon) >= lb, f'rH_forward_{rid}'
            prob += v + y_r*(ub + epsilon) <= ub, f'rH_reverse_{rid}'
            #Find the associated relevant gene for the reaction
            x = find_leading_expression(model.reactions.get_by_id(rid),gene_expression, ignore_human = ignore_human)
            #print(f'rH: Leading expression for reaction {rid} is {x}')
            c = (x - threshold_high)+1
            c_vars[rid] = c
        else: 
            prob += lb * (1 - y_f) <= v, f'rL_leftHandSide_{rid}'
            prob += v <= (1 - y_f) * ub, f'rL_rightHandSide_{rid}'
            x = find_leading_expression(model.reactions.get_by_id(rid),gene_expression, ignore_human = ignore_human)
            #print(f'rL: leading expression for reaction {rid} is {x}')
            c = x/threshold_low
            c_vars[rid] = c

    return prob, v_vars, y_vars, c_vars



def run_weightedIMat(args):
    #pd.DataFrame generation of expression data, discretization, 
    RNASeqData = iMAT_utils.generate_RNASeqDataFrame(args.model, args.expressionFile,args.geneColName, args.expressionColName, args.ignore_human)

    #Discetization on scaled expression data 
    
    [lower_threshold, upper_threshold],[lower_threshold_scaled, upper_threshold_scaled],RNASeqData = iMAT_utils.discretization(RNASeqData, args.expressionColName, args.discretization, args.quantiles)
   
    #map GPR to reaction
    penalizedReactions = []
    rewardedReactions = []
    moderateReactions = []
    active_genes = RNASeqData.index[RNASeqData['discretization'] == 1].tolist()
    nonActive_genes = RNASeqData.index[RNASeqData['discretization']==-1].tolist()
    moderate_genes = RNASeqData.index[RNASeqData['discretization']==0].tolist()

    
    for rct in args.model.reactions:
        #Remove mouse/human genes from gpr rules 
        gpr_rule = iMAT_utils.filter_gpr(rct.gpr, ignore_human=args.ignore_human)
        gpr_rule = GPR().from_string(gpr_rule)
        if not GPR(gpr_rule).eval(nonActive_genes):
            penalizedReactions.append(rct.id)
        elif not GPR(gpr_rule).eval(nonActive_genes + moderate_genes):
            moderateReactions.append(rct.id)
        elif not GPR(gpr_rule).eval(active_genes+moderate_genes+nonActive_genes):
            rewardedReactions.append(rct.id)
 
    print(f"There are {len(penalizedReactions)} penalized reaction")
    print(f"There are {len(rewardedReactions)} rewarded reaction")

    #Create PUlp Variables with individual weight (reward or penalty) for reactions that were classified
    #Either for the whole model
    prob, v_vars, y_vars, c_vars = create_expression_variable_full(args.model, rewardedReactions, penalizedReactions, RNASeqData, args.optimum ,args.epsilon,
                                                           lower_threshold_scaled, upper_threshold_scaled,o2_constraint= args.oxygenConstraint, ignore_human=args.ignore_human)
    single_weighted =True 

    
    #Constrain objective function
    terms =[]
    for rct in rewardedReactions:
        y_tot, y_f, y_r = y_vars[rct]
        c = c_vars[rct]
        terms.append((y_f+y_r)*c)
    for rct in penalizedReactions:
        y_tot, y_f, y_r = y_vars[rct]
        c = c_vars[rct]
        if single_weighted:
            terms.append(y_f)   #single weighted iMAT
        elif not single_weighted:
            terms.append((y_f+(1-y_f))*c) #double weighted iMAT

    atpObjective = False     #Add ATP hydrolysis to objective function
    if atpObjective:
        terms.append(v_vars['OF_ATP_mitoMap'])
    prob += pulp.lpSum(terms) 
                                                                                                                  
    prob.solve()

    #Collect Fluxes
    if args.output is None:
        output_utils.print_output(prob, v_vars, rewardedReactions, penalizedReactions)
        
    #Output files
    elif args.output:
        output_utils.create_outputFiles(prob,args.output, RNASeqData, args.discretization, args.expressionColName, penalizedReactions, rewardedReactions,
                                        moderateReactions, upper_threshold, lower_threshold, args.epsilon, v_vars, y_vars, c_vars, o2Constraint= args.oxygenConstraint)
        




