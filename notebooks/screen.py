import cobra
from cobra.io import (
    load_json_model,
    save_json_model,
    load_matlab_model,
    save_matlab_model,
    read_sbml_model,
    write_sbml_model,
)
from cobra.flux_analysis import flux_variability_analysis

M_xanthus = read_sbml_model("../M_xanthus_model.sbml")

exchange_id = []
for i in M_xanthus.exchanges:
    for j in i.metabolites:
        if "C" in j.formula:
            exchange_id.append(i.id)

Objective_value = {}
x = 0
for i in M_xanthus.reactions:
    if i.id in exchange_id:
        i.lower_bound = 0
        FBA = M_xanthus.optimize()
        Objective_value[i.id] = FBA.objective_value
        i.lower_bound = -1000
    
    for k in M_xanthus.reactions:
        if k.id in exchange_id:
            j.lower_bound = 0
            k.lower_bound = 0
            FBA = M_xanthus.optimize()
            Objective_value[str(j.id)+ ' + ' + str(k.id)] = FBA.objective_value
            j.lower_bound = -1000
            k.lower_bound = -1000
    x += 1
    print(x)

Objective_value