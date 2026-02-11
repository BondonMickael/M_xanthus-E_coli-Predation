#instead from command line, I want to run PulpIntegrationMetho here
from main import run_method 
from types import SimpleNamespace  # mimic argparse Namespace
import os

# iMAT without oxygen levels and IMSH model
exprColNames =['AP..cycling', 'AP..cycling.II', 'AP..cycling.III', 'BP', 'BP..cycling', 'DL.neurons', 'DL.neurons.II', 'UL..differentiating', 'UL..migrating']
celltype_names =['AP_cycling', 'AP_cycling_2', 'AP_cycling_3', 'BP', 'BP_cycling', 'DL_neurons', 'DL_neurons_2', 'UL_differentiating', 'UL_migrating']

geneColName = 'ensembl_gene_id'
exprFile = './Data/E14WT_normalized_MEAN_ensembl_transcriptomics.tsv'
mdl = './mitoMammal/Model_test_IMSH.sbml'

outputpath= ['regular_iMAT', 'weighted_iMAT', 'triple_weighted_iMAT']
for c, exprColName in enumerate(exprColNames):
    for i,method in enumerate(["iMAT", "weighted_iMAT", "triple_weighted_iMAT"]):
        output_dir = f"IntegrationOutputs/iMAT_outputs/woObjective/{outputpath[i]}/{celltype_names[c]}"

        args = SimpleNamespace(
            method=method,
            discretization="quantile",
            quantiles = None, #default quantiles will be used
            model = mdl,
            ignore_human = True, 
            epsilon = 1.0,
            optimum = 0.0,
            oxygenConstraint = None,
            expressionFile = exprFile,
            geneColName = geneColName,
            expressionColName = exprColName,
            output = output_dir
        )
        run_method(args)


# #With oxygen levels and IMSH model 
# exprColNames =['BP', 'BP..cycling', 'DL.neurons', 'DL.neurons.II', 'UL..differentiating', 'UL..migrating']
# celltype_names =['BP', 'BP_cycling', 'DL_neurons', 'DL_neurons_2', 'UL_differentiating', 'UL_migrating']
# for c, exprColName in enumerate(exprColNames):
#     geneColName = 'ensembl_gene_id'
#     exprFile = './Data/E14WT_normalized_MEAN_ensembl_transcriptomics.tsv'
#     mdl = './mitoMammal/Model_test_IMSH_glucoseImport.sbml'

#     oxygen_levels = [0.0, -1.0, -2.0, -3.0, -4.0, -5.0, -6.0, -7.0, -8.0, -9.0, -10.0, -11.0, -12.0, -13.0, 
#                     -14.0, -15.0, -16.0, -17.0, -18.0, -19.0]
#     methods = ["iMAT", "weighted_iMAT", "triple_weighted_iMAT"]
#     outputpath= ['regular_iMAT', 'weighted_iMAT', 'triple_weighted_iMAT']

#     for i, method in enumerate(methods):
#         for j, level in enumerate(oxygen_levels):
#             folder_name = f"{abs(int(level)):02d}"
#             output_dir = f"IntegrationOutputs/iMAT_outputs/woObjective/{outputpath[i]}/{celltype_names[c]}/oxygenConstraint/{folder_name}"

#             #Create folder if it doesn't exist
#             os.makedirs(output_dir, exist_ok=True)

#             args = SimpleNamespace(
#                 method=method,
#                 discretization="quantile",
#                 quantiles=None, 
#                 model=mdl,
#                 ignore_human=True,
#                 epsilon=1.0,
#                 optimum=0.0,
#                 oxygenConstraint=oxygen_levels[j],
#                 expressionFile=exprFile,
#                 geneColName=geneColName,
#                 expressionColName=exprColName,
#                 output=output_dir
#             )

#             # Run the method
#             run_method(args)

