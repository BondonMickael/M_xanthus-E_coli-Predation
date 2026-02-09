import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import cobra
from cobra.core.configuration import Configuration
from cobra.io import read_sbml_model, write_sbml_model
from cobra.flux_analysis import flux_variability_analysis
from tqdm import tqdm

from imatpy.model_utils import read_model
from imatpy.parse_gpr import gene_to_rxn_weights
from imatpy.imat import imat

# Read in the model
M_xanthus = read_model("../M_xanthus_model.sbml")

# Set the solver to glpk
Configuration().solver = "glpk"

# Read in the table
Table = pd.read_csv("../data/iMat_modified/WT_vs_Ecol_ratio1_3.csv")

# Create a pandas Series representing gene expression weights
iMatDic = {}

for i in range(len(Table.index)):  # take the genes from the experimental data
    iMatDic[Table["orgdb_old_MXAN"][i]] = int(Table["iMat"][i])

for i in M_xanthus.genes._dict:  # take the genes from the model
    if i not in iMatDic:
        iMatDic[i] = 0

iMatWeights = pd.Series(iMatDic)

# Convert the gene weights into reaction weights
reaction_weights = gene_to_rxn_weights(
    model=M_xanthus, gene_weights=iMatWeights, fill_val=0
)

# Run iMAT
imat_results = imat(
    model=M_xanthus, rxn_weights=reaction_weights, epsilon=1, threshold=0.01
)

# Print the imat objective
print(f"iMAT Objective: {imat_results.objective_value}")

# Print the imat flux distribution
print(f"iMAT Flux Distribution: \n{imat_results.fluxes}")
