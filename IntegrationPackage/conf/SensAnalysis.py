import numpy as np

# The `exprColNames` list in the code snippet contains the names of different expression columns in a
# dataset. These column names represent different types of gene expression data related to various
# cell types or biological processes. Each element in the `exprColNames` list corresponds to a
# specific type of gene expression data, such as 'AP..cycling', 'BP', 'DL.neurons', etc. These names
# are used to identify and extract specific expression data from the dataset for further analysis or
# processing.
exprColNames = ['AP..cycling', 'AP..cycling.II', 'AP..cycling.III','BP', 'BP..cycling', 'DL.neurons', 'DL.neurons.II', 'UL..differentiating', 'UL..migrating']
celltype_names = ['AP_cycling', 'AP_cycling_2', 'AP_cycling_3','BP', 'BP_cycling', 'DL_neurons', 'DL_neurons_2', 'UL_differentiating', 'UL_migrating']
methods = ["iMAT", 'triple_weighted_iMAT']
outputpaths = ['regular_iMAT', 'weighted_iMAT', 'triple_weighted_iMAT']
oxygen_levels = []
geneColName = 'ensembl_gene_id'
inputDir = ['.Data/SensitivityData']    #Directory to files with expression data
metabolicModel = ['./mitoMammal/Model_test_IMSH_glucoseImport.sbml']
discretization = ['quantile']
outputpaths = ['regular_iMAT', 'weighted_iMAT']
# Parameter Ranges
epsilon_range = np.logspace(-3, 1, 5)
lower_q_range = np.linspace(1, 75, 7)
upper_q_range = np.linspace(25, 99, 7)
oxygen_levels = None

# Grid Search

    print(f'Total parameter combinations: {len(param_grid)}')
    # Expression Files
    expr_folder = Path('./Data/SensitivityData')
    expr_files = list(expr_folder.glob("*.tsv"))

    # Generate all tasks
    params = []
    for exprFile in expr_files:
        for c, exprColName in enumerate(exprColNames):  #per celltype
            for i, method in enumerate(methods):
                for epsilon, lower_q, upper_q in param_grid:
                    params.append((method, outputpaths[i], exprColName, celltype_names[c], mdl, str(exprFile), geneColName, epsilon, lower_q, upper_q))

    # Parallelization
    num_processes = 16
    print(f"Total tasks to run: {len(params)}")

    with Pool(processes=num_processes) as pool:
        pool.map(run_single_case, params)

    print("Done!")
    