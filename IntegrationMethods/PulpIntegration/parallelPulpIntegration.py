from multiprocessing import Pool
from main import run_method
from types import SimpleNamespace
from cobra.io import read_sbml_model
import os
import numpy as np


def run_single_case(params):
    (
        method,
        outputpath,
        exprColName,
        celltype_name,
        level,
        mdl,
        exprFile,
        geneColName,
    ) = params
    if level is None:
        output_dir = f"parallelIntegrationOutputs/iMAT_outputs/woObjective/{outputpath}/{celltype_name}/oxygenConstraint/None"
    else:
        folder_name = f"{abs(int(level)):02d}"
        output_dir = f"parallelIntegrationOutputs/iMAT_outputs/woObjective/{outputpath}/{celltype_name}/oxygenConstraint/{folder_name}"
    os.makedirs(output_dir, exist_ok=True)

    args = SimpleNamespace(
        method=method,
        discretization="quantile",
        quantiles=None,
        model=mdl,
        ignore_human=True,
        epsilon=0.8,
        optimum=0.0,
        oxygenConstraint=level,
        expressionFile=exprFile,
        geneColName=geneColName,
        expressionColName=exprColName,
        output=output_dir,
    )
    run_method(args)


if __name__ == "__main__":
    exprColNames = [
        "AP..cycling",
        "AP..cycling.II",
        "AP..cycling.III",
        "BP",
        "BP..cycling",
        "DL.neurons",
        "DL.neurons.II",
        "UL..differentiating",
        "UL..migrating",
    ]
    celltype_names = [
        "AP_cycling",
        "AP_cycling_2",
        "AP_cycling_3",
        "BP",
        "BP_cycling",
        "DL_neurons",
        "DL_neurons_2",
        "UL_differentiating",
        "UL_migrating",
    ]
    methods = ["iMat", "weighted_iMAT", "triple_weighted_iMAT"]
    outputpaths = ["regular_iMAT", "weighted_iMAT", "triple_weighted_iMAT"]
    oxygen_levels = [-19.5, -20.0]
    geneColName = "ensembl_gene_id"
    exprFile = "./Data/E14WT_normalized_MEAN_ensembl_transcriptomics.tsv"
    mdl = read_sbml_model("./mitoMammal/Model_test_IMSH_glucoseImport.sbml")

    params = []
    for c, exprColName in enumerate(exprColNames):
        for i, method in enumerate(methods):
            for level in oxygen_levels:
                params.append(
                    (
                        method,
                        outputpaths[i],
                        exprColName,
                        celltype_names[c],
                        level,
                        mdl,
                        exprFile,
                        geneColName,
                    )
                )

    num_processes = 16
    print(f"Starting parallel run using {num_processes} processes")
    print(f"Total tasks to run: {len(params)}")

    with Pool(processes=num_processes) as pool:
        pool.map(run_single_case, params)

    print("Done!")
