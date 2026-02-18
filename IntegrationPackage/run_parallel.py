# NOT DONE IMPLEMENTING!

# Creates the tasks from file SensAnalysis
import conf.SensAnalysis as SensAnalysis


def run_parallel():
    model = read_sbml_model(SensAnalysis.metabolicModel)
    tasks = []  # list of dict

    # Create combination of epsilon and quantiles
    param_grid = [
        (SensAnalysis.epsilon, SensAnalysis.lower_q, SensAnalysis.upper_q)
        for epsilon, lower_q, upper_q in itertools.product(
            SensAnalysis.epsilon_range,
            SensAnalysis.lower_q_range,
            SensAnalysis.upper_q_range,
        )
        if lower_q < upper_q
    ]

    input_files = list(Path(SensAnalysis.inputDir).glob(".tsv"))
    for file in input_files:
        # create ouputpath name
        fileName = Path(file).stem
        expression_df = read_file.read_expression_file(file)
        for c, exprColName in enumerate(SensAnalysis.exprColNames):
            rnaSeq_df = utils.generate_RNASeqDf(
                model, expression_df, SensAnalysis.geneColName
            )
            for i, method in enumerate(SensAnalysis.methods):
                if SensAnalysis.level is not None:
                    for epsilon, lower_q, upper_q in param_grid:
                        tasks.append(
                            {
                                "model": model,
                                "method": method,
                                "rnaSeq_df": rnaSeq_df,
                                "discretization": SensAnalysis.discretization,
                                "quantiles": [lower_q, upper_q],
                                "epsilon": epsilon,
                                "outputpath": f"sensitivityAnalysis_output/{SensAnalysis.outputpaths[i]}/{SensAnalysis.celltype_names[c]}/{fileName}_epsilon_{epsilon}_quantiles_{lower_q}_{upper_q}.tsv",
                            }
                        )
    # Multiprocessing on tasks
    num_processes = 16
    with Pool(processes=num_processes) as pool:
        pool.map(run_single, tasks)
