import argparse 
import cli 
from pathlib import Path
from methods.IMATConfig import IMATConfig
from methods.iMAT import iMAT
from methods.weighted_iMAT import weighted_iMAT
from utils.CreateOutput import CreateOutput
from cobra.io import read_sbml_model
from utils.generate_RNASeqDf import generate_RNASeqDf
from utils.read_file import read_expression_file

def run_single(args):
        model = read_sbml_model(args.model)
        df = read_expression_file(args.expressionFile)
        expression_df = generate_RNASeqDf(model, df, args.geneColName, args.expressionColName)
        config = IMATConfig(
                expression_df = expression_df,
                discretization_method = args.discretization, 
                quantiles = args.quantiles, 
                epsilon = args.epsilon,
                metabolicModel = model
                )
        config.prepare()
        
        if args.method == 'iMAT':
            solver = iMAT(
                metabolicModel = config.metabolicModel,
                RH = config.RH,
                RM = config.RM,
                RL = config.RL,
                epsilon = config.epsilon,
                oxygenLevel = args.oxygenLevel
            )
        
        elif args.method == 'weighted_iMAT':
            solver = weighted_iMAT(
                metabolicModel = config.metabolicModel,
                RH = config.RH,
                RM = config.RM,
                RL = config.RL,
                epsilon = config.epsilon,
                oxygenLevel = args.oxygenLevel,
                gpr_mapper = config.gpr_mapper,
                lower_threshold_scaled = config.lower_threshold_scaled,
                upper_threshold_scaled = config.upper_threshold_scaled
            )
            
        solver.build_problem()
        try:
            status, fluxes, sol_y_values, sol_c_values = solver.solve()
            genOutput = CreateOutput(
                output_dir = args.output, 
                method= args.method, 
                prob= solver.prob, 
                RH = config.RH,
                RM = config.RM,
                RL = config.RL,
                flux_distribution= fluxes, 
                y_values= sol_y_values, 
                c_values= sol_c_values, 
                epsilon=config.epsilon, 
                oxygenLevel= solver.oxygenLevel,
                expression_df=config.expression_df, 
                discretization_method=config.discretization_method, 
                quantiles= config.quantiles
            )
            genOutput.create_output()
        except InterruptedError:
            #infeasible Problem
            genOutput = CreateOutput(
                output_dir = args.output, 
                method= config.discretization_method, 
                prob= solver.prob, 
                RH = config.RH,
                RM = config.RM,
                RL = config.RL,
                flux_distribution= {}, 
                y_values= {}, 
                c_values= None, 
                epsilon=config.epsilon, 
                oxygenLevel= args.oxygenLevel,
                expression_df=config.expression_df, 
                discretization_method=args.discretization, 
                quantiles= config.quantiles
            ) 
            genOutput.handle_infeasibility()


def main():
    parallelization = False #if True, reads from SensAnalysis.py
    if parallelization == False:
        # Terminal input 
        args = cli.build_parser().parse_args()
        run_single(args)

    else: 
        raise ValueError('NOT IMPLEMENTED!')
            
            
            
if __name__ == '__main__':
    main()
    
    
# python ./IntegrationPackage/main.py iMAT -m ./DataIntegration/mitoMammal/6_universal_mito_model.xml -f ./DataIntegration/Data/E14WT_normalized_MEAN_ensembl_transcriptomics.tsv -o ./Output -g ensembl_gene_id -i AP..cycle