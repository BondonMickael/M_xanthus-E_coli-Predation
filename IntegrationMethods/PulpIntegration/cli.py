"""
iMAT command line interface
AIM: different iMAT integration methods:
- standard iMAT (as in original publication or only penalize lowly expressed reactions)
- weighted iMAT (weight for all reactions based on expression value or only for low expression)

different discretization methods:
- mean based (mean +/- 0.5*sd)
- quantile based (user defined quantiles)
-
"""

import argparse


def get_parser():
    parser = argparse.ArgumentParser(description="Metabolic model integration toolkit")
    subparsers = parser.add_subparsers(
        dest="method",
        required=True,
        help="Integration method to run. Available integration methods: FBA, FVA, irreversible_iMAt_FVA, full_iMAT_FVA, iMAT, weighted_iMAT (single or double), triple_weighted_iMAT, GIMME",
    )

    # Shared arguments (except FBA)
    def add_shared_args(p):
        p.add_argument(
            "-f",
            "--expressionFile",
            type=str,
            required=True,
            help="Path tothe RNAseq file (.csv or .tsv)",
        )
        p.add_argument(
            "-g",
            "--geneColName",
            required=True,
            type=str,
            help="Name of the column with the gene id",
        )
        p.add_argument(
            "-i",
            "--expressionColName",
            required=True,
            type=str,
            help="Column name containing expression values",
        )
        p.add_argument(
            "-m",
            "--model",
            required=True,
            type=str,
            help="Path to metabolic model file",
        )
        p.add_argument(
            "-a",
            "--ignore_human",
            action="store_true",
            default=True,
            help="Ignore human gene IDs (default: True), else ignores mice genes.",
        )
        p.add_argument(
            "-o",
            "--output",
            default=None,
            type=str,
            help="Output directory path. If given it will produce following files:\n"
            "- .tsv file with the fluxes and penalty/reward\n"
            "- .txt file with the summary of the parameter",
        )

    def add_iMAT_shared_args(p):
        p.add_argument(
            "-e",
            "--epsilon",
            type=float,
            default=1.0,
            help="Threshold value to consider a reaction as active, default is 1.",
        )
        p.add_argument(
            "-d",
            "--discretization",
            choices=["mean", "quantile"],
            type=str,
            default="mean",
            help="Discretization method (default: mean)",
        )
        p.add_argument(
            "-q",
            "--quantiles",
            type=int,
            nargs=2,
            default=None,
            metavar=("qL", "qH"),
            help=(
                "Quantile thresholds for cutoff of low gene expression (qL) and high gene expression (qH). qL must be lower than qH. Default values are qL=40, qH=70."
            ),
        )
        p.add_argument(
            "-p",
            "--optimum",
            type=float,
            default=80.0,
            help="Minimum optimum percentage of objective reaction (default: 80)",
        )
        p.add_argument(
            "-x",
            "--oxygenConstraint",
            type=float,
            default=None,
            help="Flux value for oxygen exchange reaction EX_o2_e",
        )

    # weighted_iMAT
    weighted_imat_parser = subparsers.add_parser(
        "weighted_iMAT",
        help="Run weighted_iMAT integration. Highly acive reactions are further "
        "rewarded with weight values according their respective gene expression value.",
    )
    add_shared_args(weighted_imat_parser)
    add_iMAT_shared_args(weighted_imat_parser)

    # triple-weighted iMAT
    triple_weighted_imat_parser = subparsers.add_parser(
        "triple_weighted_iMAT", help="Run triple_weighted_iMAT integration."
    )
    add_shared_args(triple_weighted_imat_parser)
    add_iMAT_shared_args(triple_weighted_imat_parser)

    # regular iMAT
    imat_parser = subparsers.add_parser("iMAT", help="Run iMAT integration")
    add_shared_args(imat_parser)
    add_iMAT_shared_args(imat_parser)

    # irreversible iMAT_FVA
    irreversible_imat_fva_parser = subparsers.add_parser(
        "irreversible_iMAT_FVA",
        help="Run FVA integration using iMAT principles by constraining irreversible highly active reactions to carry flux f +/- epsilon",
    )
    add_shared_args(irreversible_imat_fva_parser)
    add_iMAT_shared_args(irreversible_imat_fva_parser)
    irreversible_imat_fva_parser.add_argument(
        "-z",
        "--fva_reactionOfInterest",
        required=True,
        type=str,
        help="Reaction ID of interest for FVA analysis",
    )

    # full iMAT_FVA
    full_imat_fva_parser = subparsers.add_parser(
        "full_iMAT_FVA",
        help="Run FVA integration using iMAT principles by constraining all highly active reactions to carry flux f +/- epsilon",
    )
    add_shared_args(full_imat_fva_parser)
    add_iMAT_shared_args(full_imat_fva_parser)
    full_imat_fva_parser.add_argument(
        "-z",
        "--fva_reactionOfInterest",
        required=True,
        type=str,
        help="Reaction ID of interest for FVA analysis",
    )

    # FVA
    fva_parser = subparsers.add_parser(
        "FVA", help="Run plain Flux Variability Analysis"
    )
    fva_parser.add_argument("--model", required=True)
    fva_parser.add_argument("-z", "--fva_reactionOfInterest", required=True)

    # GIMME
    gimme_parser = subparsers.add_parser("GIMME", help="Run GIMME integration")
    add_shared_args(gimme_parser)

    # FBA
    fba_parser = subparsers.add_parser(
        "FBA", help="Run plain Flux Balance Analysis (no integration)"
    )
    fba_parser.add_argument("--model", required=True)
    fba_parser.add_argument("--output", default="results/")

    return parser


def parse_args():
    parser = get_parser()
    args = parser.parse_args()
    return args

    # ADD OPTION TO CHOOSE LOWER QUANTILE


# python ./IntegrationMethods/PulpIntegration/main.py iMAT -g "ensembl_gene_id" -i "AP..cycling" -m ./mitoMammal/Model_test_IMSH_glucoseImport.sbml -f ./Data/E14WT_normalized_MEAN_ensembl_transcriptomics.tsv -o './iMAT_outputs/regular_iMAT_IMSH_model/AP_cycling' -d quantile -p 0.0
