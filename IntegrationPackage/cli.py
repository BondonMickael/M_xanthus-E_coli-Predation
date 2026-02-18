import argparse


def build_parser():
    parser = argparse.ArgumentParser(description="Metabolic model integration Toolkit")
    subparser = parser.add_subparsers(
        dest="method",
        required=True,
        help="Chose integration method. Following are available: iMAT, weighted_iMAT",
    )

    # iMAT
    imat_parser = subparser.add_parser("iMAT", help="Run standard iMAT")
    add_shared_args(imat_parser)
    add_iMAT_shared_args(imat_parser)

    # weighted iMAT
    weighted_imat_parser = subparser.add_parser(
        "weighted_iMAT", help="Run weighted iMAT"
    )
    add_shared_args(weighted_imat_parser)
    add_iMAT_shared_args(weighted_imat_parser)

    return parser


def add_shared_args(p):
    p.add_argument(
        "-f",
        "--expressionFile",
        type=str,
        required=True,
        help="Path to the RNAseq file (.csv or .tsv)",
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
        "-m", "--model", required=True, type=str, help="Path to cobra.Model file"
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
        "-x",
        "--oxygenLevel",
        type=float,
        default=None,
        help="Value for oxygen exchange reaction EX_o2_e",
    )
