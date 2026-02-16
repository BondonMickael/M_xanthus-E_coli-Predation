from cli import parse_args
from cobra.io import read_sbml_model


def check_discretization(args):
    if args.discretization == "mean" and args.quantiles is not None:
        raise ValueError(
            f"You do not need to specify quantiles for mean discretization."
        )
    elif args.discretization == "quantile" and args.quantiles is None:
        print(f"Quantile discretization with default quantiles {[40, 70]}.")
        args.quantiles = [40, 70]
    elif args.discretization == "quantile" and args.quantiles is not None:
        # check quantiles
        q1, q2 = args.quantiles
        if q1 >= q2:
            raise ValueError(f"qL needs to be smaller than qH. Your input [{q1,q2}.")
        if q2 > 99:
            raise ValueError(f"qH < 100. Your input: qH = {q2}")


def run_method(args):
    if args.method == "FBA":
        import methods.FBA_core as FBA_core

        FBA_core.run_classical_FBA(args, only_FBA=True)

    elif args.method == "FVA":
        print("Classical FVA is not yet implemented, Sorry!")

    elif args.method == "iMAT":
        import methods.iMAT_core as iMAT_core

        check_discretization(args)
        args.model = read_sbml_model(args.model)
        iMAT_core.run_iMAT(args)

    elif args.method == "weighted_iMAT":
        import methods.weighted_iMAT_core as weighted_iMAT_core

        check_discretization(args)
        weighted_iMAT_core.run_weightedIMat(args)

    elif args.method == "triple_weighted_iMAT":
        import methods.triple_weighted_iMAT_core as triple_weighted_iMAT_core

        check_discretization(args)
        triple_weighted_iMAT_core.run_tripleWeighted_iMAT(args)

    elif args.method == "irreversible_iMAT_FVA":
        import methods.irreversible_iMATlike_FVA_core as irreversible_iMATlike_FVA_core

        check_discretization(args)
        irreversible_iMATlike_FVA_core.run_irreversible_iMAT_FVA(args)

    elif args.method == "full_iMAT_FVA":
        import methods.full_iMATlike_FVA_core as full_iMATlike_FVA_core

        check_discretization(args)
        full_iMATlike_FVA_core.run_full_iMAT_FVA(args)

    elif args.method == "GIMME":
        print("GIMME is still in process")


def main():
    args = parse_args()
    run_method(args)


if __name__ == "__main__":
    main()
