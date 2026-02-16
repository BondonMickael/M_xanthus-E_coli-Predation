def get_infeasible_reactions(reaction_list, model, epsilon):
    """
    From list of reactions, it checks whether forcing flux through respective
    reaction make the model infeasibles

    Args
    --------
    reaction_list: list
        List of reactions to test model feasibility
    model: cobra.Model
        a cobrapy Model
    epsilon:
        value used in iMAT

    Returns
    --------
    infeasible_rct: List
        list of reactions that result in infeasible model
    """
    infeasible_rct = []
    for rct in reaction_list:
        rct = model.reactions.get_by_id(rct)
        with model:
            if rct.lower_bound >= 0:
                rct.bounds = (epsilon, rct.upper_bound)
                sol = model.optimize()
                if sol.status == "infeasible":
                    infeasible_rct.append(rct.id)
            elif rct.upper_bound <= 0:
                rct.bounds = (rct.lower_bound, -epsilon)
                sol = model.optimize()
                if sol.status == "infeasible":
                    infeasible_rct.append(rct.id)
            else:
                # reversible reactions
                rct.bounds = (epsilon, 1000)
                sol = model.optimize()
                if sol.status == "infeasible":
                    rct.bounds = (-1000, -epsilon)
                    sol = model.optimize()
                    if sol.status == "infeasible":
                        infeasible_rct.append(rct.id)
    print(
        f"Out of {len(reaction_list)} reactions, {len(infeasible_rct)} reactions are infeasible when forcing flux of +/-{epsilon}.\nInfeasible reactions:"
    )
    for i in infeasible_rct:
        print(i)
    return infeasible_rct
