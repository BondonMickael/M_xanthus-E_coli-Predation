def check_discretization(dis, quant=None):
    '''
    Checks if discretization method and variables can be used.

    :param dis: str: Discretization method chosen
    :param quant: List of Int quantile
    '''
    if dis == 'mean' and quant is not None:
        raise ValueError(f'You do not need to specify quantiles for mean discretization.')
    elif dis =='quantile' and quant is None:
        print(f'Quantile discretization with default quantiles {[40, 70]}.')
        quant =[40, 70]
    elif dis == 'quantile' and quant is not None:
        #check quantiles 
        q1, q2 = quant
        if q1 >=q2:
            raise ValueError(f'qL needs to be smaller than qH. Your input [{q1,q2}.')
        if q2 >99:
            raise ValueError(f'qH < 100. Your input: qH = {q2}')
