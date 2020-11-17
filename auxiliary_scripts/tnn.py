def get_tmask_list(paras):
    '''
    Input:  paras.tm1:[tmi,1], paras.tm2:[tmi,2]
    Output: tmask_list: [[tm1,1, tm1,2],
                         [tm2,1, tm2,2],
                         ...,
                         [tmM,1, tmM,2]]
    '''
    num_masks = len(paras.tm1)
    tmask_list = []
    for idx in range(num_masks):
        inward  = paras.tm1[idx]
        outward = paras.tm2[idx]
        mask_effeciency = [inward, outward]
        tmask_list.append(mask_effeciency)
    return tmask_list
    

def generate_new_transmissibilities_mask(tmask_list, T,):
    '''
    Input:  
            infectious - T ->|- T * outward_e -> |- T * outward_e * inward_e -> susceptible
            ----- Tij = Tmaski,1 * Tmaskj,2 * T -----
            
            effeciency_list: List of (outward effeciency(notated as 1), inward effeciency(notated as 2))
            T: Orignial transmissibility of the virus
    Output: 
            List of Tij(T_list):[[T11, T12, ..., T1M],
                                 [T21, T22, ..., T2M],
                                 ...,
                                 [TM1, TM2, ..., TMM]]
    '''
    typen = len(tmask_list)
    T_list = []
    for u in range(typen):
        u_row = []
        for v in range(typen):
            Tij = tmask_list[u][0] * tmask_list[v][1] * T
            u_row.append(Tij)
        T_list.append(u_row)
    return T_list    

def generate_new_transmissibilities_mutation(tmask_list, T, m):
    typen = len(tmask_list)
    if typen != 2:
        print("Currrently mutation only support 2 types of masks: no-mask and mask")
        assert False
    T_list = generate_new_transmissibilities_mask(tmask_list, T,)
    T1 = T_list[0][1]
    T2 = T_list[0][0]
    T3 = T_list[1][1]
    T4 = T_list[1][0]

    Q1 = T1 * (1 - m) + T2 * m
    Q2 = T3 * (1 - m) + T4 * m

    mu11 = T2 * m / Q1
    mu12 = T1 * (1 - m) / Q1
    mu22 = T3 * (1 - m) / Q2
    mu21 = T4 * m / Q2

    Q_dict = {
        "Q1": Q1,
        "Q2": Q2}
    
    mu_dict = {'mu11': mu11,
               'mu12': mu12,
               'mu22': mu22,
               'mu21': mu21,}

    return Q_dict, mu_dict