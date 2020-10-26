def get_effeciency_list(paras):
    return [[paras.tm1, paras.tm2],
            [1, 1]]
    

def generate_new_transmissibilities_mask(effeciency_list, T,):
    '''
    Input:  
            infectious - T ->|- T * inward_e -> |- T * inward_e * out_e -> susceptible
            
            effeciency_list: List of (inward effeciency(notated as 1), outward effeciency(notated as 2))
            T: Orignial transmissibility of the virus
    Output: 
            List of Tij
    '''
#     T1 = T * T_mask1 #T12
#     T2 = T * T_mask1 * T_mask2 #T11
#     T3 = T #T22
#     T4 = T * T_mask2 #T21

#     trans_dict = {'T1': T1,
#                   'T2': T2,
#                   'T3': T3,
#                   'T4': T4}
#     print(trans_dict)

#     effeciency_list = [[T_mask1, T_mask2],
#                        [1, 1]]
    typen = len(effeciency_list)
    T_list = []
    for u in range(typen):
        u_row = []
        for v in range(typen):
            Tij = effeciency_list[u][0] * effeciency_list[v][1] * T
            u_row.append(Tij)
        T_list.append(u_row)
    trans_dict = {'T1': T_list[0][1],
              'T2': T_list[0][0],
              'T3': T_list[1][1],
              'T4': T_list[1][0]}    
    print(trans_dict)
    return trans_dict    

def generate_new_transmissibilities_mutation(effeciency_list, T, m):
    trans_dict = generate_new_transmissibilities_mask(effeciency_list, T,)
    T1 = trans_dict['T1']
    T2 = trans_dict['T2']
    T3 = trans_dict['T3']
    T4 = trans_dict['T4']

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