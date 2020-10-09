import multiprocessing
import sys, os
sys.path.append(os.path.abspath("."))
from tnn import *

def get_mean_degree_list(paras):
    if paras.mdl is not None: 
        mean_degree_list = paras.mdl
    else:
        mean_degree_list = np.linspace(paras.mind, paras.maxd, paras.ns)
    return mean_degree_list

def paras_check(paras):
    model_names = {'mask', 'mutation'}
    item_names = {'es', 'pe'}
    if paras.modelname not in model_names:
        print("Error: NO model named %s!" %paras.modelname)
        assert False
    if paras.itemname not in item_names:
        print("Error: NO analyzed item named %s!" %paras.itemname)
        assert False
           
def resolve_paras(paras):
    mean_degree_list = get_mean_degree_list(paras)
    num_cores = min(paras.nc,multiprocessing.cpu_count())
    rho = 1.0 / paras.n
    k_max = 4 * paras.maxd
    q_dict, mu_dict = generate_new_transmissibilities_mutation(paras.tm1, paras.tm2, paras.T, paras.m)
    T_list  = list(generate_new_transmissibilities_mask(paras.tm1, paras.tm2, paras.T, paras.m).values())
    Q_list  = list(q_dict.values())
    mu_list = list(mu_dict.values())
    print('-------Parameter Setting-------\n', vars(paras))
    print("K_max: ", k_max)
    print("num_cores:", num_cores)
    print("mean_degree_list:", mean_degree_list)
    print('-------Parameter Setting-------\n')
    
    return num_cores, rho, k_max, T_list, Q_list, mu_list, mean_degree_list