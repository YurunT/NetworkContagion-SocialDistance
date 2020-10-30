import multiprocessing
import sys, os
import numpy as np
sys.path.append(os.path.abspath("."))
from tnn import *
from global_vars import *

def get_mean_degree_list(paras):
    if paras.mdl is not None: 
        mean_degree_list = paras.mdl
    else:
        mean_degree_list = np.linspace(paras.mind, paras.maxd, paras.ns)
    return mean_degree_list

def paras_check(paras):
    if paras.tm1 is None:
        print("------------ CMD INPUT INVALID ------------")
        print("Error: NO tm1 specified!")
        print("Please specify list of inward effeciencies!")
        print("------------ CMD INPUT INVALID ------------")
        assert False
    if paras.tm2 is None:
        print("------------ CMD INPUT INVALID ------------")
        print("Error: NO tm2 specified!")
        print("Please specify list of outward effeciencies!")
        print("------------ CMD INPUT INVALID ------------")
        assert False
    if len(paras.tm2) != len(paras.tm1):
        print("------------ CMD INPUT INVALID ------------")
        print("Error: len(tm1) %d not matching len(tm2) %d!" %(len(paras.tm1), len(paras.tm2)))
        print("Please specify correct num of in/out-ward effeciencies!")
        print("------------ CMD INPUT INVALID ------------")
        assert False
#     if len(paras.m) + 1 != len(paras.tm2):
#         print("------------ CMD INPUT INVALID ------------")
#         print("Error: len(m) should equals len(tm) - 1!")
#         print("Please check m list!")
#         print("------------ CMD INPUT INVALID ------------")
#         assert False
    if paras.modelname not in model_names:
        print("------------ CMD INPUT INVALID ------------")
        print("Error: NO model named %s!" %paras.modelname)
        print("Please select from", model_names)
        print("------------ CMD INPUT INVALID ------------")
        assert False
    if paras.itemname not in item_names:
        print("------------ CMD INPUT INVALID ------------")
        print("Error: NO analyzed item named %s!" %paras.itemname)
        print("Please select from", item_names)
        print("------------ CMD INPUT INVALID ------------")
        assert False
    if paras.change not in change_metrics_dict.values():
        print("------------ CMD INPUT INVALID ------------")
        print("Error: NO change metric as %s!" %paras.change)
        print("Please select numbers from", change_metrics_dict)
        print("------------ CMD INPUT INVALID ------------")
        assert False
        
           
def resolve_paras(paras):
    num_cores = min(paras.nc,multiprocessing.cpu_count())
    rho = 1.0 / paras.n
#     k_max = 4 * paras.maxd
    k_max = 20
    tmask_list = get_tmask_list(paras)
    T_list = generate_new_transmissibilities_mask(tmask_list, paras.T,)
    
    if paras.modelname == 'mutation':
        q_dict, mu_dict = generate_new_transmissibilities_mutation(tmask_list, paras.T, paras.m)
        Q_list  = list(q_dict.values())
        mu_list = list(mu_dict.values())    
    if paras.modelname == 'mutation':
        return rho, k_max, T_list, Q_list, mu_list
    else:
        return k_max, T_list