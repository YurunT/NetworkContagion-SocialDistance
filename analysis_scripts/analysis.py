import sys, os
import multiprocessing
from multiprocessing import Manager
from joblib import Parallel, delayed
from datetime import datetime
sys.path.append(os.path.abspath("../auxiliary_scripts/"))
sys.path.append(os.path.abspath("../theories"))
from mask import *
from mutation import *
from tnn import *
from input_module import *
from output_module import write_analysis_results

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
    
def main():
    ###### Get commandline input ######
    paras = parse_args_analysis(sys.argv[1:])
    paras_check(paras)
    num_cores, rho, k_max, T_list, Q_list, mu_list, mean_degree_list = resolve_paras(paras)

    ###### Run on multiple cores using parellel ###### 
    infection_size = Manager().dict()
    infection_size0 = Manager().dict()
    infection_size1 = Manager().dict()

    if paras.modelname == 'mask' and paras.itemname == 'es':
        Parallel(n_jobs = num_cores)(delayed(get_EpidemicSize)(mean_degree, k_max, T_list, paras, infection_size, infection_size0, infection_size1) for mean_degree in mean_degree_list)
    elif paras.modelname == 'mutation' and paras.itemname == 'es':
        Parallel(n_jobs = num_cores)(delayed(cascade_size)(mean_degree, Q_list, mu_list, k_max, rho, infection_size, infection_size0, infection_size1) for mean_degree in mean_degree_list)

    ######### Save the results for all Mean Degrees ########## 
    write_analysis_results(paras, [infection_size0, infection_size1, infection_size], paras.modelname, paras.itemname, mean_degree_list)
    print("All done!")
main()