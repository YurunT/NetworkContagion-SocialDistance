import sys, os
from multiprocessing import Manager
from joblib import Parallel, delayed
from datetime import datetime
sys.path.append(os.path.abspath("../auxiliary_scripts/"))
sys.path.append(os.path.abspath("../theories"))
from mask import *
from mutation import *
from input_module import *
from output_module import write_analysis_results
from main_aux import *
    
def main():
    ###### Get CMD input ######
    paras = parse_args(sys.argv[1:])
    paras_check(paras)
    mean_degree_list = get_mean_degree_list(paras)
    
    print('-------Parameter Setting-------\n', vars(paras))
    print("mean_degree_list:", mean_degree_list)
    print("num_cores:", paras.nc)
    print('-------Parameter Setting-------\n')
    
    ###### Run on multiple cores using parellel ###### 
    infection_size = Manager().dict()
    infection_size0 = Manager().dict()
    infection_size1 = Manager().dict()
    

    if paras.modelname   == 'mask'     and paras.itemname == 'es':
        Parallel(n_jobs  = paras.nc)(delayed(get_EpidemicSize)(mean_degree, paras, infection_size0, infection_size1, infection_size) for mean_degree in mean_degree_list)
        
    elif paras.modelname == 'mutation' and paras.itemname == 'es':
        Parallel(n_jobs  = paras.nc)(delayed(cascade_size)(mean_degree, paras, infection_size0, infection_size1, infection_size) for mean_degree in mean_degree_list)
        
    elif paras.modelname == 'mask'     and paras.itemname == 'pe':
        Parallel(n_jobs  = paras.nc)(delayed(get_ProbEmergence)(mean_degree, paras, infection_size0, infection_size1, infection_size) for mean_degree in mean_degree_list)
        
    elif paras.modelname == 'mutation' and paras.itemname == 'pe':
        Parallel(n_jobs  = paras.nc)(delayed(cascade_prob)(mean_degree, paras, infection_size0, infection_size1, infection_size) for mean_degree in mean_degree_list)

    ######### Save the results for all Mean Degrees ########## 
    write_analysis_results(paras, [infection_size0, infection_size1, infection_size], mean_degree_list)
    print("All done!")
main()