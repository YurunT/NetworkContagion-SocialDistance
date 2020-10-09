from __future__ import division
import argparse
import math
import sys, site, os
import numpy as np
from scipy.optimize import fsolve
from scipy.special import comb
from scipy.stats import poisson
import scipy.optimize
import scipy.misc
import multiprocessing
from multiprocessing import Manager
from joblib import Parallel, delayed
import json
from datetime import datetime
sys.path.append(os.path.abspath("../../auxiliary_scripts/"))
from tnn import *
from input_module import parse_args
from output_module import write_analysis_results


def obtain_val_r_1(v1, v2, t1, u_r_11, u_r_12, mean_degree, max_degree):
    val = 0

    for d_r in range(0, max_degree):
        prob_r = poisson.pmf(d_r, mean_degree)
        val += d_r*prob_r*1.0/mean_degree * ((1 - t1 + t1*u_r_11*v1 + t1*u_r_12*v2)**(d_r-1))

    return val

def obtain_val_r_2(v1, v2, t2, u_r_21, u_r_22, mean_degree, max_degree):
    val = 0

    for d_r in range(0, max_degree):
        prob_r = poisson.pmf(d_r, mean_degree)
        val += d_r*prob_r*1.0/mean_degree * ((1 - t2 + t2*u_r_22*v2 + t2*u_r_21*v1)**(d_r-1))

    return val

def equations(p, mean_degree, t_r_1, t_r_2, u_r_11, u_r_12, u_r_21, u_r_22, max_degree,):
    v1, v2 = p
    val_r_1 = obtain_val_r_1(v1, v2, t_r_1, u_r_11, u_r_12, mean_degree, max_degree)
    val_r_2 = obtain_val_r_2(v1, v2, t_r_2, u_r_21, u_r_22, mean_degree, max_degree)

    return (v1 - val_r_1, v2 - val_r_2)


def cascade_prob(mean_degree, t_r_1, t_r_2, u_r_11, u_r_12, u_r_21, u_r_22, max_degree, pe_0_list, pe_1_list, pe_list):
    h_r_1, h_r_2 = fsolve(equations, (0.01, 0.01), args=(mean_degree, t_r_1, t_r_2, u_r_11, u_r_12, u_r_21, u_r_22, max_degree,), xtol=1e-6)

    H_1 = 0
    for d_r in range(0, max_degree):
        prob_r = poisson.pmf(d_r, mean_degree)
        H_1 += prob_r*((1 - t_r_1 + t_r_1*u_r_11*h_r_1 + t_r_1*u_r_12*h_r_2)**d_r)

    H_2 = 0
    for d_r in range(0, max_degree):
        prob_r = poisson.pmf(d_r, mean_degree)
        H_2 += prob_r*((1 - t_r_2 + t_r_2*u_r_22*h_r_2 + t_r_2*u_r_21*h_r_1)**d_r)
        
    pe_0_list[mean_degree] = 1 - H_1
    pe_1_list[mean_degree] = 1 - H_2
    pe_list[mean_degree] = 2 - H_2 - H_1
    return (1 - H_1, 1 - H_2)


def main():
    ###### Paras setting ######
    paras = parse_args(sys.argv[1:])
#     numNodes = paras.n
#     rho = 1.0/numNodes

#     thrVal = paras.th
#     num_cores = min(paras.nc,multiprocessing.cpu_count())
#     # num_cores = multiprocessing.cpu_count()
#     mask_prob = paras.m
#     T_mask1 = paras.tm1
#     T_mask2 = paras.tm2
#     T = paras.T
#     degree_max = paras.maxd
#     degree_min = paras.mind
#     num_samples = paras.ns
#     change = paras.change

#     paras = dict()
#     paras['n'] = numNodes
#     paras['m'] = mask_prob
#     paras['T'] = T
#     paras['tm1'] = T_mask1
#     paras['tm2'] = T_mask2
#     paras['th'] = thrVal
#     paras['md'] = degree_max
#     paras['ns'] = num_samples

    print(vars(paras))

    mean_degree_list = np.linspace(paras.mind, paras.maxd, paras.ns)
    max_degree = 2 * paras.maxd # inf
    num_cores = min(multiprocessing.cpu_count(), paras.nc)
    q_dict, mu_dict = generate_new_transmissibilities_mutation(paras.tm1, paras.tm2, paras.T, paras.m)
    
    t1 = q_dict['Q1']
    t2 = q_dict['Q2']
    m1 = mu_dict['mu11']
    m2 = mu_dict['mu22']
    
    t_r_1 = t1
    t_r_2 = t2
    u_r_11 = m1
    u_r_12 = 1 - u_r_11
    u_r_22 = m2
    u_r_21 = 1 - u_r_22



    ###### Run on multiple cores ###### 
    pe_0_list = Manager().dict()
    pe_1_list = Manager().dict()
    pe_list = Manager().dict()
    
    Parallel(n_jobs = num_cores)(delayed(cascade_prob)(i, t_r_1, t_r_2, u_r_11, u_r_12, u_r_21, u_r_22, max_degree, pe_0_list, pe_1_list, pe_list) for i in mean_degree_list)

    ######### Save the results for all Mean Degrees ########## 
    infection_size_list = []
    infection_size_list.append(pe_0_list)
    infection_size_list.append(pe_1_list)
    infection_size_list.append(pe_list)
    
    write_analysis_results(paras, infection_size_list, 'Mutation', 'PE', mean_degree_list)
#     print("Parrell finished! Start wrting json...")

#     if change == 0:
#         change_folder = 'change_m'
#     elif change == 1:
#         change_folder = 'change_T'
#     else:
#         change_folder = 'change_tm'

#     now_finish = datetime.now() # current date and time
#     timeExp = now_finish.strftime("%m%d%H:%M")

#     ExpPath = 'Mutation_PE_Analysis_'+ change_folder +'/' + 'm' + str(mask_prob) + '_T' + "{0:.2f}".format(T) + '_tm1_' + "{0:.2f}".format(T_mask1) + '_tm2_' + "{0:.2f}".format(T_mask2) + '/' + timeExp

#     if not os.path.exists(ExpPath):
#         os.makedirs(ExpPath)
#     print("Mask Analysis results stored in: ", ExpPath)


#     setting_path = ExpPath + '/' + 'Settings'
#     if not os.path.exists(setting_path):
#         os.mkdir(setting_path)

#     res_path = ExpPath + '/' + 'Results'
#     if not os.path.exists(res_path):
#         os.mkdir(res_path) 

#     with open(setting_path + "/paras.json", "w") as fp:
#         json.dump(paras,fp) 

#     with open(res_path + "/pe_0_list.json", "w") as fp:
#         json.dump(pe_0_list.copy(),fp) 

#     with open(res_path + "/pe_1_list.json", "w") as fp:
#         json.dump(pe_1_list.copy(),fp) 

#     with open(res_path + "/pe_list.json", "w") as fp:
#         json.dump(pe_list.copy(),fp) 

    print("All done!")
main()