from __future__ import division
import sys, site, os
sys.path.append(os.path.abspath("../../auxiliary_scripts/"))
import argparse
import math
import numpy as np
from scipy import optimize 
from scipy.special import comb
from scipy.stats import poisson
import scipy.optimize
import scipy.misc
import multiprocessing
from multiprocessing import Manager
from joblib import Parallel, delayed
import json
from tnn import *
from datetime import datetime
from input_module import parse_args
from output_module import write_analysis_results

########### Mask Model PE Analysis -- Parellel ########### 
def PE(i, is_intermediate, E0, E1, T_list, m, mean_degree, max_degree):
    res = 0
    for k in range(0, max_degree):
        prob_r = poisson.pmf(k, mean_degree)
        
        if is_intermediate: # intermediate q using excess degree distribution
            pb = prob_r * k * 1.0 / mean_degree 
        else:
            pb = prob_r
            
        res += pb * PE_B(i, is_intermediate, k, E0, E1, T_list, m)
    return res

def PE_B(i, is_intermediate, k, E0, E1, T_list, m):
    res = 0
    one_minus_m = 1 - m
    
    if is_intermediate: # intermediate q, powers sum up to k - 1
        n_range = k
    else:               # generation 0, powers sum up to k
        n_range = k + 1
        
    for n in range(n_range):
        pe_bn = PE_BN(i, is_intermediate, n, k, E0, E1, T_list, m)
        res += pe_bn * comb(n_range - 1, n) * (m ** n) * (one_minus_m ** (n_range - 1 - n))
    return res

def PE_BN(i, is_intermediate, n, k, E0, E1, T_list, m):
    T1 = T_list[0]
    T2 = T_list[1]
    T3 = T_list[2]
    T4 = T_list[3]
    
    res = 0 
    
    if i == 0:
        t_mask = T2
        t_no_mask = T1
    else:
        t_mask = T4
        t_no_mask = T3
    
    one_minus_mask = 1 - t_mask
    one_minus_no_mask = 1 - t_no_mask
        
    if is_intermediate:
        k_range = k 
    else:
        k_range = k + 1
        
    for k0 in range(n + 1):
        for k1 in range(k_range - n):
            res += comb(n, k0) * (t_mask ** k0) * (one_minus_mask ** (n - k0)) *\
            comb(k_range - 1 - n, k1) * (t_no_mask ** k1) * (one_minus_no_mask ** (k_range - 1 - n - k1)) *\
            (E0 ** k0) * (E1 ** k1)
            
    return res


def PE_vec(mean_degree, is_intermediate,  T_list, m, E0, E1, max_degree):
    E0 = PE(0, is_intermediate, E0, E1, T_list, m, mean_degree, max_degree)
    E1 = PE(1, is_intermediate, E0, E1, T_list, m, mean_degree, max_degree)
    return np.array([E0, E1])

def func_root(E, mean_degree, T_list, m, max_degree):
    return PE_vec(mean_degree, True, T_list, m, E[0], E[1], max_degree) - np.array(E)

def get_ProbEmergence(mean_degree, paras, k_max, T_list, pe_list_m, pe_0_list_m, pe_1_list_m, ):
    E0, E1 = optimize.fsolve(func_root, (0.01, 0.01), args=(mean_degree, T_list, paras.m, k_max), xtol=1e-6)    
    E0, E1 = 1 - PE_vec(mean_degree, False,  T_list, paras.m, E0, E1, k_max)
    pe_list_m[mean_degree]   = paras.m * E0 + (1 - paras.m) * E1
    pe_0_list_m[mean_degree] = E0
    pe_1_list_m[mean_degree] = E1
    print(E0, E1) 

def main():
    ###### Get commadline input ######
    paras = parse_args(sys.argv[1:])
    mean_degree_list = np.linspace(paras.mind, paras.maxd, paras.ns)
    k_max = 4 * paras.maxd # inf
    num_cores = min(paras.nc,multiprocessing.cpu_count())
    T_list = list(generate_new_transmissibilities_mask(paras.tm1, paras.tm2, paras.T, paras.m).values())
    
    print('-------Parameter Setting-------\n', vars(paras))
    print("K_max: ", k_max)
    print("num_cores:", num_cores)
    print("mean_degree_list:", mean_degree_list)
    print('-------Parameter Setting-------\n')


    ###### Run on multiple cores ###### 
    pe_list_m = Manager().dict()
    pe_0_list_m = Manager().dict()
    pe_1_list_m = Manager().dict()

    Parallel(n_jobs = num_cores)(delayed(get_ProbEmergence)(mean_degree, paras, k_max, T_list, pe_list_m, pe_0_list_m, pe_1_list_m, ) for mean_degree in mean_degree_list)

    ######### Save the results for all Mean Degrees ########## 
    pe_list = [pe_0_list_m, pe_1_list_m, pe_list_m]
    write_analysis_results(paras, pe_list, 'Mask', 'PE', mean_degree_list)
    print("All done!")
main()