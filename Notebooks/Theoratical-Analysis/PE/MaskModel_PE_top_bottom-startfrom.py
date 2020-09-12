from __future__ import division
import argparse
import math
import sys, site, os
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
from datetime import datetime

"""
$time python MaskModel_PE_top_bottom-startfrom.py  -n 200000  -th 0.01 -m 0.45 -T 0.6 -tm1 0.4 -tm2 0.6 -md 2 -ns 2 -nc 40 -change 0
"""

########### Mask Model PE Analysis -- Parellel ########### 

def generate_new_transmissibilities_mask(T_mask1, T_mask2, T, m):
    T1 = T * T_mask1 
    T2 = T * T_mask1 * T_mask2
    T3 = T 
    T4 = T * T_mask2
    

    trans_dict = {'T1': T1,
                  'T2': T2,
                  'T3': T3,
                  'T4': T4}
    
    return trans_dict    

def generate_new_transmissibilities_mutation(T_mask1, T_mask2, T, m):
    trans_dict = generate_new_transmissibilities_mask(T_mask1, T_mask2, T, m)
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

    print("Q1: %.5f" %Q1)
    print("Q2: %.5f" %Q2)

    print("mu11: %.5f" %mu11)
    print("mu12: %.5f" %mu12)
    print("mu22: %.5f" %mu22)
    print("mu21: %.5f" %mu21)
    return Q_dict, mu_dict


def PE(i, is_intermediate, E0, E1, T_list, m, mean_degree):
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


def PE_vec(mean_degree, is_intermediate,  T_list, m, E0, E1):
    E0 = PE(0, is_intermediate, E0, E1, T_list, m, mean_degree)
    E1 = PE(1, is_intermediate, E0, E1, T_list, m, mean_degree)
    return np.array([E0, E1])

def func_root(E, mean_degree, T_list, m):
    return PE_vec(mean_degree, True, T_list, m, E[0], E[1]) - np.array(E)

def get_ProbEmergence(mean_degree, nodeN, T_list, m):
    E0, E1 = optimize.fsolve(func_root, (0.01, 0.01), args=(mean_degree, T_list, m), xtol=1e-6)    
    E0, E1 = 1 - PE_vec(mean_degree, False,  T_list, m, E0, E1)
    pe_list_m[mean_degree] = m * E0 + (1 - m) * E1
    pe_0_list_m[mean_degree] = E0
    pe_1_list_m[mean_degree] = E1
    print(E0, E1) 


def parse_args(args):
    parser = argparse.ArgumentParser(description = 'Parameters')
    parser.add_argument('-n', type = int, default = 200000, help='200,000 (default); the number of nodes')
    parser.add_argument('-m', type = float, default = 0.6, help='0.6 (default); the prob of wearing a mask')
    parser.add_argument('-tm1', type = float, default = 0.5, help='0.5 (default); T_mask1')
    parser.add_argument('-tm2', type = float, default = 0.5, help='0.5 (default); T_mask2')
    parser.add_argument('-T', type = float, default = 0.6, help='0.6 (default); transmissibility of the orinigal virus')
    parser.add_argument('-th', type = float, default = 0.001, help='0.001 (default); the treshold to consider a component giant')
    parser.add_argument('-nc', type = int, default = 40, help='number of Cores')
    parser.add_argument('-md', type = int, default = 10, help='[0, max_degree]')
    parser.add_argument('-ns', type = int, default = 50, help='Num of sample within [0, max_degree]')
    parser.add_argument('-change', type = int, default = 0, help='Change m(0), T(1), tm(2)')
    return parser.parse_args(args)

###### Paras setting ######
paras = parse_args(sys.argv[1:])
numNodes = paras.n
rho = 1.0/numNodes

thrVal = paras.th
num_cores = min(paras.nc,multiprocessing.cpu_count())
mask_prob = paras.m
T_mask1 = paras.tm1
T_mask2 = paras.tm2
T = paras.T
degree_max = paras.md
num_samples = paras.ns
change = paras.change


mean_degree_list = np.linspace(0, degree_max, num_samples)
max_degree = 2 * degree_max # inf
T_list = list(generate_new_transmissibilities_mask(T_mask1, T_mask2, T, mask_prob).values())


###### Run on multiple cores ###### 
pe_list_m = Manager().dict()
pe_0_list_m = Manager().dict()
pe_1_list_m = Manager().dict()


Parallel(n_jobs = num_cores)(delayed(get_ProbEmergence)(mean_degree, numNodes, T_list, mask_prob) for mean_degree in mean_degree_list)

######### Save the results for all Mean Degrees ########## 
print("Parrell finished! Start wrting json...")
print(pe_list_m)
print(pe_0_list_m)
print(pe_1_list_m)

if change == 0:
    change_folder = 'change_m'
elif change == 1:
    change_folder = 'change_T'
else:
    change_folder = 'change_tm'
    
now_finish = datetime.now() # current date and time
timeExp = now_finish.strftime("%m%d%H:%M")

ExpPath = 'Mask_PE_Analysis_'+ change_folder +'/' + 'm' + str(mask_prob) + '_T' + "{0:.2f}".format(T) + '_tm1_' + "{0:.2f}".format(T_mask1) + '_tm2_' + "{0:.2f}".format(T_mask2) + '/' + timeExp


if not os.path.exists(ExpPath):
    os.makedirs(ExpPath)
print("Mask Analysis results stored in: ", ExpPath)


setting_path = ExpPath + '/' + 'Settings'
if not os.path.exists(setting_path):
    os.mkdir(setting_path)

res_path = ExpPath + '/' + 'Results'
if not os.path.exists(res_path):
    os.mkdir(res_path) 
    
paras = dict()
paras['n'] = numNodes
paras['th'] = thrVal
paras['tm1'] = T_mask1
paras['tm2'] = T_mask2
paras['m'] = mask_prob
paras['T'] = T
paras['md'] = degree_max
paras['ns'] = num_samples

with open(setting_path + "/paras.json", "w") as fp:
    json.dump(paras,fp) 


with open(res_path + "/pe_0_list_m.json", "w") as fp:
    json.dump(pe_0_list_m.copy(),fp) 
    
with open(res_path + "/pe_1_list_m.json", "w") as fp:
    json.dump(pe_1_list_m.copy(),fp) 
    
with open(res_path + "/pe_list_m.json", "w") as fp:
    json.dump(pe_list_m.copy(),fp) 
    
print("All done!")

