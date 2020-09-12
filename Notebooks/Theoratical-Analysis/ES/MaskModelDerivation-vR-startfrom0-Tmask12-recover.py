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
$time python MaskModelDerivation-vR-startfrom0-Tmask12.py  -n 200000  -th 0.01 -m 0.45 -T 0.6 -tm1 0.4 -tm2 0.6 -md 2 -ns 2 -nc 40 -change 0
"""

########### Mask Model ES Analysis -- Parellel ########### 

def generate_new_transmissibilities_mask(T_mask1, T_mask2, T, m):
    T2 = T * T_mask1 
    T1 = T * T_mask1 * T_mask2
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

def P_A_given_R(i, T_list, k0, k1):
    one_minus_T1 = 1 - T_list[0]
    one_minus_T2 = 1 - T_list[1]
    one_minus_T3 = 1 - T_list[2]
    one_minus_T4 = 1 - T_list[3]
    
    if i == 0:
        res = 1 - (one_minus_T2 ** k0) * (one_minus_T4 ** k1)
    else:
        res = 1 - (one_minus_T1 ** k0) * (one_minus_T3 ** k1)
    return res

def P_A_given_B_N(i, is_intermediate, k, n, T_list, A_0, A_1):
    one_minus_A0 = 1 - A_0
    one_minus_A1 = 1 - A_1
    p_abn = 0
    
    if is_intermediate:
        k1_range = k - n
    else:
        k1_range = k - n + 1
        
    for k0 in range(0, n + 1): # Start from 0! Or the curves will be totally fucked up!
        for k1 in range(0, k1_range):
            p_a_given_r = P_A_given_R(i, T_list, k0, k1)
            p_abn += p_a_given_r * \
                     comb(n, k0) * \
                     comb(k1_range - 1, k1) * \
                     (A_0 ** k0) * \
                     (A_1 ** k1) * \
                     (one_minus_A0 ** (n - k0)) * \
                     (one_minus_A1 ** (k1_range - 1 - k1))
    return p_abn

def P_A_given_B(i, is_intermediate, k, T_list, A_0, A_1, m):
    p_ab = 0
    
    if is_intermediate:
        n_range = k
    else:
        n_range = k + 1
        
    for n in range(n_range):
        p_abn = P_A_given_B_N(i, is_intermediate, k, n, T_list, A_0, A_1)
        p_ab += p_abn * comb(n_range - 1, n) * \
                (m ** n) * \
                ((1 - m) ** (n_range - 1 - n))
    return p_ab

def P_A(i, is_intermediate, mean_degree, nodeN, T_list, m, A_0, A_1):
    pa_L = 0

    for k in range(0, max_degree):
        p_k = poisson.pmf(k, mean_degree)
        
        if is_intermediate:
            if k == 0: continue #### Sep.5 Doesn't influence the final results
            p_b = k * p_k * 1.0 / mean_degree
        else:
            p_b = p_k
        
        p_ab = P_A_given_B(i, is_intermediate, k, T_list, A_0, A_1, m)
        pa_L += p_ab * p_b
        
    return rho + (1 - rho)*pa_L

def p_A_vec(mean_degree, nodeN, T_list, m, A_0, A_1, is_intermediate):
    A0 = P_A(0, is_intermediate, mean_degree, nodeN, T_list, m, A_0, A_1)
    A1 = P_A(1, is_intermediate, mean_degree, nodeN, T_list, m, A_0, A_1)
    return [A0, A1]

def func_fix(A, mean_degree, nodeN, T_list, m):
    return np.array(p_A_vec(mean_degree, nodeN, T_list, m, A[0], A[1]))

def func_root(A, mean_degree, nodeN, T_list, m,):
    return np.array(p_A_vec(mean_degree, nodeN, T_list, m, A[0], A[1], True)) - np.array(A)

"""
The only change from the orignial jupyter
"""
def get_EpidemicSize(mean_degree, nodeN, T_list, m):
    init_A = (0.9, 0.9)
         
    A_0_1_root = optimize.fsolve(func_root, init_A, args=(mean_degree, nodeN, T_list, m))
    
    pa_L_0 = 0
    pa_L_1 = 0
    pa_L = 0
    pa_L_0, pa_L_1 = p_A_vec(mean_degree, nodeN, T_list, m, A_0_1_root[0], A_0_1_root[1], False)
    
    infection_size[mean_degree] = m * pa_L_0 + (1 - m) * pa_L_1
    infection_size0[mean_degree] = pa_L_0
    infection_size1[mean_degree] = pa_L_1
    
    print(pa_L_0, pa_L_1)


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
infection_size = Manager().dict()
infection_size0 = Manager().dict()
infection_size1 = Manager().dict()

Parallel(n_jobs = num_cores)(delayed(get_EpidemicSize)(mean_degree, numNodes, T_list, mask_prob) for mean_degree in mean_degree_list)


######### Save the results for all Mean Degrees ########## 
print("Parrell finished! Start wrting json...")
print(infection_size)
print(infection_size0)
print(infection_size1)

if change == 0:
    change_folder = 'change_m'
elif change == 1:
    change_folder = 'change_T'
else:
    change_folder = 'change_tm'
    
now_finish = datetime.now() # current date and time
timeExp = now_finish.strftime("%m%d%H:%M")

ExpPath = 'Mask_ES_Analysis_'+ change_folder +'/' + 'm' + str(mask_prob) + '_T' + "{0:.2f}".format(T) + '_tm1_' + "{0:.2f}".format(T_mask1) + '_tm2_' + "{0:.2f}".format(T_mask2) + '/' + timeExp


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


with open(res_path + "/infection_size.json", "w") as fp:
    json.dump(infection_size.copy(),fp) 
    
with open(res_path + "/infection_size0.json", "w") as fp:
    json.dump(infection_size0.copy(),fp) 
with open(res_path + "/infection_size1.json", "w") as fp:
    json.dump(infection_size1.copy(),fp) 
    
print("All done!")