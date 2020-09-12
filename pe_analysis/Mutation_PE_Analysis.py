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

"""
$time python anaEvProb-Parellel.py  -n 200000  -th 0.01 -m 0.45 -T 0.6 -tm1 0.4 -tm2 0.6 -md 2 -ns 2 -nc 40 -change 0
"""
print('Running Mutation_PE_Analysis.py ...')

########### Mutation Model PE Analysis -- Parellel ########### 
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

#     print("Q1: %.5f" %Q1)
#     print("Q2: %.5f" %Q2)

#     print("mu11: %.5f" %mu11)
#     print("mu12: %.5f" %mu12)
#     print("mu22: %.5f" %mu22)
#     print("mu21: %.5f" %mu21)
    return Q_dict, mu_dict


def obtain_val_r_1(v1, v2, t1, mean_degree):
    val = 0

    for d_r in range(0, max_degree):
        prob_r = poisson.pmf(d_r, mean_degree)
        val += d_r*prob_r*1.0/mean_degree * ((1 - t1 + t1*u_r_11*v1 + t1*u_r_12*v2)**(d_r-1))

    return val

def obtain_val_r_2(v1, v2, t2, mean_degree):
    val = 0

    for d_r in range(0, max_degree):
        prob_r = poisson.pmf(d_r, mean_degree)
        val += d_r*prob_r*1.0/mean_degree * ((1 - t2 + t2*u_r_22*v2 + t2*u_r_21*v1)**(d_r-1))

    return val

def equations(p, mean_degree):
    v1, v2 = p
    val_r_1 = obtain_val_r_1(v1, v2, t_r_1, mean_degree)
    val_r_2 = obtain_val_r_2(v1, v2, t_r_2, mean_degree)

    return (v1 - val_r_1, v2 - val_r_2)


def cascade_prob(mean_degree):
    h_r_1, h_r_2 = fsolve(equations, (0.01, 0.01), args=(mean_degree), xtol=1e-6)

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
# num_cores = multiprocessing.cpu_count()
mask_prob = paras.m
T_mask1 = paras.tm1
T_mask2 = paras.tm2
T = paras.T
degree_max = paras.md
num_samples = paras.ns
change = paras.change

paras = dict()
paras['n'] = numNodes
paras['th'] = thrVal
paras['tm1'] = T_mask1
paras['tm2'] = T_mask2
paras['m'] = mask_prob
paras['T'] = T
paras['md'] = degree_max
paras['ns'] = num_samples

print(paras)

mean_degree_list = np.linspace(0, degree_max, num_samples)
max_degree = 2 * degree_max # inf
q_dict, mu_dict = generate_new_transmissibilities_mutation(T_mask1, T_mask2, T, mask_prob)

t1 = q_dict['Q1']
t2 = q_dict['Q2']
m1 = mu_dict['mu11']
m2 = mu_dict['mu22']

###### Run on multiple cores ###### 
pe_0_list = Manager().dict()
pe_1_list = Manager().dict()
pe_list = Manager().dict()

numExp = len(mean_degree_list)

t_r_1 = t1
t_r_2 = t2
u_r_11 = m1
u_r_12 = 1 - u_r_11
u_r_22 = m2
u_r_21 = 1 - u_r_22


Parallel(n_jobs = num_cores)(delayed(cascade_prob)(i) for i in mean_degree_list)

######### Save the results for all Mean Degrees ########## 
print("Parrell finished! Start wrting json...")

if change == 0:
    change_folder = 'change_m'
elif change == 1:
    change_folder = 'change_T'
else:
    change_folder = 'change_tm'

now_finish = datetime.now() # current date and time
timeExp = now_finish.strftime("%m%d%H:%M")
    
ExpPath = 'Mutation_PE_Analysis_'+ change_folder +'/' + 'm' + str(mask_prob) + '_T' + "{0:.2f}".format(T) + '_tm1_' + "{0:.2f}".format(T_mask1) + '_tm2_' + "{0:.2f}".format(T_mask2) + '/' + timeExp

if not os.path.exists(ExpPath):
    os.makedirs(ExpPath)
print("Mask Analysis results stored in: ", ExpPath)


setting_path = ExpPath + '/' + 'Settings'
if not os.path.exists(setting_path):
    os.mkdir(setting_path)

res_path = ExpPath + '/' + 'Results'
if not os.path.exists(res_path):
    os.mkdir(res_path) 

with open(setting_path + "/paras.json", "w") as fp:
    json.dump(paras,fp) 

with open(res_path + "/pe_0_list.json", "w") as fp:
    json.dump(pe_0_list.copy(),fp) 
    
with open(res_path + "/pe_1_list.json", "w") as fp:
    json.dump(pe_1_list.copy(),fp) 
    
with open(res_path + "/pe_list.json", "w") as fp:
    json.dump(pe_list.copy(),fp) 
    
print("All done!")