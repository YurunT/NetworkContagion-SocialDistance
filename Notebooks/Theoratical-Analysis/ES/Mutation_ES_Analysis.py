from __future__ import division
import argparse
import math
import sys, site, os

site.addsitedir('/afs/ece.cmu.edu/usr/reletreb/dep/lib/python2.7/site-packages')


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
# import ray
from datetime import datetime

# ray.init()

def generate_new_transmissibilities_mask(T_mask1, T_mask2, T, m):
    T1 = T * T_mask1 
    T2 = T * T_mask1 * T_mask2
    T3 = T 
    T4 = T * T_mask2

#     T1 = T * T_mask1 * T_mask2 * m
#     T2 = T * T_mask1 * (1 - m)
#     T3 = T * (1 - m)
#     T4 = T * T_mask2 * m
    

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


# numNodes = 100000;
# rho = 1.0/numNodes

# max_degree = 20  #### P.S. 1 I am using 22, I stop when the p(degree) < 1/node_num ####

def obtain_val_r_1(v1, v2, t1, lambda_r):
    val = 0

    for d_r in range(0, max_degree):
        if d_r == 0: continue

        prob_r = poisson.pmf(d_r, lambda_r)
        tmp_val = 0

        for k1 in range(0, d_r):
            for k2 in range(0, d_r - k1):
                if k1 == 0 and k2 == 0: continue

                extra_term = 0

                for x in range(1,k1+1):
                    for y in range(1,k2+1):
                        extra_term += comb(k1, x) * comb(k2, y) * t1**x * t2**y *\
                        (1-t1)**(k1-x) * (1-t2)**(k2-y) *\
                        ((x*u_r_11+y*u_r_21)/(x+y))

                a1 = (1-t1)**k1
                a2 = (1-t2)**k2
                tmp_val += comb(d_r - 1, k1) * comb(d_r - 1 - k1, k2) *(v1**k1) * (v2 ** k2) * \
                ((1 - v1 - v2) ** (d_r - 1 - k1 - k2)) *\
                (a2*(1-a1)*u_r_11 + a1*(1-a2)*u_r_21  + extra_term )

        val += d_r*prob_r*1.0/lambda_r * tmp_val

    return rho + (1 - rho)*val  #### Q1: Why need to this? Shouldn't return val as q_{l + 1, i}? ####

def obtain_val_r_2(v1, v2, t2, lambda_r):
    val = 0

    for d_r in range(0, max_degree):
        if d_r == 0: continue
        prob_r = poisson.pmf(d_r, lambda_r)
        tmp_val = 0


        for k1 in range(0, d_r):
            for k2 in range(0, d_r - k1):
                if k1 == 0 and k2 == 0: continue

                extra_term = 0

                for x in range(1,k1+1):
                    for y in range(1,k2+1):
                        extra_term += comb(k1, x) * comb(k2, y) * t1**x * \
                        t2**y * (1-t1)**(k1-x) * (1-t2)**(k2-y) *\
                        ((x*u_r_12+y*u_r_22)/(x+y))

                a1 = (1-t1)**k1
                a2 = (1-t2)**k2
                tmp_val += comb(d_r - 1, k1) * comb(d_r - 1 - k1, k2) *(v1**k1) * (v2 ** k2) * \
                ((1 - v1 - v2) ** (d_r - 1 - k1 - k2)) *\
                (a2*(1-a1)*u_r_12 + a1*(1-a2)*u_r_22  + extra_term )

        val += d_r*prob_r*1.0/lambda_r * tmp_val

    return rho + (1 - rho)*val #### Q1: Why need to this? Shouldn't return val as q_{l + 1, i}? ####

def equations(p, lambda_r):
    v1, v2 = p
    val_r_1 = obtain_val_r_1(v1, v2, t1, lambda_r)
    val_r_2 = obtain_val_r_2(v1, v2, t2, lambda_r)

    return (v1 - val_r_1, v2 - val_r_2) #### f(vl, v2) = v1, v2 ####

# @ray.remote
def cascade_size(lambda_r):
    h_r_1, h_r_2 = fsolve(equations, (0.9, 0.9), args=(lambda_r), xtol=1e-10)

    H = 0
    H1 = 0
    H2 = 0
    for d_r in range(0, max_degree):
        prob_r = poisson.pmf(d_r, lambda_r)
        tmp_val1 = 0
        tmp_val2 = 0

        for k1 in range(0, d_r + 1):
            for k2 in range(0, d_r + 1 - k1):
                extra_term_1 = 0
                extra_term_2 = 0
                for x in range(1,k1+1):
                    for y in range(1,k2+1):
                        extra_term_1 += comb(k1, x) * comb(k2, y) * t1**x *\
                        t2**y * (1-t1)**(k1-x) * (1-t2)**(k2-y) *\
                        ((x*u_r_11+y*u_r_21)/(x+y))
                        extra_term_2 += comb(k1, x) * comb(k2, y) * t1**x *\
                        t2**y * (1-t1)**(k1-x) * (1-t2)**(k2-y) *\
                        ((x*u_r_12+y*u_r_22)/(x+y))

                a1 = (1-t1)**k1
                a2 = (1-t2)**k2
                tmp_val1 += comb(d_r , k1) * comb(d_r - k1, k2) *(h_r_1**k1) * (h_r_2 ** k2) * \
                ((1 - h_r_1 - h_r_2) ** (d_r - k1 - k2)) *\
                (a2*(1-a1)*u_r_11 + a1*(1-a2)*u_r_21  + extra_term_1 )

                tmp_val2 += comb(d_r, k1) * comb(d_r - k1, k2) *(h_r_1**k1) * (h_r_2 ** k2) * \
                ((1 - h_r_1 - h_r_2) ** (d_r - k1 - k2)) *\
                (a2*(1-a1)*u_r_12 + a1*(1-a2)*u_r_22  + extra_term_2 )

        H += prob_r * (tmp_val1 + tmp_val2) #### Q2: Need to add all of them together? ####
        H1 += prob_r * tmp_val1
        H2 += prob_r * tmp_val2

    print(lambda_r, H, H1, H2)
    infection_size_mu[lambda_r] = H
    infection_size0_mu[lambda_r] = H1
    infection_size1_mu[lambda_r] = H2
    
    return (lambda_r, H, H1, H2)

def parse_args(args):
    parser = argparse.ArgumentParser(description = 'Parameters')
    parser.add_argument('-n', type = int, default = 200000, help='200,000 (default); the number of nodes')
    parser.add_argument('-e', type = int, default = 50, help='50 (default); the number of experiments')
    parser.add_argument('-m', type = float, default = 0.6, help='0.6 (default); the prob of wearing a mask')
    parser.add_argument('-tm1', type = float, default = 0.5, help='0.5 (default); T_mask1')
    parser.add_argument('-tm2', type = float, default = 0.5, help='0.5 (default); T_mask2')
    parser.add_argument('-T', type = float, default = 0.6, help='0.6 (default); transmissibility of the orinigal virus')
    parser.add_argument('-th', type = float, default = 0.001, help='0.001 (default); the treshold to consider a component giant')
    parser.add_argument('-nc', type = int, default = 12, help='number of Cores')
    parser.add_argument('-md', type = int, default = 10, help='[0, max_degree]')
    parser.add_argument('-ns', type = int, default = 50, help='Num of sample within [0, max_degree]')
    parser.add_argument('-change', type = int, default = 0, help='Change m(0), T(1), tm(2)')
    return parser.parse_args(args)

paras = parse_args(sys.argv[1:])
# t1 = paras.t1
# t2 = paras.t2
# m1 = paras.m1
# m2 = paras.m2
# mean_degree_list = paras.m






numNodes = paras.n
rho = 1.0/numNodes

numExp = paras.e
thrVal = paras.th
# num_cores = min(paras.numCores,multiprocessing.cpu_count())
num_cores = multiprocessing.cpu_count()
mask_prob = paras.m
T_mask1 = paras.tm1
T_mask2 = paras.tm2
T = paras.T
degree_max = paras.md
num_samples = paras.ns
change = paras.change

mean_degree_list = np.linspace(0, degree_max, num_samples)
max_degree = 2 * degree_max # inf


q_dict, mu_dict = generate_new_transmissibilities_mutation(T_mask1, T_mask2, T, mask_prob)

t1 = q_dict['Q1']
t2 = q_dict['Q2']
m1 = mu_dict['mu11']
m2 = mu_dict['mu22']
print("Q_dict:", q_dict)
print("mu_dict:", mu_dict)

######
infection_size_mu = Manager().dict()
infection_size0_mu = Manager().dict()
infection_size1_mu = Manager().dict()
# infection_size_mu = []
# infection_size0_mu = []
# infection_size1_mu = []

num_cores = multiprocessing.cpu_count()
numExp = len(mean_degree_list)



# results_ids = []
# for lambda_r in mean_degree_list:
#     u_r_11 = m1
#     u_r_12 = 1-u_r_11
#     u_r_22 = m2
#     u_r_21 = 1-u_r_22
    
#     results_ids.append(cascade_size.remote(lambda_r))
#     lambda_r, H, H1, H2, h_r_1, h_r_2 = cascade_size.remote(lambda_r)

    
#     lambda_r, H, H1, H2 = cascade_size(lambda_r)
#     infection_size_mu.append(H) 
#     infection_size0_mu.append(H1)
#     infection_size1_mu.append(H2)
# ray.get(infection_size_mu)
# ray.get(infection_size0_mu)
# ray.get(infection_size1_mu)



# results = ray.get(results_ids)
# print(results)
u_r_11 = m1
u_r_12 = 1-u_r_11
u_r_22 = m2
u_r_21 = 1-u_r_22
Parallel(n_jobs = num_cores)(delayed(cascade_size)(lambda_r) for lambda_r in mean_degree_list)

print("Parrell finished! Start wrting json...")

# for idx, mean_degree in enumerate(mean_degree_list):
#     infection_size_mu.append(results[idx][1])
#     infection_size0_mu.append(results[idx][2])
#     infection_size1_mu.append(results[idx][3])
    
'''
With ray
[0.0, 0.31999068069847]
[0.0, 0.18461000809527114]
[0.0, 0.13538067260319883]
'''

'''
NO ray
[0.0, 2.4007115114959465e-07]
[0.0, 1.3879791630781522e-07]
[0.0, 1.0127323484177944e-07]
'''
######### Save the results for all Mean Degrees ########## 
print("Parellel finished! Start wrting json...")
print(infection_size_mu)
print(infection_size0_mu)
print(infection_size1_mu)

if change == 0:
    change_folder = 'change_m'
elif change == 1:
    change_folder = 'change_T'
else:
    change_folder = 'change_tm'
    
now_finish = datetime.now() # current date and time
timeExp = now_finish.strftime("%m%d%H:%M")

ExpPath = 'MutationParellel_ES_Analysis_'+ change_folder +'/' + 'm' + str(mask_prob) + '_T' + "{0:.2f}".format(T) + '_tm1_' + "{0:.2f}".format(T_mask1) + '_tm2_' + "{0:.2f}".format(T_mask2) + '/' + timeExp


if not os.path.exists(ExpPath):
    os.makedirs(ExpPath)
print("MutationRay Analysis results stored in: ", ExpPath)


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


with open(res_path + "/infection_size_mu.json", "w") as fp:
    json.dump(infection_size_mu.copy(),fp) 
    
with open(res_path + "/infection_size0_mu.json", "w") as fp:
    json.dump(infection_size0_mu.copy(),fp) 
with open(res_path + "/infection_size1_mu.json", "w") as fp:
    json.dump(infection_size1_mu.copy(),fp) 
    
print("All done!")

# with open("results.json", "w") as fp:
#     json.dump(results.copy(),fp) 

# with open("infection_size_mu.json", "w") as fp:
#     json.dump(infection_size_mu.copy(),fp) 
    
# with open("infection_size0_mu.json", "w") as fp:
#     json.dump(infection_size0_mu.copy(),fp) 
# with open("infection_size1_mu.json", "w") as fp:
#     json.dump(infection_size1_mu.copy(),fp) 

# np.save("infection_size_mu.npy", list(infection_size_mu.values()))
# np.save("infection_size0_mu.npy", list(infection_size0_mu.values()))
# np.save("infection_size1_mu.npy", list(infection_size1_mu.values()))