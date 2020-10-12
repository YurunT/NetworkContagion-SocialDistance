from __future__ import division
import math
import sys, site, os
import numpy as np
from scipy import optimize 
from scipy.special import comb
from scipy.stats import poisson
import scipy.optimize
import scipy.misc
from datetime import datetime

########### Mask Model ES Analysis -- Parellel ########### 
def generate_degree_list(mean_degree, nodeN):
    degree_max = nodeN - 1
    p_k = dict() # k: degree, v: prob
    for degree in range(1 , degree_max + 1):
        p_k[degree] = poisson.pmf(degree, mean_degree)
        if p_k[degree] < 10 ** (- math.log10(nodeN) - 1): # Stop when the prob < 1/node_num (like a thr)
            break
    return p_k, degree

def P_A_given_R(i, T_list, k0, k1):
    one_minus_T1 = 1 - T_list[0]
    one_minus_T2 = 1 - T_list[1]
    one_minus_T3 = 1 - T_list[2]
    one_minus_T4 = 1 - T_list[3]
    if i == 0:
        res = 1 - (one_minus_T2 ** k0) * (one_minus_T4 ** k1)
    else:
        res = 1 - (one_minus_T1 ** k0) * (one_minus_T3 ** k1)
    assert res >= 0, "P_A_given_R should be greater than 0"
    assert res <= 1, "P_A_given_R should be less than 1"
    return res

def P_A_given_B_N(i, k, n, T_list, A_0, A_1):
    one_minus_A0 = 1 - A_0
    one_minus_A1 = 1 - A_1
    
    p_abn = 0
    
    for k0 in range(n + 1):
        for k1 in range(k - n):
            p_a_given_r = P_A_given_R(i, T_list, k0, k1)
            p_abn += p_a_given_r * \
                     comb(n, k0) * \
                     comb(k - 1 - n, k1) * \
                     (A_0 ** k0) * \
                     (A_1 ** k1) * \
                     (one_minus_A0 ** (n - k0)) * \
                     (one_minus_A1 ** (k - 1 - n - k1))
    return p_abn

def P_A_given_B(i, k, T_list, A_0, A_1, m):
    p_ab = 0
    for n in range(k):
        p_abn = P_A_given_B_N(i, k, n, T_list, A_0, A_1)
        p_ab += p_abn * comb(k - 1, n) * \
                (m ** n) * \
                ((1 - m) ** (k - 1 - n))
    return p_ab

def P_A(i, mean_degree, nodeN, T_list, m, A_0, A_1, k_max):
#     P_k_dict, k_max0 = generate_degree_list(mean_degree, nodeN)
    
    pa_L = 0
    for k in range(1, k_max):
#         if k not in P_k_dict.keys():
#             p_k = 0
#         else:
#             p_k = P_k_dict[k]
        p_k = poisson.pmf(k, mean_degree)
        p_b = k * p_k / mean_degree
        p_ab = P_A_given_B(i, k, T_list, A_0, A_1, m)
        pa_L += p_ab * p_b
    return pa_L

def p_A_vec(mean_degree, nodeN, T_list, m, A_0, A_1, k_max):
    A0 = P_A(0, mean_degree, nodeN, T_list, m, A_0, A_1, k_max)
    A1 = P_A(1, mean_degree, nodeN, T_list, m, A_0, A_1, k_max)
    return [A0, A1]

def func_fix(A, mean_degree, nodeN, T_list, m):
    return np.array(p_A_vec(mean_degree, nodeN, T_list, m, A[0], A[1]))

def func_root(A, mean_degree, nodeN, T_list, m, k_max):
    return np.array(p_A_vec(mean_degree, nodeN, T_list, m, A[0], A[1], k_max)) - np.array(A)

def get_EpidemicSize(mean_degree, k_max, T_list, paras, infection_size, infection_size0, infection_size1):
    '''
    S
    '''
    init_A = (0.9, 0.9)
    m = paras.m

    A_0_1_root = optimize.fsolve(func_root, init_A, args=(mean_degree, paras.n, T_list, m, k_max))
    pa_L_0 = 0
    pa_L_1 = 0
    pa_L = 0
    
    for k in range(0, k_max): # k now can start from 0, because for the last level the k could be 0 [0, kmax)
        p_k = poisson.pmf(k, mean_degree) 
        p_ab_0 = 0
        p_ab_1 = 0
        
        for n in range(k + 1): # n can be k
       
            p_abn_0 = 0
            p_abn_1 = 0
            
            for k0 in range(n + 1): # [0, n] 
                for k1 in range(k - n + 1): # k1 can reach k - n
                    p_a_given_r_0 = P_A_given_R(0, T_list, k0, k1)
                    p_a_given_r_1 = P_A_given_R(1, T_list, k0, k1)
                    
                    p_abn_0 += p_a_given_r_0 * \
                             comb(n, k0) * \
                             comb(k - n, k1) * \
                             (A_0_1_root[0] ** k0) * \
                             (A_0_1_root[1] ** k1) * \
                             ((1 - A_0_1_root[0]) ** (n - k0)) * \
                             ((1 - A_0_1_root[1]) ** (k - n - k1))
                    
                    p_abn_1 += p_a_given_r_1 * \
                             comb(n, k0) * \
                             comb(k - n, k1) * \
                             (A_0_1_root[0] ** k0) * \
                             (A_0_1_root[1] ** k1) * \
                             ((1 - A_0_1_root[0]) ** (n - k0)) * \
                             ((1 - A_0_1_root[1]) ** (k - n - k1))
                    
            p_ab_0 += p_abn_0 * comb(k, n) * \
                    (m ** n) * \
                    ((1 - m) ** (k - n))
            p_ab_1 += p_abn_1 * comb(k, n) * \
                    (m ** n) * \
                    ((1 - m) ** (k - n))
            
            
        pa_L_0 += p_k * p_ab_0
        pa_L_1 += p_k * p_ab_1
        pa_L += p_k * (m * p_ab_0 + (1 - m) * p_ab_1) 
    
    print(mean_degree, pa_L, pa_L_0, pa_L_1)
    infection_size0[mean_degree] = pa_L_0
    infection_size1[mean_degree] = pa_L_1
    infection_size[mean_degree]  = pa_L
    
    return pa_L_0, pa_L_1, pa_L, 1 - np.array(A_0_1_root)


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

def func_root_pe(E, mean_degree, T_list, m, max_degree):
    return PE_vec(mean_degree, True, T_list, m, E[0], E[1], max_degree) - np.array(E)

def get_ProbEmergence(mean_degree, paras, k_max, T_list, pe_list_m, pe_0_list_m, pe_1_list_m, ):
    E0, E1 = optimize.fsolve(func_root_pe, (0.01, 0.01), args=(mean_degree, T_list, paras.m, k_max), xtol=1e-6)    
    E0, E1 = 1 - PE_vec(mean_degree, False,  T_list, paras.m, E0, E1, k_max)
    pe_list_m[mean_degree]   = paras.m * E0 + (1 - paras.m) * E1
    pe_0_list_m[mean_degree] = E0
    pe_1_list_m[mean_degree] = E1
    print(mean_degree, E0, E1) 