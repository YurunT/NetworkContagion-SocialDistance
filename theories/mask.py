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
sys.path.append(os.path.abspath("../auxiliary_scripts/"))
from main_aux import *


########### Mask Model ES Analysis -- Parellel ########### 
def P_A_given_R(i, T_list, k0, k1):
    one_minus_T1 = 1 - T_list[0][1]
    one_minus_T2 = 1 - T_list[0][0]
    one_minus_T3 = 1 - T_list[1][1]
    one_minus_T4 = 1 - T_list[1][0]
    if i == 0:
        res = 1 - (one_minus_T2 ** k0) * (one_minus_T4 ** k1)
    else:
        res = 1 - (one_minus_T1 ** k0) * (one_minus_T3 ** k1)
    assert res >= 0, "P_A_given_R should be greater than 0"
    assert res <= 1, "P_A_given_R should be less than 1"
    return res

def P_A_given_B_N(i, is_intermediate, k, n, T_list, A_0, A_1):
    one_minus_A0 = 1 - A_0
    one_minus_A1 = 1 - A_1
    
    p_abn = 0
    
    if is_intermediate:
        k1_range = k - n
    else:
        k1_range = k + 1 - n
    
    for k0 in range(n + 1):
        for k1 in range(k1_range):
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

def P_A(i, is_intermediate, mean_degree, T_list, m, A, k_max):    
    pa_L = 0
    for k in range(1, k_max):
        prob_r = poisson.pmf(k, mean_degree)
        if is_intermediate: # intermediate q using excess degree distribution
            pb = prob_r * k * 1.0 / mean_degree 
        else:
            pb = prob_r
        p_ab = P_A_given_B(i, is_intermediate, k, T_list, A[0], A[1], m)
        pa_L += p_ab * pb
    return pa_L

def pA_vec(mean_degree, is_intermediate, T_list, m, A, k_max, num_mask_types):
    P_A_list = []
    for i in range(num_mask_types):
        P_A_list.append(P_A(i, is_intermediate, mean_degree, T_list, m, A, k_max))
    return P_A_list


def func_root(A, mean_degree, T_list, m, k_max, num_mask_types):
    return np.array(pA_vec(mean_degree, True, T_list, m, A, k_max, num_mask_types)) - np.array(A)

def get_EpidemicSize(mean_degree, paras, infection_sizes):
    '''
    S
    '''    
    infection_size0 = infection_sizes[0]
    infection_size1 = infection_sizes[1]
    infection_size  = infection_sizes['ttl']
    
    k_max, T_list = resolve_paras(paras)
    num_mask_types = len(T_list)
    init_A = np.ones(num_mask_types) * 0.9
    m = paras.m

    A_root = optimize.fsolve(func_root, init_A, args=(mean_degree, T_list, m, k_max, num_mask_types))
    print()
    P_A_list = pA_vec(mean_degree, False,  T_list, paras.m, A_root, k_max, num_mask_types)
    A = P_A_list[0] * paras.m + P_A_list[1] * (1 - paras.m)
    print(mean_degree, P_A_list,)
    infection_size0[mean_degree] = P_A_list[0]
    infection_size1[mean_degree] = P_A_list[1]
    infection_size[mean_degree]  = A
    return P_A_list[0], P_A_list[1], A

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
    T1 = T_list[0][1]
    T2 = T_list[0][0]
    T3 = T_list[1][1]
    T4 = T_list[1][0]
    
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

def get_ProbEmergence(mean_degree, paras, pe_0_list_m, pe_1_list_m, pe_list_m):
    k_max, T_list = resolve_paras(paras)
    E0, E1 = optimize.fsolve(func_root_pe, (0.01, 0.01), args=(mean_degree, T_list, paras.m, k_max), xtol=1e-6)    
    E0, E1 = 1 - PE_vec(mean_degree, False,  T_list, paras.m, E0, E1, k_max)
    pe_list_m[mean_degree]   = paras.m * E0 + (1 - paras.m) * E1
    pe_0_list_m[mean_degree] = E0
    pe_1_list_m[mean_degree] = E1
    print(mean_degree, E0, E1)