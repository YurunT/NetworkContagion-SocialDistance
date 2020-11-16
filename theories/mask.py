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
from theory_aux import * 

########### Mask Model ES Analysis -- Parellel ########### 
def get_par(i, T_list, end, idx, vec_I):
    mid_res = 1
    if idx <= end:
        mid_res *= (1 - T_list[idx][i]) ** vec_I[idx] * get_par(i, T_list, end, idx + 1, vec_I)
        return mid_res
    else:
        return 1
    
def P_A_given_R(i, T_list, vec_I, end):    
    return 1 - get_par(i, T_list, end, 0, vec_I)

def get_p_abn(i, T_list, A, end, idx, vec_N, vec_I,):
    ####### !!!!! PE has no par
    one_minus_A = 1 - np.array(A)
    p_abn = 0
    vec_N = vec_N.astype(int)
    
    if idx <= end:
        for ki in range(vec_N[idx] + 1):
            vec_I[idx] = ki
            p_abn += comb(vec_N[idx], ki) * (A[idx] ** ki) * (one_minus_A[idx] ** (vec_N[idx] - ki)) * get_p_abn(i, T_list, A, end, idx + 1, vec_N, vec_I,)
        return p_abn
    else:
        p_a_given_r = P_A_given_R(i, T_list, vec_I, end)
        
        return p_a_given_r

def condition_on_active_nodes(i, item, k, vec_N, T_list, A, num_mask_types,):
    vec_I = np.zeros(num_mask_types)
    if item == 'es':
        p_abn = get_p_abn(i, T_list, A, num_mask_types - 1, 0, vec_N, vec_I,)
    return p_abn

def recursion_N(i, item, is_intermediate, k, T_list, A, end, idx, m, n_i_range, vec_N, num_mask_types):
    '''
    Recursion through all the possible values for the vec N = (N1, N2, ..., NM)
    
    Input:  end: The last index of vec N, end = M - 1. M refers to the M in the derivation notebook.
            idx: Loop over all the possible values for N_idx in vector N, 0 <= idx <= end (M - 1)
            m  : List of mask type possibibities.
            n_i_range: Biggest value for N_m to have
    '''
    pab = 0
    if idx < end: # not the last ele in the N vec
        for n_i in range(n_i_range):
            vec_N[idx] = n_i
            pab += (m[idx] ** n_i) * recursion_N(i, item, is_intermediate, k, T_list, A, end, idx + 1, m, n_i_range - n_i, vec_N, num_mask_types)
        return pab
    else:
        n_end = n_i_range - 1
        vec_N[idx] = n_end
        pab = (m[idx] ** n_end)   
        multinomial_coeffecient = np.exp(get_log_multinomial_coeffecient(vec_N))
        p_abn = condition_on_active_nodes(i, item, k, vec_N, T_list, A, num_mask_types)
        pab *= (multinomial_coeffecient * p_abn)
        return pab


def condition_on_mask_types(i, item, is_intermediate, k, T_list, A, m, num_mask_types):
    p_ab = 0
    if is_intermediate:
        n_range = k
    else:
        n_range = k + 1
    
    vec_N = np.zeros(num_mask_types) # N = (N1, N2, ..., NM) [N_0, ..., N_M-1] [0, 0]
    p_ab = recursion_N(i, item, is_intermediate, k, T_list, A, num_mask_types - 1, 0, m, n_range, vec_N, num_mask_types) 
    return p_ab

def condition_on_neighbors(i, item, is_intermediate, mean_degree, T_list, m, vec, k_max, num_mask_types):    
    pa_L = 0
    for k in range(0, k_max):
        prob_r = poisson.pmf(k, mean_degree)
        
        if is_intermediate: # intermediate q using excess degree distribution
            pb = div(prob_r * k, mean_degree)
        else:
            pb = prob_r
            
        if item == 'es':
            p_ab = condition_on_mask_types(i, item, is_intermediate, k, T_list, vec, m, num_mask_types)
        else:
            p_ab = PE_B(i, is_intermediate, k, vec, T_list, m)
            
        pa_L += p_ab * pb
    return pa_L

def get_var_vec(item, mean_degree, is_intermediate, T_list, m, vec, k_max, num_mask_types):
    P_A_list = []
    
    if item not in item_names or item == 'sim':
        print("get_var_vec: item value of %s incorrect!" %item)
        print("Possible options are(except sim):", item_names)
        assert False
        
    for i in range(num_mask_types):
        P_A_list.append(condition_on_neighbors(i, item, is_intermediate, mean_degree, T_list, m, vec, k_max, num_mask_types))


    return np.array(P_A_list)


def func_root(A, mean_degree, T_list, m, k_max, num_mask_types):
    return np.array(get_var_vec('es', mean_degree, True, T_list, m, A, k_max, num_mask_types)) - np.array(A)

def get_EpidemicSize(mean_degree, paras, infection_sizes):
    '''
    S
    '''    
    k_max, T_list = resolve_paras(paras)
    num_mask_types = len(T_list)
    init_A = np.ones(num_mask_types) * 0.9
    m = paras.m

    A_root = optimize.fsolve(func_root, init_A, args=(mean_degree, T_list, m, k_max, num_mask_types))
    P_A_list = get_var_vec('es', mean_degree, False,  T_list, paras.m, A_root, k_max, num_mask_types)
    A = np.dot(P_A_list, m)
    print(mean_degree, P_A_list,)
    
    for i in range(num_mask_types):
        infection_sizes[i][mean_degree] = P_A_list[i]
    infection_sizes['ttl'][mean_degree] = A
    return P_A_list, A

########### Mask Model PE Analysis -- Parellel ########### 
def PE_B(i, is_intermediate, k, E_list, T_list, m):
    res = 0
    
    if is_intermediate: # intermediate q, powers sum up to k - 1
        n_range = k
    else:               # generation 0, powers sum up to k
        n_range = k + 1
        
    for n in range(n_range):
        pe_bn = PE_BN(i, is_intermediate, n, k, E_list, T_list, m)
        res += pe_bn * comb(n_range - 1, n) * (m[0] ** n) * (m[1] ** (n_range - 1 - n))
    return res

def PE_BN(i, is_intermediate, n, k, E_list, T_list, m):
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
            (E_list[0] ** k0) * (E_list[1] ** k1)
            
    return res

def func_root_pe(E_list, mean_degree, T_list, m, k_max, num_mask_types):
    return get_var_vec('pe', mean_degree, True, T_list, m, E_list, k_max, num_mask_types) - np.array(E_list)

def get_ProbEmergence(mean_degree, paras, pes):
    k_max, T_list = resolve_paras(paras)
    num_mask_types = len(T_list)
    init_E = np.ones(num_mask_types) * 0.1
    E_list = optimize.fsolve(func_root_pe, init_E, args=(mean_degree, T_list, paras.m, k_max, num_mask_types), xtol=1e-6)    
    E_list = 1 - get_var_vec('pe', mean_degree, False, T_list, paras.m, E_list, k_max, num_mask_types)
    for i in range(num_mask_types):
        pes[i][mean_degree] = E_list[i]
    pes['ttl'][mean_degree]   = np.dot(E_list, paras.m)
    print(mean_degree, E_list)