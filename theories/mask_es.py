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
    pa_L = 0
    for k in range(1, k_max):
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

def get_EpidemicSize(mean_degree, paras, infection_size0, infection_size1, infection_size):
    '''
    S
    '''    
    k_max, T_list = resolve_paras(paras)
    
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