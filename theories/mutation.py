from __future__ import division
import math
import sys, site, os
import numpy as np
from scipy.optimize import fsolve
from scipy.special import comb
from scipy.stats import poisson
import scipy.optimize
import scipy.misc
from datetime import datetime
sys.path.append(os.path.abspath("../auxiliary_scripts/"))
from main_aux import *

def get_t_mu(Q_list, mu_list):
    t1 = Q_list[0]
    t2 = Q_list[1]
    u_r_11 = mu_list[0]
    u_r_12 = mu_list[1]
    u_r_22 = mu_list[2]
    u_r_21 = mu_list[3]
    
    return t1, t2, u_r_11, u_r_12, u_r_22, u_r_21

########### Mutation Model ES Analysis -- Parellel ########### 
def obtain_val_r_1(v1, v2, t1, t2, u_r_11, u_r_12, u_r_21, u_r_22, lambda_r, max_degree, rho):
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

    return rho + (1 - rho)*val  

def obtain_val_r_2(v1, v2, t1, t2, u_r_11, u_r_12, u_r_21, u_r_22, lambda_r, max_degree, rho):
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

    return rho + (1 - rho)*val 

def equations(p, lambda_r, t1, t2, u_r_11, u_r_12, u_r_21, u_r_22, max_degree, rho):
    v1, v2 = p
    val_r_1 = obtain_val_r_1(v1, v2, t1, t2, u_r_11, u_r_12, u_r_21, u_r_22, lambda_r, max_degree, rho)
    val_r_2 = obtain_val_r_2(v1, v2, t1, t2, u_r_11, u_r_12, u_r_21, u_r_22, lambda_r, max_degree, rho)

    return (v1 - val_r_1, v2 - val_r_2) 

def cascade_size(lambda_r, paras, infection_size0_mu, infection_size1_mu, infection_size_mu):
    rho, max_degree, T_list, Q_list, mu_list = resolve_paras(paras)
    t1, t2, u_r_11, u_r_12, u_r_22, u_r_21 = get_t_mu(Q_list, mu_list) 
    
    h_r_1, h_r_2 = fsolve(equations, (0.9, 0.9), args=(lambda_r, t1, t2, u_r_11, u_r_12, u_r_21, u_r_22, max_degree, rho, ), xtol=1e-10)

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

        H += prob_r * (tmp_val1 + tmp_val2) 
        H1 += prob_r * tmp_val1
        H2 += prob_r * tmp_val2

    print(lambda_r, H, H1, H2)
    infection_size_mu[lambda_r] = H
    infection_size0_mu[lambda_r] = H1
    infection_size1_mu[lambda_r] = H2
    
    return (lambda_r, H, H1, H2)

########### Mutation Model PE Analysis -- Parellel ########### 

def obtain_val_r_1_pe(v1, v2, t1, u_r_11, u_r_12, mean_degree, max_degree):
    val = 0

    for d_r in range(0, max_degree):
        prob_r = poisson.pmf(d_r, mean_degree)
        val += d_r*prob_r*1.0/mean_degree * ((1 - t1 + t1*u_r_11*v1 + t1*u_r_12*v2)**(d_r-1))

    return val

def obtain_val_r_2_pe(v1, v2, t2, u_r_21, u_r_22, mean_degree, max_degree):
    val = 0

    for d_r in range(0, max_degree):
        prob_r = poisson.pmf(d_r, mean_degree)
        val += d_r*prob_r*1.0/mean_degree * ((1 - t2 + t2*u_r_22*v2 + t2*u_r_21*v1)**(d_r-1))

    return val

def equations_pe(p, mean_degree, t_r_1, t_r_2, u_r_11, u_r_12, u_r_21, u_r_22, max_degree,):
    v1, v2 = p
    val_r_1 = obtain_val_r_1_pe(v1, v2, t_r_1, u_r_11, u_r_12, mean_degree, max_degree)
    val_r_2 = obtain_val_r_2_pe(v1, v2, t_r_2, u_r_21, u_r_22, mean_degree, max_degree)

    return (v1 - val_r_1, v2 - val_r_2)


def cascade_prob(mean_degree, paras, pe_0_list, pe_1_list, pe_list):
    rho, max_degree, T_list, Q_list, mu_list = resolve_paras(paras)
    t_r_1, t_r_2, u_r_11, u_r_12, u_r_22, u_r_21 = get_t_mu(Q_list, mu_list) 

    h_r_1, h_r_2 = fsolve(equations_pe, (0.01, 0.01), args=(mean_degree, t_r_1, t_r_2, u_r_11, u_r_12, u_r_21, u_r_22, max_degree,), xtol=1e-6)

    H_1 = 0
    for d_r in range(0, max_degree):
        prob_r = poisson.pmf(d_r, mean_degree)
        H_1 += prob_r*((1 - t_r_1 + t_r_1*u_r_11*h_r_1 + t_r_1*u_r_12*h_r_2)**d_r)

    H_2 = 0
    for d_r in range(0, max_degree):
        prob_r = poisson.pmf(d_r, mean_degree)
        H_2 += prob_r*((1 - t_r_2 + t_r_2*u_r_22*h_r_2 + t_r_2*u_r_21*h_r_1)**d_r)
        
    print(2 - H_2 - H_1, 1 - H_1, 1 - H_2)
    pe_0_list[mean_degree] = 1 - H_1
    pe_1_list[mean_degree] = 1 - H_2
    pe_list[mean_degree] = 2 - H_2 - H_1
    return (1 - H_1, 1 - H_2)