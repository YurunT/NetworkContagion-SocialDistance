import re
from os import listdir
from os.path import isfile, join
import json
from operator import mul    # or mul=lambda x,y:x*y
from fractions import Fraction
from functools import reduce
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.stats import poisson
from scipy.special import comb
from scipy import optimize 
import multiprocessing, time
from multiprocessing import Manager
from collections import defaultdict
import collections
import sys, os
sys.path.append(os.path.abspath("auxiliary_scripts/"))

'''
Collections of data loading functions
'''


def div(x, y):
    if y == 0:
        return 0
    else:
        return x*1.0/y
    
def separate_number_chars(s):
    
    res = re.split('([-+]?\d+\.\d+)|([-+]?\d+)', s.strip())
    res_f = [r.strip() for r in res if r is not None and r.strip() != '']
    return res_f

def load_ES_analysis_results(change_metric, paras_setting, model_names, exp_time=''):
    '''
    Load ES Analysis Results
    Mask: parellel
    Mutation: parellel
    '''
    
    change_metrics_names = ['m', 'T', 'tm']

    infection_size = dict()
    for model in model_names:
        infection_size[model] = dict()
        print("Model ", model)
        print("change ", change_metric)
        print(paras_setting)

                
        infection_size[model][change_metric] = dict()
        change_metric_path = 'es_analysis/' + model + '_ES_Analysis_change_' + change_metric
        print(change_metric_path)
        paras_settings = [f for f in listdir(change_metric_path) if 'tm' in f]   
        infection_size[model][change_metric][paras_setting] = dict()

        if exp_time == '':
            exp_times = [f for f in listdir(change_metric_path + '/' + paras_setting) if f != '.ipynb_checkpoints']
            exp_time = max(exp_times) # Get the latest results
            print(exp_time)
            

        json_path = change_metric_path + '/' + paras_setting + '/' + exp_time + '/' + 'Results/'
        json_path_list = []

        json_path_list.append(json_path + 'infection_size0.json')
        json_path_list.append(json_path + 'infection_size1.json')
        json_path_list.append(json_path + 'infection_size.json')

        for i, json_path_i in enumerate(json_path_list):
            with open(json_path_i) as json_file:
                infection_size_mask = json.load(json_file)

            infection_size[model][change_metric][paras_setting][i] = \
            list(collections.OrderedDict(sorted(infection_size_mask.items())).values())

            infection_size[model][change_metric][paras_setting][i] = \
            [float(i) for i in infection_size[model][change_metric][paras_setting][i]]
    return infection_size  


            
def plot_mask_mutation_theory_change(change_metric, mean_degree_list, infection_size, model_names):    
    for paras_setting in infection_size[model_names[0]][change_metric].keys():
        
        paras = paras_setting.split('_')
        m = float(separate_number_chars(paras[0])[1])
        T = float(separate_number_chars(paras[1])[1])
        tm1 = float(paras[3])
        tm2 = float(paras[5])
        print("m: %.2f, T:%.2f, tm1:%.2f, tm2:%.2f" %(m,T,tm1,tm2))

        infection_size_mask = infection_size[model_names[0]][change_metric][paras_setting]
        infection_size_mu   = infection_size[model_names[1]][change_metric][paras_setting]

    
        plt.figure()
#         plt.plot(mean_degree_list, np.array(infection_size_mask[2]), 'r.')
        plt.plot(mean_degree_list, np.array(infection_size_mask[0]) * m + np.array(infection_size_mask[1]) * (1 - m), 'r.')
        plt.plot(mean_degree_list, np.array(infection_size_mask[0]) * m, 'g.')
        plt.plot(mean_degree_list, np.array(infection_size_mask[1]) * (1 - m), 'b.')

#         plt.plot(mean_degree_list, np.array(infection_size_mask[2]) * m, 'g.')
#         plt.plot(mean_degree_list, np.array(infection_size_mask[0]) * (1 - m), 'b.')
# #         plt.plot(mean_degree_list, np.array(infection_size_mask[1]) , 'r.')

        plt.plot(mean_degree_list, np.array(infection_size_mu[2]), 'rx')
        plt.plot(mean_degree_list, np.array(infection_size_mu[0]) , 'gx')
        plt.plot(mean_degree_list, np.array(infection_size_mu[1]) , 'bx')

        plt.legend(["Avg        (Mask model)", "Mask      (Mask model)", "No mask(Mask model)", 
                    "Avg        (Mutation model)", "Strain-1 (Mutation model)", "Strain-2 (Mutation model)"])
        plt.xlabel("Mean Degree")
        plt.ylabel("Infection Fraction")
        title = "Theoretical Epidemic Size\nChange " + change_metric + "\nm:%.2f, T:%.2f, tm1:%.2f, tm2:%.2f"%(m, T, tm1, tm2)
        plt.title(title)         
        
def load_happyfeet_ES_results(base_happy_path, metric='', para_setting='',
                              exp_setting='', sim_model_name='', start_strain=0, from_sirius = True):    
    base_happy_path = data_path + base_happy_path 
    if from_sirius:
        happy_path = base_happy_path + sim_model_name + 'Results_change_' + \
                     metric + '/' + para_setting + '/' +  exp_setting
        exp_times = [f for f in listdir(happy_path) 
                 if f != '.ipynb_checkpoints']
        exp_time = max(exp_times) # Compare the latest results
        print(exp_time)

        if start_strain == 0:
            happy_path = happy_path + '/' + exp_time + '/' + 'ss1/'
            print("Data from ss1")
        else :
            print("Data from ss1")
            happy_path = happy_path + '/' + exp_time + '/' + 'ss2/'

    else:
        happy_path = base_happy_path
    
    mean_dgree_files = [f for f in listdir(happy_path)]
    
    infection_size0 = defaultdict(list)
    infection_size1 = defaultdict(list)
    cp_raw_results      = defaultdict(list)
    
    
    for idx, mdf in enumerate(mean_dgree_files):
        md = float(separate_number_chars(mdf)[1])
        cpfiles = [f for f in listdir(happy_path + mdf)]
        infection_size0_cp_list = []
        infection_size1_cp_list = []
        cp_results_list = []
        
        for cp_idx, cp in enumerate(cpfiles):
            
            if idx == 0 and cp_idx == 0:
                with open (happy_path + mdf + '/' + cp + '/Settings/paras.json') as jf:
                    paras = json.load(jf)
                    
            np_path = happy_path + mdf + '/' + cp + '/Results/'
            
            infection_size0_cp_list.append(np.load(np_path + 'infSt1.npy'))
            infection_size1_cp_list.append(np.load(np_path + 'infSt2.npy'))  
            
            try:
                with open(np_path + 'results.json', 'r') as fp:# Load raw results
                    exp_results = json.load(fp)
                cp_results_list.append(exp_results)
            except:
                print("cp %s doesn't have raw results" %cp)
            
        infection_size0[md] = np.array(infection_size0_cp_list).mean()
        infection_size1[md] = np.array(infection_size1_cp_list).mean()
        cp_raw_results[md] = cp_results_list
    
    infection_size0_list = list(collections.OrderedDict(sorted(infection_size0.items())).values())
    infection_size1_list = list(collections.OrderedDict(sorted(infection_size1.items())).values())
    return infection_size0_list, infection_size1_list, paras, cp_raw_results     

def load_happyfeet_PE_results(base_happy_path, metric='', para_setting='', 
                              exp_setting='', sim_model_name='', start_strain=0, from_sirius = True):    
    base_happy_path = data_path + base_happy_path 
    if from_sirius:
        happy_path = base_happy_path + sim_model_name + 'Results_change_' + \
                     metric + '/' + para_setting + '/' +  exp_setting
        exp_times = [f for f in listdir(happy_path) 
                 if f != '.ipynb_checkpoints']
        exp_time = max(exp_times) # Compare the latest results
        print(exp_time)

        if start_strain == 0:
            happy_path = happy_path + '/' + exp_time + '/' + 'ss1/'
        else:
            happy_path = happy_path + '/' + exp_time + '/' + 'ss2/'
    else:
        happy_path = base_happy_path
    
    mean_dgree_files = [f for f in listdir(happy_path)]
    infection_size0 = defaultdict(list)

    
    for idx, mdf in enumerate(mean_dgree_files):
        md = float(separate_number_chars(mdf)[1])
        cpfiles = [f for f in listdir(happy_path + mdf)]
        infection_size0_cp_list = []
        
        for cp_idx, cp in enumerate(cpfiles):
            
            if idx == 0 and cp_idx == 0:
                with open (happy_path + mdf + '/' + cp + '/Settings/paras.json') as jf:
                    paras = json.load(jf)
                    
            np_path = happy_path + mdf + '/' + cp + '/Results/'
#             print(np_path)
            infection_size0_cp_list.append(np.load(np_path + 'Prob_Emergence.npy'))
        
        infection_size0[md] = np.array(infection_size0_cp_list).mean()

    
    infection_size0_list = list(collections.OrderedDict(sorted(infection_size0.items())).values())

    return infection_size0_list, paras 

# def load_old_path_sim_results(metric, exp_setting):
#     sim_mask_path = '../../../Mask2Results_change_' + metric + '/'
#     exp_times = [f for f in listdir(sim_mask_path) 
#                      if f != '.ipynb_checkpoints' and exp_setting in f]
#     exp_time = max(exp_times) # Compare the latest results
#     np_path = sim_mask_path + exp_time + "/Results/start-mask/"
#     if0 = np.load(np_path + "infSt1.npy")
#     if1 = np.load(np_path + "infSt2.npy")
    
#     with open (sim_mask_path + exp_time + "/Settings/paras.json") as js:
#         paras_sim = json.load(js)
#     return if0, if1, paras_sim

def load_ES_sim_results(change_metric, paras_setting,):
    sim_mask_path = data_path + 'simulation/Mask_Simulation_change_' + change_metric + '/' + paras_setting + '/'
    exp_times = [f for f in listdir(sim_mask_path) 
                     if f != '.ipynb_checkpoints']
    exp_time = max(exp_times) # Compare the latest results
    print(exp_time)
    np_path0 = sim_mask_path + exp_time + "/Results/start-mask/"

    if0 = np.load(np_path0 + "infSt1.npy")
    if1 = np.load(np_path0 + "infSt2.npy")
    
    with open (sim_mask_path + exp_time + "/Settings/paras.json") as js:
        paras_sim = json.load(js)
        
    return if0, if1, paras_sim

def load_PE_sim_results(change_metric, paras_setting,exp_time = ''):
    sim_mask_path = data_path + 'simulation/Mask_Simulation_change_' + change_metric + '/' + paras_setting + '/'
    if exp_time == '':
        exp_times = [f for f in listdir(sim_mask_path) 
                         if f != '.ipynb_checkpoints']
        exp_time = max(exp_times) # Compare the latest results
        print(exp_time)
    np_path0 = sim_mask_path + exp_time + "/Results/start-mask/"
    np_path1 = sim_mask_path + exp_time + "/Results/start-nomask/"

    
    pe0 = np.load(np_path0 + "Prob_Emergence.npy")
    pe1 = np.load(np_path1 + "Prob_Emergence.npy")
    
    with open (sim_mask_path + exp_time + "/Settings/paras.json") as js:
        paras_sim = json.load(js)
        
    return pe0, pe1, paras_sim



def div(x, y):
    if y == 0:
        return 0
    else:
        return x*1.0/y
    
def resolve_rawdata(raw_results, check_point, thrVal):
    '''
    Load raw ES data from Ray
    raw_results from only 1 strain
    '''
    Prob_Emergence = []
    AvgValidSize = []
#     AvgSize = []
#     StdValidSize = []
    infSt1 = []
    infSt2 = []

    ordered_raw = list(collections.OrderedDict(sorted(raw_results.items())).values())
    
    for idx, res in enumerate(ordered_raw):
        num_cps = len(res) 

        ttlEpidemicsSize = 0
        numEpidemics_1 = 0
        numEpidemics_2 = 0
        numEpidemics = 0
        Epidemics = []
        EpidemicsPerSt = [0,0,0]
        fractionDic = dict()
        infectedPerStDic = dict()
        ttlFrac = 0


        for cp in range(num_cps):
            result = res[cp]
            for ii in range(check_point):
                fractionDic[ii] = result[ii][0] ###
                infectedPerStDic[ii] = result[ii][1] ###
                if fractionDic[ii] >= thrVal:
                    numEpidemics += 1
                    ttlEpidemicsSize += fractionDic[ii]
                    Epidemics.append(fractionDic[ii])
                    EpidemicsPerSt[0] += infectedPerStDic[ii][0]
                    EpidemicsPerSt[1] += infectedPerStDic[ii][1]

                ttlFrac += fractionDic[ii]

        Prob_Emergence.append(numEpidemics*1.0/(check_point * num_cps))
        AvgValidSize.append(div(ttlEpidemicsSize*1.0, numEpidemics))
#         AvgSize.append(ttlFrac*1.0/check_point)
#         StdValidSize.append(np.std(Epidemics))
        infSt1.append(div(EpidemicsPerSt[0],numEpidemics))
        infSt2.append(div(EpidemicsPerSt[1],numEpidemics))
        
    return np.array(infSt1), np.array(infSt2), np.array(AvgValidSize), np.array(Prob_Emergence), sorted(raw_results.keys())
    
    
def load_happyfeet_raw_results(base_happy_path, metric='', para_setting='',
                              exp_setting='', sim_model_name='', start_strain=0, exp_time='', from_sirius = True):    
    print('load_happyfeet_raw_results:', para_setting)
    base_happy_path = data_path + base_happy_path 
    if from_sirius:
        happy_path = base_happy_path + sim_model_name + 'Results_change_' + \
                     metric + '/' + para_setting + '/' +  exp_setting
        
        if exp_time == '':
            exp_times = [f for f in listdir(happy_path) 
                     if f != '.ipynb_checkpoints' ]
            exp_time = max(exp_times) # Compare the latest results
            
        print(exp_time)

        if start_strain == 0:
            happy_path = happy_path + '/' + exp_time + '/' + 'ss1/'
            print("Data from ss1")
        else :
            print("Data from ss2")
            happy_path = happy_path + '/' + exp_time + '/' + 'ss2/'

    else:
        happy_path = base_happy_path
    
    mean_dgree_files = [f for f in listdir(happy_path)]
    
#     infection_size0 = defaultdict(list)
#     infection_size1 = defaultdict(list)
    cp_raw_results = defaultdict(list)
    
    
    for idx, mdf in enumerate(mean_dgree_files):
        md = float(separate_number_chars(mdf)[1])
        cpfiles = [f for f in listdir(happy_path + mdf)]
#         infection_size0_cp_list = []
#         infection_size1_cp_list = []
        cp_results_list = []
        
        for cp_idx, cp in enumerate(cpfiles):
            
            if idx == 0 and cp_idx == 0:
                with open (happy_path + mdf + '/' + cp + '/Settings/paras.json') as jf:
                    paras = json.load(jf)
                    
            np_path = happy_path + mdf + '/' + cp + '/Results/'
            
#             infection_size0_cp_list.append(np.load(np_path + 'infSt1.npy'))
#             infection_size1_cp_list.append(np.load(np_path + 'infSt2.npy'))  
            
            try:
                with open(np_path + 'results.json', 'r') as fp:# Load raw results
                    exp_results = json.load(fp)
                cp_results_list.append(exp_results)
            except:
                print("cp %s doesn't have raw results" %cp)
            
#         infection_size0[md] = np.array(infection_size0_cp_list).mean()
#         infection_size1[md] = np.array(infection_size1_cp_list).mean()
        cp_raw_results[md] = cp_results_list
    
#     infection_size0_list = list(collections.OrderedDict(sorted(infection_size0.items())).values())
#     infection_size1_list = list(collections.OrderedDict(sorted(infection_size1.items())).values())
    return cp_raw_results, paras,    


    
def load_all_valid_results(metric='m', 
                           para_setting='m0.45_T0.60_tm1_0.30_tm2_0.70', 
                           base_happy_path = 'happyfeet/', sim_model_names = ['Mask2', ], thrVal=0.05):
    '''
    Not very useful
    '''
    base_happy_path = data_path + base_happy_path 
    sim_model_name = sim_model_names[0]
    
    happy_path = base_happy_path + sim_model_name + 'Results_change_' + \
                     metric + '/' + para_setting + '/'
    print(happy_path)
    exp_setting_files = [f for f in listdir(happy_path) if f != '.ipynb_checkpoints']

    es0 = []
    es1 = []
    es  = []
    pe0 = []
    pe1 = []
    total_exp_es = 0
    total_exp_pe = 0
    for exp_setting in exp_setting_files:
        print(exp_setting)
        exp_num = int(separate_number_chars(exp_setting)[3])
        exp_time_files = [f for f in listdir(happy_path + exp_setting + '/') if f != '.ipynb_checkpoints']
        for exp_time in exp_time_files:
            total_exp_es += exp_num * 2
            total_exp_pe += exp_num 
            print(exp_time)
            
            raw_ss1, paras = load_happyfeet_raw_results(base_happy_path, metric, 
                                                        para_setting, 
                                                        exp_setting, 
                                                        sim_model_name, 
                                                        exp_time=exp_time,
                                                        start_strain=0)
            raw_ss2, paras = load_happyfeet_raw_results(base_happy_path, metric, 
                                                        para_setting, 
                                                        exp_setting, 
                                                        sim_model_name, 
                                                        exp_time=exp_time,
                                                        start_strain=1)
            check_point = paras['e']
            m = paras['m']
            es0_s1, es1_s1, es_s1, pe0_s1 = resolve_rawdata(raw_ss1, check_point, thrVal)
            es0_s2, es1_s2, es_s2, pe1_s2 = resolve_rawdata(raw_ss2, check_point, thrVal)
            
            
            es0.append((es0_s1  + es0_s2)/2)
            es1.append((es1_s1  + es1_s2)/2)
            es.append((es_s1 + es_s2)/2)
            pe0.append(pe0_s1)
            pe1.append(pe1_s2)
    
    es0 = np.array(es0)
    es1 = np.array(es1)
    es  = np.array(es)
    pe0 = np.array(pe0)
    pe1 = np.array(pe1)
    
    
    
    return np.mean(es0, axis = 0), np.mean(es1, axis = 0), np.mean(es, axis = 0),np.mean(pe0, axis = 0), np.mean(pe1, axis = 0), total_exp_es, total_exp_pe
        
        
    

        
def load_sirius_rawdata(change_metric, paras_setting,exp_time=''):
    sim_mask_path = data_path + 'simulation/Mask_Simulation_change_' + change_metric + '/' + paras_setting + '/'
#     sim_mask_path = 'simulation/Mask_Simulation_change_' + change_metric + '/' + paras_setting + '/'
    if exp_time == '':
        exp_times = [f for f in listdir(sim_mask_path) 
                         if f != '.ipynb_checkpoints']
        exp_time = max(exp_times) # Compare the latest results
    
    print(exp_time)
    np_path0 = sim_mask_path + exp_time + "/Results/start-mask/"
    np_path1 = sim_mask_path + exp_time + "/Results/start-nomask/"
    
    with open(np_path0 + 'raw_data.json') as jf:
        raw_0 = json.load(jf)
    with open(np_path1 + 'raw_data.json') as jf:
        raw_1 = json.load(jf)
    
#     pe0 = np.load(np_path0 + "Prob_Emergence.npy")
#     pe1 = np.load(np_path1 + "Prob_Emergence.npy")
    
    with open (sim_mask_path + exp_time + "/Settings/paras.json") as js:
        paras_sim = json.load(js)
        
    return raw_0, raw_1, paras_sim    

def load_dm5_change_m_ES(thrVal0 = 0.1, thrVal1 = 0.1):

    m_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    metric_analysis = 'm'
#     para_setting = 'm0.1_T0.60_tm1_0.30_tm2_0.70'
    sim_model_names = ['Mask2', ]
    model_names_es = ['Mask', ]
    model_names_pe = ['Mask', 'Mutation']
    base_happy_path = 'happyfeet/'
    sim_metric = 'm'
    exp_setting = 'n5000000_totalexp200'
    sim_model_name = sim_model_names[0]
    
    es0_theory_list = []
    es1_theory_list = []
    es_theory_list  = []
    pe0_theory_list = []
    pe1_theory_list = []
    pe_theory_list  = []
    
    
    es0_sim_list = []
    es1_sim_list = []
    es_sim_list  = []
    pe0_sim_list = []
    pe1_sim_list = []
    
    other_para_setting = 'T0.60_tm1_0.30_tm2_0.70'
    
    for m in m_list:
        
        para_setting = 'm' + str(m) + '_' + other_para_setting
        print(para_setting)
        
        raw_ss1, paras = load_happyfeet_raw_results(base_happy_path, sim_metric, 
                                    para_setting, 
                                    exp_setting, 
                                    sim_model_name, 
                                    0)
        raw_ss2, paras = load_happyfeet_raw_results(base_happy_path, sim_metric, 
                                    para_setting, 
                                    exp_setting, 
                                    sim_model_name, 
                                    1)
        check_point = paras['e']
        
        es0, es1, es, pe0, mlss1 = resolve_rawdata(raw_ss1, check_point, thrVal0)
        es0_, es1_, es_, pe1, mlss2 = resolve_rawdata(raw_ss2, check_point, thrVal1)
        
        print(mlss1)
        print(mlss2)
        es0_sim_list.append(es0)
        es1_sim_list.append(es1)
        es_sim_list.append(np.array(es0) * m + np.array(es1) * (1 - m))
        pe0_sim_list.append(pe0)
        pe1_sim_list.append(pe1)
        
        

        
        infection_size = load_ES_analysis_results(metric_analysis, para_setting, model_names_es)
        probability_e  = load_PE_analysis_results(metric_analysis, para_setting, model_names_pe)
        
        es0_theory = np.array(infection_size[model_names[0]][metric_analysis][para_setting][0]) * m
        es1_theory = np.array(infection_size[model_names[0]][metric_analysis][para_setting][1]) * (1 - m)
        
        pe0_theory = np.array(probability_e[model_names[0]][metric_analysis][para_setting][0]) 
        pe1_theory = np.array(probability_e[model_names[0]][metric_analysis][para_setting][1]) 

        
        
        es0_theory_list.append(es0_theory)
        es1_theory_list.append(es1_theory)
        es_theory_list.append(es0_theory + es1_theory)
        
        pe0_theory_list.append(pe0_theory)
        pe1_theory_list.append(pe1_theory)
        pe_theory_list.append(pe0_theory * m + pe1_theory * (1 - m))
#     print('pe0_theory_list:', pe0_theory_list)
#     print('pe1_theory_list:', pe1_theory_list)
#     print(pe_theory_list)
    
#     print(pe0_sim_list)
#     print(pe1_sim_list)
    
        
    print(es0_theory_list)
        
    plt.figure()
    plt.plot(m_list, es0_theory_list, 'g--', )
    plt.plot(m_list, es1_theory_list, 'b--', )
    plt.plot(m_list, es_theory_list, 'r--', )
    
    plt.scatter(m_list, es0_sim_list, facecolors='none', edgecolors='g')
    plt.scatter(m_list, es1_sim_list, facecolors='none', edgecolors='b')
    plt.scatter(m_list, np.array(es0_sim_list) + np.array(es1_sim_list),  facecolors='none', edgecolors='r')
    
    plt.legend([
            "Theory(mask)", 
            "Theory(no mask)",
            "Theory(total)",

           "Sim(mask)", 
            "Sim(no mask)",
            "Sim(total)"
    ],)
    
    plt.xlabel('m')
    plt.ylabel('frac')
    
    title = "Mask model Epidemic Size for different m"  + "\n" + other_para_setting + "\n" + exp_setting + '\nmean_degree = 5'
    # title = "Mutation model Epidemic Size"  + "\n" + 'T1: %.2f, T2: %.2f'%(q_dict['Q1'], q_dict['Q2']) + "\n" 
    plt.title(title)
    
#     plt.figure()
#     plt.plot(m_list, pe0_theory_list, 'g--', )
#     plt.plot(m_list, pe1_theory_list, 'b--', )
#     plt.plot(m_list, pe_theory_list, 'r--', )
    
#     plt.scatter(m_list, np.array(pe0_sim_list) , facecolors='none', edgecolors='g')
#     plt.scatter(m_list, np.array(pe1_sim_list) , facecolors='none', edgecolors='b')
# #     plt.scatter(m_list, np.array(pe0_sim_list) * m + np.array(pe1_sim_list) * (1 - m),  facecolors='none', edgecolors='r')
    
    
    
#     plt.legend([
#             "Theory(mask)", 
#             "Theory(no mask)",
#             "Theory(total)",

#            "Sim(mask)", 
#             "Sim(no mask)",
#             "Sim(total)"
#     ],)
    
#     plt.xlabel('m')
#     plt.ylabel('Probability')
    
#     title = "Mask model Probability of Emergence for different m"  + "\n" + other_para_setting + "\n" + exp_setting + '\nmean_degree = 5'
#     # title = "Mutation model Epidemic Size"  + "\n" + 'T1: %.2f, T2: %.2f'%(q_dict['Q1'], q_dict['Q2']) + "\n" 
#     plt.title(title)

def get_md5_pe(thrVal0 = 0.05):
    metric_analysis = 'm'
    sim_model_names = ['Mask2', ]
    base_happy_path = 'simulation/'
    sim_metric = 'm'
    exp_setting = 'n5000_totalexp1000'
    sim_model_name = sim_model_names[0]
    model_names = ['Mask',]
# para_setting = 'm0.1_T0.60_tm1_0.30_tm2_0.70'
    tran_theory_list = []
    tran_sim_list = []
    
    tran_theory_list0 = []
    tran_sim_list0 = []
    
    tran_theory_list1 = []
    tran_sim_list1 = []
    m_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    for m in m_list:
        # '''Load  raw simulation results for ES(new path)'''

        para_setting = 'm'+ str(m) + '_T0.60_tm1_0.30_tm2_0.70'


        raw_ss1_0925, paras = load_happyfeet_raw_results(base_happy_path, sim_metric, 
                                            para_setting, 
                                            exp_setting, 
                                            sim_model_name, 
                                            start_strain=0, )

        raw_ss2_0925, paras = load_happyfeet_raw_results(base_happy_path, sim_metric, 
                                            para_setting, 
                                            exp_setting, 
                                            sim_model_name, 
                                            start_strain=1, )
        prob_emergence = load_PE_analysis_results(metric_analysis, para_setting, model_names)
        pe0_theory = np.array(prob_emergence[model_names[0]][metric_analysis][para_setting][0]) 
        pe1_theory = np.array(prob_emergence[model_names[0]][metric_analysis][para_setting][1]) 
        pe_theory  = pe0_theory * m + pe1_theory * (1 - m)
        
        check_point = paras['e']
#         check_point = 50
        es0_s1_0923, es1_s1_0923, es_s1_0923, pe0, mlist0 = resolve_rawdata(raw_ss1_0925, check_point, thrVal0)
        es0_s1_0923, es1_s1_0923, es_s1_0923, pe1, mlist1 = resolve_rawdata(raw_ss2_0925, check_point, thrVal0)
        
        pe0_sim = pe0[25]
        pe1_sim = pe1[25]
        pe_sim = pe0 * m + pe1 * (1 - m)
        
        
#         pe_sim = np.array(pe0) * m + np.array(pe1) * (1 - m)
#         pe_sim_right_shift = list(pe_sim[1:])
#         pe_sim_right_shift.append(1)
#         pe_sim_right_shift = np.array(pe_sim_right_shift)
#         pe_sim_res = abs(pe_sim_right_shift - pe_sim)
# #         if m != 0.9:
#         tran_theory = mean_degree_list[np.where(pe_sim_res > 10 ** (-4))[0][0]]
#         tran_sim = mean_degree_list[np.where(pe_sim_res > 0)[0][0]]
        tran_theory_list.append(pe_theory[25])
#         tran_sim_list.append(pe_sim)
        
        tran_theory_list0.append(pe0_theory[25])
        tran_sim_list0.append(pe0_sim)
        
        tran_theory_list1.append(pe1_theory[25])
        tran_sim_list1.append(pe1_sim)
    
    plt.figure()
    plt.plot(m_list, tran_theory_list0, 'g--')
    plt.plot(m_list, tran_theory_list1, 'b--')
    plt.plot(m_list, np.array(tran_theory_list0) * m + np.array(tran_theory_list1) * (1 - m), 'r--')
    
    
    plt.scatter(m_list, tran_sim_list0 , facecolors='none', edgecolors='g')
    plt.scatter(m_list, tran_sim_list1 , facecolors='none', edgecolors='b')
    plt.scatter(m_list, np.array(tran_sim_list0) * m + np.array(tran_sim_list1)* (1 - m), facecolors='none', edgecolors='r')
    
    plt.xlabel('m')
    plt.ylabel('Probability')
    plt.title('Probability of Emergence for Mask Model\nMean degree = 5')
    plt.legend([
                'Theory(mask)',
                'Theory(no-mask)',
                'Theory(total)',
                'Simulation(mask)',
                'Simulation(no-mask)',
                'Simulation(total)',
    ])
    
def get_transition_points(thrVal0 = 0.15):
    metric_analysis = 'm'
    sim_model_names = ['Mask2', ]
    base_happy_path = 'simulation/'
    sim_metric = 'm'
    exp_setting = 'n5000_totalexp1000'
    sim_model_name = sim_model_names[0]
    model_names = ['Mask',]
# para_setting = 'm0.1_T0.60_tm1_0.30_tm2_0.70'
    tran_theory_list = []
    tran_sim_list = []
    m_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    for m in m_list:
        # '''Load  raw simulation results for ES(new path)'''

        para_setting = 'm'+ str(m) + '_T0.60_tm1_0.30_tm2_0.70'


        raw_ss1_0925, paras = load_happyfeet_raw_results(base_happy_path, sim_metric, 
                                            para_setting, 
                                            exp_setting, 
                                            sim_model_name, 
                                            start_strain=0, )

        raw_ss2_0925, paras = load_happyfeet_raw_results(base_happy_path, sim_metric, 
                                            para_setting, 
                                            exp_setting, 
                                            sim_model_name, 
                                            start_strain=1, )
        prob_emergence = load_PE_analysis_results(metric_analysis, para_setting, model_names)
        
        check_point = paras['e']
#         check_point = 50
        es0_s1_0923, es1_s1_0923, es_s1_0923, pe0, mlist0 = resolve_rawdata(raw_ss1_0925, check_point, thrVal0)
        es0_s1_0923, es1_s1_0923, es_s1_0923, pe1, mlist1 = resolve_rawdata(raw_ss2_0925, check_point, thrVal0)
        
        pe_theory = np.array(prob_emergence[model_names[0]][metric_analysis][para_setting][0]) * m + np.array(prob_emergence[model_names[0]][metric_analysis][para_setting][1]) * (1 - m)
        pe_sim = np.array(pe0) * m + np.array(pe1) * (1 - m)
        pe_sim_right_shift = list(pe_sim[1:])
        pe_sim_right_shift.append(1)
        pe_sim_right_shift = np.array(pe_sim_right_shift)
        pe_sim_res = abs(pe_sim_right_shift - pe_sim)
#         if m != 0.9:
        tran_theory = mean_degree_list[np.where(pe_sim_res > 10 ** (-4))[0][0]]
        tran_sim = mean_degree_list[np.where(pe_sim_res > 0)[0][0]]
        tran_theory_list.append(tran_theory)
        tran_sim_list.append(tran_sim)
    
    plt.figure()
    plt.plot(m_list, tran_theory_list, 'r+')
#     plt.plot(m_list, tran_sim_list, 'cx')
    plt.scatter(m_list, tran_sim_list , facecolors='none', edgecolors='c')
    plt.xlabel('m')
    plt.ylabel('Mean degree')
    plt.title('Transition point of Probability for Mask model')
    plt.legend([
                'Theory',
                'Simulation'])

def resolve_sirius_raw(raw_data, thrVal):
    
    Prob_Emergence = []
    AvgValidSize = []
    
    for mean_degree, md_res in raw_data.items():
        
#         print("Mean degree:", mean_degree)
        
        ttlEpidemicsSize = 0
        numEpidemics_1 = 0
        numEpidemics_2 = 0
        numEpidemics = 0
        Epidemics = []
        EpidemicsPerSt = [0,0,0]
        fractionDic = dict()
        infectedPerStDic = dict()
        ttlFrac = 0
        
        fractionDic = md_res
        numExp = len(fractionDic.keys())
        
        for ii in fractionDic.keys():
            if fractionDic[ii] >= thrVal:
                numEpidemics += 1
                ttlEpidemicsSize += fractionDic[ii]
                Epidemics.append(fractionDic[ii])


            ttlFrac += fractionDic[ii]

        if len(Epidemics) == 0:
            Epidemics.append(0)


        ######### Record the results for this Mean Degree ##########    
        Prob_Emergence.append(numEpidemics*1.0/(numExp))
        AvgValidSize.append(div(ttlEpidemicsSize*1.0, numEpidemics))


    return np.array(AvgValidSize), np.array(Prob_Emergence), sorted(raw_data.keys())


degree_min = 0
degree_max = 10 # para md in simulation script
interval_num = 50 # para ns in simulation script
mean_degree_list = np.linspace(degree_min, degree_max, interval_num)

data_path = '/mnt/hdd-storage/ytian/ns/'

def load_dm5_change_ES(thrVal0 = 0.05, thrVal1 = 0.05, metric_analysis = 'T', which_tm=1, 
                       exp_setting = 'n5000000_totalexp200'):

    m_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
#     m_list = [ 0.8, 0.9]
# 'n5000000_totalexp200'
    

    sim_model_names = ['Mask2', ]
    model_names = ['Mask',]
    base_happy_path = 'simulation/'
    sim_metric = metric_analysis
    
    sim_model_name = sim_model_names[0]
    
    es0_theory_list = []
    es1_theory_list = []
    es_theory_list  = []
    
    es0_sim_list = []
    es1_sim_list = []
    es_sim_list  = []
    
    
    pe0_theory_list = []
    pe1_theory_list = []
    pe_theory_list  = []
    
    pe0_sim_list = []
    pe1_sim_list = []
    pe_sim_list  = []
    
    
    for m in m_list:

        print('\n\nchange ' + metric_analysis + '=%.2f'%m )
        
        if metric_analysis == 'm':
            para_setting = ('m' + str(m) + '_' + 'T0.60_tm1_0.30_tm2_0.70')
            
        elif metric_analysis == 'tm':
            if which_tm == 1:
                para_setting = 'm0.45_T0.60_tm1_%.2f_tm2_0.70'%m
                
            else:
                para_setting = 'm0.45_T0.60_tm1_0.30_tm2_%.2f'%m

        elif metric_analysis == 'T':
            para_setting = ('m0.45_T%.2f_tm1_0.30_tm2_0.70'%m)

        else:
            assert False
            print("No such metric!!")
#         print(para_setting)
        
        raw_ss1, paras = load_happyfeet_raw_results(base_happy_path, sim_metric, 
                                    para_setting, 
                                    exp_setting, 
                                    sim_model_name, 
#                                     exp_time = '092801:13',
                                    start_strain=0)
        raw_ss2, paras = load_happyfeet_raw_results(base_happy_path, sim_metric, 
                                    para_setting, 
                                    exp_setting, 
                                    sim_model_name, 
#                                     exp_time = '092801:13',
                                    start_strain=1)
        check_point = paras['e']
        print(paras)
        
        es0, es1, es, pe0, ml = resolve_rawdata(raw_ss1, check_point, thrVal0)
        es0_, es1_, es_, pe1, ml = resolve_rawdata(raw_ss2, check_point, thrVal1)
        
        pe = pe0 * paras['m'] + pe1 * (1 - paras['m'])
        
        print('pe', pe)
        
#         if m < 0.7:
#             es0_sim_list.append(es0[25])
#             es1_sim_list.append(es1[25])
#             es_sim_list.append(es[25])

#             pe0_sim_list.append(pe0[25])
#             pe1_sim_list.append(pe1[25])
#             pe_sim_list.append(pe[25])
#         else:
        es0_sim_list.append(es0)
        es1_sim_list.append(es1)
        es_sim_list.append(es)

        pe0_sim_list.append(pe0)
        pe1_sim_list.append(pe1)
        pe_sim_list.append(pe)
        
        

        print("\nLoad analysis!")
        infection_size = load_ES_analysis_results(metric_analysis, para_setting, model_names)
        prob_emergence = load_PE_analysis_results(metric_analysis, para_setting, model_names, )
        print("Load analysis done!")
        
        if m == 0.9:
            es0_theory = np.array(infection_size[model_names[0]][metric_analysis][para_setting][0]) * 0.45
            es1_theory = np.array(infection_size[model_names[0]][metric_analysis][para_setting][1]) * (1 - 0.45)
            print(es0_theory)
        else:
            es0_theory = np.array(infection_size[model_names[0]][metric_analysis][para_setting][0]) * 0.45
            es1_theory = np.array(infection_size[model_names[0]][metric_analysis][para_setting][1]) * (1 - 0.45)
            print(es0_theory)
        

        es0_theory_list.append(es0_theory)
        es1_theory_list.append(es1_theory)
        es_theory_list.append(es0_theory + es1_theory)

        pe0_theory = np.array(prob_emergence[model_names[0]][metric_analysis][para_setting][0][25]) 
        pe1_theory = np.array(prob_emergence[model_names[0]][metric_analysis][para_setting][1][25]) 
        pe_theory  = pe0_theory * paras['m'] + pe1_theory * (1 - paras['m'])

        pe0_theory_list.append(pe0_theory)
        pe1_theory_list.append(pe1_theory)
        pe_theory_list.append(pe_theory)
      
            
        
#     pe0_sim_list[0] = pe0_theory_list[0]
#     pe1_sim_list[0] = pe1_theory_list[0]
#     pe_sim_list[0] = pe_theory_list[0]
        
#     pe0_sim_list[-1] = pe0_theory_list[-1]
#     pe1_sim_list[-1] = pe1_theory_list[-1]
#     pe_sim_list[-1] = pe_theory_list[-1]
      
    if metric_analysis == 'tm':
        if which_tm == 1:
            
            metric_analysis = 'Tmask1'
        else:
            metric_analysis = 'Tmask2'
        
    plt.figure(figsize=(7, 5))
    plt.plot(m_list, es0_theory_list, 'g--', )
    plt.plot(m_list, es1_theory_list, 'b--', )
    plt.plot(m_list, es_theory_list, 'r--', )
    
    plt.scatter(m_list, es0_sim_list, facecolors='none', edgecolors='g')
    plt.scatter(m_list, es1_sim_list, facecolors='none', edgecolors='b')
    plt.scatter(m_list, np.array(es0_sim_list) + np.array(es1_sim_list),  facecolors='none', edgecolors='r')
    
    plt.legend([
            "Theory(mask)", 
            "Theory(no mask)",
            "Theory(total)",

           "Simulation(mask)", 
            "Simulation(no mask)",
            "Simulation(total)"
    ],)
    
    plt.xlabel(metric_analysis)
    plt.ylabel('Fraction')
    
    title = "Epidemic Size for Mask model\nMean degree = 5" 
    # title = "Mutation model Epidemic Size"  + "\n" + 'T1: %.2f, T2: %.2f'%(q_dict['Q1'], q_dict['Q2']) + "\n" 
    plt.title(title) 
    
    
    plt.figure(figsize=(7, 5))
    plt.plot(m_list, pe0_theory_list, 'g--', )
    plt.plot(m_list, pe1_theory_list, 'b--', )
    plt.plot(m_list, pe_theory_list, 'r--', )
    
    plt.scatter(m_list, pe0_sim_list, facecolors='none', edgecolors='g')
    plt.scatter(m_list, pe1_sim_list, facecolors='none', edgecolors='b')
    plt.scatter(m_list, pe_sim_list,  facecolors='none', edgecolors='r')
    
    plt.legend([
            "Theory(mask)", 
            "Theory(no mask)",
            "Theory(total)",

           "Simulation(mask)", 
            "Simulation(no mask)",
            "Simulation(total)"
    ],)
    
    
            
    plt.xlabel(metric_analysis)
    plt.ylabel('Fraction')
    
    title = "Probability for Mask model\nMean degree = 5" 
    # title = "Mutation model Epidemic Size"  + "\n" + 'T1: %.2f, T2: %.2f'%(q_dict['Q1'], q_dict['Q2']) + "\n" 
    plt.title(title) 