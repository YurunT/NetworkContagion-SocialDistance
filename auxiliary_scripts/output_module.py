import json
import matplotlib.pyplot as plt
from collections import defaultdict 
import numpy as np
import os
from datetime import datetime
import time
import sys, os
sys.path.append(os.path.abspath("."))
from global_vars import *

def div(x, y):
    if y == 0:
        return 0
    else:
        return x*1.0/y
    
def get_change_name(change):
    if change == 0:
        change_folder = 'change_m'
    elif change == 1:
        change_folder = 'change_T'
    else:
        change_folder = 'change_tm'
    return change_folder

def generate_path(path_name):
    if not os.path.exists(path_name):
        os.makedirs(path_name)

def json_save(file_name, file):
    with open(file_name, 'w') as fp:
        json.dump(file, fp)
        
def get_para_setting_str(paras):
    para_setting_str = 'm' + str(paras.m) + '_T' + "{0:.2f}".format(paras.T) + '_tm1_' + "{0:.2f}".format(paras.tm1) + '_tm2_' + "{0:.2f}".format(paras.tm2)
    return para_setting_str

def get_common_path(paras):
    change_folder = get_change_name(paras.change)
    common_path = paras.modelname + '/' + paras.itemname + '/' + change_folder + '/' + paras.msg + '/' + get_para_setting_str(paras) + '/'
    return common_path

def get_setting_path(paras, time_exp):
    setting_path = base_path + 'simulation/' + get_common_path(paras) + 'n' + str(paras.n) + '_ttle' + str(paras.e) + '/' + time_exp + '/Settings/' 
    print("Setting_path:", setting_path)
    return setting_path

def get_exp_path(paras, cp, mean_degree, time_exp, start_strain):
    exp_path = base_path + 'simulation/' + get_common_path(paras) + 'n' + str(paras.n) + '_ttle' + str(paras.e) + '/' + time_exp +'/ss'+ str(start_strain) + '/meandegree'+ str(mean_degree) +'/' +'cp' + str(cp)
    print("Experiment path:", exp_path)
    return exp_path

def get_analysis_path(paras, time_analysis,):
    analysis_path = base_path + 'analysis/'   + get_common_path(paras) + time_analysis
    print("Analysis path:", analysis_path)
    return analysis_path
    
def write_analysis_results(paras, infection_size_list, mean_degree_list):
    ''' Save the results for anaysis.
        Analysis code are accelarated by parellel
    '''

    print("Analysis finished! Start wrting json...")
   
    infection_size0 = infection_size_list[0]
    infection_size1 = infection_size_list[1]
    infection_size = infection_size_list[2]
    
    print("from write_analysis_results's infection_size:", infection_size)

    


    ######### Generate paths ########## 
    time_analysis = datetime.now().strftime("%m%d%H:%M")
    analysis_path = get_analysis_path(paras, time_analysis,)
    res_path = analysis_path + '/' + 'Results'
    setting_path = analysis_path + '/' + 'Settings'

    generate_path(analysis_path)
    generate_path(setting_path)
    generate_path(res_path)

    ######### Save results ########## 
    json_save(setting_path + "/paras.json",vars(paras))
    json_save(res_path + "/total.json",    infection_size.copy())
    json_save(res_path + "/withmask.json", infection_size0.copy())
    json_save(res_path + "/nomask.json",   infection_size1.copy())
    np.save(setting_path + '/mean_degree_list.npy', mean_degree_list)
    
    
def write_cp_raw_results(results, start_strain, mean_degree, cp, time_exp, start_time, paras,):
    ''' Save the checkponint raw results for simulation.
        Sim codes are accelarated by Ray.
    '''
    ######### Generate paths ########## 
    ExpPath = get_exp_path(paras, cp, mean_degree, time_exp, start_strain)
    generate_path(ExpPath)

    ######### Save results ########## 
    json_save(ExpPath + '/results.json', results)
    
    ######### Time setting ########## 
    now_finish = datetime.now() # current date and time
    timeExp = now_finish.strftime("%m%d%H:%M")
    print("checkpoint %d Done! for exp %s" %(cp, timeExp))
    print("--- %.2s seconds ---" % (time.time() - start_time))

def write_exp_settings(time_exp, paras, mean_degree_list):
    setting_path = get_setting_path(paras, time_exp)
    generate_path(setting_path)
    json_save(setting_path + 'paras.json', vars(paras))
    np.save(setting_path + 'mean_degree_list.npy', np.array(mean_degree_list))