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
    
def process_sim_res(results, paras,):
#     Prob_Emergence = 
#     AvgValidSize = defaultdict(list)
#     AvgSize = defaultdict(list)
#     StdValidSize = defaultdict(list)
#     infSt1 = defaultdict(list)
#     infSt2 = defaultdict(list)
    
    ttlEpidemicsSize = 0
    numEpidemics_1 = 0
    numEpidemics_2 = 0
    numEpidemics = 0
    Epidemics = []
    EpidemicsPerSt = [0,0,0]
    fractionDic = dict()
    infectedPerStDic = dict()
    ttlFrac = 0

    for ii in range(paras.cp):
        fractionDic[ii] = results[ii][0] ###
        infectedPerStDic[ii] = results[ii][1] ###
        if fractionDic[ii] >= paras.th:
            numEpidemics += 1
            ttlEpidemicsSize += fractionDic[ii]
            Epidemics.append(fractionDic[ii])
            EpidemicsPerSt[0] += infectedPerStDic[ii][0]
            EpidemicsPerSt[1] += infectedPerStDic[ii][1]

        ttlFrac += fractionDic[ii]

    if len(Epidemics) == 0:
        Epidemics.append(0)

    ######### Record the results for this Mean Degree ##########    
    Prob_Emergence = (numEpidemics*1.0/(paras.cp))
    AvgValidSize = (div(ttlEpidemicsSize*1.0, numEpidemics))
    AvgSize = ttlFrac*1.0/paras.cp
    StdValidSize = (np.std(Epidemics))
    infSt1 = (div(EpidemicsPerSt[0],numEpidemics))
    infSt2 = (div(EpidemicsPerSt[1],numEpidemics))
    
    return Prob_Emergence, AvgValidSize, AvgSize, StdValidSize, infSt1, infSt2
    
    
def write_cp_raw_results(results, start_strain, mean_degree, cp, time_exp, start_time, paras,):
    ''' Save the checkponint raw results for simulation.
        Analysis code are accelarated by Ray.
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

def write_exp_settings(time_exp, paras, mean_degree_list,):
    setting_path = get_setting_path(paras, time_exp)
    generate_path(setting_path)
    json_save(setting_path + 'paras.json', vars(paras))
    np.save(setting_path + 'mean_degree_list.npy', np.array(mean_degree_list)) 
    return 
    
def draw_figures(mean_degree_list, Prob_Emergence, AvgValidSize, AvgSize, ExpPath, m):
    """
    Added by Yurun
    """
    figure_path = ExpPath + '/' + 'Figures'
    if not os.path.exists(figure_path):
        os.mkdir(figure_path)
    
    ### Probability of Emergence for 2 strains ###    
    plt.figure()
    
    plt.plot(mean_degree_list, np.array(Prob_Emergence[1]) * m + np.array(Prob_Emergence[2]) * (1 - m), 'yo')
    plt.plot(mean_degree_list, Prob_Emergence[1], 'bo')
    plt.plot(mean_degree_list, Prob_Emergence[2], 'ro')
    
    plt.legend(["Avg", "Mask", "No mask"])
    plt.xlabel("Mean Degree")
    plt.ylabel("Prob of Emergence")
    title = "Simulated Probability of Emergence for Mask Model(separate)"
    plt.title(title)
    plt.savefig(figure_path + '/' + title.replace(" ", "_") + '.png')
    
    ### Probability of Emergence ###    
    plt.figure()
    plt.plot(mean_degree_list, np.array(Prob_Emergence[1]) * m + np.array(Prob_Emergence[2]) * (1 - m), 'go')
    plt.xlabel("Mean Degree")
    plt.ylabel("Prob of Emergence")
    title = "Simulated Probability of Emergence for Mask Model"
    plt.title(title)
    plt.savefig(figure_path + '/' + title.replace(" ", "_") + '.png')
    
    
    '''
    Start from strain-1
    '''
    ### Epidemic Size for 2 strains ###
    plt.figure()
    plt.plot(mean_degree_list, AvgValidSize[1], 'yo')
    plt.plot(mean_degree_list, infSt1[1], 'bo')
    plt.plot(mean_degree_list, infSt2[1], 'ro')
    plt.xlabel("Mean Degree")
    plt.ylabel("Epidemic Size")
    plt.legend(["Avg", "Mask", "No mask"])
    title = "Simulated Epidemic Size for Mask Model(separate)\n Seed wears mask"
    plt.title(title)
    plt.savefig(figure_path + '/' + title.replace(" ", "_") + '.png')
    
    ### Epidemic Size  ###
    plt.figure()
    plt.plot(mean_degree_list, AvgValidSize[1], 'go')
    plt.xlabel("Mean Degree")
    plt.ylabel("Epidemic Size")
    title = "Simulated Epidemic Size for Mask Model\n Seed wears mask"
    plt.title(title)
    plt.savefig(figure_path + '/' + title.replace(" ", "_") + '.png')
   
    
    ### Infected Frac ###
    plt.figure()
    plt.plot(mean_degree_list, AvgSize[1], 'go')
    plt.xlabel("Mean Degree")
    plt.ylabel("Infected Frac")
    title = "Simulated Infected Frac for Mask Model\n Seed wears mask"
    plt.title(title)
    plt.savefig(figure_path + '/' + title.replace(" ", "_") + '.png')
    
    '''
    Start from strain-2
    '''
    ### Epidemic Size for 2 strains ###
    plt.figure()
    plt.plot(mean_degree_list, AvgValidSize[2], 'yo')
    plt.plot(mean_degree_list, infSt1[2], 'bo')
    plt.plot(mean_degree_list, infSt2[2], 'ro')
    plt.xlabel("Mean Degree")
    plt.ylabel("Epidemic Size")
    plt.legend(["Avg", "Mask", "No mask"])
    title = "Simulated Epidemic Size for Mask Model(separate)\n Seed no mask"
    plt.title(title)
    plt.savefig(figure_path + '/' + title.replace(" ", "_") + '.png')
    
    ### Epidemic Size  ###
    plt.figure()
    plt.plot(mean_degree_list, infSt1[2], 'go')
    plt.xlabel("Mean Degree")
    plt.ylabel("Epidemic Size")
    title = "Simulated Epidemic Size for Mask Model\n Seed no mask"
    plt.title(title)
    plt.savefig(figure_path + '/' + title.replace(" ", "_") + '.png')
   
    
    ### Infected Frac ###
    plt.figure()
    plt.plot(mean_degree_list, AvgSize[2], 'go')
    plt.xlabel("Mean Degree")
    plt.ylabel("Infected Frac")
    title = "Simulated Infected Frac for Mask Model\n Seed no mask"
    plt.title(title)
    plt.savefig(figure_path + '/' + title.replace(" ", "_") + '.png')