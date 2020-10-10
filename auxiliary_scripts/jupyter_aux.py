import re
import sys, os
from os import listdir
import collections
sys.path.append(os.path.abspath("."))
from global_vars import *
from output_module import *
from para_settings import *
from main_aux import *
from collections import defaultdict

def test_import_aux():
    print("This is jupyter aux")
    
def separate_number_chars(s):
    res = re.split('([-+]?\d+\.\d+)|([-+]?\d+)', s.strip())
    res_f = [r.strip() for r in res if r is not None and r.strip() != '']
    return res_f
    
def json_load(json_path):
    with open(json_path) as json_file:
        item = json.load(json_file)
    return item

def get_time(exp_time, path):
    if exp_time == '':
        exp_times = [f for f in listdir(path) if f != '.ipynb_checkpoints']
        exp_time = max(exp_times) # Get the latest results
        print(exp_time)
    return exp_time

def get_mdl(path):
    mean_degrees = [separate_number_chars(f)[1] for f in listdir(path) if f != '.ipynb_checkpoints']
    return mean_degrees
    
def get_ordered_values_by_key(ordered_dict):
    ordered_values_list = list(collections.OrderedDict(sorted(ordered_dict.items())).values())
    return ordered_values_list

def get_parasobj(m=0.45, T=0.6, tm1=0.3, tm2=0.7, msg='test', modelname='mask', itemname='es', change_metric='m', n=50, e=10, checkpoint='5'):
    paras = para_setting()
    paras.m = m
    paras.T = T
    paras.tm1 = tm1
    paras.tm2 = tm2
    paras.msg = msg
    paras.modelname = modelname
    paras.itemname = itemname
    paras.change = change_metrics_dict[change_metric]
    paras.n = n
    paras.e = e
    paras.cp = checkpoint
    paras_check(paras)
    return paras

def load_sim_settings(paras, time_exp):
    setting_path = get_setting_path(paras, time_exp)
    paras_arg = json_load(setting_path + 'paras.json')
    return paras_arg

def load_sim_raw_results(m=0.45, T=0.6, tm1=0.3, tm2=0.7, msg='test', modelname='mask', change_metric='m', n=50, e=10, checkpoint=5, time_exp='',):
    '''Load Simulation Results'''
  
    # Prepare paras obj
    infection_size = dict()
    itemname='es' # es and pe are using the same script for simulation
    paras = get_parasobj(m, T, tm1, tm2, msg, modelname, itemname, change_metric, n, e, checkpoint)
    paras.print_paras()
    print("time_exp: ", time_exp)
    
    # Prepare paths
    path = base_path + 'simulation/' + get_common_path(paras) + 'n' + str(paras.n) + '_ttle' + str(paras.e) + '/' 
    time_exp = get_time(time_exp, path)
    mdl_path = path + time_exp +'/ss'+ str(1) + '/'
    mean_degree_list = get_mdl(mdl_path)
    
    raw = defaultdict(dict)
    for start_strain in [1, 2]:
        for mean_degree in mean_degree_list:
            raw[start_strain][mean_degree] = dict()
            for cp in range(1, int(paras.e/paras.cp) + 1): 
                exp_path = get_exp_path(paras, cp, mean_degree, time_exp, start_strain)
                raw[start_strain][mean_degree][cp] = json_load(exp_path + '/results.json') # results has paras.cp exp results
    
    paras_arg = load_sim_settings(paras, time_exp)
    
    return raw, paras_arg

def process_sim_res(results, checkpoint, thr):    
    ttlEpidemicsSize = 0
    numEpidemics_1 = 0
    numEpidemics_2 = 0
    numEpidemics = 0
    Epidemics = []
    EpidemicsPerSt = [0,0,0]
    fractionDic = dict()
    infectedPerStDic = dict()
    ttlFrac = 0

    for ii in range(checkpoint):
        fractionDic[ii] = results[ii][0] ###
        infectedPerStDic[ii] = results[ii][1] ###
        if fractionDic[ii] >= thr:
            numEpidemics += 1
            ttlEpidemicsSize += fractionDic[ii]
            Epidemics.append(fractionDic[ii])
            EpidemicsPerSt[0] += infectedPerStDic[ii][0]
            EpidemicsPerSt[1] += infectedPerStDic[ii][1]

        ttlFrac += fractionDic[ii]

    if len(Epidemics) == 0:
        Epidemics.append(0)

    ######### Record the results for this Mean Degree ##########    
    Prob_Emergence = (numEpidemics*1.0/(checkpoint))
    AvgValidSize = (div(ttlEpidemicsSize*1.0, numEpidemics))
    AvgSize = ttlFrac*1.0/checkpoint
    StdValidSize = (np.std(Epidemics))
    infSt1 = (div(EpidemicsPerSt[0],numEpidemics))
    infSt2 = (div(EpidemicsPerSt[1],numEpidemics))
    return Prob_Emergence, AvgValidSize, infSt1, infSt2

def process_raw(raw, paras, thr,):
    processed_res_strains = []
    for start_strain, results in raw.items():
        processed_res = defaultdict(dict)
        for mean_degree, result in results.items():
            Prob_Emergence_list = []
            AvgValidSize_list = []
            infSt1_list = []
            infSt2_list = []
            
            for cp, cp_res in result.items():
                Prob_Emergence, AvgValidSize, infSt1, infSt2 = process_sim_res(cp_res, paras['cp'], thr=thr)
                Prob_Emergence_list.append(Prob_Emergence)
                AvgValidSize_list.append(AvgValidSize)
                infSt1_list.append(infSt1)
                infSt2_list.append(infSt2)
                
            processed_res['pe'][mean_degree] = np.array(Prob_Emergence_list).mean()
            processed_res['es'][mean_degree] = np.array(AvgValidSize_list).mean()
            processed_res['es0'][mean_degree] = np.array(infSt1_list).mean()
            processed_res['es1'][mean_degree] = np.array(infSt2_list).mean()
        
        ordered_res = dict()
        ordered_res['pe']  = get_ordered_values_by_key(processed_res['pe'])
        ordered_res['es']  = get_ordered_values_by_key(processed_res['es'])
        ordered_res['es0'] = get_ordered_values_by_key(processed_res['es0'])
        ordered_res['es1'] = get_ordered_values_by_key(processed_res['es1'])
        
        processed_res_strains.append(ordered_res)
    
    return processed_res_strains
        
                
    
def load_analysis_results(m=0.45, T=0.6, tm1=0.3, tm2=0.7, msg='test', modelname='mask', itemname='pe', change_metric='m', time_analysis=''):
    '''Load Analysis Results'''

    # Prepare paras obj
    infection_size = dict()
    paras = get_parasobj(m, T, tm1, tm2, msg, modelname, itemname, change_metric)
    paras.print_paras()
    print("time_analysis: ", time_analysis)
    
    # Prepare paths
    path = base_path + 'analysis/' + get_common_path(paras)
    time_analysis = get_time(time_analysis, path)
    analysis_path = get_analysis_path(paras, time_analysis,)
    res_path = analysis_path + '/' + 'Results'
    setting_path = analysis_path + '/' + 'Settings'
    
    # Load results
    paras_json = json_load(setting_path + "/paras.json")
    total = json_load(res_path + "/total.json",)
    withmask = json_load(res_path + "/withmask.json")
    nomask = json_load(res_path + "/nomask.json")
    mean_degree_list = np.load(setting_path + '/mean_degree_list.npy')
    return get_ordered_values_by_key(total), get_ordered_values_by_key(withmask), get_ordered_values_by_key(nomask), paras_json, mean_degree_list