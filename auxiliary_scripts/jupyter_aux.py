import sys, os
from os import listdir
sys.path.append(os.path.abspath("."))
from global_vars import *
from output_module import *
from para_settings import *
from main_aux import *

def test_import_aux():
    print("This is jupyter aux")
    
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
    
def load_analysis_results(m=0.45, T=0.6, tm1=0.3, tm2=0.7, msg='test', modelname='mask', itemname='pe', change_metric='m', time_analysis=''):
    '''Load Analysis Results'''

    infection_size = dict()
    
    paras = para_setting()
    paras.m = m
    paras.T = T
    paras.tm1 = tm1
    paras.tm2 = tm2
    paras.msg = msg
    paras.modelname = modelname
    paras.itemname = itemname
    paras.change = change_metrics_dict[change_metric]
    paras_check(paras)
    paras.print_paras()
    print("time_analysis: ", time_analysis)
    
    path = base_path + 'analysis/' + get_common_path(paras)
    time_analysis = get_time(time_analysis, path)
    analysis_path = get_analysis_path(paras, time_analysis,)

    res_path = analysis_path + '/' + 'Results'
    setting_path = analysis_path + '/' + 'Settings'
    
    paras_json = json_load(setting_path + "/paras.json")
    total = json_load(res_path + "/total.json",)
    withmask = json_load(res_path + "/withmask.json")
    nomask = (res_path + "/nomask.json")
    mean_degree_list = np.load(setting_path + '/mean_degree_list.npy')
    return total, withmask, nomask, paras_json, mean_degree_list