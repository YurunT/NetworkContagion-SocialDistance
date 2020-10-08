import json
import argparse
import collections
import random
import sys, pdb
import site, os
import matplotlib.pyplot as plt
import numpy as np
import igraph as ig
from datetime import datetime
import time
from joblib import Parallel, delayed
import multiprocessing, time
from multiprocessing import Manager
from collections import defaultdict 
import ray

ray.init()

# Modified as least as possible from Rashad's code

# Last update: sep 14th: remove all the global variable

# This works as good as the "Parellel" version

def create_network(mean_degree, num_nodes):
    degree_sequence = np.random.poisson(mean_degree, num_nodes)
    while (np.sum(degree_sequence) % 2 !=0):
        degree_sequence = np.random.poisson(mean_degree, num_nodes)
    return ig.Graph.Degree_Sequence(list(degree_sequence))

def div(x, y):
    if y == 0:
        return 0
    else:
        return x*1.0/y
    
@ray.remote
def runExp(i, mean_degree, num_nodes, T_list, mask_prob, start_strain): #### should add start_strain
    network = create_network(mean_degree, num_nodes)
    size, size1, size2 = evolution(network, T_list, mask_prob, start_strain)
#     fractionDic[i] = div(size,num_nodes)
#     infectedPerStDic[i] = [div(size1,num_nodes), div(size2,num_nodes)]
    return div(size,num_nodes), [div(size1,num_nodes), div(size2,num_nodes)]

def infected_rule(infected_neighbors_dict, T_list, susceptible_nodes, num_strain, mask_prob):
    """
    Changes:
    1. Change paramenter mutation_prob to mask_prob
    2. Add one line in line 48 to 61
    3. Comment line 66 to 72
    4. Add line 73
    """
    new_nodes_list = [set(), set()]
    if len(infected_neighbors_dict.keys()) != 0:
        for node in infected_neighbors_dict:
            trial_list = infected_neighbors_dict[node] # All of his parents
            random.shuffle(trial_list)
            for neighbor in trial_list:
                
                #### Change 1: see if he wears a mask ####
                if random.random() < mask_prob: # susceptible wears a mask
                    strain_type_idx = 0
                    if neighbor == 0: # if neighbor also wears a mask
                        T = T_list[1]
                    else:
                        T = T_list[3]
                else:
                    strain_type_idx = 1
                    if neighbor == 0: # if neighbor also wears a mask
                        T = T_list[0]
                    else:
                        T = T_list[2]

                if random.random() < T:
                    susceptible_nodes.remove(node)
                    new_nodes_list[strain_type_idx].add(node)
                    
                    break
    return new_nodes_list

def evolution(g, t_list, mask_prob, start_strain):
    """
    Changes:
    1. Change paramenter mutation_prob to mask_prob
    """
    g.simplify()
    node_set = set(g.vs.indices)
    num_nodes = len(node_set)
    if start_strain == 1: # Select one node who wears a mask to be the seed(active)
        strain_set_1 = set([int(np.random.randint(0, num_nodes - 1))])
        strain_set_2 = set()
    elif start_strain == 2:
        strain_set_1 = set()
        strain_set_2 = set([int(np.random.randint(0, num_nodes - 1))])
    else:
        exit()

    strain_list = [strain_set_1, strain_set_2] # strain_set_1: active and wear masks; 
    num_strain = len(strain_list) # Total # nodes who are active

    susceptible_nodes = node_set
    
    for strain_set in strain_list:
        susceptible_nodes = susceptible_nodes.difference(strain_set)
    new_nodes_list = [strain_set_1, strain_set_2] # level L - 1

    while(sum([len(new_nodes) for new_nodes in new_nodes_list])):
        neighbor_dict = collections.defaultdict(list) # susceptible nodes in level L, its parents are in the list

        for strain_type, strain_set in enumerate(new_nodes_list): # string type == 0: wear mask
            strain_neighbors_list = []
            for node in strain_set:
                strain_neighbors_list += g.neighbors(node)
            if len(strain_neighbors_list) == 0: continue
            for node in strain_neighbors_list:
                if node not in susceptible_nodes: continue
                neighbor_dict[node].append(strain_type) 
                
        new_nodes_list = infected_rule(neighbor_dict, t_list, susceptible_nodes, num_strain, mask_prob) # Get next level

        strain_list = [strain_list[s_idx].union(s) for s_idx, s in enumerate(new_nodes_list)]
    num_infected = sum([len(s) for s in strain_list])
    num_infected1, num_infected2 = map(len, strain_list)
    return num_infected, num_infected1, num_infected2


"""
Everything changed/added by Yurun
"""
def parse_args(args):
    parser = argparse.ArgumentParser(description = 'Parameters')
    parser.add_argument('-n', type = int, default = 200000, help='200,000 (default); the number of nodes')
    parser.add_argument('-e', type = int, default = 50, help='50 (default); the number of experiments')
    parser.add_argument('-m', type = float, default = 0.6, help='0.6 (default); the prob of wearing a mask')
    parser.add_argument('-tm1', type = float, default = 0.5, help='0.5 (default); T_mask1')
    parser.add_argument('-tm2', type = float, default = 0.5, help='0.5 (default); T_mask2')
    parser.add_argument('-T', type = float, default = 0.6, help='0.6 (default); transmissibility of the orinigal virus')
    parser.add_argument('-th', type = float, default = 0.001, help='0.001 (default); the treshold to consider a component giant')
    parser.add_argument('-nc', type = int, default = 40, help='number of Cores')
    parser.add_argument('-md', type = int, default = 10, help='[0, max_degree]')
    parser.add_argument('-ns', type = int, default = 50, help='Num of sample within [0, max_degree]')
    parser.add_argument('-cp', type = int, default = 500, help='Num of exps to save a checkpoint')
    parser.add_argument('-change', type = int, default = 0, help='Change m(0), T(1), tm(2)')
    return parser.parse_args(args)


def generate_new_transmissibilities_mask(T_mask1, T_mask2, T, m):
    T1 = T * T_mask1
    T2 = T * T_mask1 * T_mask2
    T3 = T
    T4 = T * T_mask2

    trans_dict = {'T1': T1,
                  'T2': T2,
                  'T3': T3,
                  'T4': T4}

    print("T1: %.5f" %T1)
    print("T2: %.5f" %T2)
    print("T3: %.5f" %T3)
    print("T4: %.5f" %T4)
    
    return trans_dict    


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

    
def write_results(results, start_strain, mean_degree, cp, timeExp, check_point, thrVal, change, mask_prob, T, T_mask1, T_mask2, num_nodes, numExp, degree_max, num_samples, mean_degree_list, T_list, start_time):
#     high_level_length = len(mean_degree_list) * check_point ###
#     for start_strain in [1, 2]:
#         loop_1_length = (start_strain - 1) * high_level_length ###
#           for idx, mean_degree in enumerate(mean_degree_list):
#     loop_2_length = idx * check_point ###
    Prob_Emergence = defaultdict(list)
    AvgValidSize = defaultdict(list)
    AvgSize = defaultdict(list)
    StdValidSize = defaultdict(list)
    infSt1 = defaultdict(list)
    infSt2 = defaultdict(list)
    
    ttlEpidemicsSize = 0
    numEpidemics_1 = 0
    numEpidemics_2 = 0
    numEpidemics = 0
    Epidemics = []
    EpidemicsPerSt = [0,0,0]
    fractionDic = dict()
    infectedPerStDic = dict()
    ttlFrac = 0

    

    for ii in range(check_point):
        fractionDic[ii] = results[ii][0] ###
        infectedPerStDic[ii] = results[ii][1] ###
        if fractionDic[ii] >= thrVal:
            numEpidemics += 1
            ttlEpidemicsSize += fractionDic[ii]
            Epidemics.append(fractionDic[ii])
            EpidemicsPerSt[0] += infectedPerStDic[ii][0]
            EpidemicsPerSt[1] += infectedPerStDic[ii][1]

        ttlFrac += fractionDic[ii]

    if len(Epidemics) == 0:
        Epidemics.append(0)


    ######### Record the results for this Mean Degree ##########    
    Prob_Emergence[start_strain].append(numEpidemics*1.0/(check_point))
    AvgValidSize[start_strain].append(div(ttlEpidemicsSize*1.0, numEpidemics))
    AvgSize[start_strain].append(ttlFrac*1.0/check_point)
    StdValidSize[start_strain].append(np.std(Epidemics))
    infSt1[start_strain].append(div(EpidemicsPerSt[0],numEpidemics))
    infSt2[start_strain].append(div(EpidemicsPerSt[1],numEpidemics))

    ######### Save the results for all Mean Degrees ########## 
    if change == 0:
        change_folder = 'change_m'
    elif change == 1:
        change_folder = 'change_T'
    else:
        change_folder = 'change_tm'
        
    ExpPath = 'Mask2Results_'+ change_folder +'/' + 'm' + str(mask_prob) + '_T' + "{0:.2f}".format(T) + '_tm1_' + "{0:.2f}".format(T_mask1) + '_tm2_' + "{0:.2f}".format(T_mask2) + '/'  + 'n' + str(num_nodes) + '_totalexp' + str(numExp) + '/' + timeExp +'/ss'+ str(start_strain) + '/meandegree'+ str(mean_degree) +'/e'+str(check_point) +'_cp' + str(cp)

    if not os.path.exists(ExpPath):
        os.makedirs(ExpPath)
    print("Experiment results stored in: ", ExpPath)

#     draw_figures(mean_degree_list, Prob_Emergence, AvgValidSize, AvgSize, ExpPath, mask_prob)


    setting_path = ExpPath + '/' + 'Settings'
    if not os.path.exists(setting_path):
        os.mkdir(setting_path)

    res_path = ExpPath + '/' + 'Results'
    if not os.path.exists(res_path):
        os.mkdir(res_path) 

#     res_paths1 = res_path+ '/start-mask'
#     if not os.path.exists(res_paths1):
#         os.mkdir(res_paths1)

#     res_paths2 = res_path+ '/start-nomask'
#     if not os.path.exists(res_paths2):
#         os.mkdir(res_paths2)

    ### Experiment Parameters ###
    paras = dict()
    paras['e'] = check_point
    paras['n'] = num_nodes
    paras['th'] = thrVal
    paras['tm1'] = T_mask1
    paras['tm2'] = T_mask2
    paras['m'] = mask_prob
    paras['T'] = T
    paras['md'] = degree_max
    paras['ns'] = num_samples
    paras['meandegree'] = mean_degree
    paras['start_strain '] = start_strain
    paras['check_point'] = cp
    
    




    with open(setting_path + '/paras.json', 'w') as fp:
        json.dump(paras, fp)

    ### Degree list ###
    np.save(setting_path + '/mean_degree_list.npy', np.array(mean_degree_list)) 

    ### Transmissibilites and mutation probs for Mutation Model ###
    np.save(setting_path + '/trans_dict_mu.npy', np.array(T_list)) 


    ### Results start from mask ###
    # Processed data #
    np.save(res_path + '/Prob_Emergence.npy', np.array(Prob_Emergence[start_strain])) 
    np.save(res_path + '/AvgValidSize.npy', np.array(AvgValidSize[start_strain])) 
    np.save(res_path + '/StdValidSize.npy', np.array(StdValidSize[start_strain])) 
    np.save(res_path + '/infSt1.npy', np.array(infSt1[start_strain])) 
    np.save(res_path + '/infSt2.npy', np.array(infSt2[start_strain])) 
    
    # Raw data #
    with open(res_path + '/results.json', 'w') as fp:
        json.dump(results, fp)

    
    now_finish = datetime.now() # current date and time
    timeExp = now_finish.strftime("%m%d%H:%M")
    print("checkpoint %d Done! at: %s" %(cp, timeExp))
    print("--- %.2s seconds ---" % (time.time() - start_time))



def main():
    ########### Paras & Path preparation ###########
    paras = parse_args(sys.argv[1:])

    num_nodes = paras.n
    numExp = paras.e
    thrVal = paras.th
    num_cores = min(paras.nc,multiprocessing.cpu_count())
    # num_cores = multiprocessing.cpu_count()
    mask_prob = paras.m
    T_mask1 = paras.tm1
    T_mask2 = paras.tm2
    T = paras.T
    degree_max = paras.md
    num_samples = paras.ns
    check_point = paras.cp
    change = paras.change

    mean_degree_list = np.linspace(0, degree_max, num_samples)


    print("Node number:", num_nodes)
    print("Exp number:", numExp)
    print("Num of cores used: ", num_cores)
    print("mean_degree_list: [0, %d], num = %d" %(degree_max, num_samples))



    trans_dict = generate_new_transmissibilities_mask(T_mask1, T_mask2, T, mask_prob)

    t1 = trans_dict['T1']
    t2 = trans_dict['T2']
    t3 = trans_dict['T3']
    t4 = trans_dict['T4']

    T_list = [t1, t2, t3, t4]


    ############ Start Exp ############
    now = datetime.now() # current date and time
    start_time = time.time()
    timeExp = now.strftime("%m%d%H:%M")
    print("Exp start at:" + timeExp)


    for start_strain in [1, 2]:
        for mean_degree in mean_degree_list:
            for cp in range(1, int(numExp/check_point) + 1): 
                results_ids = []
                for i in range(check_point):
                    results_ids.append(runExp.remote(i, mean_degree, num_nodes, T_list, mask_prob, start_strain,))  
                results = ray.get(results_ids)
                print(len(results))
                write_results(results, start_strain, mean_degree, cp, timeExp, check_point, thrVal, change, mask_prob, T, T_mask1, T_mask2, num_nodes, numExp, degree_max, num_samples, mean_degree_list, T_list, start_time)


    now_finish = datetime.now() # current date and time
    timeExp = now_finish.strftime("%m%d%H:%M")
    print("Done! at:" + timeExp)
    print("--- %.2s seconds ---" % (time.time() - start_time))
    
    
main()
