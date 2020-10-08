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
from joblib import Parallel, delayed
import multiprocessing, time
from multiprocessing import Manager
from collections import defaultdict 

# Modified as least as possible from Rashad's code

# Sep 15 copied from Mask_simulation.py

# Modified it into the wrapped version(no global). Sep 17th: THIS WORKS FOR BOTH ES AND PE!!!!!!

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

def runExp(i, mean_degree,num_nodes, T_list, mask_prob, start_strain, fractionDic, infectedPerStDic):
    network = create_network(mean_degree, num_nodes)
    size, size1, size2 = evolution(network, T_list, mask_prob, start_strain)
    fractionDic[i] = div(size,num_nodes)
    infectedPerStDic[i] = [div(size1,num_nodes), div(size2,num_nodes)]

def infected_rule(infected_neighbors_dict, T_list, susceptible_nodes, num_strain, mask_prob, mask_status):
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

                if neighbor == 0 and mask_status[node] == 0: # both wear a mask
                    T = T_list[1]
                elif neighbor == 0 and mask_status[node] == 1: # parent wears mask, node does not
                    T = T_list[0]
                elif neighbor == 1 and mask_status[node] == 0: # parent no mask, node mask
                    T = T_list[3]
                elif neighbor == 1 and mask_status[node] == 1: # parent no mask, node no mask
                    T = T_list[2]
                else:
                    print('Error in checking mask combinations')
                    assert(False)

                strain_type_idx = mask_status[node]
                if random.random() < T:
                    susceptible_nodes.remove(node)
                    new_nodes_list[strain_type_idx].add(node)
                    
                    break
    return new_nodes_list, susceptible_nodes

def evolution(g, T_list, mask_prob, start_strain):
    """
    Changes:
    1. Change paramenter mutation_prob to mask_prob
    """
    g.simplify()
    node_set = set(g.vs.indices)
    num_nodes = len(node_set)
    mask_status = []
    for i in range(num_nodes):
        coin = 0 if random.random() < mask_prob else 1
        mask_status += [coin] # 0 means mask, 1 means no mask
    
    seed = int(np.random.random_integers(0, num_nodes - 1))
    if start_strain == 1: # Select one node who wears a mask to be the seed(active)
        strain_set_1 = set([seed])
        strain_set_2 = set()
#         mask_status[seed] == 0
    elif start_strain == 2:
        strain_set_1 = set()
        strain_set_2 = set([seed])
#         mask_status[seed] == 1
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
                
        new_nodes_list, susceptible_nodes = infected_rule(neighbor_dict, T_list, susceptible_nodes, num_strain, mask_prob, mask_status) # Get next level

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
    parser.add_argument('-maxd', type = int, default = 10, help='[min_degree, max_degree]')
    parser.add_argument('-mind', type = int, default = 10, help='[min_degree, max_degree]')
    parser.add_argument('-ns', type = int, default = 50, help='Num of sample within [0, max_degree]')
    parser.add_argument('-change', type = int, default = 0, help='Change m(0), T(1), tm(2)')
    return parser.parse_args(args)


   

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
    degree_max = paras.maxd
    degree_min = paras.mind
    num_samples = paras.ns
    change = paras.change

    paras = dict()
    paras['n'] = num_nodes
    paras['e'] = numExp
    paras['m'] = mask_prob
    paras['T'] = T
    paras['tm1'] = T_mask1
    paras['tm2'] = T_mask2
    paras['th'] = thrVal
    paras['md'] = degree_max
    paras['ns'] = num_samples

    print("Node number:", num_nodes)
    print("Exp number:", numExp)
    print(paras)

    mean_degree_list = np.linspace(degree_min, degree_max, num_samples)
#     mean_degree_list = [2.44897959,  2.65306122,  2.85714286,
#         3.06122449,  3.26530612,  3.46938776,  3.67346939,  3.87755102,
#         4.08163265,  4.28571429,  4.48979592,  4.69387755,  4.89795918,
#         5.10204082,  5.30612245,  5.51020408,  5.71428571,  5.91836735,
#         6.12244898,  6.32653061,  6.53061224,  6.73469388,  6.93877551,
#         7.14285714,  7.34693878,  7.55102041,  7.75510204,  7.95918367,
#         8.16326531,  8.36734694,  8.57142857,  8.7755102 ,  8.97959184,
#         9.18367347,  9.3877551 ,  9.59183673,  9.79591837, 10.        ]
    mean_degree_list = [2.6, 2.64, 2.68, 2.72, 2.76, 2.8, 2.84, 2.88, 2.92, 2.96, 3.0]
    print("Mean degree list:",mean_degree_list )

    trans_dict = generate_new_transmissibilities_mask(T_mask1, T_mask2, T, mask_prob)
    t1 = trans_dict['T1']
    t2 = trans_dict['T2']
    t3 = trans_dict['T3']
    t4 = trans_dict['T4']

    T_list = [t1, t2, t3, t4]

    ############ Start Exp ############
    Prob_Emergence = defaultdict(list)
    AvgValidSize = defaultdict(list)
    AvgSize = defaultdict(list)
    StdValidSize = defaultdict(list)
    infSt1 = defaultdict(list)
    infSt2 = defaultdict(list)

    now = datetime.now() # current date and time
    timeExp = now.strftime("%m%d%H:%M")
    print("Exp start at:" + timeExp)

    raw_data = defaultdict(dict)
    
    for start_strain in [1, 2]:
        
        for mean_degree in mean_degree_list:
            print("mean_degree:", mean_degree)
            a = time.time()
            ttlEpidemicsSize = 0
            numEpidemics_1 = 0
            numEpidemics_2 = 0
            numEpidemics = 0
            Epidemics = []
            EpidemicsPerSt = [0,0,0]
            fractionDic = Manager().dict()
            infectedPerStDic = Manager().dict()
            ttlFrac = 0

            Parallel(n_jobs = num_cores)(delayed(runExp)(i, mean_degree,num_nodes, T_list, mask_prob, start_strain, fractionDic, infectedPerStDic) for i in range(numExp))

            raw_data[start_strain][mean_degree] = fractionDic.copy()
            for ii in range(numExp):
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
            Prob_Emergence[start_strain].append(numEpidemics*1.0/(numExp))
            AvgValidSize[start_strain].append(div(ttlEpidemicsSize*1.0, numEpidemics))
            AvgSize[start_strain].append(ttlFrac*1.0/numExp)
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

    # ExpPath = '../Mask2Results_'+ change_folder +'/' + timeExp + '_n' + str(num_nodes) + '_e' + str(numExp)
    ExpPath =  '/mnt/hdd-storage/ytian/ns/simulation/' + 'Mask_Simulation_'+ change_folder +'/' + 'm' + str(mask_prob) + '_T' + "{0:.2f}".format(T) + '_tm1_' + "{0:.2f}".format(T_mask1) + '_tm2_' + "{0:.2f}".format(T_mask2) + '/' + timeExp

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

    res_paths1 = res_path+ '/start-mask'
    if not os.path.exists(res_paths1):
        os.mkdir(res_paths1)

    res_paths2 = res_path+ '/start-nomask'
    if not os.path.exists(res_paths2):
        os.mkdir(res_paths2)

    ### Experiment Parameters ###
    with open(setting_path + '/paras.json', 'w') as fp:
        json.dump(paras, fp)

    ### Degree list ###
    np.save(setting_path + '/mean_degree_list.npy', np.array(mean_degree_list)) 

    ### Transmissibilites and mutation probs for Mutation Model ###
    np.save(setting_path + '/trans_dict_mu.npy', np.array(T_list)) 


    ### Results start from mask ###
    np.save(res_paths1 + '/Prob_Emergence.npy', np.array(Prob_Emergence[1])) 
    np.save(res_paths1 + '/AvgValidSize.npy', np.array(AvgValidSize[1])) 
    np.save(res_paths1 + '/StdValidSize.npy', np.array(StdValidSize[1])) 
    np.save(res_paths1 + '/infSt1.npy', np.array(infSt1[1])) 
    np.save(res_paths1 + '/infSt2.npy', np.array(infSt2[1])) 
    with open(res_paths1 + '/raw_data.json', 'w') as jf:
        json.dump(raw_data[1], jf)

    ### Results start from no mask ###
    np.save(res_paths2 + '/Prob_Emergence.npy', np.array(Prob_Emergence[2])) 
    np.save(res_paths2 + '/AvgValidSize.npy', np.array(AvgValidSize[2])) 
    np.save(res_paths2 + '/StdValidSize.npy', np.array(StdValidSize[2])) 
    np.save(res_paths2 + '/infSt1.npy', np.array(infSt1[2])) 
    np.save(res_paths2 + '/infSt2.npy', np.array(infSt2[2])) 
    with open(res_paths2 + '/raw_data.json', 'w') as jf:
        json.dump(raw_data[2], jf)

    now = datetime.now() # current date and time
    timeExp = now.strftime("%m%d%H:%M")
    print("Done! at:" + timeExp)
main()