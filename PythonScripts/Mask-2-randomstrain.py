import json
import argparse
import collections
import random
import sys, pdb
import site, os
# site.addsitedir('/afs/ece.cmu.edu/usr/reletreb/dep/lib/python2.7/site-packages')
import matplotlib.pyplot as plt
import numpy as np
import igraph as ig
from datetime import datetime
from joblib import Parallel, delayed
import multiprocessing, time
from multiprocessing import Manager
from collections import defaultdict 

# Modified as least as possible from Rashad's code

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

def runExp(i, mean_degree,num_nodes, T_list, mask_prob):
    ### Decide start_strain ###
    if random.random() < mask_prob: # seed wears a mask
        start_strain = 1
    else:
        start_strain = 2
        
    network = create_network(mean_degree, num_nodes)
    size, size1, size2 = evolution(network, T_list, mask_prob)
    fractionDic[i] = div(size,num_nodes)
    infectedPerStDic[i] = [div(size1,num_nodes), div(size2,num_nodes)]
    startStrainDic[i] = start_strain

def infected_rule(infected_neighbors_dict, t_list, susceptible_nodes, num_strain, mask_prob):
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

#                     strain_type_idx = neighbor % num_strain
#                     next_strain_type_idx = (neighbor + 1) % num_strain
#                     if random.random() < mutation_prob[neighbor]:
#                         new_nodes_list[strain_type_idx].add(node)
#                     else:
#                         new_nodes_list[next_strain_type_idx].add(node)
                    new_nodes_list[strain_type_idx].add(node)
                    
                    break
    return new_nodes_list

def evolution(g, t_list, mask_prob):
    """
    Changes:
    1. Change paramenter mutation_prob to mask_prob
    """
    g.simplify()
    node_set = set(g.vs.indices)
    num_nodes = len(node_set)
    if start_strain == 1: # Select one node who wears a mask to be the seed(active)
        strain_set_1 = set([int(np.random.random_integers(0, num_nodes - 1))])
        strain_set_2 = set()
    elif start_strain == 2:
        strain_set_1 = set()
        strain_set_2 = set([int(np.random.random_integers(0, num_nodes - 1))])
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
    parser.add_argument('-tm', type = float, default = 0.5, help='0.5 (default); T_mask')
    parser.add_argument('-T', type = float, default = 0.6, help='0.6 (default); transmissibility of the orinigal virus')
    parser.add_argument('-th', type = float, default = 0.001, help='0.001 (default); the treshold to consider a component giant')
    parser.add_argument('-numCores', type = int, default = 12, help='number of Cores')
    return parser.parse_args(args)


def generate_new_transmissibilities_mask(T_mask, T, m):
    """
    Added by Yurun
    """
    roundN = 5 # Round T to roundN digits
    T1 = round(T * T_mask * T_mask * m, roundN)
    T2 = round(T * T_mask * (1 - m), roundN)
    T3 = round(T * (1 - m), roundN)
    T4 = round(T * T_mask * m , roundN)

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
    
    plt.plot(mean_degree_list, Prob_Emergence['avg'], 'yo')
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
    plt.plot(mean_degree_list, Prob_Emergence['avg'], 'go')
    plt.xlabel("Mean Degree")
    plt.ylabel("Prob of Emergence")
    title = "Simulated Probability of Emergence for Mask Model"
    plt.title(title)
    plt.savefig(figure_path + '/' + title.replace(" ", "_") + '.png')
    
    
    ### Epidemic Size for 2 strains ###
    plt.figure()
    plt.plot(mean_degree_list, AvgValidSize, 'yo')
    plt.plot(mean_degree_list, infSt1, 'bo')
    plt.plot(mean_degree_list, infSt2, 'ro')
    plt.xlabel("Mean Degree")
    plt.ylabel("Epidemic Size")
    plt.legend(["Avg", "Mask", "No mask"])
    title = "Simulated Epidemic Size for Mask Model(separate)"
    plt.title(title)
    plt.savefig(figure_path + '/' + title.replace(" ", "_") + '.png')
    
    ### Epidemic Size  ###
    plt.figure()
    plt.plot(mean_degree_list, AvgValidSize, 'go')
    plt.xlabel("Mean Degree")
    plt.ylabel("Epidemic Size")
    title = "Simulated Epidemic Size for Mask Model"
    plt.title(title)
    plt.savefig(figure_path + '/' + title.replace(" ", "_") + '.png')
   
    
    ### Infected Frac ###
    plt.figure()
    plt.plot(mean_degree_list, AvgSize, 'go')
    plt.xlabel("Mean Degree")
    plt.ylabel("Infected Frac")
    title = "Simulated Infected Frac for Mask Model"
    plt.title(title)
    plt.savefig(figure_path + '/' + title.replace(" ", "_") + '.png')
    

    



########### Paras & Path preparation ###########
paras = parse_args(sys.argv[1:])
mean_degree_list = np.linspace(0, 10, 50)


num_nodes = paras.n
numExp = paras.e
thrVal = paras.th
num_cores = min(paras.numCores,multiprocessing.cpu_count())
mask_prob = paras.m
T_mask = paras.tm
T = paras.T
start_strain = 1
# start_strain = 1


print("Node number:", num_nodes)
print("Exp number:", numExp)

trans_dict = generate_new_transmissibilities_mask(T_mask, T, mask_prob)
t1 = trans_dict['T1']
t2 = trans_dict['T2']
t3 = trans_dict['T3']
t4 = trans_dict['T4']

T_list = [t1, t2, t3, t4]




############ Start Exp ############
Prob_Emergence = defaultdict(list)
AvgValidSize = []
AvgSize = []
StdValidSize = []
infSt1 = []
infSt2 = []

now = datetime.now() # current date and time
timeExp = now.strftime("%m%d%H:%M")
print("Exp start at:" + timeExp)




for mean_degree in mean_degree_list:
    a = time.time()
    ttlEpidemicsSize = 0
    numEpidemics_1 = 0
    numEpidemics_2 = 0
    numEpidemics = 0
    Epidemics = []
    EpidemicsPerSt = [0,0,0]
    fractionDic = Manager().dict()
    infectedPerStDic = Manager().dict()
    startStrainDic = Manager().dict() # Added
    ttlFrac = 0


    Parallel(n_jobs = num_cores)(delayed(runExp)(i, mean_degree,num_nodes, T_list, mask_prob) for i in range(numExp))


    for ii in range(numExp):
        if fractionDic[ii] >= thrVal:
            numEpidemics += 1
            '''Added'''
            if startStrainDic[ii] == 1: 
                numEpidemics_1 += 1
            else:
                numEpidemics_2 += 1

            ttlEpidemicsSize += fractionDic[ii]
            Epidemics.append(fractionDic[ii])
            EpidemicsPerSt[0] += infectedPerStDic[ii][0]
            EpidemicsPerSt[1] += infectedPerStDic[ii][1]

        ttlFrac += fractionDic[ii]

    if len(Epidemics) == 0:
        Epidemics.append(0)


    ######### Record the results for this Mean Degree ##########    
    Prob_Emergence['avg'].append(div(numEpidemics, numExp))
    Prob_Emergence[1].append(div(numEpidemics_1, numExp))
    Prob_Emergence[2].append(div(numEpidemics_2, numExp))

    AvgValidSize.append(div(ttlEpidemicsSize*1.0, numEpidemics))
    AvgSize.append(ttlFrac*1.0/numExp)
    StdValidSize.append(np.std(Epidemics))
    infSt1.append(div(EpidemicsPerSt[0],numEpidemics))
    infSt2.append(div(EpidemicsPerSt[1],numEpidemics))

######### Save the results for all Mean Degrees ########## 
ExpPath = '../Mask2Results/' + timeExp + '_n' + str(num_nodes) + '_e' + str(numExp)

if not os.path.exists(ExpPath):
    os.makedirs(ExpPath)
print("Experiment results stored in: ", ExpPath)

draw_figures(mean_degree_list, Prob_Emergence, AvgValidSize, AvgSize, ExpPath, mask_prob)



setting_path = ExpPath + '/' + 'Settings'
if not os.path.exists(setting_path):
    os.mkdir(setting_path)
    
res_path = ExpPath + '/' + 'Results'
if not os.path.exists(res_path):
    os.mkdir(res_path) 



### Experiment Parameters ###
paras = dict()
paras['ExpN'] = numExp
paras['nodeN'] = num_nodes
paras['threshold'] = thrVal

with open(setting_path + '/paras.json', 'w') as fp:
    json.dump(paras, fp)
    
### Degree list ###
np.save(setting_path + '/mean_degree_list.npy', np.array(mean_degree_list)) 
    
### Transmissibilites and mutation probs for Mutation Model ###
np.save(setting_path + '/trans_dict_mu.npy', np.array(T_list)) 


### Results start from mask ###
np.save(res_path + '/Prob_Emergence_avg.npy', np.array(Prob_Emergence['avg'])) 
np.save(res_path + '/Prob_Emergence_start1.npy', np.array(Prob_Emergence[1])) 
np.save(res_path + '/Prob_Emergence_start2.npy', np.array(Prob_Emergence[2])) 
np.save(res_path + '/AvgValidSize.npy', np.array(AvgValidSize)) 
np.save(res_path + '/StdValidSize.npy', np.array(StdValidSize)) 
np.save(res_path + '/infSt1.npy', np.array(infSt1)) 
np.save(res_path + '/infSt2.npy', np.array(infSt2)) 


now = datetime.now() # current date and time
timeExp = now.strftime("%m%d%H:%M")
print("Done! at:" + timeExp)