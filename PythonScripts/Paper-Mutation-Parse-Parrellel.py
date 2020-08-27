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
from collections import defaultdict 
import multiprocessing, time
from multiprocessing import Manager

"""
usage:

python Paper-Mutation-Parse-Parrellel.py -n 10 -e 10 -t 0.01 -m 0.6 -T 0.6 -tm 0.5
"""

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

def runExp(i, mean_degree,num_nodes, T_list, mutation_probability):
    network = create_network(mean_degree, num_nodes)
    size, size1, size2 = evolution(network, T_list, mutation_probability)
    fractionDic[i] = div(size,num_nodes)
    infectedPerStDic[i] = [div(size1,num_nodes), div(size2,num_nodes)]

def infected_rule(infected_neighbors_dict, t_list, susceptible_nodes, num_strain, mutation_prob):
    new_nodes_list = [set(), set()]
    if len(infected_neighbors_dict.keys()) != 0:
        for node in infected_neighbors_dict:
            trial_list = infected_neighbors_dict[node]
            random.shuffle(trial_list)
            for neighbor in trial_list:
                if random.random() < t_list[neighbor]:
                    susceptible_nodes.remove(node)

                    strain_type_idx = neighbor % num_strain
                    next_strain_type_idx = (neighbor + 1) % num_strain
                    if random.random() < mutation_prob[neighbor]:
                        new_nodes_list[strain_type_idx].add(node)
                    else:
                        new_nodes_list[next_strain_type_idx].add(node)
                    break
    return new_nodes_list

def evolution(g, t_list, mutation_prob):
    g.simplify()
    node_set = set(g.vs.indices)
    num_nodes = len(node_set)
    if start_strain == 1:
        strain_set_1 = set([int(np.random.random_integers(0, num_nodes - 1))])
        strain_set_2 = set()
    elif start_strain == 2:
        strain_set_1 = set()
        strain_set_2 = set([int(np.random.random_integers(0, num_nodes - 1))])
    else:
        exit()
    strain_list = [strain_set_1, strain_set_2]
    num_strain = len(strain_list)

    susceptible_nodes = node_set
    for strain_set in strain_list:
        susceptible_nodes = susceptible_nodes.difference(strain_set)
    new_nodes_list = [strain_set_1, strain_set_2]

    while(sum([len(new_nodes) for new_nodes in new_nodes_list])):
        
        neighbor_dict = collections.defaultdict(list)

        for strain_type, strain_set in enumerate(new_nodes_list):
            strain_neighbors_list = []
            for node in strain_set:
                strain_neighbors_list += g.neighbors(node)
            if len(strain_neighbors_list) == 0: continue
                
            for node in strain_neighbors_list:
                if node not in susceptible_nodes: continue
                neighbor_dict[node].append(strain_type)
        new_nodes_list = infected_rule(neighbor_dict, t_list, susceptible_nodes, num_strain, mutation_prob)

        strain_list = [strain_list[s_idx].union(s) for s_idx, s in enumerate(new_nodes_list)]
    num_infected = sum([len(s) for s in strain_list])
    num_infected1, num_infected2 = map(len, strain_list)
    return num_infected, num_infected1, num_infected2

# def parse_args(args):
#     parser = argparse.ArgumentParser(description = 'Parameters')
#     parser.add_argument('-m', type = float, nargs = '+', default = np.arange(1, 10.1, 0.1), help='np.linspace(0.001, 7, 50) (default); list of mean degree: you can type 1 3 5')
#     parser.add_argument('-n', type = int, default = 200000, help='10,000 (default); the number of nodes')
#     parser.add_argument('-e', type = int, default = 50, help='100 (default); the number of experiments')
#     parser.add_argument('-t1', type = float, default = 0.2, help='0.5 (default); the transmissibility of strain-1')
#     parser.add_argument('-t2', type = float, default = 0.5, help='0.5 (default); the transmissibility of strain-2')
#     parser.add_argument('-m1', type = float, default = 0.9, help='0.5 (default); the mutation probability from 1 to 1')
#     parser.add_argument('-m2', type = float, default = 1.0, help='0.5 (default); the mutation probability from 2 to 2')
#     parser.add_argument('-thrVal', type = float, default = 0.005, help='0.001 (default); the treshold to consider a component giant')
#     parser.add_argument('-numCores', type = int, default = 12, help='number of Cores')
#     parser.add_argument('-logName', default = 'logfile', help='The name of the log file')
#     parser.add_argument('-i', type = int, default = 1, help='1 (default); starting from type-i node')
#     return parser.parse_args(args)

def parse_args(args):
    '''
    Added by Yurun
    '''
    parser = argparse.ArgumentParser(description = 'Parameters')
    parser.add_argument('-n', type = int, default = 200000, help='200,000 (default); the number of nodes')
    parser.add_argument('-e', type = int, default = 50, help='50 (default); the number of experiments')
    parser.add_argument('-m', type = float, default = 0.6, help='0.6 (default); the prob of wearing a mask')
    parser.add_argument('-tm', type = float, default = 0.5, help='0.5 (default); T_mask')
    parser.add_argument('-T', type = float, default = 0.6, help='0.6 (default); transmissibility of the orinigal virus')
    parser.add_argument('-th', type = float, default = 0.001, help='0.001 (default); the treshold to consider a component giant')
    parser.add_argument('-numCores', type = int, default = 12, help='number of Cores')
    return parser.parse_args(args)

def generate_new_transmissibilities_mutation(T_mask, T, m):
    '''
    Added by Yurun
    '''
    roundN = 5 # Round T to roundN digits
    T1 = round(T * T_mask * T_mask * m, roundN)
    T2 = round(T * T_mask * (1 - m), roundN)
    T3 = round(T * (1 - m), roundN)
    T4 = round(T * T_mask * m , roundN)

    Q1 = T1 * (1 - m) + T2 * m
    Q2 = T3 * (1 - m) + T4 * m

    mu11 = T2 * m / Q1
    mu12 = T1 * (1 - m) / Q1
    mu22 = T3 * (1 - m) / Q2
    mu21 = T4 * m / Q2

    trans_dict = {'Q1': Q1,
                  'Q2': Q2}
    
    mu_dict = {'mu11': mu11,
               'mu12': mu12,
               'mu22': mu22,
               'mu21': mu21, }

    print("Q1: %.5f" %Q1)
    print("Q2: %.5f" %Q2)

    print("mu11: %.5f" %mu11)
    print("mu12: %.5f" %mu12)
    print("mu22: %.5f" %mu22)
    print("mu21: %.5f" %mu21)
    return trans_dict, mu_dict


def draw_figures(mean_degree_list, Prob_Emergence, AvgValidSize, AvgSize, ExpPath, infSt1, infSt2, m):
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
    
    plt.legend(["Avg", "Strain-1", "Strain-2"])
    plt.xlabel("Mean Degree")
    plt.ylabel("Prob of Emergence")
    title = "Simulated Probability of Emergence for Mutation Model(separate)"
    plt.title(title)
    plt.savefig(figure_path + '/' + title.replace(" ", "_") + '.png')
    
    ### Probability of Emergence ###    
    plt.figure()
    plt.plot(mean_degree_list, np.array(Prob_Emergence[1]) * m + np.array(Prob_Emergence[2]) * (1 - m), 'go')
    plt.xlabel("Mean Degree")
    plt.ylabel("Prob of Emergence")
    title = "Simulated Probability of Emergence for Mutation Model"
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
    plt.legend(["Avg", "Strain-1", "Strain-2"])
    title = "Simulated Epidemic Size for Mutation Model(separate)\n Seed strain-1"
    plt.title(title)
    plt.savefig(figure_path + '/' + title.replace(" ", "_") + '.png')
    
    ### Epidemic Size  ###
    plt.figure()
    plt.plot(mean_degree_list, AvgValidSize[1], 'go')
    plt.xlabel("Mean Degree")
    plt.ylabel("Epidemic Size")
    title = "Simulated Epidemic Size for Mutation Model\n Seed strain-1"
    plt.title(title)
    plt.savefig(figure_path + '/' + title.replace(" ", "_") + '.png')
   
    
    ### Infected Frac ###
    plt.figure()
    plt.plot(mean_degree_list, AvgSize[1], 'go')
    plt.xlabel("Mean Degree")
    plt.ylabel("Infected Frac")
    title = "Simulated Infected Frac for Mutation Model\n Seed strain-1"
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
    plt.legend(["Avg", "Strain-1", "Strain-2"])
    title = "Simulated Epidemic Size for Mutation Model(separate)\n Seed strain-2"
    plt.title(title)
    plt.savefig(figure_path + '/' + title.replace(" ", "_") + '.png')
    
    ### Epidemic Size  ###
    plt.figure()
    plt.plot(mean_degree_list, infSt1[2], 'go')
    plt.xlabel("Mean Degree")
    plt.ylabel("Epidemic Size")
    title = "Simulated Epidemic Size for Mutation Model\n Seed strain-2"
    plt.title(title)
    plt.savefig(figure_path + '/' + title.replace(" ", "_") + '.png')
   
    
    ### Infected Frac ###
    plt.figure()
    plt.plot(mean_degree_list, AvgSize[2], 'go')
    plt.xlabel("Mean Degree")
    plt.ylabel("Infected Frac")
    title = "Simulated Infected Frac for Mutation Model\n Seed strain-2"
    plt.title(title)
    plt.savefig(figure_path + '/' + title.replace(" ", "_") + '.png')
    
# paras = parse_args(sys.argv[1:])
# mean_degree_list = paras.m
# t1 = paras.t1
# t2 = paras.t2
# m1 = paras.m1
# m2 = paras.m2
# num_nodes = paras.n
# numExp = paras.e
# start_strain = paras.i
# num_cores = min(paras.numCores,multiprocessing.cpu_count())
# thrVal = paras.thrVal

########### Paras & Path preparation ###########
paras = parse_args(sys.argv[1:])
mean_degree_list = np.linspace(0, 10, 50)

nodeN = paras.n
ExpN = paras.e
threshold = paras.th
num_cores = min(paras.numCores,multiprocessing.cpu_count())
mask_prob = paras.m
T_mask = paras.tm
T = paras.T

trans_dict = generate_new_transmissibilities_mutation(T_mask, T, mask_prob)
t1 = trans_dict[0]['Q1'] # Q1
t2 = trans_dict[0]['Q2'] # Q2


m1 = trans_dict[1]['mu11'] # mu11
m2 = trans_dict[1]['mu22'] ## mu22

num_nodes = nodeN
numExp = ExpN


num_cores = min(2,multiprocessing.cpu_count())
thrVal = 0.005

print("Node number:", num_nodes)
print("Exp number:", numExp)

T_list = [t1, t2]
mutation_probability = [m1, m2]
# ff = open("log1"+'Det','w+')
# f = open("log1", 'w+')


Prob_Emergence = defaultdict(list)
AvgValidSize = defaultdict(list)
AvgSize = defaultdict(list)
StdValidSize = defaultdict(list)
infSt1 = defaultdict(list)
infSt2 = defaultdict(list)

now = datetime.now() # current date and time
timeExp = now.strftime("%m%d%H:%M")
print("Exp start at:" + timeExp)


for start_strain in [1, 2]:
    for mean_degree in mean_degree_list:
        a = time.time()
        ttlEpidemicsSize = 0
        numEpidemics = 0
        Epidemics = []
        EpidemicsPerSt = [0,0,0]
        fractionDic = Manager().dict()
        infectedPerStDic = Manager().dict()
        ttlFrac = 0

        Parallel(n_jobs = num_cores)(delayed(runExp)(i, mean_degree,num_nodes, T_list, mutation_probability) 
                                     for i in range(numExp))

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
ExpPath = '../MutationResults/' + timeExp + '_n' + str(num_nodes) + '_e' + str(numExp)

if not os.path.exists(ExpPath):
    os.makedirs(ExpPath)
    
print("Experiment results stored in: ", ExpPath)

draw_figures(mean_degree_list, Prob_Emergence, AvgValidSize, AvgSize, ExpPath, infSt1, infSt2, mask_prob)

setting_path = ExpPath + '/' + 'Settings'
if not os.path.exists(setting_path):
    os.mkdir(setting_path)
 
res_path = ExpPath + '/' + 'Results'
if not os.path.exists(res_path):
    os.mkdir(res_path) 

res_paths1 = res_path+ '/start-s1'
if not os.path.exists(res_paths1):
    os.mkdir(res_paths1)
    
res_paths2 = res_path+ '/start-s2'
if not os.path.exists(res_paths2):
    os.mkdir(res_paths2)

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
np.save(setting_path + '/mu_dict.npy', np.array(mutation_probability))


### Results start from strain-1 ###
np.save(res_paths1 + '/Prob_Emergence.npy', np.array(Prob_Emergence[1])) 
np.save(res_paths1 + '/AvgValidSize.npy', np.array(AvgValidSize[1])) 
np.save(res_paths1 + '/StdValidSize.npy', np.array(StdValidSize[1])) 
np.save(res_paths1 + '/infSt1.npy', np.array(infSt1[1])) 
np.save(res_paths1 + '/infSt2.npy', np.array(infSt2[1])) 

### Results start from strain-2 ###
np.save(res_paths2 + '/Prob_Emergence.npy', np.array(Prob_Emergence[2])) 
np.save(res_paths2 + '/AvgValidSize.npy', np.array(AvgValidSize[2])) 
np.save(res_paths2 + '/StdValidSize.npy', np.array(StdValidSize[2])) 
np.save(res_paths2 + '/infSt1.npy', np.array(infSt1[2])) 
np.save(res_paths2 + '/infSt2.npy', np.array(infSt2[2])) 
    

now = datetime.now() # current date and time
timeExp = now.strftime("%m%d%H:%M")
print("Done! at:" + timeExp)