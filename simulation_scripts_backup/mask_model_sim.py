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
    network = create_network(mean_degree, num_nodes)
    size = evolution(network, T_list, mask_prob)
    fractionDic[i] = div(size,num_nodes)

def infected_rule(infected_neighbors_dict, t_list, susceptible_nodes, mask_prob, mask_status):

    new_nodes_list = set()
    if len(infected_neighbors_dict.keys()) != 0:
        for node in infected_neighbors_dict:
            parentTypes = infected_neighbors_dict[node] # bunch of zeros and ones
            #print(infected_neighbors_dict[node])
            for maskType in parentTypes:
                if maskType == 1 and mask_status[node] == 1: # both wear a mask
                    T = T_list[1]
                elif maskType == 1 and mask_status[node] == 0: # parent wears mask, node does not
                    T = T_list[0]
                elif maskType == 0 and mask_status[node] == 1: # parent no mask, node mask
                    T = T_list[3]
                elif maskType == 0 and mask_status[node] == 0: # parent no mask, node no mask
                    T = T_list[2]
                else:
                    print('Error in checking mask combinations')
                    assert(False)

                if random.random() < T:
                    susceptible_nodes.remove(node)
                    new_nodes_list.add(node)
                    break


    return new_nodes_list, susceptible_nodes

def evolution(g, t_list, mask_prob):

    g.simplify()
    node_set = set(g.vs.indices)
    num_nodes = len(node_set)
    mask_status = []
    for i in range(num_nodes):
        coin = 1 if random.random() < mask_prob else 0
        mask_status += [coin] # 1 means mask, 0 means no mask
    strain_set = set([int(np.random.random_integers(0, num_nodes - 1))])

    susceptible_nodes = node_set

    susceptible_nodes = susceptible_nodes.difference(strain_set)
    new_nodes_list = strain_set # level L - 1

    while(len(new_nodes_list)):
        neighbor_dict = collections.defaultdict(list) # susceptible nodes in level L, its parents are in the list

        #for strain_type, strain_set in enumerate(new_nodes_list): # string type == 0: wear mask
        strain_neighbors_list = []
        for node in new_nodes_list:
            #strain_neighbors_list += g.neighbors(node)
            for node2 in g.neighbors(node):
                if node2 not in susceptible_nodes: continue
                neighbor_dict[node2].append(mask_status[node])
                
        new_nodes_list, susceptible_nodes = infected_rule(neighbor_dict, t_list, susceptible_nodes, mask_prob, mask_status) # Get next level

        strain_set = strain_set.union(new_nodes_list)

    num_infected = len(strain_set)
    return num_infected


"""
Everything changed/added by Yurun
"""
# def parse_args(args):
#     parser = argparse.ArgumentParser(description = 'Parameters')
#     parser.add_argument('-n', type = int, default = 200000, help='200,000 (default); the number of nodes')
#     parser.add_argument('-e', type = int, default = 50, help='50 (default); the number of experiments')
#     parser.add_argument('-m', type = float, default = 0.6, help='0.6 (default); the prob of wearing a mask')
#     parser.add_argument('-tm', type = float, default = 0.5, help='0.5 (default); T_mask')
#     parser.add_argument('-T', type = float, default = 0.6, help='0.6 (default); transmissibility of the orinigal virus')
#     parser.add_argument('-th', type = float, default = 0.001, help='0.001 (default); the treshold to consider a component giant')
#     parser.add_argument('-numCores', type = int, default = 12, help='number of Cores')
#     return parser.parse_args(args)

def parse_args(args):
    parser = argparse.ArgumentParser(description = 'Parameters')
    parser.add_argument('-n', type = int, default = 100000, help='200,000 (default); the number of nodes')
    parser.add_argument('-e', type = int, default = 50, help='50 (default); the number of experiments')
    parser.add_argument('-m', type = float, default = 0.45, help='0.6 (default); the prob of wearing a mask')
    parser.add_argument('-tm1', type = float, default = 0.5, help='0.5 (default); T_mask1')
    parser.add_argument('-tm2', type = float, default = 0.5, help='0.5 (default); T_mask2')
    parser.add_argument('-T', type = float, default = 0.6, help='0.6 (default); transmissibility of the orinigal virus')
    parser.add_argument('-th', type = float, default = 0.001, help='0.001 (default); the treshold to consider a component giant')
    parser.add_argument('-nc', type = int, default = 4, help='number of Cores')
    parser.add_argument('-md', type = int, default = 10, help='[0, max_degree]')
    parser.add_argument('-ns', type = int, default = 50, help='Num of sample within [0, max_degree]')
    parser.add_argument('-change', type = int, default = 0, help='Change m(0), T(1), tm(2)')
    return parser.parse_args(args)

    
def generate_new_transmissibilities_mask(T_mask1, T_mask2, T, m):
    T1 = T * T_mask1
    T2 = T * T_mask1 * T_mask2
    T3 = T
    T4 = T * T_mask2

    T2 = 0.126
    T1 = 0.18
    T4 = 0.42
    T3 = 0.6

    trans_dict = {'T1': T1,
                  'T2': T2,
                  'T3': T3,
                  'T4': T4}

    print("T1: %.5f" %T1)
    print("T2: %.5f" %T2)
    print("T3: %.5f" %T3)
    print("T4: %.5f" %T4)
    
    return trans_dict    


########### Paras & Path preparation ###########
paras = parse_args(sys.argv[1:])
mean_degree_list = [2.2, 2.4, 2.5, 2.6, 2.64, 2.68, 2.72, 2.76, 2.8, 2.84, 2.88, 2.92, 2.96, 3.0] #np.linspace(0, 10, 50)
print("mead degree list:", mean_degree_list)
# mean_degree = 5



num_nodes = paras.n
numExp = paras.e
thrVal = paras.th
num_cores = min(paras.nc,multiprocessing.cpu_count())
# num_cores = multiprocessing.cpu_count()
mask_prob = paras.m
T_mask1 = paras.tm1
T_mask2 = paras.tm2
T = paras.T
start_strain = 1
degree_max = paras.md
num_samples = paras.ns
change = paras.change
# start_strain = 1


print("Node number:", num_nodes)
print("Exp number:", numExp)

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

# pList = [0.0, 0.1, 0.2, 0.3, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9]
# pList = [0.45]/

pe_list = []
fractionDic_list = []
for mean_degree in mean_degree_list:
    print('Running simulations for mean_degree = ' + str(mean_degree))
    
    fractionDic = [0 for i in range(numExp)]

    # running experiments now
    for i in range(numExp):
        #print('Running simulation ' + str(i))
        runExp(i, mean_degree,num_nodes, T_list, mask_prob)
    fractionDic_list.append(fractionDic)
    # output average epidemic sizes
    av_epidemic_size = 0
    num_epidemics = 0
    for i in range(numExp):
        if fractionDic[i] >= thrVal:
            av_epidemic_size += fractionDic[i]
            num_epidemics += 1
    av_epidemic_size  = div(av_epidemic_size,num_epidemics)
    av_prob_emergence = div(num_epidemics, numExp)
    pe_list.append(av_prob_emergence)
    print('Average size of the epidemic (conditioned on emergence) is ' + str(av_epidemic_size))
    print('Average probability of emergence is ' + str(av_prob_emergence))
np.save('m_%.2f_T%.2f_tm1_%.2f_tm2_%.2f_ani_thr.npy'%(mask_prob, T, T_mask1, T_mask2), np.array(pe_list))
np.save('m_%.2f_T%.2f_tm1_%.2f_tm2_%.2f_ani_thr_raw.npy'%(mask_prob, T, T_mask1, T_mask2), np.array(fractionDic_list))