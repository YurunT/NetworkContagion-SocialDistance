import random
import sys, pdb
import site, os
import collections
import numpy as np
import igraph as ig
from datetime import datetime
import time
from joblib import Parallel, delayed
import multiprocessing
from multiprocessing import Manager
import ray
sys.path.append(os.path.abspath("../auxiliary_scripts/"))
from input_module import *
from output_module import *
from tnn import generate_new_transmissibilities_mask
from main_aux import *

# Modified as least as possible from Rashad's code
ray.init()

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

def infected_rule(infected_neighbors_dict, T_list, susceptible_nodes, mask_prob, mask_status):

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

@ray.remote
def runExp(i, mean_degree, num_nodes, T_list, mask_prob):
    network = create_network(mean_degree, num_nodes)
    size = evolution(network, T_list, mask_prob)
    return div(size,num_nodes), [0, 0]

    
def main():
    ########### Get commandline input ###########
    paras = parse_args(sys.argv[1:])
    paras_check(paras)
    num_cores, rho, k_max, T_list, Q_list, mu_list, mean_degree_list = resolve_paras(paras)

    ############ Start Exp ############
    now = datetime.now() # current date and time
    start_time = time.time()
    time_exp = now.strftime("%m%d%H:%M")
    print("-------Exp start at:" + time_exp + '-------')

    for start_strain in [1, 2]:
        for mean_degree in mean_degree_list:
            for cp in range(1, int(paras.e/paras.cp) + 1): 
                results_ids = []
                for i in range(paras.cp):
                    results_ids.append(runExp.remote(i, mean_degree, paras.n, T_list, paras.m))  
                results = ray.get(results_ids)
                write_results(results, start_strain, mean_degree, cp, time_exp, mean_degree_list, T_list, start_time, paras,)


    now_finish = datetime.now() # current date and time
    print("All Done! for:" + time_exp)
    print("--- %.2s seconds in total ---" % (time.time() - start_time))
main()