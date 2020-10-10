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


ray.init()

def div(x, y):
    if y == 0:
        return 0
    else:
        return x*1.0/y

def create_network(mean_degree, num_nodes):
    degree_sequence = np.random.poisson(mean_degree, num_nodes)
    while (np.sum(degree_sequence) % 2 !=0):
        degree_sequence = np.random.poisson(mean_degree, num_nodes)
    return ig.Graph.Degree_Sequence(list(degree_sequence))

def infected_rule(infected_neighbors_dict, T_list, susceptible_nodes, num_strain, mask_prob, mask_status):
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
        
    seed = int(np.random.randint(0, num_nodes - 1))   
    if start_strain == 1: # Select one node who wears a mask to be the seed(active) 
        while mask_status[seed] != 0:
            seed = int(np.random.randint(0, num_nodes - 1))
        strain_set_1 = set([seed])
        strain_set_2 = set()
    elif start_strain == 2:
        while mask_status[seed] != 1:
            seed = int(np.random.randint(0, num_nodes - 1))
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
                
        new_nodes_list, susceptible_nodes = infected_rule(neighbor_dict, T_list, susceptible_nodes, num_strain, mask_prob, mask_status) # Get next level

        strain_list = [strain_list[s_idx].union(s) for s_idx, s in enumerate(new_nodes_list)]
    num_infected = sum([len(s) for s in strain_list])
    num_infected1, num_infected2 = map(len, strain_list)
    return num_infected, num_infected1, num_infected2

@ray.remote
def runExp(i, mean_degree, num_nodes, T_list, mask_prob, start_strain): #### should add start_strain
    network = create_network(mean_degree, num_nodes)
    size, size1, size2 = evolution(network, T_list, mask_prob, start_strain)
    return div(size,num_nodes), [div(size1,num_nodes), div(size2,num_nodes)]

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
            for cp in range(1, int(paras.e/paras.cp) + 1): # cp order
                results_ids = []
                for i in range(paras.cp):
                    results_ids.append(runExp.remote(i, mean_degree, paras.n, T_list, paras.m, start_strain,))  
                results = ray.get(results_ids)
                write_cp_raw_results(results, start_strain, mean_degree, cp, time_exp, start_time, paras,)
    
    write_exp_settings(time_exp, paras,)


    now_finish = datetime.now() # current date and time
    print("All Done! for:" + time_exp)
    print("--- %.2s seconds in total ---" % (time.time() - start_time))
main()