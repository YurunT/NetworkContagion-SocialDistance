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
    new_nodes_list = [set() for i in range(num_strain)]
    
    if len(infected_neighbors_dict.keys()) != 0:
        
        for node in infected_neighbors_dict:
            trial_list = infected_neighbors_dict[node] # All of his parents
            random.shuffle(trial_list)
            
            for neighbor in trial_list:
                T = T_list[neighbor][mask_status[node]]
                strain_type_idx = mask_status[node]

                if random.random() < T:
                    susceptible_nodes.remove(node)
                    new_nodes_list[strain_type_idx].add(node)
                    break
    return new_nodes_list, susceptible_nodes

def get_accumulated_mask_probs(mask_probs,):
    '''
    Input: 
    mask_probs: A list of mask probabilities, where sum(mask_probs) = 1. 
    Output: 
    accumulated_mask_probs: An length 1 interval with each segment having the length of each ele in mask prob
    e.g. 
    mask_probs = [0.3, 0.2, 0.1, 0.4]
    accumulated_mask_probs = [0, 0.3, 0.5, 0.6, 1]
    '''
    accumulated_mask_probs = [0]
    
    for idx, m in enumerate(mask_probs):
        accumulated_mask_probs.append(m + accumulated_mask_probs[idx])
        
    return accumulated_mask_probs

def get_node_status(accumulated_mask_probs):
    '''
    Input:  accumulated_mask_probs
    Output: int, corresponding to the idx of the mask_prob, 
            the maks wearing type of a single node.
    '''
    roll_dice = random.random() # [0.0, 1.0)
    distance = roll_dice - np.array(accumulated_mask_probs)
    
    mask_type = 0
    for idx, dis in enumerate(distance):
        if dis >= 0 and distance[idx + 1] < 0:
            mask_type = idx
            break
            
    return mask_type

def get_mask_status(mask_probs, num_nodes):
    '''
    Input:  
    mask_prob: A list of mask probabilities(len(mask_prob = num_mask_types)
    num_nodes: Graph size
    Output: 
    mask_status: A list of mask wearing states for each node in the graph.
    '''
    accumulated_mask_probs = get_accumulated_mask_probs(mask_probs,)
    mask_status = []
    for i in range(num_nodes):
        mask_type = get_node_status(accumulated_mask_probs)
        mask_status.append(mask_type)
        
    return mask_status

def get_seed(start_strain, num_nodes, mask_status, num_strain):
    if start_strain < 0 or start_strain >= num_strain:
        print("sim-ray.py: start_strain out of index!")
        assert False
        
    seed = int(np.random.randint(0, num_nodes - 1))   
    strain_list = [set() for i in range(num_strain)]
    while mask_status[seed] != start_strain:
        seed = int(np.random.randint(0, num_nodes - 1))
    strain_list[start_strain] = set([seed])
    
    return seed, strain_list
    

def evolution(g, T_list, mask_prob, start_strain):
    g.simplify()
    node_set = set(g.vs.indices)
    num_nodes = len(node_set)
    num_strain = len(T_list) # num_mask_types
    mask_status = get_mask_status(mask_prob, num_nodes, )
    seed, strain_list = get_seed(start_strain, num_nodes, mask_status, num_strain) 
    susceptible_nodes = node_set
    
    for strain_set in strain_list:
        susceptible_nodes = susceptible_nodes.difference(strain_set)
    new_nodes_list = strain_list # level L - 1

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
    num_infected_list = map(len, strain_list)
    return num_infected, num_infected_list

@ray.remote
def runExp(i, mean_degree, num_nodes, T_list, mask_prob, start_strain): #### should add start_strain
    network = create_network(mean_degree, num_nodes)
    size, size_list = evolution(network, T_list, mask_prob, start_strain)
    return div(size,num_nodes), [div(sizei, num_nodes) for sizei in size_list]

def main():
    ########### Get commandline input ###########
    paras = parse_args(sys.argv[1:])
    paras_check(paras)
    mean_degree_list = get_mean_degree_list(paras)
    k_max, T_list = resolve_paras(paras)
    num_mask_types = len(T_list)
    
    ############ Start Exp ############
    now = datetime.now() # current date and time
    start_time = time.time()
    time_exp = now.strftime("%m%d%H:%M")
    print("-------Exp start at:" + time_exp + '-------')

    for start_strain in range(num_mask_types):
        for mean_degree in mean_degree_list:
            for cp in range(1, int(paras.e/paras.cp) + 1): # cp order
                results_ids = []
                for i in range(paras.cp):
                    results_ids.append(runExp.remote(i, mean_degree, paras.n, T_list, paras.m, start_strain,))  
                results = ray.get(results_ids)
                write_cp_raw_results(results, start_strain, mean_degree, cp, time_exp, start_time, paras,)
    
    write_exp_settings(time_exp, paras, mean_degree_list)

    now_finish = datetime.now() # current date and time
    print("All Done! for:" + time_exp)
    print("--- %.2s seconds in total ---" % (time.time() - start_time))
main()