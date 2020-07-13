#!/usr/bin/env python
# coding: utf-8

# In[4]:


import os
import time
import json
import random
import threading
import numpy as np
from os import path
import networkx as nx
from datetime import datetime
import matplotlib.pyplot as plt
from networkx.utils import powerlaw_sequence
from networkx.algorithms.dag import dag_longest_path


# # Function 1.Generate graph by Configuration Model

# Generate a degree sequence with **nodeN** nodes, and **meanDegree**
# 
# 

# In[10]:


def generate_graph(nodeN, mean_degree):
    # Generate degree sequence
    sequence = np.random.poisson(mean_degree, nodeN)
    while (np.sum(sequence) % 2 !=0):
        sequence = np.random.poisson(mean_degree, nodeN)
    # sequence = nx.random_powerlaw_tree_sequence(nodeN, tries=5000)

    # Generate Graph according to the degree sequence
    G = nx.configuration_model(sequence)

    # Remove parallel edges
    G = nx.Graph(G)

    # Remove self-loops
    G.remove_edges_from(nx.selfloop_edges(G))
    return G    


# # Function 2. Decide nodes' mask wearing states
# 
# ## P(A person wears a mask) = m

# Generate the mask wearing states of each node.

# In[11]:


def init_mask(G, mask_prob):
    # A list of 1 and 0 indicating mask wearing or not
    # 1 means wear mask, 0 means not wearing a mask
    masks = np.random.binomial(1, mask_prob, nodeN)

    # Node idx
    nodes = np.linspace(0, nodeN - 1, nodeN, dtype = int)
    # Dict of node attributes
    mask_dict = dict(zip(nodes, masks))

    # Set nodes attributes
    nx.set_node_attributes(G, mask_dict, 'mask')    
    return G, mask_dict


# # Function 3. Init nodes infection states to all 0 (not infected)

# In[12]:


def init_infected(G):
    # Init all nodes to be healthy
    infected = np.zeros(nodeN, dtype = int)

    # Node idx
    nodes = np.linspace(0, nodeN - 1, nodeN, dtype = int)

    # Dict of node attributes
    infected_dict = dict(zip(nodes, infected))

    # Set nodes attributes
    nx.set_node_attributes(G, infected_dict, 'infected')
    
    return G, infected_dict   


# # Function 4. Generate the BFS tree structures for each components

# In[13]:


def generate_BFS_tree(G, mask_dict, infected_dict):
    # Get all the connected components
    components = list(nx.connected_components(G))

    # Roots stores the randomly selected root for each components
    roots = []
    for component in components:
        if len(component) < 100: # ingnore the small components
            continue
        roots.append(random.choice(list(component)))
        
    # Convert the components to a BFS tree-like structure, with randomly selected roots
    # Trees stores all the tree-structured components
    Trees = []
    for root in roots:
        T = nx.bfs_tree(G, source = root)
        nx.set_node_attributes(T, mask_dict, 'mask')
        nx.set_node_attributes(T, infected_dict, 'infected')
        Trees.append(T)
    
    return roots, Trees


# # Function 5. Starting the infection process with 1 virus strain 
# 
# 
# ## Function 5.1 Mask wearing 
# **Notice**
# Infection starts from root to the leaves
# 
# | Infectious    |  Susceptible   | Transmissibillity     | Notation |
# | :------------- | :----------: | :----------- | :----------- |
# | 1             | 0              | T * T_mask^2 * m      | T_1 |
# | 1             | 1              | T * T_mask * (1 - m)  | T_2 |
# | 0             | 0              | T * (1 - m)           | T_3 |
# | 0             | 1              | T * T_mask * m        | T_4 |

# In[14]:


def generate_new_transmissibilities_mask(T_mask, T, m):
    roundN = 5 # Round T to roundN digits
    T1 = round(T * T_mask * T_mask * m, roundN)
    T2 = round(T * T_mask * (1 - m), roundN)
    T3 = round(T * (1 - m), roundN)
    T4 = round(T * T_mask * m , roundN)

    trans_dict = {'T1': T1,
                  'T2': T2,
                  'T3': T3,
                  'T4': T4}

#     print("T1: %.5f" %T1)
#     print("T2: %.5f" %T2)
#     print("T3: %.5f" %T3)
#     print("T4: %.5f" %T4)
    
    return trans_dict    


# ## Function 5.2 Mutation

# In[15]:


def generate_new_transmissibilities_mutation(T_mask, T, m):
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

#     print("Q1: %.5f" %Q1)
#     print("Q2: %.5f" %Q2)

#     print("mu11: %.5f" %mu11)
#     print("mu12: %.5f" %mu12)
#     print("mu22: %.5f" %mu22)
#     print("mu21: %.5f" %mu21)
    return trans_dict, mu_dict


# # Function 6: Traverse from the root nodes using bfs search

# ## Function 6.1 infection with mask

# In[16]:


def start_infection_mask(Trees, roots, trans_dict, infected_dict):
    for idx, tree in enumerate(Trees): 
#         print("TREE No.", idx)

        edge_attr = dict()
        root = roots[idx]
        dfs_edges = list(nx.bfs_edges(tree, source = root))
        total_depth = nx.dag_longest_path_length(tree)

        nx.set_node_attributes(tree, infected_dict, 'infected')
        nx.set_node_attributes(tree, {root: 1}, 'infected') # Let root node be infected by nature


        for depth in range(1, total_depth + 1): # Transmitted level by level
#             print('LEVEL %d' % depth)

            if depth == 1: # Get only this level's node pairs
                dfs_edges = list(nx.dfs_edges(tree, source = root, depth_limit = depth))

            else: 

                dfs_edges = set(nx.dfs_edges(tree, source = root, depth_limit = depth)) -                             set(nx.dfs_edges(tree, source = root, depth_limit = depth - 1))


            for father, son in dfs_edges: # Check each node pairs in this level
#                 print("(%d, %d), node %d is_infected = %d" %(father, son, father, tree.nodes[father]['infected'] ))

                if tree.nodes[father]['infected'] == 1: 

                    # Decide which transmissibility
                    if tree.nodes[father]['mask'] == 1 and tree.nodes[son]['mask'] == 0:
                        T_edge = 'T1'
                    elif tree.nodes[father]['mask'] == 1 and tree.nodes[son]['mask'] == 1:
                        T_edge = 'T2'
                    elif tree.nodes[father]['mask'] == 0 and tree.nodes[son]['mask'] == 0:
                        T_edge = 'T3'
                    else:
                        T_edge = 'T4'

                    edge_attr[(father,son)] = {'T': T_edge}


                    # Set the 'Transmissibility'edge attrs
                    nx.set_edge_attributes(tree, edge_attr)


                    # Decide if the susceptible is infected
                    is_infected = int(random.random() < trans_dict[T_edge])

#                     if is_infected:
#                         print("node %d is infected with %s" %(son, T_edge))
#                     else:
#                         print("node %d is not infected" %(son))


                    # Set the 'infected' node attr accordingly
                    nx.set_node_attributes(tree, {son: is_infected}, 'infected')    


# ## Function 6.2 infection with mutation

# In[17]:


def start_infection_mutation(Trees, roots, trans_dict, infected_dict, mu_dict):

    for idx, tree in enumerate(Trees): 
#         print("TREE No.", idx)

        edge_attr = dict()
        root = roots[idx]
        dfs_edges = list(nx.bfs_edges(tree, source = root))
        total_depth = nx.dag_longest_path_length(tree)

        nx.set_node_attributes(tree, infected_dict, 'infected')
        nx.set_node_attributes(tree, {root: 1}, 'infected') # Let root node be infected by nature


        for depth in range(1, total_depth + 1): # Transmitted level by level
#             print('LEVEL %d' % depth)

            if depth == 1: # Get only this level's node pairs
                dfs_edges = list(nx.dfs_edges(tree, source = root, depth_limit = depth))

            else: 

                dfs_edges = set(nx.dfs_edges(tree, source = root, depth_limit = depth)) -                             set(nx.dfs_edges(tree, source = root, depth_limit = depth - 1))


            for father, son in dfs_edges: # Check each node pairs in this level
#                 print("(%d, %d), node %d is_infected = %d" %(father, son, father, tree.nodes[father]['infected'] ))

                if tree.nodes[father]['infected'] == 1: 



                    # Decide which transmissibility
                    if tree.nodes[father]['mask'] == 1:
                        T_edge = 'Q1'
                    else:
                        T_edge = 'Q2'
                    # Set the 'Transmissibility'edge attrs
                    edge_attr[(father,son)] = {'T': T_edge}
                    nx.set_edge_attributes(tree, edge_attr)



                    # Decide if the susceptible is infected
                    is_infected = int(random.random() < trans_dict[T_edge])
                    # Set the 'infected' node attr accordingly
                    nx.set_node_attributes(tree, {son: is_infected}, 'infected')

#                     if is_infected:
#                         print("node %d is infected with %s" %(son, T_edge))



                    # Decide if the susceptible mutates
                    if tree.nodes[son]['mask'] == 1: # mu11 or mu12

                        # Decide if the susceptible is mutated
                        is_mutated = int(random.random() > mu_dict['mu11'])
                        # Set the 'mask' node attr accordingly
                        if is_mutated: # mask: 1 -> 0
                            nx.set_node_attributes(tree, {son: 0}, 'mask')
#                             print("node %d is mutated to not wearking mask!" %(son))

                    else:

                        # Decide if the susceptible is mutated
                        is_mutated = int(random.random() > mu_dict['mu22'])
                        # Set the 'mask' node attr accordingly
                        if is_mutated: # mask: 0 -> 1
                            nx.set_node_attributes(tree, {son: 1}, 'mask')
#                             print("node %d is mutated to wearking mask!" %(son))
                        



def cal_EpdSize_Trees(Trees):
    E_S = []
    for idx, tree  in enumerate(Trees):
        res = np.array(list(tree.nodes.data('infected')))
#         print(res)
        es = sum(res[:,1]) / res.shape[0]
        E_S.append(es)
#         print("Tree %d Epidemic Size: %.3f" %(idx, es))
    avg_ES = div(sum(E_S), len(E_S))
#     print('Avg Epidemic Size:', avg_ES)
    return avg_ES    


# In[20]:


# threshold = 0.05 # threshold to check if a process becomes a epidemic
# mean_degree = 2
# nodeN = 20
# mask_prob = 0.6
# T_mask = 0.1
# T = 0.8


# # Experiments Auxiliary Functions

# In[23]:


def div(x, y):
    if y == 0:
        return 0
    else:
        return x*1.0/y


# In[24]:


def draw_figures(mean_degree_list, Prob_Emergence, AvgValidSize, AvgSize, ExpPath):
    figure_path = ExpPath + '/' + 'Figures'
    if not os.path.exists(figure_path):
        print("make path ", figure_path)
        os.mkdir(figure_path)
    
    ### Probability of Emergence ###
    plt.figure()
    plt.plot(mean_degree_list, Prob_Emergence['mask'], 'bo')
    plt.xlabel("Mean Degree")
    plt.ylabel("Prob of Emergence")
    title = "Probability of Emergence for Mask Model"
    plt.title(title)
    plt.savefig(figure_path + '/' + title.replace(" ", "_") + '.png')
    
    plt.figure()
    plt.plot(mean_degree_list, Prob_Emergence['mutation'], 'ro')
    plt.xlabel("Mean Degree")
    plt.ylabel("Prob of Emergence")
    title = "Probability of Emergence for Mutation Model"
    plt.title(title)
    plt.savefig(figure_path + '/' + title.replace(" ", "_") + '.png')
    
    ### Epidemic Size ###
    plt.figure()
    plt.plot(mean_degree_list, AvgValidSize['mask'], 'bo')
    plt.xlabel("Mean Degree")
    plt.ylabel("Epidemic Size")
    title = "Epidemic Size for Mask Model"
    plt.title(title)
    plt.savefig(figure_path + '/' + title.replace(" ", "_") + '.png')
    
    plt.figure()
    plt.plot(mean_degree_list, AvgValidSize['mutation'], 'ro')
    plt.xlabel("Mean Degree")
    plt.ylabel("Epidemic Size")
    title = "Epidemic Size for Mutation Model"
    plt.title(title)
    plt.savefig(figure_path + '/' + title.replace(" ", "_") + '.png')
    
    ### Infected Frac ###
    plt.figure()
    plt.plot(mean_degree_list, AvgSize['mask'], 'bo')
    plt.xlabel("Mean Degree")
    plt.ylabel("Infected Frac")
    title = "Infected Frac for Mask Model"
    plt.title(title)
    plt.savefig(figure_path + '/' + title.replace(" ", "_") + '.png')
    
    plt.figure()
    plt.plot(mean_degree_list, AvgSize['mutation'], 'ro')
    plt.xlabel("Mean Degree")
    plt.ylabel("Infected Frac")
    title = "Infected Frac for Mutation Model"
    plt.title(title)
    plt.savefig(figure_path + '/' + title.replace(" ", "_") + '.png')
    


# In[49]:


def runExp(i, mean_degree, nodeN, trans_dict_mask, trans_dict_mu, mu_dict):
    '''
    Arg: 
    mean_degree: Random graph's mean degree
    nodeN:       Node number
    mask_prob:   P(a person wearing a mask)
    T_mask:      The impact of wearing a mask
    T:           The transmissibility of a real strain
    
    '''
    
    G = generate_graph(nodeN, mean_degree)
    G, mask_dict = init_mask(G, mask_prob)
    G, infected_dict = init_infected(G)
    roots, Trees = generate_BFS_tree(G, mask_dict, infected_dict)
    Trees2 = Trees.copy()
    
    # Transmission with mask
    start_infection_mask(Trees, roots, trans_dict_mask, infected_dict)
    mask_size = cal_EpdSize_Trees(Trees)       
    
    # Transmission with mutation
    start_infection_mutation(Trees2, roots, trans_dict_mu, infected_dict, mu_dict)
    mutation_size = cal_EpdSize_Trees(Trees2) 
    
    # Count Results
    if mask_size > threshold:
        Epidemics['mask'][i] = mask_size
        numEpidemics['mask'] += 1
        ttlEpidemicsSize['mask'] += mask_size

    if mutation_size > threshold:
        Epidemics['mutation'][i] = mask_size
        numEpidemics['mutation'] += 1
        ttlEpidemicsSize['mutation'] += mutation_size


    ttlFrac['mask'] += mask_size
    ttlFrac['mutation'] += mutation_size
    
    
    


# In[53]:


#!/usr/bin/python3

import threading
import time


class myThread (threading.Thread):
    def __init__(self, i, mean_degree, nodeN, trans_dict_mask, trans_dict_mu, mu_dict):
        threading.Thread.__init__(self)
        
        self.i = i
        self.mean_degree = mean_degree
        self.nodeN = nodeN
        self.trans_dict_mask = trans_dict_mask
        self.trans_dict_mu = trans_dict_mu
        self.mu_dict = mu_dict

    def run(self):
#         print("Starting Thread Exp-%d, meanDeg-%.2f" % (self.i, self.mean_degree))
        runExp(self.i, self.mean_degree, self.nodeN, self.trans_dict_mask, self.trans_dict_mu, self.mu_dict)
#         print("Exiting Thread Exp-%d, meanDeg-%.2f" % (self.i, self.mean_degree))

    
    
# # Start new Threads
# thread1.start()
# thread2.start()
# thread1.join()
# thread2.join()
# print ("Exiting Main Thread")


# # Exp 1: Check mean degree imapct

# ## Parameters setting

# In[54]:
def parse_args(args):
    parser = argparse.ArgumentParser(description = 'Parameters')
#     parser.add_argument('-mean', type = float, nargs = '+', default = np.arange(1, 10.1, 0.1), help='np.linspace(0.001, 7, 50) (default); list of mean degree: you can type 1 3 5')
    parser.add_argument('-n', type = int, default = 200000, help='10,000 (default); the number of nodes')
    parser.add_argument('-e', type = int, default = 50, help='50 (default); the number of experiments')
    parser.add_argument('-mask', type = float, default = 0.6, help='0.6 (default); the prob of wearing a mask')
    parser.add_argument('-tm', type = float, default = 0.1, help='0.1 (default); T_mask')
    parser.add_argument('-T', type = float, default = 0.1, help='0.1 (default); T_mask')
    parser.add_argument('-t', type = float, default = 0.005, help='0.001 (default); the treshold to consider a component giant')
    parser.add_argument('-numCores', type = int, default = 12, help='number of Cores')
    return parser.parse_args(args)

paras = parse_args(sys.argv[1:])
mean_degree_list = np.linspace(0, 10, 50)

nodeN = paras.n
ExpN = paras.e
threshold = paras.t
num_cores = min(paras.numCores,multiprocessing.cpu_count())
mask_prob = paras.mask


T_mask = paras.tm
T = paras.T


# ## Start Experiment

# In[ ]:


########### Experiment Time & Path Setup ###########
now = datetime.now() # current date and time
timeExp = now.strftime("%m%d%H:%M")
ExpPath = 'Results/' + timeExp + '_n' + str(nodeN) + '_e' + str(ExpN)

if not os.path.exists(ExpPath):
    print("make path ", ExpPath)
    os.mkdir(ExpPath)


######### Experiment Variables Declararion #########
Prob_Emergence = dict()
Prob_Emergence['mask'] = []
Prob_Emergence['mutation'] = []

AvgValidSize = dict()
AvgValidSize['mask'] = []
AvgValidSize['mutation'] = []

AvgSize = dict()
AvgSize['mask'] = []
AvgSize['mutation'] = []

StdValidSize = dict()
StdValidSize['mask'] = []
StdValidSize['mutation'] = []

################# Run Experiments #################

# Caculate the transmissibilities for Mask Model, i.e. T1 to T4
trans_dict_mask = generate_new_transmissibilities_mask(T_mask, T, mask_prob) 
# Caulculate the transmissibilities and mutation probs for Mutation Model
trans_dict_mu, mu_dict = generate_new_transmissibilities_mutation(T_mask, T, mask_prob) 


for mean_degree in mean_degree_list:
    print("MeanDegree:", mean_degree)
    a = time.time()
    
    # The epidemic size list for each mean degree value
    Epidemics = dict()
    Epidemics['mask'] = dict()
    Epidemics['mutation'] = dict()
    # The count of successful epidemics
    numEpidemics = dict()
    numEpidemics['mask'] = 0
    numEpidemics['mutation'] = 0
    # Total Epidemics Size
    ttlEpidemicsSize = dict()
    ttlEpidemicsSize['mask'] = 0
    ttlEpidemicsSize['mutation'] = 0
    # Total frac of infected people, not only epidemics
    ttlFrac = dict()
    ttlFrac['mask'] = 0
    ttlFrac['mutation'] = 0
    
    threads = []
    
    for i in range(ExpN): # Each degree run ExpN times
        threads.append(myThread(i, mean_degree, nodeN, trans_dict_mask, trans_dict_mu, mu_dict))
        threads[-1].start()
    
    for i in range(ExpN):
        threads[i].join()
    
    print("All Experiments for mean_degree %.2f ends!" %mean_degree)

        

        
        
    ######### Record the results for this Mean Degree ##########    
    timedeltta = time.time()-a
       
    Prob_Emergence['mask'].append(numEpidemics['mask']*1.0/(ExpN))
    AvgValidSize['mask'].append(div(ttlEpidemicsSize['mask']*1.0, numEpidemics['mask']))
    AvgSize['mask'].append(ttlFrac['mask']*1.0/ExpN)
    StdValidSize['mask'].append(np.std(list(Epidemics['mask'].values())))
    
    Prob_Emergence['mutation'].append(numEpidemics['mutation']*1.0/(ExpN))
    AvgValidSize['mutation'].append(div(ttlEpidemicsSize['mutation']*1.0, numEpidemics['mutation']))
    AvgSize['mutation'].append(ttlFrac['mutation']*1.0/ExpN)
    StdValidSize['mutation'].append(np.std(list(Epidemics['mutation'].values())))
    

######### Save the results for all Mean Degrees ########## 

draw_figures(mean_degree_list, Prob_Emergence, AvgValidSize, AvgSize, ExpPath)

setting_path = ExpPath + '/' + 'Settings'
if not os.path.exists(setting_path):
    os.mkdir(setting_path)
    
res_paths = ExpPath + '/' + 'Results'
if not os.path.exists(res_paths):
    os.mkdir(res_paths)

### Experiment Parameters ###
paras = dict()
paras['ExpN'] = ExpN
paras['nodeN'] = nodeN
paras['mask_prob'] = mask_prob
paras['T_mask'] = T_mask
paras['T'] = T
paras['threshold'] = threshold
with open(setting_path + '/paras.json', 'w') as fp:
    json.dump(paras, fp)
    
### Degree list ###
np.save(setting_path + '/mean_degree_list.npy', np.array(mean_degree_list)) 

### Transmissibilites for Mask Model ###
with open(setting_path + '/trans_dict_mask.json', 'w') as fp:
    json.dump(trans_dict_mask, fp)
    
### Transmissibilites and mutation probs for Mutation Model ###
with open(setting_path + '/trans_dict_mu.json', 'w') as fp:
    json.dump(trans_dict_mu, fp)
    
with open(setting_path + '/mu_dict.json', 'w') as fp:
    json.dump(mu_dict, fp)

### Results ###
with open(res_paths + '/Prob_Emergence.json', 'w') as fp:
    json.dump(Prob_Emergence, fp)

with open(res_paths + '/AvgValidSize.json', 'w') as fp:
    json.dump(AvgValidSize, fp)

with open(res_paths + '/StdValidSize.json', 'w') as fp:
    json.dump(StdValidSize, fp)
    
with open(res_paths + '/AvgSize.json', 'w') as fp:
    json.dump(AvgSize, fp)
    


# In[ ]:




