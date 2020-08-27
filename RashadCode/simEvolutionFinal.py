import argparse
import collections
import random
import sys, pdb
import site, os
site.addsitedir('/afs/ece.cmu.edu/usr/reletreb/dep/lib/python2.7/site-packages')

import numpy as np
import igraph as ig
from joblib import Parallel, delayed
import multiprocessing, time
from multiprocessing import Manager



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

def parse_args(args):
    parser = argparse.ArgumentParser(description = 'Parameters')
    parser.add_argument('-m', type = float, nargs = '+', default = np.arange(1, 10.1, 0.1), help='np.linspace(0.001, 7, 50) (default); list of mean degree: you can type 1 3 5')
    parser.add_argument('-n', type = int, default = 200000, help='10,000 (default); the number of nodes')
    parser.add_argument('-e', type = int, default = 50, help='100 (default); the number of experiments')
    parser.add_argument('-t1', type = float, default = 0.2, help='0.5 (default); the transmissibility of strain-1')
    parser.add_argument('-t2', type = float, default = 0.5, help='0.5 (default); the transmissibility of strain-2')
    parser.add_argument('-m1', type = float, default = 0.9, help='0.5 (default); the mutation probability from 1 to 1')
    parser.add_argument('-m2', type = float, default = 1.0, help='0.5 (default); the mutation probability from 2 to 2')
    parser.add_argument('-thrVal', type = float, default = 0.005, help='0.001 (default); the treshold to consider a component giant')
    parser.add_argument('-numCores', type = int, default = 12, help='number of Cores')
    parser.add_argument('-logName', default = 'logfile', help='The name of the log file')
    parser.add_argument('-i', type = int, default = 1, help='1 (default); starting from type-i node')
    return parser.parse_args(args)

paras = parse_args(sys.argv[1:])
mean_degree_list = paras.m
t1 = paras.t1
t2 = paras.t2
m1 = paras.m1
m2 = paras.m2
num_nodes = paras.n
numExp = paras.e
start_strain = paras.i
num_cores = min(paras.numCores,multiprocessing.cpu_count())
thrVal = paras.thrVal

T_list = [t1, t2]
mutation_probability = [m1, m2]
ff = open(paras.logName+'Det','w+')
f = open(paras.logName, 'w+')

for mean_degree in mean_degree_list:
    a = time.time()
    ttlEpidemicsSize = 0
    numEpidemics = 0
    Epidemics = []
    EpidemicsPerSt = [0,0,0]
    fractionDic = Manager().dict()
    infectedPerStDic = Manager().dict()
    ttlFrac = 0

    Parallel(n_jobs = num_cores)(delayed(runExp)(i, mean_degree,num_nodes, T_list, mutation_probability) for i in xrange(numExp))

    for ii in xrange(numExp):
        resultsFrac = 'meanDeg: {0} Size: {1} infSt1: {2} infSt2: {3}\n'.format(mean_degree, fractionDic[ii], infectedPerStDic[ii][0],infectedPerStDic[ii][1] )

        if fractionDic[ii] >= thrVal:
            ff.write(resultsFrac)
            ff.flush()
            numEpidemics += 1
            ttlEpidemicsSize += fractionDic[ii]
            Epidemics.append(fractionDic[ii])
            EpidemicsPerSt[0] += infectedPerStDic[ii][0]
            EpidemicsPerSt[1] += infectedPerStDic[ii][1]

        ttlFrac += fractionDic[ii]

    if len(Epidemics) == 0:
        Epidemics.append(0)


    results = 'numExp: {0} Threshold: {1} n: {2} meanDeg: {3} Prob: {4} AvgValidSize: {5} StdValidSize: {6} infSt1: {7} infSt2: {8} AvgSize: {9} T: {10} Mu: {11} Time: {12} \n'.format(numExp, thrVal , num_nodes, mean_degree, numEpidemics*1.0/(numExp), div(ttlEpidemicsSize*1.0,numEpidemics), np.std(Epidemics) , div(EpidemicsPerSt[0],numEpidemics),div(EpidemicsPerSt[1],numEpidemics), ttlFrac*1.0/numExp,' '.join(map(str, T_list)), ' '.join(map(str, mutation_probability)), time.time()-a )

    print results
    f.write(results)
    f.flush()
