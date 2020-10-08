import json
import matplotlib.pyplot as plt
from collections import defaultdict 
import numpy as np
import os
from datetime import datetime
import time

def div(x, y):
    if y == 0:
        return 0
    else:
        return x*1.0/y
    
def write_results(results, start_strain, mean_degree, cp, timeExp, check_point, thrVal, change, mask_prob, T, T_mask1, T_mask2, num_nodes, numExp, degree_max, num_samples, mean_degree_list, T_list, start_time):
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
        
    ExpPath = '/mnt/hdd-storage/ytian/ns/simulation/' + 'Mask2Results_'+ change_folder +'/' + 'm' + str(mask_prob) + '_T' + "{0:.2f}".format(T) + '_tm1_' + "{0:.2f}".format(T_mask1) + '_tm2_' + "{0:.2f}".format(T_mask2) + '/'  + 'n' + str(num_nodes) + '_totalexp' + str(numExp) + '/' + timeExp +'/ss'+ str(start_strain) + '/meandegree'+ str(mean_degree) +'/e'+str(check_point) +'_cp' + str(cp)

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