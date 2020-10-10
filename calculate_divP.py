"""
READ ME 

INPUT:  path_to_wdir
        sample_name
        identity
        qs
        list_of_lambdas
        cc_size

OUTPUT: 
divP identity for all the qs and lambda values 
divP sim q= 0,1,inf for all different lambda values

USAGE: python3 path_to_wdir sample_name identity n_number qs list_of_lambdas



Takes in as input already calculated distance matrices. 


It can calculate 2 types of diverisites: with and without taking similarity
into account. 

Both are calculated explicitly with the new script. 

"""
from sys import argv
from collections import defaultdict
import pandas as pd 
import numpy as np
import time
import math
import itertools
import multiprocessing as mp
import gc


##########################################################
# 1 Parse commandline 
##########################################################

def parse_commandline():
    path_to_wdir = argv[1]  
    sample_name = argv[2] 
    identity = argv[3]
    #n = int(float(argv[4]))     # Number of random samples for 2<=q<infinity
    qs = argv[4]              # List of qs
    if 'inf' in qs:
        qs = qs.strip('[]').split(',')
        qs = [q for q in qs if q!= 'inf']
        qs = list(map(int, qs)) + [11]
        #qs = qs + ['inf']
    else:
        qs = list(map(int, qs.strip('[]').split(',')))
    list_of_lambdas = argv[5] # List of Lambdas
    list_of_lambdas = list(map(float, list_of_lambdas.strip('[]').split(',')))
    cc_size = int(argv[6]) # chunk size 
    return path_to_wdir, sample_name, identity, qs, list_of_lambdas, cc_size


def read_in_filtered_sample(path_to_wdir, sample_name):
    with open(path_to_wdir + sample_name + "/"+ sample_name + "_filtered.tsv", "r") as f: 
        lines = [line.strip("\n").split("\t") for line in f if not line.startswith("#")]
    cdr3_freq_subsample_dict = defaultdict()
    cdr3_order_subsample_dict = defaultdict()
    clones = list()
    for (line, enumerator) in zip(lines, range(len(lines))): 
        cdr3_freq_subsample_dict[line[0]]=float(line[1])
        cdr3_order_subsample_dict[line[0]]= enumerator
        clones.append(line[0])
    return cdr3_freq_subsample_dict, cdr3_order_subsample_dict, clones


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def read_dist_from_chunk(chunk_number, cdr31_file_indx, cdr32_indx, path_to_wdir, sample_name):
    with open(path_to_wdir + sample_name + "/" + chunk_number + "_" + sample_name + "_filtered_dmat.tsv", "r") as fl: 
        distance_list_cdr3_i = [line.strip("\n").split("\t") for line in fl if not line.startswith("#")]
    #
    if fl.closed == False:
        print(fl.closed)
    dist_cdr3 = distance_list_cdr3_i[cdr31_file_indx][cdr32_indx]
    return float(dist_cdr3)

def read_dist_ij_parallel(cdr3_tuple_list, path_to_wdir, sample_name):
    dist_ij = [read_dist_from_chunk(chunk_number, cdr31_file_indx, cdr32_indx, path_to_wdir, sample_name) for (chunk_number, cdr31_file_indx, cdr32_indx) in cdr3_tuple_list]
    return dist_ij

def calculate_Zpi_parallel(path_to_wdir,sample_name, chunk_number, cdr31_indx, lambdas_list, frequencies):
    # read in distance line
    with open(path_to_wdir + sample_name + "/" + chunk_number + "_" + sample_name + "_filtered_dmat.tsv", "r") as fl: 
        lines = fl.readlines()
        line = lines[cdr31_indx]
        del lines 
    distance_list_cdr3_i = line.strip("\n").split("\t")
    distance_list_cdr3_i = np.array(distance_list_cdr3_i)
    distance_list_cdr3_i = distance_list_cdr3_i.astype(np.float)
    Zpi_dict_i = defaultdict()
    for lambda_k in lambdas_list:
        print(lambda_k)
        Zp_i = distance_list_cdr3_i*lambda_k
        Zp_i = np.exp(-Zp_i)
        Zp_i = np.multiply(Zp_i, frequencies)
        Zpi_dict_i[lambda_k] = (np.sum(Zp_i))
        print(Zpi_dict_i[lambda_k])
    return Zpi_dict_i





##########################################################
# 2 Diversity Functions 
##########################################################

#---------------------Identity---------------------------#

def diversity_q_identity(freqs,q):
    """
    Returns the q-diversity of a list of frequencies
    """
    diversity = 0.0
    #
    #if (q != 1.0) and (q != 'inf'):
    if (q != 1.0) and (q != 11):
        for freq in freqs:
            if freq != 0.0:
                diversity += freq**q
        diversity = diversity**(1.0/(1.0-q))
    #
    elif q == 1.0:
        for freq in freqs:
            if freq != 0.0:
                diversity -= freq*math.log(freq)       
        diversity = math.exp(diversity)
    #elif q == 'inf':
    elif q == 11:
        diversity = 1.0 /max(freqs)
    #
    return diversity

def make_diversity_w_identity_dict(cdr3_freq_subsample_dict, list_of_qs):
    """
    Dictionary with keys [q, 'identity'], diversity calculated with the identity matrix as values.  
    """
    diversity_w_identity = dict()
    for qu in list_of_qs: 
        diversity_w_identity[qu, 'identity'] = diversity_q_identity(list(cdr3_freq_subsample_dict.values()), qu)
    #
    return diversity_w_identity

#---------------------Similarity--------------------------#

def diversity_q_similarity(list_Z, q, number_of_random_samples):
    #if (q != 1) and (q != 'inf') and (q!= 0):
    if (q != 1) and (q != 11) and (q!= 0):    
        Z = np.array(list_Z[0])
        for ii in (range(q-1)[1:]):
            Z = Z*np.array(list_Z[ii])
        assert float(number_of_random_samples) == len(Z)
        diversity = (sum(Z)/len(Z))**(1/(1-float(q)))
    #   
    return diversity


def make_diversity_w_similarity_dict(cdr3_freq_subsample_dict, cdr3_order_subsample_dict, clones, path_to_wdir, sample_name , list_of_qs, lambdas_list, cc_size):
    """
    z = np.array([vals])
    z = z**(1./(1.-q))
    
    gives a speed up
    
    """
    # read in Zpi in a parallel fashion 
    Zp_i_dict = defaultdict(list)# (Zp)_i is stored in Zp[lambda]
    frequencies = [cdr3_freq_subsample_dict[key] for key in clones]
    frequencies = np.array(frequencies)
    cdr3_i_indx = [cdr3_order_subsample_dict[key] for key in clones]
    start_overall_time = time.time()
    cores = 8
    Que = list()
    pool = mp.Pool(processes = cores, maxtasksperchild = 1000)
    for cdr3_i in clones:
        chunk_number = str(cdr3_order_subsample_dict[cdr3_i]//cc_size)
        cdr31_file_indx = cdr3_order_subsample_dict[cdr3_i]%cc_size
        process = pool.apply_async(calculate_Zpi_parallel, args =(path_to_wdir,sample_name, chunk_number,cdr31_file_indx, lambdas_list, frequencies))
        Que.append(process)
    start_time = time.time()
    for oo in range(len(Que)):
        function_output = Que[oo].get()
        for lambda_k in lambdas_list: 
            Zp_i_dict[lambda_k].append(function_output[lambda_k])
        del function_output
    pool.terminate()
    del Que
    gc.collect()
    print("Used overall_time: --- %s seconds ---" % (time.time() - start_overall_time))
    print("Used loop time: --- %s seconds ---" % (time.time() - start_time))
    diversity_w_similarity = dict()
    for lambda_k in lambdas_list: 
        for qu in list_of_qs: 
            if qu not in [11,0,1]:
                print(qu)
                Zpi_power_q = np.power(np.array(Zp_i_dict[lambda_k]), (qu-1))
                diversity_w_similarity[qu, lambda_k] = np.power(np.sum(freq_i*Zpiq for (freq_i, Zpiq) in zip(frequencies, Zpi_power_q)), (1/(1-qu)))
                #diversity_w_similarity[qu, lambda_k] = np.power(np.sum(freq_i*np.power(Zp_i,(qu-1))  for (freq_i,Zp_i) in zip(frequencies, Zp_i_dict[lambda_k])),(1/1-qu))
            elif qu == 0: 
                diversity_w_similarity[0, lambda_k] = np.sum(freq_i/Zp_i for (freq_i,Zp_i) in zip(frequencies, Zp_i_dict[lambda_k]))
            elif qu == 1:
                ln_D_q1 = -np.sum(freq_i*np.log(Zp_i) for (freq_i,Zp_i) in zip(frequencies, Zp_i_dict[lambda_k]))
                diversity_w_similarity[1, lambda_k] = np.exp(ln_D_q1)
            elif qu == 11:
                diversity_w_similarity[11, lambda_k] = 1.0/max(Zp_i_dict[lambda_k])
    return diversity_w_similarity


def write_out_table(sample_name, diversity_w_similarity, identity):
    """
    identity= 'ID'/'SIM'
    """
    with open(sample_name + '_' + identity + ".tsv", 'w') as f:
        # write out argv options:
        command_line = ['# ', 'folder_path: ',  argv[1], 'sample_name: ', argv[2], 'similarity: ', argv[3], 'qs: ', argv[4], 'lambdas: ', argv[5]]
        f.writelines('\t'.join(command_line) + '\n')
        #
        q_list = sorted(list(set([k1 for (k1,k2) in diversity_w_similarity])))
        lambdas_list = sorted(list(set([k2 for (k1,k2) in diversity_w_similarity])))
        #
        header = ["q"] + [str(lambda_k) for lambda_k in lambdas_list]
        f.writelines('\t'.join(header) + '\n')
        row = list()
        #
        for q in q_list:
            row = [q] + [diversity_w_similarity[q, lambda_k] for lambda_k in lambdas_list]
            row = map(str, row)
            f.writelines('\t'.join(row) + '\n')


##########################################################
# 3 Executable part
##########################################################


if __name__ == '__main__':
    start_overall_time = time.time()
    path_to_wdir, sample_name, identity, qs, lambdas_list, cc_size = parse_commandline()    
    cdr3_freq_dict, cdr3_order_dict, clones = read_in_filtered_sample(path_to_wdir, sample_name)
    if identity == 'ID':
        diversity_dict = make_diversity_w_identity_dict(cdr3_freq_dict, qs)
    #
    else:
        diversity_dict = make_diversity_w_similarity_dict(cdr3_freq_dict, cdr3_order_dict, clones, path_to_wdir, sample_name , qs, lambdas_list, cc_size)
        
    #
    write_out_table(sample_name, diversity_dict, identity)
    print("Used time: --- %s seconds ---" % (time.time() - start_overall_time))

