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
import os
from subprocess import Popen, PIPE





##########################################################
# 1 Parse commandline 
##########################################################

def parse_commandline():
    path_to_wdir = argv[1]  
    randomization_option = argv[2]
    return path_to_wdir, randomization_option


def read_in_filtered_sample_shuffle_freqs(path_to_wdir, sample_name):
    with open(path_to_wdir + sample_name + "/" + sample_name + "_filtered.tsv", "r") as f: 
        lines = [line.strip("\n").split("\t") for line in f if not line.startswith("#")]
    clones = list()
    frequencies = list()
    for (line, enumerator) in zip(lines, range(len(lines))): 
        clones.append(line[0])
        frequencies.append(float(line[1]))
    np.random.shuffle(frequencies)
    cdr3_freq_subsample_dict = {cdr3:freq for (cdr3,freq) in zip(clones, frequencies)}
    return cdr3_freq_subsample_dict, clones




def normalize(d, target=1.0):
    """
    Returns a transformed dictionary with the values normalized
    """
    raw = sum(d.values())
    factor = target/raw
    return {key:value*factor for (key,value) in d.items()}



def read_in_filtered_sample_random_freq(path_to_wdir, sample_name):
    with open(path_to_wdir + sample_name + "/"+ sample_name + "_filtered.tsv", "r") as f: 
        lines = [line.strip("\n").split("\t") for line in f if not line.startswith("#")]
    clones = list()
    frequencies = list()
    clones = [line[0] for line in lines]
    frequencies = np.random.random(size=len(clones))
    cdr3_freq_subsample_dict = {cdr3:freq for (cdr3,freq) in zip(clones, frequencies)}
    cdr3_freq_subsample_dict = normalize(cdr3_freq_subsample_dict, target=1.0)
    return cdr3_freq_subsample_dict, clones








def write_out_filtered_sample(directory_path, sample_name, cdr3_freq):
    """
    Writes out dictionary from the subsampled filtered dataframe. 
    Each row: cdr3 /t freq
    """
    with open(directory_path + "/" + sample_name + "_filtered.tsv", 'w') as f:
        for cdr3,freq in sorted(cdr3_freq.items()): # here added sorted
            f.writelines('\t'.join([cdr3,str(freq)]) + '\n')
    del cdr3_freq




##########################################################
# 3 Executable part
##########################################################


if __name__ == '__main__':

    # READ IN PATH TO FILE.
    path_to_wdir, randomization_option = parse_commandline()    
    # READ IN NUMBER OF FILES
    sample_folders = [f for f in os.listdir(path_to_wdir) if os.path.isdir(os.path.join(path_to_wdir, f))]
    if randomization_option == "shuffle": 
        for sample in sample_folders: 
            # read in filtered sample and shuffle it: 
            cdr3_freq_dict_shuffled, clones = read_in_filtered_sample_shuffle_freqs(path_to_wdir, sample)
            # write out shuffled filtered sample 
            write_out_filtered_sample(path_to_wdir + sample, sample, cdr3_freq_dict_shuffled)
    if randomization_option == "randomize": 
        for sample in sample_folders: 
            # read in filtered sample and assign random freqs:
            cdr3_freq_dict_random, clones = read_in_filtered_sample_random_freq(path_to_wdir, sample)
            # write out random freq sample: 
            write_out_filtered_sample(path_to_wdir + sample, sample, cdr3_freq_dict_random)
