"""
INPUT:  filtered_"sample".tsv file containing In frame cdr3 seq and freq in each row. Freq are normalized. 

OUTPUT: cc_"sample_name"_distance_matrix.tsv file distances of 50 cdr3i against all cdr3j, 

USAGE: python3 calculate_dmat.py filtered_"sample".tsv chunk_size chunk_num distance_option

OPTIONS: 

The distance can be calculated using different metrics therefore the options are: 

"AF_euclidean" "Levenshtein" "BLOSUM45_score_dist"

Atchley factor distance:  calculated using Atchley factors for each amino acid, 
as described in : # Atchley et al. "Solving the protein sequence metric problem"
					2005, PNAS, vol 102, p 6395-6400
It is the euclidean distance between average 5 AFs for each CDR3. 


Levenshtein distance: string metric: distance between two words is the minimum number of 
single-character edits (insertions, deletions or substitutions) required to change one word into the other




"""

from sys import argv 
from collections import defaultdict
import pandas as pd 
import numpy as np
import os 				
import multiprocessing as mp
import time
import gc

from Bio.Align import substitution_matrices
from Bio import Align

##########################################################
# 1.1. Parse commandline 
##########################################################

def parse_commandline():
	"""
	commandline call: python3 calculate_dmat.py input_filename distance_option 
	
		"""
	filtered_input_filename = argv[1] #path to file
	chunk_size = int(float(argv[2]))
	chunk_num = int(float(argv[3]))
	distance_option = argv[4] #"AF_euclidean" "Levenshtein"n "BLOSUM45_score_dist"
	return filtered_input_filename, chunk_size, chunk_num, distance_option


##########################################################
# 1.2. Distance Functions 
##########################################################
#OPTION Levenshtein
def levenshteinDistance(s1, s2):
	if len(s1) > len(s2):
		s1, s2 = s2, s1
	#
	distances = range(len(s1) + 1)
	for i2, c2 in enumerate(s2):
		distances_ = [i2+1]
		for i1, c1 in enumerate(s1):
			if c1 == c2:
				distances_.append(distances[i1])
			else:
				distances_.append(1 + min((distances[i1], distances[i1 + 1], distances_[-1])))
		distances = distances_
	return distances[-1]

##########################################################

#OPTION "AF_euclidean"
def make_AA_to_Atchley_dict():
	"""
	REFERENCE: 
	Atchley et al. "Solving the protein sequence metric problem", 2005, PNAS, vol. 102, no. 18, pp: 6395-6400
	
	Atchley_factor_data order of the rows follows: 
	list_of_AAs=['A','C', 'D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','Z']
	each element of the list follows the format: ['AF1 AF2 AF3 AF4 AF5', ...]

	# AA_to_Atchley dict  has (letter symbol of AminoAcid  as sting, AF index as integer) as keys, and the corresponding AF as float as values
	"""
	Atchley_factor_data = ['-0.591 -1.302 -0.733 1.570 -0.146', 
	'-1.343 0.465 -0.862 -1.020 -0.255',
	'1.050 0.302 -3.656 -0.259 -3.242',
	'1.357 -1.453 1.477 0.113 -0.837',
	'-1.006 -0.590 1.891 -0.397 0.412',
	'-0.384 1.652 1.330 1.045 2.064',
	'0.336 -0.417 -1.673 -1.474 -0.078',
	'-1.239 -0.547 2.131 0.393 0.816',
	'1.831 -0.561 0.533 -0.277 1.648',
	'-1.019 -0.987 -1.505 1.266 -0.912',
	'-0.663 -1.524 2.219 -1.005 1.212',
	'0.945 0.828 1.299 -0.169 0.993',
	'0.189 2.081 -1.628 0.421 -1.392',
	'0.931 -0.179 -3.005 -0.503 -1.853',
	'1.538 -0.055 1.502 0.440 2.897',
	'-0.228 1.339 -4.760 0.670 -2.647',
	'-0.032 0.326 2.213 0.908 1.313',
	'-1.337 -0.279 -0.544 1.242 -1.262',
	'-0.595 0.009 0.672 -2.128 -0.184',
	'0.260 0.830 3.097 -0.838 1.512',
	'0.000 0.000 0.000 0.000 0.000']
	#
	list_of_AAs=['A','C', 'D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','Z']
	#
	AA_to_Atchley = dict()
	for (AA, row) in zip(list_of_AAs, Atchley_factor_data): 
		for (ii, entry) in enumerate(row.split(" "), start =1):
			AA_to_Atchley[AA, ii] = float(entry)
	#
	return AA_to_Atchley


def Atchley_euclidean_dist(s1,s2, AA_to_Atch = make_AA_to_Atchley_dict(), AF_list = [1,2,3,4,5]):
	"""
	Returns the distance calculated as euclidean distance of average AFs for each CDR3
	E.g. each AF is calculated as average AF value for the CDR3 resulting in a 5-tuple 
	corresponding to each AF. For 2 CDR3s the distance between them is calculated as 
	euleadian distance between the two 5-tuples.
	"""
	s1_list = np.array([sum(AA_to_Atch[AA, AF] for AA in list(s1))/float(len(s1)) for AF in AF_list])
	s2_list = np.array([sum(AA_to_Atch[AA, AF] for AA in list(s2))/float(len(s2)) for AF in AF_list])
	distance = np.sqrt(np.sum((np.array(s1_list)-np.array(s2_list))**2))
	del s1_list
	del s2_list
	return distance


def BLOSUM45_score_dist(s1,s2):
	aligner = Align.PairwiseAligner()
	aligner.open_gap_score = -10
	aligner.substitution_matrix = substitution_matrices.load("BLOSUM45")
	aligner.mode = "global"
	score_s12 = aligner.score(s1,s2)
	score11 = aligner.score(s1,s1)
	score22 = aligner.score(s2,s2)
	distance = 1- score_s12/max(score11,score22)
	return distance


##########################################################

def choose_distance_function(dist_option): 
	if dist_option == "AF_euclidean": 
		dist_func = Atchley_euclidean_dist
	elif dist_option == "Levenshtein":
		dist_func = levenshteinDistance
	elif dist_option == "BLOSUM45_score_dist":
		dist_func = BLOSUM45_score_dist
	return dist_func

############################################################
# 1.2. Extract sequence list  
##########################################################

def read_in_seq_freq_list(filtered_input_filename):
	"""
	Reads in the filtered sample line by line; 
	appends each cdr3 to a list 
	"""
	with open(filtered_input_filename,'r') as f:
		lines = [line.strip('\n').split('\t') for line in f if not line.startswith('#')]
	#
	clones = list()
	for line in lines: 
		clones.append(line[0])
	#    
	return clones
	

def dist_i_list_parallel(cdr3i, cdr3s_list, distance_function_used):
	dist_i = np.array([distance_function(cdr3i, cdr3j) for cdr3j in cdr3s_list])
	return dist_i


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def write_out_dist_i_chunk(chunk_enumerator, sample_name, distance_line):
	"""
	Adds a distance row for each cdr3_i in the format X X X d(cdr3_i, cdr3_j)
	"""
	with open(chunk_enumerator + "_" + sample_name + "_dmat.tsv", 'a+b') as f:
		np.savetxt(f, distance_line, delimiter = "", newline = "\t")
		last_tab_in_line_size_in_bytes = -1
		f.seek(last_tab_in_line_size_in_bytes, 2)
		f.truncate()
		f.write(b"\n")
	del distance_line




##########################################################
# 2.1. Executable
##########################################################

if __name__ == '__main__':
	start_overall_time = time.time()
	filtered_input, chunk_size, chunk_num, distance_option  = parse_commandline()
	sample_name = os.path.splitext(os.path.basename(filtered_input))[0]
	distance_function = choose_distance_function(distance_option)
	cdr3_list = read_in_seq_freq_list(filtered_input)
	cdr3_chunk = list(chunks(cdr3_list, chunk_size))[chunk_num]
	# make sure you calculate self score if any scoring matrix is used: 
	start_parallel_time = time.time()
	#
	#this doesn't work on the cluster because it counts all the cores on a node 
	#doesn't limit to the number of cores requested in the job script
	#cores = mp.cpu_count() 
	cores = 8
	Que = list()
	if os.path.isfile(str(chunk_num) + "_" + sample_name + "_dmat.tsv"):
		os.remove(str(chunk_num) + "_" + sample_name + "_dmat.tsv")
	#	
	pool = mp.Pool(processes = cores, maxtasksperchild = 1000)
	for ii, cdr3_i in enumerate(cdr3_chunk):
		print(ii, cdr3_i)
		process = pool.apply_async(dist_i_list_parallel, (cdr3_i, cdr3_list, distance_function))
		Que.append(process)
		#del process
	for oo in range(len(Que)):
		function_output = Que[oo].get()
		write_out_dist_i_chunk(str(chunk_num), sample_name, function_output)
		del function_output
	pool.terminate()
	del Que
	gc.collect()
	print("parallel calculating time")
	print("--- %s seconds ---" % (time.time() - start_parallel_time))
	print("overall used time")
	print("--- %s seconds ---" % (time.time()- start_overall_time))
