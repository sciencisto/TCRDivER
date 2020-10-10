"""
READ ME 

INPUT:  folder of .tsv T cell sequencing files from Adaptive Biotechnologies. 
OUTPUT: Folder with subfoders for each sample. 
		In each subfolder: 
			filtered_"sample_name".tsv file containing In frame cdr3 seq and freq in each row. Freq are normalized. 
 			(# seqs/ chunk size) job scripts calling calculate_dmat.py

USAGE: 
python3 filter_df_and_make_jobs.py 
									path_to_files_folder
									subsample_option kk_subsample_size filter_singlets
									path_to_scripts_folder chunk_size 
									path_to_working_dir distance_option

OPTIONS: 

Subsampling options: None, highest, lowest, random 

Subsample size: None, or int

distance option: "AF_euclidean" "Levenshtein" "BLOSUM45_score_dist"

"""

from sys import argv 
from collections import defaultdict
import pandas as pd 
import numpy as np
import itertools
import heapq            #
import os
from os import listdir
from os.path import isfile, join
from subprocess import Popen, PIPE
import time 
import gc
from collections import Counter

##########################################################
# 1.1. Parse commandline 
##########################################################

def parse_commandline():
	"""
	ARGUMENTS: input_filepath, subsample_option, kk, filter_singlets, path_to_job_sc, chunk_size, path_to_working_dir, distance_option
	subsample_option: subsampling method
		if subsample_option = None doesn't subsample
		if subsample_option = 'highest' subsamples kk highest clones based on frequency
		if subsample_option = 'lowest' subsamples kk lowest clones based on frequency
		if subsample_option = 'random' subsamples kk random clones based on frequency distribution
		"""
	input_filepath = argv[1] #path to files
	#sample_name = os.path.splitext(os.path.basename(input_filename))[0]
	subsample_option = argv[2] # subsample_option 
	if subsample_option == "None":
		subsample = None
	if argv[3] == "None":       # Subsample size: None/number  
		kk = None
	else: 
		kk = int(float(argv[3]))
	filter_singlets = argv[4]
	path_to_scripts_folder = argv[5]
	chunk_size = int(float(argv[6]))
	path_to_working_dir = argv[7]
	distance_option = argv[8]
	return input_filepath, subsample_option, kk, filter_singlets, path_to_scripts_folder, chunk_size, path_to_working_dir, distance_option


############################################################
# 1.2. Extract sequence list  
##########################################################

def make_dataframe_from_tsv(input_filename):
	"""
	Returns a pandas dataframe from Adaptive Biotechnologies sequencing .tsv file
	"""
	adaptive_bio_df = pd.read_csv(input_filename, sep = '\t')
	return adaptive_bio_df


def choose_column_names(df):
	if 'amino_acid' in list(df.columns.values):
		type_of_seq = 'amino_acid'
		type_of_freq = 'productive_frequency'
	else: 
		type_of_seq = 'aminoAcid'
		type_of_freq = 'count (reads)'  
	type_of_frame = 'In'
	if type_of_freq not in list(df.columns.values): 
		type_of_freq = 'count (templates/reads)'
	return type_of_seq, type_of_freq, type_of_frame


def normalize(d, target=1.0):
	"""
	Returns a transformed dictionary with the values normalized
	"""
	raw = sum(d.values())
	factor = target/raw
	return {key:value*factor for (key,value) in d.items()}



def make_cdr3_freq(adaptive_bio_df, type_of_seq, type_of_freq, type_of_frame, subsample_choice, kk_choice, filter_singlets):
	"""
	INPUT: 
	type_of_seq  string, can take values  "rearrangement" for nucleic acid or "amino_acid" in the new files, or "aminoAcid" in the old files 
	type_of_frame string,  can take values "In", "Stop", "Out"
	type_of_freq string, can take values: "frequency", "productive_frequency"
	adaptive_bio_df pandasDataFrame from Adaptive BioTech to be filtered based on type of seq, frame and freq. 
	subsample_choice None or string: "highest", "lowest", "random"
	kk int number of seqs to be sampled 

	OUTPUT: 
	Returns : Subsampled dictionary of CDR3 seq as key, and corresponding productive freqs as values, normalized from 0 to 1.
	"""
	cdr3_freq = defaultdict(float)
	# 
	if type_of_seq == 'amino_acid':
		df = adaptive_bio_df[(adaptive_bio_df.frame_type == type_of_frame)]
		if filter_singlets == "filter_singlets":
			df = df[(df.templates != 1)]
	elif type_of_seq == 'aminoAcid':
		df = adaptive_bio_df[(adaptive_bio_df.sequenceStatus == type_of_frame)]
	del adaptive_bio_df
	for index, row in df.iterrows():
		cdr3_freq[row[type_of_seq]] += row[type_of_freq]
	if filter_singlets == "filter_singlets":
		cdr3_freq = {key:value for (key, value) in sorted(cdr3_freq.items()) if value != 1.0} #here added sorted
	mean_subsample_before_normalization = np.array(list(cdr3_freq.values())).mean()
	var_subsample_before_normalization = np.array(list(cdr3_freq.values())).var()
	cdr3_freq = normalize(cdr3_freq, target=1.0)
	del df
	frequencies = []
	clones = []
	if subsample_choice and kk_choice:
		if subsample_choice == 'highest': 
			kk_keys_sampled = heapq.nlargest(kk_choice, cdr3_freq, key=cdr3_freq.get)
			for key in sorted(kk_keys_sampled):						# here added sorted
				frequencies = np.append(frequencies, cdr3_freq[key])
				clones.append(key)
		elif subsample_choice == 'lowest':
			kk_keys_sampled = heapq.nsmallest(kk_choice, cdr3_freq, key=cdr3_freq.get)
			for key in sorted(kk_keys_sampled):						# here added sorted
				frequencies = np.append(frequencies, cdr3_freq[key])
				clones.append(key)
		elif subsample_choice == 'random':
			np.random.seed(seed=1)
			kk_keys_sampled = np.array(np.random.choice(list(sorted(cdr3_freq.keys())), size=kk_choice, replace=True, p=list(sorted(cdr3_freq.values()))))
			cdr3_freq_subsample = dict(Counter(kk_keys_sampled))
			try:
				assert sum(cdr3_freq_subsample.values()) == kk_choice
			except AssertionError:
				print("Not the same subsample size and sample size choice")
		if subsample_choice != "random":
			try: 
				assert abs(1 - sum(frequencies)) < min(frequencies)
			except AssertionError:
				print("Sum of frequencies deviates from 1 by more than the min frequency in sample!")
			cdr3_freq_subsample = dict()
			for (clone, freq) in zip(clones, frequencies):
				cdr3_freq_subsample[clone] = freq
		cdr3_freq_subsample = normalize(cdr3_freq_subsample, target=1.0)
		cdr3_freq = cdr3_freq_subsample
	else: 
		try: 
			assert abs(1 - sum(cdr3_freq.values())) < min(cdr3_freq.values())
		except AssertionError:
			print("Sum of frequencies deviates from 1 by more than the min frequency in sample!")
	return cdr3_freq, mean_subsample_before_normalization, var_subsample_before_normalization


def write_out_filtered_sample(directory_path, sample_name, cdr3_freq):
	"""
	Writes out dictionary from the subsampled filtered dataframe. 
	Each row: cdr3 /t freq
	"""
	with open(directory_path + "/" + sample_name + "_filtered.tsv", 'w') as f:
		for cdr3,freq in sorted(cdr3_freq.items()): # here added sorted
			f.writelines('\t'.join([cdr3,str(freq)]) + '\n')
	del cdr3_freq

def write_out_mean_and_var(directory_path, sample_name, mean_subsample_before_normalization, var_subsample_before_normalization):
	with open(directory_path + "/" + sample_name + "_mean_var.tsv", 'w') as f: 
		f.writelines('\t'.join(["mean","variance"]) + '\n')
		f.writelines('\t'.join([str(mean_subsample_before_normalization),str(var_subsample_before_normalization)]) + '\n')
	del mean_subsample_before_normalization, var_subsample_before_normalization

	
def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


##########################################################
# 2.1. Executable
##########################################################

if __name__ == '__main__':
	start_overall_time = time.time()
	# READ IN PATH TO FILE.
	filepath, subsample, kk, filter_singlets,script_folder_path, cc_size, w_dir_path, distance_option = parse_commandline()
	# READ IN NUMBER OF FILES
	files = [f for f in listdir(filepath) if isfile(join(filepath, f))]
	files = [f for f in files if f.endswith('tsv')]
	files.sort() 
	make_big_folder = Popen(['mkdir', w_dir_path], stdout = PIPE, stderr = PIPE) # commandline options in one script
	make_big_folder.wait() # wait until it makes the folder 
	for file_name in files:
		start_file_time = time.time()
		sample_name = os.path.splitext(os.path.basename(file_name))[0]
		sample_dir_path = "".join([w_dir_path, sample_name])
		print(w_dir_path) 
		# MAKE A FOLDER
		make_w_dir = Popen(['mkdir', sample_dir_path], stdout=PIPE, stderr=PIPE)
		make_w_dir.wait()
		# DO THE FILTERING 
		input_filename = "".join([filepath,file_name])
		adaptive_df = make_dataframe_from_tsv(input_filename)
		type_of_seq, type_of_freq, type_of_frame = choose_column_names(adaptive_df)
		cdr3_freq_dict, mean_subsample_before_normalization, var_subsample_before_normalization = make_cdr3_freq(adaptive_df, type_of_seq, type_of_freq, type_of_frame, subsample, kk, filter_singlets)
		del adaptive_df
		# WRITE OUT FILTERED DICT
		write_out_filtered_sample(sample_dir_path, sample_name, cdr3_freq_dict) # writes out a dictionary of selected freqs and cdr3s 
		write_out_mean_and_var(sample_dir_path, sample_name, mean_subsample_before_normalization, var_subsample_before_normalization)
		# CHUNK OUT AND MAKE THE JOB FILES
		number_of_chunks = len(cdr3_freq_dict.keys()) // cc_size
		if (len(cdr3_freq_dict.keys()) % cc_size) != 0:
			number_of_chunks += 1
		del cdr3_freq_dict
		for cc in range(number_of_chunks):
			print(cc)
			job_filename = sample_dir_path + "/job_chunk_"+str(cc)+"_calc_dmat.sh" 
			#make chunk job line 
			chunk_job_line = "python3 " + script_folder_path + "calculate_dmat.py " + sample_name +"_filtered.tsv "+ str(cc_size)+ " " + str(cc)+ " " + distance_option
			#append to job script
			job_sc_name = script_folder_path + "templates/job_parallel_chunk_dmat_template.txt"
			make_job_chunk_script = Popen(['cp', job_sc_name, job_filename], stdout = PIPE, stderr = PIPE) # commandline options in one script
			make_job_chunk_script.wait() # wait until it copies the file 
			with open(job_filename, 'a') as f_job:
				f_job.write(chunk_job_line)
			#
		print("file_time", file_name)
		print( time.time()- start_file_time)
	print("overall time")
	print(time.time() - start_overall_time)
	# EXIT

