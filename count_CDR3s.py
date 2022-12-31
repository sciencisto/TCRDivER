"""
Counting script for all the samples: 

Read the folder
Take in all the samples in the folder
Count the number of CDR3s in a dictionary 
Write out a file for each sample with the CDR3 and the count 

USAGE: 

python3 count_CDR3s.py /path_to_wdir/ 


"""

from sys import argv 
import pandas as pd 
import numpy as np      #
import os
from os import listdir
from os.path import isfile, join
from collections import Counter


##########################################################
# 1.1. Parse commandline 
##########################################################

def parse_commandline():
	"""
	ARGUMENTS: path_to_wdir
	"""
	wdir_path = argv[1] #path to files
	return wdir_path


def make_dict_from_file(wdir_path, fname):
	with open(wdir_path+fname, "r") as f: 
		lines = [line.strip("\n") for line in f.readlines()]
	cdr3_count_dict = dict(Counter(lines))
	return cdr3_count_dict


def write_out_CDR3_tsv_file(wdir_path, fname, cdr3_count_dict):
	fname = fname.strip(".txt")
	with open(wdir_path+fname+".tsv", "w") as f: 
		f.writelines( "aminoAcid"+"\t"+ "count (reads)" + "\t" "sequenceStatus"+"\n")
		for cdr3,count in sorted(cdr3_count_dict.items()): # here added sorted
			f.writelines('\t'.join([cdr3,str(count),"In"]) + '\n')

##########################################################
# 2.1. Executable
##########################################################

if __name__ == '__main__':
	# READ IN WORKING DIR PATH.
	wdir_path = parse_commandline()
	# READ IN NUMBER OF FILES
	files = [f for f in listdir(wdir_path) if isfile(join(wdir_path, f))]
	files = [f for f in files if f.endswith('cdr3.txt')]
	files.sort()
	for file_name in files:
		# MAKE THE DICT 
		cdr3_count_dict = make_dict_from_file(wdir_path, file_name)
		# WRITE OUT FILE 
		write_out_CDR3_tsv_file(wdir_path, file_name, cdr3_count_dict)
	# EXIT





