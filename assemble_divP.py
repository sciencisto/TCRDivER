"""
2019_11_08MV
Assemble divP by taking the SIM and ID files for each sample. 
Output: assembled divP_file


"""

from sys import argv
from collections import defaultdict

def parse_commandline():
	"""
	ARGUMENTS: 	path_to_file
				 sample_name 

		"""
	path_to_working_dir = argv[1] #path to files
	#sample_name = os.path.splitext(os.path.basename(input_filename))[0]
	sample_name = argv[2] # subsample_option 
	return path_to_working_dir, sample_name


def read_in_SIM_calc(path_to_working_dir, sample_name): 
	with open(path_to_working_dir + sample_name + "/" + sample_name + "_SIM.tsv", "r") as fl: 
		lines = [line.strip("\n").split("\t") for line in fl.readlines()]
	divP_SIM = defaultdict()
	list_of_qs = list()
	for line in lines[1:]:
		divP_SIM[line[0]] = line[1:]
		list_of_qs.append(line[0])
	return divP_SIM, list_of_qs


def read_in_ID_calc(path_to_working_dir, sample_name):
	with open(path_to_working_dir + sample_name + "/" + sample_name + "_ID.tsv", "r") as fl: 
		lines = [line.strip("\n").split("\t") for line in fl.readlines()]
	divP_ID = defaultdict()
	for line in lines[1:]:
		divP_ID[line[0]] = line[1:]
	return divP_ID

def write_divP(path_to_working_dir, sample_name, divP_SIM, divP_ID, list_of_qs):
	divP = defaultdict(list)
	for key in divP_SIM.keys():
		divP[key] = divP_SIM[key]+divP_ID[key]
	with open(path_to_working_dir + sample_name +"/" + sample_name +  "_divP" + ".tsv", 'w') as f:
		for qu in list_of_qs:
			row = [qu] + divP[qu]
			f.writelines('\t'.join(row) + '\n')


###### EXECUTABLE ####### 

if __name__ == '__main__':
	path_to_working_dir, sample_name = parse_commandline()
	diversity_SIM, qs = read_in_SIM_calc(path_to_working_dir, sample_name)
	diversity_ID = read_in_ID_calc(path_to_working_dir, sample_name)
	write_divP(path_to_working_dir, sample_name, diversity_SIM, diversity_ID, qs)



