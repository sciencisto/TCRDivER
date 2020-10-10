#########################################################################
# ImmuneDivER algorithm
#########################################################################

Thank you for the interest and using ImmuneDivER algorithm! 

We are constantly improving the code and if you wish to contribute or
have any questions feel free to contact. 

Milena Vujovic (milvu@dtu.dk)


The software runs in a Python3.6 environment. We have tested in on UNIX 
based systems. If you have any issues running the software, again, 
feel free to contact us- we are more than happy to help!

The algorithm is split up into several parts. Each script is designed 
to perform a step of the algorithm. Since these calculations can be 
computationally expensive, we have designed it like this so that one runs 
into the least ammount of problems. 

If you are running any of the scripts locally, the usage is: 

python3.6 [python_script.py] [arguments]

IMPORTANT: The code is optimised for working remotely on a cluster system.
We have run the scripts on the DTU university cluster, and have included 
the bash scripts used to make and run these calculations in an automated 
fashion. 

If you are good in bash you will be able to augment the bash scripts to suit 
your needs. If not - shoot us a line --> we might be able to help :) 



The algorithm consits of several parts: 
#########################################################################
I. Filter sample and make distance matrix chunk jobs: 
#########################################################################

	1. A folder containing all sample files is read in. 

	2. A new working directory is created in a specified path. 

	3. For each file a new subdirectory is created in the working directory. 

	4. Each sample is read in as a list of CDR3s with their respective 
	frequencies or counts. Preferably Adaptive Biotechnologies format. 
	It takes in only the "In" frame sequences of CDR3s.

	5. Each sample is filtered for the "In" frame sequences. Depending on the 
	option selected the counts of "In" frame CDR3s can stay the same (no 
	subsampling), or a subsampling technique can be employed: 
		- highest: takes in a number of highest expanded CDR3s
		- lowest: takes in a number of lowest expanded CDR3s. 
		- random: takes in a number of random CDR3s. 

	6. The subsample is then normalized so the frequencies of the subsample 
	sum up to 1. And bash files containing the command to run individual 
	distance matrix chunk calculatons are created. In this step distance measures
	are chosen such as BLOSUM45, Levenstein distance and AF euclidean distance.  
	When calculated all the distance matrix chunks can be concatenaded into the 
	full size distance matrix containing all the pairwise comparisons between 
	CDR3s in the subsample. 


##############################
USAGE:
##############################

RUN: Bsub < filter_df_and_make_chunk_jobs.sh

Opens the folder with all the samples. Makes a new working directory with sample subdirectories. Calls python 
script to filter based on filtering criteria and makes chunk distance matrix jobs. Inside the bash script 
filter_df_and_make_chunk_jobs.sh you can change arguments for the python script as follows: 

Python script command: 

python3 filter_df_and_make_jobs.py 
									- sample_data_dir 
									- subsampling option: None, highest, lowest, random 
									- subsample size: None or integer
									- path_to_scripts_folder: this is the templates folder 
									- chunk_size: integer 
									- new_working_dir_name
									- distance option: "AF_euclidean" "Levenshtein" "BLOSUM45_score_dist"

e.g.

python3 filter_df_and_make_jobs.py 	- "/work3/milvu/Data/2018_Formenti_Demaria_Lung/" 
									- "random" 
									- "50000" 
									- "None" 
									- "/work3/milvu/divP_parallel/" 
									- "100" 
									- "/work3/milvu/2018_Formenti_Demaria_Lung_50000_random_dist_random_freq/" 
									- "BLOSUM45_score_dist"

#########################################################################
II. Run all distance matrix calculations. 
#########################################################################

	In order to run all the distance matrix calculations it is enough to run the bash scripts. 
	Here is a neat bach loop to do it for one sample: 

	for D in *;do if [ -d "${D}" ]; then  cd "$D"; for jb in *sh; do bsub < $jb;done;fi;cd ../; done

	This will submit all the "calculate chunk job scripts" within a sample folder. If you want to run all 
	samples at once do another level in the loop. 


##############################
USAGE:
##############################

Navigate to the working directory and RUN:  

for D in *;do if [ -d "${D}" ]; then  cd "$D"; for jb in *sh; do bsub < $jb;done;fi;cd ../; done



#########################################################################
III. Make and run true diversity calculations. 
#########################################################################

Now that you've calculated the distance matrix, you can proceed with calculating the tru divesity with 
varying values of lambda and q. A bash script takes in the whole working directory and makes naive and 
similarity scaled diversity jobs. The outputs of these two scripts are two .tsv files containing a 
dataframe with the values of naive diversity with varying q and similarity scaled diversity with varying
q and lambda values denoted [sample_name]_ID.tsv and [sample_name]_SIM.tsv, respectively. 


Inside the make_n_run_divP_calcs.sh bash script you can set individual parameters for the diversity calculations. Specify
the wdir folder and the chunk size. NOTE: MAke sure that the chunk size is the same as the size that you have specified 
for the distance matrix chunk calculation. 

Arguments to change: 
	# specify the working directory created by the filter_df_and_make_chunk_jobs.sh
	w_dir_path="/work3/milvu/Chain_OVAexpt_data1_B45_50000_random_subsmpl_random_freqs/"

	#specify chunk size: same as filter_df_and_make_chunk_jobs.sh value 
	chunk_size="100"

	# specify the path to the python script: 
	line_beginning_ID="python3 /work3/milvu/divP_parallel/calculate_divP.py "$w_dir_path" "
	
	# specify the q values and the lambda values you want calculated for the naive diversity calculation 
	line_ending_ID=" "ID" "[0,1,2,3,4,5,6,inf]" "[0.1,0.2,0.25,0.3,0.4,0.5,0.75,1.0,1.5,2.0,4.0,8.0,16.0,32.0,64.0]" "

	# specify the q values and the lambda values you want calculated for the similarity scaled diversity calculation 
	line_beginning_SIM_2="python3 /work3/milvu/divP_parallel/calculate_divP.py \$TMPDIR "
	line_ending_SIM_2=" "SIM" "[0,1,2,3,4,5,6,inf]" "[0.1,0.2,0.25,0.3,0.4,0.5,0.75,1.0,1.5,2.0,4.0,8.0,16.0,32.0,64.0]" "



##############################
USAGE: 
##############################

RUN: bsub < make_n_run_divP_calcs.sh

These jobs execute:  

python3 /zhome/cf/1/129355/divP_parallel/calculate_divP.py  - path_to_wdir
															- sample_name 
															- similarity_option: "ID"/"SIM"
															- qs:  "[0,1,2,3,4,5,6,inf]" 
															- lambdas: "[0.1,0.2,0.25,0.3,0.4,0.5,0.75,1.0,1.5,2.0,4.0,8.0,16.0,32.0,64.0]" 
															- chunk_size: "100"


e.g. 
python3 /work3/milvu/divP_parallel/calculate_divP.py 	- "$w_dir_path" 
														- "$sample_name" 
														- "ID" or "SIM"
														- "[0,1,2,3,4,5,6,inf]" 
														- "[0.1,0.2,0.25,0.3,0.4,0.5,0.75,1.0,1.5,2.0,4.0,8.0,16.0,32.0,64.0]" 
														- 100

#########################################################################
IV. Assemble complete diversity calculation 
#########################################################################

One the naive and similarity scaled diversity calculations are completed each sample folder will contain two files: 
[sample_name]_ID.tsv and [sample_name]_SIM.tsv calculations. These are now assembled into an overall file of diversity 
using the assemble_divP_cluster.sh bash script. 

It runs the python script assemble_divP.py which takes the two files [sample_name]_ID.tsv and [sample_name]_SIM.tsv 
as input and outputs a complete diversity dataframe named [sample_name]_divP.tsv. 


You need to specify the working folder where all samples with the calculated diversity files are and the path to the 
python script assemble_divP.py. 

# specify working directory path - same as folder that is created in the filtering step (I.)
w_dir_path="/work3/milvu/2018_Formenti_Demaria_Lung_50000_random_freqs/"


for D in "$w_dir_path"*;
do if [ -d "${D}" ];
    # concatenate and delete chunk dmat 
    then  cd "$D";
    sample_name=${D##*/};
    filename="$sample_name""_divP.tsv";
    job_divP_filename="job_divP_"$sample_name".sh";
    # specify path to assemble_divP.py script
    file_line="python3 /work3/milvu/divP_parallel/assemble_divP.py "$w_dir_path" "$sample_name;
    cat /work3/milvu/divP_parallel/templates/job_template_divP.txt > $job_divP_filename;
    echo $file_line >> $job_divP_filename;
    bsub < $job_divP_filename;
fi;
cd ../;
done

##############################
USAGE: 
##############################

RUN: bsub < assemble_divP_cluster.sh

The job executes: 

python3 [path_to_script]/assemble_divP.py   working_directory_path
											sample name


e.g.
python3 work3/milvu/divP_parallel/assemble_divP.py 
		"$w_dir_path" 
		"$sample_name"



























