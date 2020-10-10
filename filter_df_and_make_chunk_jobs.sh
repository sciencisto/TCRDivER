#!/bin/sh 
### General options 
### -- specify queue -- 
#BSUB -q hpc
### -- set the job Name -- 
#BSUB -J filter_and_chunk
### -- ask for number of cores (default: 1) -- 
#BSUB -n 1 
### -- specify that the cores must be on the same host -- 
#BSUB -R "span[hosts=1]"
### -- specify that we need 2GB of memory per core/slot -- 
#BSUB -R "rusage[mem=2GB]"
### -- specify that we want the job to get killed if it exceeds 3 GB per core/slot -- 
#BSUB -M 2GB
### -- set walltime limit: hh:mm -- 
#BSUB -W 1:00 
### -- set the email address -- 
# please uncomment the following line and put in your e-mail address,
# if you want to receive e-mail notifications on a non-default address
#BSUB -u milvu@dtu.dk
### -- send notification at start -- 
#BSUB -B 
### -- send notification at completion -- 
#BSUB -N 
### -- Specify the output and error file. %J is the job-id -- 
### -- -o and -e mean append, -oo and -eo mean overwrite -- 
#BSUB -o Output_%J.out 
#BSUB -e Error_%J.err 
# here follow the commands you want to execute
module load python3/3.5.4 ;
module load pandas/0.20.3-python-3.5.4 ;
python3 filter_df_and_make_jobs.py "/work3/milvu/Data/CFA_OVA_p277/Cinelli_Chain_OVAexpt_data1/" "random" "50000" "None" "/work3/milvu/divP_parallel/" "100" "/work3/milvu/Chain_OVAexpt_data1_B45_50000_random_subsmpl22/" "BLOSUM45_score_dist"
