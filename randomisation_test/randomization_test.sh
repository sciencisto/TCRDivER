#!/bin/sh 
### General options 
### -- specify queue -- 
#BSUB -q hpc
### -- set the job Name -- 
#BSUB -J randomization_calc 
### -- ask for number of cores (default: 1) -- 
#BSUB -n 1
### -- specify that the cores must be on the same host -- 
#BSUB -R "span[hosts=1]"
### -- specify that we need 2GB of memory per core/slot -- 
#BSUB -R "rusage[mem=2GB]"
### -- specify that we want the job to get killed if it exceeds 3 GB per core/slot -- 
#BSUB -M 3GB
### -- set walltime limit: hh:mm -- 
#BSUB -W 10:00 
### -- set the email address -- 
# please uncomment the following line and put in your e-mail address,
# if you want to receive e-mail notifications on a non-default address
#BSUB -u 
### -- send notification at start -- 
#BSUB -B 
### -- send notification at completion -- 
#BSUB -N 
### -- Specify the output and error file. %J is the job-id -- 
### -- -o and -e mean append, -oo and -eo mean overwrite -- 
#BSUB -o Output_%J.out 
#BSUB -e Error_%J.err 
# -- run in the current working (submission) directory --
if test X$BSUB_ENVIRONMENT = XBSUB_BATCH; then cd $BSUB_O_WORKDIR; fi
# here follow the commands you want to execute
module load python3/3.6.2 ;
module load pandas/0.20.3-python-3.6.2 ;
python3 /work3/milvu/divP_parallel/randomize_frequencies_filtered_sample.py "/work3/milvu/2018_Formenti_Demaria_Lung_50000_random_freqs/" "randomize" 