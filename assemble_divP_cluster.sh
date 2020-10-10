#!/bin/sh 
### General options 
### -- specify queue -- 
#BSUB -q hpc
### -- set the job Name -- 
#BSUB -J assemble_divP 
### -- ask for number of cores (default: 1) -- 
#BSUB -n 1 
### -- specify that the cores must be on the same host -- 
#BSUB -R "span[hosts=1]"
### -- specify that we need 2GB of memory per core/slot -- 
#BSUB -R "rusage[mem=1GB]"
### -- specify that we want the job to get killed if it exceeds 3 GB per core/sl
#BSUB -M 1GB
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
# This script does the following: 
#			concatenate all the chunk dmat files in the working dir
# 			delete the individual chunk dmat files 
# 			make divP jobs
#			run divP jobs	
#

##############################
# specify working dir
##############################

w_dir_path="/work3/milvu/Chain_OVAexpt_data2_20000_highest_subsmpl/"

for D in "$w_dir_path"*;
do if [ -d "${D}" ]; 
    # concatenate and delete chunk dmat 
    then  cd "$D"; 
    sample_name=${D##*/};
    filename="$sample_name""_divP.tsv";
    job_divP_filename="job_divP_"$sample_name".sh";
    file_line="python3 /work3/milvu/divP_parallel/assemble_divP.py "$w_dir_path" "$sample_name;
    cat /work3/milvu/divP_parallel/templates/job_template_divP.txt > $job_divP_filename;
    echo $file_line >> $job_divP_filename;
    bsub < $job_divP_filename;
fi;
cd ../; 
done
