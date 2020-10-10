#!/bin/sh 
### General options 
### -- specify queue -- 
#BSUB -q hpc
### -- set the job Name -- 
#BSUB -J divP_calcs
### -- ask for number of cores (default: 1) -- 
#BSUB -n 1 
### -- specify that the cores must be on the same host -- 
#BSUB -R "span[hosts=1]"
### -- specify that we need 2GB of memory per core/slot -- 
#BSUB -R "rusage[mem=1GB]"
### -- specify that we want the job to get killed if it exceeds 3 GB per core/slot -- 
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
# 			delete the individual Error, Output and job_chunk files
# 			make divP jobs
#			run divP jobs	
#

##############################
# specify working dir
##############################

w_dir_path="/work3/milvu/Chain_OVAexpt_data1_B45_50000_random_subsmpl_random_freqs/"

chunk_size="100"

line_beginning_ID="python3 /work3/milvu/divP_parallel/calculate_divP.py "$w_dir_path" "
line_ending_ID=" "ID" "[0,1,2,3,4,5,6,inf]" "[0.1,0.2,0.25,0.3,0.4,0.5,0.75,1.0,1.5,2.0,4.0,8.0,16.0,32.0,64.0]" "

line_beginning_SIM_2="python3 /work3/milvu/divP_parallel/calculate_divP.py \$TMPDIR "
line_ending_SIM_2=" "SIM" "[0,1,2,3,4,5,6,inf]" "[0.1,0.2,0.25,0.3,0.4,0.5,0.75,1.0,1.5,2.0,4.0,8.0,16.0,32.0,64.0]" "



for D in "$w_dir_path"*;
do if [ -d "${D}" ]; 
    then  cd "$D"; 
    # delete Error/Output/jobchunk files 
    rm Error* Output* job_chunk*sh
    sample_name=${D##*/};
    filename="$sample_name""_filtered.tsv";
    # make ID jobs
    job_ID_filename="job_ID_"$sample_name".sh";
    line_ID=$line_beginning_ID$sample_name""$line_ending_ID$chunk_size;
    cat /work3/milvu/divP_parallel/templates/job_template_ID.txt > $job_ID_filename;
    echo $line_ID >> $job_ID_filename;
    # make SIM jobs
    job_SIM_filename="job_SIM_"$sample_name".sh";
    line_SIM_1="rsync -av "$w_dir_path""$sample_name" \$TMPDIR"
    line_SIM_2="# do the processing"
    line_SIM_3=$line_beginning_SIM_2$sample_name""$line_ending_SIM_2$chunk_size;
    line_SIM_4="#delete our "copy" from /tmp"
    line_SIM_5="rm -rfv \$TMPDIR"
    cat /work3/milvu/divP_parallel/templates/job_template_SIM.txt > $job_SIM_filename;
    echo $line_SIM_1 >> $job_SIM_filename;
    echo $line_SIM_2 >> $job_SIM_filename;
    echo $line_SIM_3 >> $job_SIM_filename;
    echo $line_SIM_4 >> $job_SIM_filename;
    echo $line_SIM_5 >> $job_SIM_filename;
    #run all divP jobs
    bsub < $job_ID_filename;
    bsub < $job_SIM_filename;
fi; 
cd ../; 
done


