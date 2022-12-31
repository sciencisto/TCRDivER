#!/bin/sh 
### General options 
### -- specify queue -- 
#BSUB -q hpc
### -- set the job Name -- 
#BSUB -J TCRDivER
### -- ask for number of cores (default: 1) -- 
#BSUB -n 1 
### -- specify that the cores must be on the same host -- 
#BSUB -R "span[hosts=1]"
### -- specify that we need 2GB of memory per core/slot -- 
#BSUB -R "rusage[mem=2GB]"
### -- specify that we want the job to get killed if it exceeds 3 GB per core/slot -- 
#BSUB -M 10GB
### -- set walltime limit: hh:mm -- 
#BSUB -W 24:00 
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
# here follow the commands you want to execute


##############################
# Read in config file 
##############################

. ./TCRDivER.config 

##############################
# load modules
##############################

module load python3/3.6.2 ;
module load pandas/0.20.3-python-3.6.2 ;
module load matplotlib/2.0.2-python-3.6.2 ;

###################################################################
# Execute TCRDivER
###################################################################

# Pre process samples if needed:

if [ "$count" = "count" ]
then  
  # Count CDR3s from txt files: 
  python3 count_CDR3s.py $data_dir_path
fi

wait 
echo "COUNT DONE"
#####################################

python3 filter_df_and_make_jobs.py $data_dir_path $filtering_mode $subsample_size  $filter_singlets  "$TCRDivER_dir" $chunk_size $w_dir_path  $distance

wait
echo " FILTER DONE"
# Run all distance matrix chunk jobs
#####################################

for D in "$w_dir_path"*;
do if [ -d "${D}" ];
    then  cd "$D";
    for jb in *sh;
    do bsub < $jb;
    done;
fi;
cd ../;
done

while [[ $(bstat | awk '$4 ~ /_dmat_job/ { print }') ]]; do wait ; done
echo " DISTANCE MATRIX CALCULATION DONE"
## Make and run divP calculations
######################################

# Create commands to insert into job templates

line_beginning_ID="python3 "$TCRDivER_dir"calculate_divP.py "$w_dir_path" "
line_ending_ID=" 'ID' "$q_values" "$lambda_values" "
line_beginning_SIM_2="python3 "$TCRDivER_dir"calculate_divP.py "$w_dir_path" "
line_ending_SIM_2=" 'SIM' "$q_values" "$lambda_values" "

# Assemble and run jobs


for D in "$w_dir_path"*;
do if [ -d "${D}" ]; 
    then  cd "$D"; 
    # delete Error/Output/jobchunk files 
    rm "${D}/"Error* "${D}/"Output* "${D}/"job_chunk*sh
    sample_name=${D##*/};
    filename="$sample_name""_filtered.tsv";
    # make ID jobs
    job_ID_filename="job_ID_"$sample_name".sh";
    line_ID=$line_beginning_ID$sample_name""$line_ending_ID$chunk_size;
    cat $TCRDivER_dir""templates/job_template_ID.txt > $job_ID_filename;
    echo $line_ID >> $job_ID_filename;
    # make SIM jobs
    job_SIM_filename="job_SIM_"$sample_name".sh";
    line_SIM_1="rsync -av "$w_dir_path$sample_name" \$TMPDIR"
    line_SIM_2="# do the processing"
    line_SIM_3=$line_beginning_SIM_2$sample_name$line_ending_SIM_2$chunk_size;
    line_SIM_4="#delete our copy from /tmp"
    line_SIM_5="rm -rfv \$TMPDIR"
    cat $TCRDivER_dir""templates/job_template_SIM.txt > $job_SIM_filename;
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

while [[ $(bstat | awk '$4 ~ /divP/ { print }') ]]; do wait ; done
echo "DIVERSITY CALCULATIONS DONE"

## Assemble diversity profiles
######################################
#
for D in "$w_dir_path"*;
do if [ -d "${D}" ]; 
    # concatenate and delete chunk dmat 
    then  cd "$D"; 
    ## delete Error/Output/jobchunk files 
    rm "${D}/"Error* "${D}/"Output* 
    sample_name=${D##*/};
    filename="$sample_name""_divP.tsv";
    job_divP_filename="job_divP_"$sample_name".sh";
    file_line="python3 "$TCRDivER_dir"assemble_divP.py "$w_dir_path" "$sample_name;
    cat $TCRDivER_dir""templates/job_template_divP.txt > $job_divP_filename;
    echo $file_line >> $job_divP_filename;
    bsub < $job_divP_filename;
fi;
cd ../; 
done
##
#

echo "ASSEMBLE DIVERSITY DONE"
