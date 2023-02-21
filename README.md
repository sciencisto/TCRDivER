#########################################################################
# TCRDivER algorithm
#########################################################################

Thank you for the interest and using TCRDivER algorithm! 

We are constantly improving the code and if you wish to contribute or
have any questions feel free to contact. 

Milena Vujovic (milvu@dtu.dk; sciencistom@gmail.com)
Joseph Kaplinsky (joseph.kaplinsky@ludwig.ox.ac.uk)

The software runs in a Python3.6 environment. We have tested in on UNIX 
based systems. If you have any issues running the software, again, 
feel free to contact us- we are more than happy to help!

The algorithm is split up into several parts. Each script is designed 
to perform a step of the algorithm. These scripts individual scripts 
are combined in an overall bash script named TCRDiVER_run.sh. 

All the bash scripts contain a preamble that enable them to run 
on our university cluster. Please keep in mind that you might need to 
change the preamble if the cluster you are accessing has different 
specifications or queing system. The queing system command needs to be replaced 
in the TCRDivER_run.sh script. The default is "bsub <". Another example could be "qsub ". 


If you are good in bash you will be able to augment the bash scripts 
to suit your needs. If not, shoot us a line --> we might be able to help :) 

We envision our tool to be run on cluster systems since the calculations 
can be computationaly expensive. If you wish to run them locally, we advise you change 
the TCRDivER_run.sh script to suit your local needs or run the scripts individually. 

For individual usage: 

python3.6 [python_script.py] [arguments]


The input can be one of three: 

1. Each sample is a .txt file with a list of CDR3s. These CDR3s are all the possible CDR3s and we expect them to be non-unique.

2. Each sample is in the ImmunoSeq Adaptive Biotechnologies sequencing .tsv file version1:  The columns needed to run TCRDivER are: aminoAcid, count (templates/reads), sequenceStatus.

3. Each sample is in the ImmunoSeq Adaptive Biotechnologies sequencing .tsv file version2:  The columns needed to run TCRDivER are: amino_acid, productive_frequency, frame_type

The frame_type or sequenceStatus columns in input options 2 and 3 are Adaptive Biotechnologies sequening marking that designated that the sequence read was "In" frame and not "Out" of frame or containing a "Stop" codon. 

Considerations around using input templates: 
	
* Templates of these inputs are given in the folder Input_templates. If populating the columns manually, please keep in mind that in option 3, the sum of the individual productive_frequencies of the CDR3s should be equal to 1. TCRDivER by default works only with "In" frame sequences. Therefore if you are populating the input template 2 and 3 please remember to populate the column frame_type or sequenceStatus with the string "In". 
	
* If choosing option 1 as input - TCRDivER.config file count variable has to be changed to: count="count"

* Please note that the input templates contain only columns needed to run TCRDivER. If you obtain the files from the ImmunoSeq Adaptive Biotechnologies platform there will be more columns in the files. 




# Algorithm outline


#########################################################################
## I. Filter sample and make distance matrix chunk jobs: 
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



#########################################################################
## II. Run all distance matrix calculations. 
#########################################################################


#########################################################################
## III. Make and run true diversity calculations. 
#########################################################################

Now that you've calculated the distance matrix, you can proceed with calculating the tru divesity with 
varying values of lambda and q. A bash script takes in the whole working directory and makes naive and 
similarity scaled diversity jobs. The outputs of these two scripts are two .tsv files containing a 
dataframe with the values of naive diversity with varying q and similarity scaled diversity with varying
q and lambda values denoted [sample_name]_ID.tsv and [sample_name]_SIM.tsv, respectively. 


#########################################################################
## IV. Assemble complete diversity calculation 
#########################################################################

One the naive and similarity scaled diversity calculations are completed each sample folder will contain two files: 
[sample_name]_ID.tsv and [sample_name]_SIM.tsv calculations. These are now assembled into an overall file of diversity 
using the assemble_divP_cluster.sh bash script. 

It runs the python script assemble_divP.py which takes the two files [sample_name]_ID.tsv and [sample_name]_SIM.tsv 
as input and outputs a complete diversity dataframe named [sample_name]_divP.tsv. 





























