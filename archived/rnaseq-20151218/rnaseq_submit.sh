#!/bin/bash


#==================================================================================================
# Created on: 2015-12-18
# Usage: ./rnase_submit.sh
# Author: javier.quilez@crg.eu
# Goal: Pipeline for analysis of RNA-seq data
#==================================================================================================



#==================================================================================================
# CONFIGURATION VARIABLES AND PATHS
#==================================================================================================

# Variables and paths
pipeline_name=rnaseq
pipeline_version=20151218
PIPELINE=$HOME/pipelines/$pipeline_name-$pipeline_version
JOB_CMD=$PIPELINE/job_cmd
JOB_OUT=$PIPELINE/job_out 
#config=$PIPELINE/$pipeline_name.config
config=$1
pipeline_file=$PIPELINE/$pipeline_name.sh



#==================================================================================================
# CONFIGURATION VARIABLES AND PATHS
#==================================================================================================

# Get variables from configuration file
if ! [[ -n "$config" ]]; then
	"configuration file with analysis parameters does not exist at $config ! Exiting..."
	exit 
else
	samples=`cat $config | grep '=' | grep samples | sed 's/[\t]//g' | sed 's/;.*//g' | cut -f2 -d"="`
	pipeline_run_mode=`cat $config | grep pipeline_run_mode | sed 's/[ \t]//g' | cut -f2 -d"="`
	queue=`cat $config | grep queue | sed 's/[ \t]//g' | sed 's/;.*//g' | cut -f2 -d"="`
	memory=`cat $config | grep memory | sed 's/[ \t]//g' | sed 's/;.*//g' | cut -f2 -d"="`
	max_time=`cat $config | grep max_time | sed 's/[ \t]//g' | sed 's/;.*//g' | cut -f2 -d"="`
	slots=`cat $config | grep slots | sed 's/[ \t]//g' | sed 's/;.*//g' | cut -f2 -d"="`
	submit_to_cluster=`cat $config | grep submit_to_cluster | sed 's/[ \t]//g' | sed 's/;.*//g' | cut -f2 -d"="`
fi

# Run pipeline for each sample
for s in $samples; do

	# Build job: parameters
	submitted_on=`date +"%Y_%m_%d"`
	job_name=${s}_${submitted_on}_${pipeline_run_mode}_${pipeline_name}-${pipeline_version}
	job_file=$JOB_CMD/$job_name.sh
	m_out=$JOB_OUT
	echo "#!/bin/bash
	#$ -N $job_name
	#$ -q $queue
	#$ -l virtual_free=$memory
	#$ -l h_rt=$max_time
	#$ -o $m_out/${job_name}_\$JOB_ID.out
	#$ -e $m_out/${job_name}_\$JOB_ID.err
	#$ -j y
	#$ -M javier.quilez@crg.eu
	#$ -m abe
	#$ -pe smp $slots" > $job_file
	sed -i 's/^\t//g' $job_file

	# Add date of submission
	echo -e "\nsubmitted_on=$submitted_on" >> $job_file
	# Add pipeline version 
	echo "pipeline_version=$pipeline_version" >> $job_file
	# Add sample ID 
	echo "sample_id=$s" >> $job_file
	# Add parameters from the configuration file
	cat $config | grep '=' | grep -v samples | sed 's/[ \t]//g' | sed 's/;.*//g' | sed '/^$/d' >> $job_file
	# Add path to pipeline files and configuration file
	echo "PIPELINE=$PIPELINE" >> $job_file
	echo "config=$config" >> $job_file

	# Add content of the pipeline script to the submission script
	cat $pipeline_file >> $job_file

	# Submit
	chmod a+x $job_file
	if [[ $submit_to_cluster == "yes" ]]; then
		qsub < $job_file
	else
		chmod a+x $job_file
		$job_file
	fi

done
