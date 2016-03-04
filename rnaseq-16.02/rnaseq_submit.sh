#!/bin/bash


#==================================================================================================
# Created on: 2016-02-22
# Usage: ./rnaseq_submit.sh
# Author: javier.quilez@crg.eu
# Goal: Pipeline for analysis of RNA-seq data
#==================================================================================================



#==================================================================================================
# CONFIGURATION VARIABLES AND PATHS
#==================================================================================================

# Variables and paths
pipeline_name=rnaseq
pipeline_version=16.02
PIPELINE=/users/GR/mb/jquilez/pipelines/$pipeline_name-$pipeline_version
config=$1
pipeline_file=$PIPELINE/$pipeline_name.sh



#==================================================================================================
# CODE EXECUTION
#==================================================================================================

# Get variables from configuration file
if ! [[ -n "$config" ]]; then
	"configuration file with analysis parameters does not exist at $config ! Exiting..."
	exit 
else
	samples=`cat $config | grep '=' | grep -v 'control_bam' | grep samples | sed 's/[\t]//g' | sed 's/;.*//g' | cut -f2 -d"="`
	pipeline_run_mode=`cat $config | grep pipeline_run_mode | sed 's/[ \t]//g' | cut -f2 -d"="`
	io_mode=`cat $config | grep "^io_mode" | sed 's/[ \t]//g' | sed 's/;.*//g' | cut -f2 -d"="`
	CUSTOM_OUT=`cat $config | grep "CUSTOM_OUT" | sed 's/[ \t]//g' | sed 's/;.*//g' | cut -f2 -d"="`
	queue=`cat $config | grep queue | sed 's/[ \t]//g' | sed 's/;.*//g' | cut -f2 -d"="`
	memory=`cat $config | grep memory | sed 's/[ \t]//g' | sed 's/;.*//g' | cut -f2 -d"="`
	max_time=`cat $config | grep max_time | sed 's/[ \t]//g' | sed 's/;.*//g' | cut -f2 -d"="`
	slots=`cat $config | grep slots | sed 's/[ \t]//g' | sed 's/;.*//g' | cut -f2 -d"="`
	submit_to_cluster=`cat $config | grep submit_to_cluster | sed 's/[ \t]//g' | sed 's/;.*//g' | cut -f2 -d"="`
	email=`cat $config | grep email | sed 's/[ \t]//g' | sed 's/;.*//g' | cut -f2 -d"="`	
fi

if [[ $io_mode == "standard" ]]; then
	JOB_CMD=$PIPELINE/job_cmd
	JOB_OUT=$PIPELINE/job_out 
elif [[ $io_mode == "custom" ]]; then
	JOB_CMD=$CUSTOM_OUT/job_cmd
	JOB_OUT=$CUSTOM_OUT/job_out
	mkdir -p $JOB_CMD
	mkdir -p $JOB_OUT
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
	#$ -M $email
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
