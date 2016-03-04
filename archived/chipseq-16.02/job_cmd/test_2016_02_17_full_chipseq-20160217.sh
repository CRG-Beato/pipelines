#!/bin/bash
#$ -N test_2016_02_17_full_chipseq-20160217
#$ -q long-sl65
#$ -l virtual_free=40G
#$ -l h_rt=10:00:00
#$ -o /users/GR/mb/jquilez/pipelines/chipseq-20160217/job_out/test_2016_02_17_full_chipseq-20160217_$JOB_ID.out
#$ -e /users/GR/mb/jquilez/pipelines/chipseq-20160217/job_out/test_2016_02_17_full_chipseq-20160217_$JOB_ID.err
#$ -j y
#$ -M javier.quilez@crg.eu
#$ -m abe
#$ -pe smp 8

submitted_on=2016_02_17
pipeline_version=20160217
sample_id=test
data_type=chipseq
pipeline_run_mode=full
io_mode=standard
CUSTOM_IN=/users/GR/mb/jquilez/pipelines/rnaseq-20160119/test
read1_fname=test_read1.fastq.gz
read2_fname=test_read2.fastq.gz
CUSTOM_OUT=/users/GR/mb/jquilez/misc
submit_to_cluster=no
queue=long-sl65
memory=40G
max_time=10:00:00
slots=8
email=javier.quilez@crg.eu
sequencing_type=SE
seedMismatches=2
palindromeClipThreshold=30
simpleClipThreshold=12
leading=3
trailing=3
minAdapterLength=1
keepBothReads=true
minQual=3
targetLength=40
strictness=0.999
minLength=36
PIPELINE=/users/GR/mb/jquilez/pipelines/chipseq-20160217
config=pipelines/chipseq-20160217/chipseq.config
time_start=$(date +"%s")




#==================================================================================================
# PATHS
#==================================================================================================


# (1) Directories

# pipeline scripts
SCRIPTS=$PIPELINE/scripts

# Primary output directory
if [[ $io_mode == "custom" ]]; then
	SAMPLE=$CUSTOM_OUT/$sample_id
else
	SAMPLE=/users/GR/mb/jquilez/data/$data_type/samples/$sample_id
fi

# Logs
LOGS=$SAMPLE/logs/$version

# Trim reads
PROCESSED=$SAMPLE/fastqs_processed/trimmomatic
SE_READS=$PROCESSED/single_end_reads
ADAPTERS=/software/mb/el6.3/Trimmomatic-0.33/adapters


# (2) Directories

# input FASTQ
if [[ $io_mode == "custom" ]]; then
	ifq1=$CUSTOM_IN/$read1_fname
else
	ifq1=/users/GR/mb/jquilez/data/$data_type/raw/*/${sample_id}*read1.fastq.gz
fi


# =================================================================================================
# CODE EXECUTION
# =================================================================================================

main() {

	echo
	if [[ $pipeline_run_mode == 'full' ]]; then
		trim_reads_single_end
		#map_star
		#make_profiles
		#quantification_htseq
		#quantification_kallisto
		#quantification_featurecounts
	fi
	echo

}



#==================================================================================================
# FUNCTIONS DEFINITIONS
#==================================================================================================


# =================================================================================================
# Pipeline progress functions
# =================================================================================================

# Outputs a message about the task being done
message_info() {
	step_name=$1
	message=$2
	echo -e "INFO \t`date +"%Y-%m-%d %T"` \t[$step_name] \t$message"
}

# Outputs a message about the error found and exits
message_error() {
	step_name=$1
	message=$2
	echo -e "ERROR \t`date +"%Y-%m-%d %T"` \t[$step_name] \t$message"
	exit	
}

# Outputs a warning message about the task being done
message_warn() {
	step_name=$1
	message=$2
	echo -e "WARN \t`date +"%Y-%m-%d %T"` \t[$step_name] \t$message"
}

# Outputs a message with the time in seconds spent in a given step
message_time_step() {
	step_name=$1
	field_name="TIME_${step_name^^}"
	time0=$2
	time1=$(date +"%s")
	length=$(($time1-$time0))
	echo -e "TIME \t`date +"%Y-%m-%d %T"` \t[$step_name] \tstep time for completion (seconds) = $length"
	echo
}

# Outputs the total time in seconds for the pipeline to run
message_time_pipeline() {
	field_name="TIME_COMPLETE_PIPELINE"
	time0=$time_start
	time1=$(date +"%s")
	length=$(($time1-$time0))
	echo -e "TIME \t`date +"%Y-%m-%d %T"` \t[pipeline] \ttotal time for completion (seconds) = $length"
	echo
}


# =================================================================================================
# Trim adapter and low-quality ends (single-end data)
# =================================================================================================

trim_reads_single_end() {

	step="trim_reads"
	time0=$(date +"%s")

	# Check that FASTQ file exists
	if [ -f $ifq1 ]; then
		mkdir -p $SAMPLE
		mkdir -p $SE_READS
		mkdir -p $LOGS
		log_trim_reads=$SAMPLE/logs/${sample_id}_${step}.log
	else
		message_error $step "$ifq1 not found. Exiting..."
	fi

	# adapter trimming: the trimmomatic program directory contains a folder with the adapter sequences for
	# the Illumina sequencers in use. 'TruSeq3-PE.fa' is used, which contains the adapter sequences for the HiSeq 
	message_info $step "trimming adapter sequences for HiSeq, NextSeq or HiSeq"
	message_info $step "trimming low-quality reads ends using trimmomatic's recommended practices"
	seqs=$ADAPTERS/TruSeq3-$sequencing_type.fa
	trimmed1=$PAIRED/${sample_id}_read1.fastq.gz
	$trimmomatic $sequencing_type \
 					$ifq1 \
 					$trimmed1 \
 					ILLUMINACLIP:$seqs:$seedMismatches:$palindromeClipThreshold:$simpleClipThreshold:$minAdapterLength:$keepBothReads \
 					LEADING:$leading \
 					TRAILING:$trailing \
 					MAXINFO:$targetLength:$strictness \
 					MINLEN:$minLength >$log_trim_reads 2>&1

	message_time_step $step $time0

}




main


