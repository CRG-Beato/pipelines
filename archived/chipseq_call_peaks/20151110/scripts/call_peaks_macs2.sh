#!/bin/bash


#==================================================================================================
# Created on: 2015-11-11
# Usage: ./call_peaks_macs2.sh
# Author: Javier Quilez (GitHub: jaquol)
# Goal: identify ChIP-seq peaks using MACS2
#==================================================================================================


#=================== (REVIEW before running script!) ==============================================

# Samples
samples="T_0_ER_MCF_7_11730_CGATGT \
			T_0_PR_MCF_7_11728_ACAGTG \
			T_30_ER_R5020_11731_CAGATC \
			T_30_PR_R5020_11729_GTGAAA"
samples="T_30_ER_R5020_11731_CAGATC"

# Variable parameters
pipeline_name=chipseq_call_peaks
pipeline_version=20151110
process=peak_calling
mapper_name=bowtie
peak_caller_name=macs2
verbose_level=3

# Paths
PIPELINE=$HOME/pipelines/$pipeline_name/$pipeline_version
JOB_CMD=$PIPELINE/job_cmd 
JOB_OUT=$PIPELINE/job_out
peak_caller=`which $peak_caller_name`

# CRG cluster parameters
queue=short-sl65
memory=40G
max_time=06:00:00
slots=1

#==================================================================================================

for s in $samples; do 

	# Sample directory
	SAMPLE=$HOME/data/chipseq/samples/$s
	BAMS=$SAMPLE/alignments/$mapper_name
	PEAKS=$SAMPLE/peaks/$peak_caller_name
	mkdir -p $PEAKS

	# Build job: parameters
	job_name=${pipeline_name}_${pipeline_version}_${process}_${s}
	job_file=$JOB_CMD/$job_name.sh
	m_out=$JOB_OUT
	echo "#!/bin/bash
	#$ -N $job_name
	#$ -q $queue
	#$ -l virtual_free=$memory
	#$ -l h_rt=$max_time
	#$ -o $m_out/${job_name}_\$JOB_ID.out
	#$ -e $m_out/${job_name}_\$JOB_ID.err
	#$ -M javier.quilez@crg.eu
	#$ -m ae
	#$ -pe smp $slots" > $job_file
	sed -i 's/^\t//g' $job_file

	# Build job: commands
	ibam=$BAMS/${s}_sorted.bam 
	job_cmd="$peak_caller callpeak \
						-t $ibam \
						-n $s \
						--outdir $PEAKS \
						-g hs \
						-B \
						--call-summits \
						--verbose $verbose_level"
	echo $job_cmd >> $job_file

	# Submit job
	chmod a+x $job_file
	qsub < $job_file

done	