#!/bin/bash


#==================================================================================================
# Created on: 2015-11-03
# Usage: ./mapping_bowtie.sh
# Author: Javier Quilez (GitHub: jaquol)
# Goal: map ChIP-seq data to the reference genome with the bowtie mapped
#==================================================================================================

#==================================================================================================
# Notes:
# + as of 2015-11-03, the `which samtools` version of samtools fails to sort the BAM of some samples
# Quique mentioned he had a similar problem and recommended me to use the '/software/mb/el6.3/samtools-1.2/samtools'
# version
#==================================================================================================



# =================== (REVIEW before running script!) ====================================

# Samples
# samples="T_0_ER_MCF_7_11730_CGATGT \
# 			T_0_PR_MCF_7_11728_ACAGTG \
# 			T_30_ER_R5020_11731_CAGATC \
# 			T_30_PR_R5020_11729_GTGAAA"
samples="T_30_ER_R5020_11731_CAGATC"

# Variable parameters
pipeline_name=chipseq_call_peaks
pipeline_version=20151110
process=mapping
mapper=bowtie
max_num_alignments=1

# Paths
PIPELINE=$HOME/pipelines/$pipeline_name/$pipeline_version
JOB_CMD=$PIPELINE/job_cmd 
JOB_OUT=$PIPELINE/job_out
FASTQ=$HOME/data/chipseq/raw/2015-11-03
pigz=`which pigz`
bowtie=`which $mapper`
samtools=/software/mb/el6.3/samtools-1.2/samtools
index=/users/GR/mb/dsoronellas/tracks/human_genome_19/bowtie_index_t47d/hg19-noMY-noMMTV

# CRG cluster parameters
queue=short-sl65
memory=40G
max_time=06:00:00
slots=8

# ========================================================================================

for s in $samples; do

	# Sample directory
	SAMPLE=$HOME/data/chipseq/samples/$s
	BAMS=$SAMPLE/alignments/$mapper
	mkdir -p $SAMPLE
	mkdir -p $BAMS

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
	ifq=$FASTQ/$s*fastq.gz
	obam=$BAMS/$s.bam
	obam_sorted=$BAMS/${s}_sorted.bam
	job_cmd="$pigz -dc -p $slots $ifq | $bowtie -p $slots -m $max_num_alignments -S $index - | $samtools view -@ $slots -Sb - > $obam"
	echo $job_cmd >> $job_file
	job_cmd="$samtools sort -o $obam_sorted -T $BAMS/tmp $obam"
	echo $job_cmd >> $job_file

	# Submit job
	chmod a+x $job_file
	#qsub < $job_file

done