#!/bin/bash
#$ -N rnaseq_20151218_full_2015_12_18_test
#$ -q short-sl65
#$ -l virtual_free=10G
#$ -l h_rt=06:00:00
#$ -o /users/GR/mb/jquilez/pipelines/rnaseq-20151218/job_out/rnaseq_20151218_full_2015_12_18_test_$JOB_ID.out
#$ -e /users/GR/mb/jquilez/pipelines/rnaseq-20151218/job_out/rnaseq_20151218_full_2015_12_18_test_$JOB_ID.err
#$ -j y
#$ -V
#$ -M javier.quilez@crg.eu
#$ -m abe
#$ -pe smp 1

submitted_on=2015_12_18
pipeline_version=20151218
sample_id=test
data_type=rnaseq
pipeline_run_mode=full
submit_to_cluster=no
queue=short-sl65
memory=10G
max_time=06:00:00
slots=1
sequencing_type=PE
trimmomatic_seedMismatches=2
trimmomatic_palindromeClipThreshold=30
trimmomatic_simpleClipThreshold=12
trimmomatic_minAdapterLength=1
trimmomatic_keepBothReads=true
trimmomatic_minQual=3
trimmomatic_targetLength=40
trimmomatic_strictness=0.5
trimmomatic_minLength=36
species=homo_sapiens
version=hg19
read_length=50
PIPELINE=/users/GR/mb/jquilez/pipelines/rnaseq-20151218
config=/users/GR/mb/jquilez/pipelines/rnaseq-20151218/rnaseq.config
time_start=$(date +"%s")




#==================================================================================================
# PATHS
#==================================================================================================

# Directories
DATA=$HOME/data
SCRIPTS=$PIPELINE/scripts
SAMPLE=$HOME/data/$data_type/samples/$sample_id
LOGS=$SAMPLE/logs
PROCESSED=$SAMPLE/fastqs_processed/trimmomatic
PAIRED=$PROCESSED/paired_reads
UNPAIRED=$PROCESSED/unpaired_reads
ADAPTERS=/software/mb/el6.3/Trimmomatic-0.33/adapters
GENOME_DIR=$HOME/assemblies/$species/$version/star_genome_index/read_length_${read_length}bp
STAR=$SAMPLE/alignments/star

# Files 
ifq1=$HOME/data/$data_type/raw/*/${sample_id}_read1.fastq.gz
ifq2=$HOME/data/$data_type/raw/*/${sample_id}_read2.fastq.gz
trimmomatic=`which trimmomatic`
star=`which STAR`


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
# Trim adapter and low-quality ends
# =================================================================================================

trim_reads() {

	step="trim_reads"
	time0=$(date +"%s")

	# Check that FASTQ file exists
	if [ -f $ifq1 ] && [ -f $ifq2 ]; then
		mkdir -p $SAMPLE
		mkdir -p $PAIRED
		mkdir -p $UNPAIRED
		mkdir -p $LOGS
		log_trim_reads=$LOGS/${sample_id}_${step}.log
	else
		message_error $step "$ifq1 and/or $ifq2 not found. Exiting..."
	fi

	# adapter trimming: the trimmomatic program directory contains a folder with the adapter sequences for
	# the Illumina sequencers in use. 'TruSeq3-PE.fa' is used, which contains the adapter sequences for the HiSeq 
	message_info $step "trimming adapter sequences for HiSeq, NextSeq or HiSeq"
	message_info $step "trimming low-quality reads ends using trimmomatic's recommended practices"
	seqs=$ADAPTERS/TruSeq3-$sequencing_type.fa
	paired1=$PAIRED/${sample_id}_read1.fastq.gz
	paired2=$PAIRED/${sample_id}_read2.fastq.gz
	unpaired1=$UNPAIRED/${sample_id}_read1.fastq.gz
	unpaired2=$UNPAIRED/${sample_id}_read2.fastq.gz
	$trimmomatic $sequencing_type \
 					$ifq1 $ifq2 \
 					$paired1 $unpaired1 \
 					$paired2 $unpaired2 \
 					ILLUMINACLIP:$seqs:$trimmomatic_seedMismatches:$trimmomatic_palindromeClipThreshold:$trimmomatic_simpleClipThreshold:$trimmomatic_minAdapterLength:$trimmomatic_keepBothReads \
 					MINLEN:$trimmomatic_minLength >$log_trim_reads 2>&1
	message_info $step "unpaired reads are deleted"
	rm -fr $UNPAIRED

	message_time_step $step $time0

}


# =================================================================================================
# Map read with STAR
# =================================================================================================

map_star() {

	step="map_star"
	time0=$(date +"%s")

	# map paired-end reads with STAR
	message_info $step "mapping trimmed paired-end reads with STAR"
	message_info $step "using ENCODE standard options for long RNA-seq pipeline"
	# Mapping options adjusted following ENCODE standard options for long RNA-seq found in the 
	# STAR manual https://github.com/alexdobin/STAR/tree/master/doc (note that although this links to the STAR 2.4 manual,
	# the options used here are also available for STAR 2.3, the version used here)
	# --outFilterType BySJout = reduces the number of ”spurious” junctions
	# --outFilterMultimapNmax 20 = max number of multiple alignments allowed for a read: if exceeded, the read is considered unmapped
	# --alignSJoverhangMin 8 = minimum overhang for unannotated junctions
	# --alignSJDBoverhangMin 1 = minimum overhang for annotated junctions
	# --outFilterMismatchNmax 999 = maximum number of mismatches per pair, large number switches off this filter
	# --outFilterMismatchNoverLmax 0.04 = max number of mismatches per pair relative to read length: for 2x50b, 
	# max number of mis- matches is 0.04*100=4 for the pair (default is 0.3 so we are much more restrictive)
	# --alignIntronMin 20: ...
	# --alignIntronMax 1000000: maximum intron length
	# --outSAMtype BAM SortedByCoordinate = output sorted by coordinate
	paired1=$PAIRED/${sample_id}_read1.fastq.gz
	paired2=$PAIRED/${sample_id}_read2.fastq.gz
	log_map_star=$LOGS/${sample_id}_${step}.log
	mkdir -p $STAR
	TMP_DIR=$STAR/tmp
	$star \
		    --genomeDir $GENOME_DIR/ \
			--genomeLoad NoSharedMemory \
			--runThreadN $slots \
			--outFilterType "BySJout" \
			--outFilterMultimapNmax 20 \
			--alignSJoverhangMin 8 \
			--alignSJDBoverhangMin 1 \
			--outFilterMismatchNmax 999 \
			--outFilterMismatchNoverLmax 0.04 \
			--alignIntronMin 20 \
			--alignIntronMax 1000000 \
			--alignMatesGapMax 1000000 \
			--readFilesIn $paired1 $paired2 \
			--outSAMtype BAM SortedByCoordinate \
			--outTmpDir $TMP_DIR \
			--outFileNamePrefix $STAR/$sample_id. \
			--readFilesCommand zcat >$log_map_star 2>&1
	rm -fr $STAR/tmp*

	message_time_step $step $time0

}


# =================================================================================================
# MAIN CODE
# =================================================================================================

echo
if [[ $pipeline_run_mode == 'full' ]]; then
	trim_reads
	map_star
#	filter_alignments
#	make_profiles
elif [[ $pipeline_run_mode == 'trim_reads' ]]; then trim_reads
elif [[ $pipeline_run_mode == 'map_star' ]]; then map_star
#elif [[ $pipeline_run_mode == 'mapping_tophat' ]]; then mapping_tophat
#elif [[ $pipeline_run_mode == 'filter_alignments' ]]; then filter_alignments
#elif [[ $pipeline_run_mode == 'make_profiles' ]]; then make_profiles

fi
echo
