#!/bin/bash
#$ -N netseq_20151125_mapping_tophat_2015_12_02_NETseq.T0-2
#$ -q long-sl65
#$ -l virtual_free=50G
#$ -l h_rt=48:00:00
#$ -o /users/GR/mb/jquilez/pipelines/netseq-20151125/job_out/netseq_20151125_mapping_tophat_2015_12_02_NETseq.T0-2_$JOB_ID.out
#$ -e /users/GR/mb/jquilez/pipelines/netseq-20151125/job_out/netseq_20151125_mapping_tophat_2015_12_02_NETseq.T0-2_$JOB_ID.err
#$ -j y
#$ -V
#$ -M javier.quilez@crg.eu
#$ -m abe
#$ -pe smp 8

submitted_on=2015_12_02
pipeline_version=20151125
sample_id=NETseq.T0-2
pipeline_run_mode=mapping_tophat
submit_to_cluster=yes
queue=long-sl65
memory=50G
max_time=48:00:00
slots=8
library_type=fr-firststrand
species=homo_sapiens
version=hg19
read_length=50
PIPELINE=/users/GR/mb/jquilez/pipelines/netseq-20151125
config=/users/GR/mb/jquilez/pipelines/netseq-20151125/netseq.config
time_start=$(date +"%s")




#==================================================================================================
# PATHS
#==================================================================================================

# Directories
DATA=$HOME/data
SCRIPTS=$PIPELINE/scripts
SAMPLE=$HOME/data/netseq/samples/$sample_id
LOGS=$SAMPLE/logs
PROCESSED=$SAMPLE/fastqs_processed/cutadapt
TOPHAT=$SAMPLE/alignments/tophat
STAR=$SAMPLE/alignments/star

# Files 
fq1=$HOME/data/netseq/raw/*/${sample_id}.read1.fastq.gz
fq2=$HOME/data/netseq/raw/*/${sample_id}.read2.fastq.gz
python=`which python`
cutadapt=`which cutadapt`
tophat=`which tophat`
# + as of 2015-11-03, the `which samtools` version of samtools fails to sort the BAM of some samples
# Quique mentioned he had a similar problem and recommended me to use the '/software/mb/el6.3/samtools-1.2/samtools'
# version
samtools=/software/mb/el6.3/samtools-1.2/samtools
index=/users/GR/mb/dsoronellas/tracks/human_genome_19/bowtie2_index_t47d/hg19-noMY-noMMTV
transcriptomeIndex=/users/GR/mb/dsoronellas/tracks/human_genome_19/gencode_annotation/transcriptome_index_tophat2/
gtf=/users/GR/mb/dsoronellas/tracks/human_genome_19/gencode_annotation/transcriptome_index_tophat2/gencode.v18.annotation-sorted.gtf
tophat=`which tophat`
star=`which STAR`
GENOME_DIR=$HOME/assemblies/$species/$version/star_genome_index/read_length_${read_length}bp

# Specifies which python to use
WORKON_HOME=/software/mb/el6.3/python/envs/.virtualenvs;
source /usr/bin/virtualenvwrapper.sh
export PYTHONZ_ROOT=/software/mb/el6.3/python
[[ -s /software/mb/el6.3/python/etc/bashrc ]] && source /software/mb/el6.3/python/etc/bashrc
workon cpython279

# Required for tophat to find bowtie
export PATH=/software/mb/bin/bowtie2:$PATH



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
# Trim adapter
# =================================================================================================

# From the following directory I understand that a first step is trimming the Truseq Small RNA 
# preparation kit adapter from the paired-end sequences:
# ls /users/GR/mb/dsoronellas/NET-seq/002_TRIMMED_FASTQS

trim_adapter() {

	step="trim_adapter"
	time0=$(date +"%s")

	# Check that FASTQ file exists
	if [ -f $fq1 ] && [ -f $fq2 ]; then
		mkdir -p $SAMPLE
		mkdir -p $PROCESSED
		mkdir -p $LOGS
		log_trim_adapter=$LOGS/${sample_id}_${step}.log
	else
		message_error $step "$fq1 and/or $fq2 not found. Exiting..."
	fi

	# remove Truseq small RNA preparation kit adapter with cutadapt
	# I use parameters as defined by Dani
	# --error-rate = Maximum allowed error rate (no. of errors divided by the length of the matching region)
	# --minimum-length = discard reads that are shorter than a given lenght after trimming
	# --time-n = trim N's on ends of reads
	# --quality-cutoff = Trim low-quality bases from 5' and/or 3' ends of reads
    #                    before adapter removal. If one value is given, only
    #                    the 3' end is trimmed. If two comma-separated cutoffs
    #                    are given, the 5' end is trimmed with the first
    #                    cutoff, the 3' end with the second. The algorithm is
    #                    the same as the one used by BWA (see documentation).
    # -a = sequence of the adapter
    # -A = 3' adapter to be removed from the second read in a pair
    # --output, --paired-output = output files for read1 and read2, respectively
	message_info $step "remove Truseq small RNA preparation kit adapter with cutadapt"
	ofq1=$PROCESSED/${sample_id}_read1.fastq.gz
	ofq2=$PROCESSED/${sample_id}_read2.fastq.gz
	$python $cutadapt \
			--error-rate=0.05 \
			--minimum-length=10 \
			--trim-n \
			--quality-cutoff=10 \
	        -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC \
	        -A GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT \
	        --output=$ofq1 \
	        --paired-output=$ofq2 \
	        $fq1 $fq2 > $log_trim_adapter

	message_time_step $step $time0

}


# =================================================================================================
# Mapping with tophat
# =================================================================================================

# From the following directory I understand that the next step is mapping with TopHat:
# ls /users/GR/mb/dsoronellas/NET-seq/003_MAPPED_SAMPLES

mapping_tophat() {

	step="mapping_tophat"
	time0=$(date +"%s")

	# map paired-end reads with tophat
	# In Dani's script some of the tophat parameters values as passed as script arguments so I 
	# a priori I cannot which values were used. On 2015-11-26 I asked him to specify but until he
	# replies I will use tophat's default values
	# Tophat parameters:
	# --min-anchor = TopHat will report junctions spanned by reads with at least this many bases on each side of the junction
	# --max-multihits = Instructs TopHat to allow up to this many alignments to the reference for a given read
	# Nojima et al use 1 as value so it makes sense Dani used it as well
	# library-type = ...
	# --GTF = gene annotation file. It makes sense this file is within the same directory as $transcriptomeIndex
	# indeed, there is a GTF so I guess it is the one 
	# --no-coverage-search = Disables the coverage based search for junctions
	# --segment-length = Each read is cut up into segments, each at least this long. These segments are mapped independently.
	# I am not sure what Dani is trying to do to set the value of this parameter. Therefore I will used the default 25
	ifq1=$PROCESSED/${sample_id}_read1.fastq.gz
	ifq2=$PROCESSED/${sample_id}_read2.fastq.gz
	log_mapping=$LOGS/${sample_id}_${step}.log
	mkdir -p $TOPHAT
	message_info $step "map paired-end reads with tophat"
	$tophat -p $slots \
	        --output-dir $TOPHAT \
	        --transcriptome-index $transcriptomeIndex \
            --min-anchor 5 \
            --max-multihits 1 \
            --library-type $library_type \
            --GTF $gtf \
            --b2-very-sensitive \
            --no-coverage-search \
            --segment-length 25 \
			$index \
			$ifq1 $ifq2 2>$log_mapping

	message_time_step $step $time0

}


# =================================================================================================
# Mapping
# =================================================================================================

# I found that mapping with tophat (see above) is very slow so I will use STAR instead

mapping() {

	step="mapping"
	time0=$(date +"%s")

	# map paired-end reads with STAR
	message_info $step "map paired-end reads with star"
	# Mapping options adjusted following ENCODE standard options for RNA-seq found in the 
	# STAR manual https://github.com/alexdobin/STAR/tree/master/doc (note that although this links to the STAR 2.4 manual,
	# the options used here are also available for STAR 2.3, the version used here)
	# --outFilterType BySJout = reduces the number of ”spurious” junctions
	# --outFilterMismatchNmax 999 = maximum number of mismatches per pair, large number switches off this filter
	# --outFilterMismatchNoverLmax 0.04 = max number of mismatches per pair relative to read length: for 1x50b, 
	# max number of mis- matches is 0.04*5=2 for the read (default is 0.3 so we are much more restrictive)
	# --alignIntronMax 1000000: maximum intron length (no default)
	# --outSAMtype BAM SortedByCoordinate = output sorted by coordinate
	# Mapping options adjusted following Dani's criteria in the original pipeline using tophat
	# (these prevail in case of conflict with ENCODE's values)
	# --outFilterMultimapNmax 1 = max number of multiple alignments allowed for a read: if exceeded, the read is considered unmapped
	# --alignSJoverhangMin 5 = minimum overhang for unannotated junctions
	# --alignSJDBoverhangMin 5 = minimum overhang for annotated junctions
	ifq1=$PROCESSED/${sample_id}_read1.fastq.gz
	ifq2=$PROCESSED/${sample_id}_read2.fastq.gz
	log_mapping=$LOGS/${sample_id}_${step}.log
	mkdir -p $STAR
	TMP_DIR=$STAR/tmp
	$star \
		    --genomeDir $GENOME_DIR/ \
			--genomeLoad NoSharedMemory \
			--runThreadN $slots \
			--outFilterType "BySJout" \
			--outFilterMultimapNmax 1 \
			--alignSJoverhangMin 8 \
			--alignSJDBoverhangMin 1 \
			--outFilterMismatchNmax 999 \
			--outFilterMismatchNoverLmax 0.04 \
			--alignIntronMin 20 \
			--alignIntronMax 1000000 \
			--readFilesIn $ifq1 $ifq2 \
			--outSAMtype BAM SortedByCoordinate \
			--outTmpDir $TMP_DIR \
			--outFileNamePrefix $STAR/$sample_id. \
			--readFilesCommand zcat 2>$log_mapping
	rm -fr $STAR/tmp*

	message_time_step $step $time0

}


# =================================================================================================
# MAIN CODE
# =================================================================================================

echo
if [[ $pipeline_run_mode == 'full' ]]; then
	trim_adapter
	mapping
elif [[ $pipeline_run_mode == 'trim_adapter' ]]; then trim_adapter
elif [[ $pipeline_run_mode == 'mapping' ]]; then mapping
elif [[ $pipeline_run_mode == 'mapping_tophat' ]]; then mapping_tophat
fi
echo

