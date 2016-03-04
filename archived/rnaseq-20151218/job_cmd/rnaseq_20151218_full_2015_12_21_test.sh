#!/bin/bash
#$ -N rnaseq_20151218_full_2015_12_21_test
#$ -q short-sl65
#$ -l virtual_free=10G
#$ -l h_rt=06:00:00
#$ -o /users/GR/mb/jquilez/pipelines/rnaseq-20151218/job_out/rnaseq_20151218_full_2015_12_21_test_$JOB_ID.out
#$ -e /users/GR/mb/jquilez/pipelines/rnaseq-20151218/job_out/rnaseq_20151218_full_2015_12_21_test_$JOB_ID.err
#$ -j y
#$ -V
#$ -M javier.quilez@crg.eu
#$ -m abe
#$ -pe smp 1

submitted_on=2015_12_21
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
n_bootstraps=5
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
# Trim reads
PROCESSED=$SAMPLE/fastqs_processed/trimmomatic
PAIRED=$PROCESSED/paired_reads
UNPAIRED=$PROCESSED/unpaired_reads
ADAPTERS=/software/mb/el6.3/Trimmomatic-0.33/adapters
# Mapping/alignment
GENOME_DIR=$HOME/assemblies/$species/$version/star_genome_index/read_length_${read_length}bp
STAR=$SAMPLE/alignments/star
KALLISTO_ALIGN=$SAMPLE/alignments/kallisto
# reads per transcript quantification
KALLISTO_QUANT=$SAMPLE/quantifications/kallisto
# Read profiles
PROFILES=$SAMPLE/profiles


# Files 
ifq1=$HOME/data/$data_type/raw/*/${sample_id}_read1.fastq.gz
ifq2=$HOME/data/$data_type/raw/*/${sample_id}_read2.fastq.gz
trimmomatic=`which trimmomatic`
star=`which STAR`
# + as of 2015-11-03, the `which samtools` version of samtools fails to sort the BAM of some samples
# Quique mentioned he had a similar problem and recommended me to use the '/software/mb/el6.3/samtools-1.2/samtools'
# version
samtools=/software/mb/el6.3/samtools-1.2/samtools
bamToBed=`which bamToBed`
bed2bb=`which bedToBigBed`
chrom_sizes=$HOME/assemblies/$species/$version/ucsc/hg19.chrom.sizes.autosomes.chrX
kallisto_index=$HOME/assemblies/$species/$version/kallisto_index/kallisto_homo_sapiens_hg19_mrna.index
kallisto=`which kallisto`



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
	message_info $step "alignments are in $STAR"

	message_time_step $step $time0

}


# =================================================================================================
# Pseudo alignment and transcript quantification with kallisto
# =================================================================================================

quantification_kallisto() {

	step="quantification_kallisto"
	time0=$(date +"%s")

	# perform pseudoalignment and quantify abundances of transcripts
	# --bias = Perform sequence based bias correction
	# --pseudobam = Output pseudoalignments in SAM format to stdout, which is converted to BAM thorugh piping to samtools
	# --plaintext = Output plaintext instead of HDF5 for the quantifications
	message_info $step "performing pseudoalignment and quantifying abundances of transcripts"
	mkdir -p $KALLISTO_ALIGN
	mkdir -p $KALLISTO_QUANT/bootstraps
	paired1=$PAIRED/${sample_id}_read1.fastq.gz
	paired2=$PAIRED/${sample_id}_read2.fastq.gz
	obam=$KALLISTO_ALIGN/${sample_id}_pseudoalignment.bam
	log_quantification_kallisto=$LOGS/${sample_id}_${step}.log
	$kallisto quant -i $kallisto_index \
					-o $KALLISTO_QUANT \
					--bias \
					--pseudobam \
					-t $slots \
					-b $n_bootstraps \
					$paired1 $paired2 | $samtools view -Sb - > $obam
	$kallisto h5dump -o $KALLISTO_QUANT $KALLISTO_QUANT/abundance.h5
	mv $KALLISTO_QUANT/bs_abundance* $KALLISTO_QUANT/bootstraps

	message_info $step "pseudoalignments are in $KALLISTO_ALIGN"
	message_info $step "quantifications are in $KALLISTO_QUANT"

	message_time_step $step $time0

}


# =================================================================================================
# Make RNA-seq read profiles
# =================================================================================================

make_profiles() {

	step="make_profiles"
	time0=$(date +"%s")

	mkdir -p $PROFILES
	log_make_profiles=$LOGS/${sample_id}_${step}_.log

	if [ -d $STAR ]; then
		
		# make read profiles from STAR alignments		
		# filter valid pairs (for paired-end data only) and convert BAM to BED
		# valid pairs are retained with samtools's -f 0x2 command
		message_info $step "make read profiles from STAR alignments"
		message_info $step "Filtering valid pairs (for paired-end data only) and convert BAM to BED"
		ibam=$STAR/${sample_id}.Aligned.sortedByCoord.out.bam
		obed=$PROFILES/${sample_id}.bed
		if [[ $sequencing_type == "PE" ]]; then
			$samtools view -bf 0x2 $ibam | $bamToBed -i stdin > $obed
		else
			$bamToBed -i $ibam > $obed
		fi

		# convert BED to bigWig (more suitable for UCSC Genome Browser uploads)
		message_info $step "converting BED to bigWig (more suitable for UCSC Genome Browser uploads)"
		obb=$PROFILES/${sample_id}.bb
		$bed2bb $obed $chrom_sizes $obb >$log_make_profiles 2>&1

		# delete intermediate BED
		message_info $step "delete intermediate BED"
		rm $obed

	fi

	message_time_step $step $time0

}




# =================================================================================================
# MAIN CODE
# =================================================================================================

echo
if [[ $pipeline_run_mode == 'full' ]]; then
	#trim_reads
	#map_star
	#quantification_kallisto
	make_profiles
elif [[ $pipeline_run_mode == 'trim_reads' ]]; then trim_reads
elif [[ $pipeline_run_mode == 'map_star' ]]; then map_star
elif [[ $pipeline_run_mode == 'quantification_kallisto' ]]; then quantification_kallisto
elif [[ $pipeline_run_mode == 'make_profiles' ]]; then make_profiles
#elif [[ $pipeline_run_mode == 'filter_alignments' ]]; then filter_alignments
#elif [[ $pipeline_run_mode == 'make_profiles' ]]; then make_profiles

fi
echo
