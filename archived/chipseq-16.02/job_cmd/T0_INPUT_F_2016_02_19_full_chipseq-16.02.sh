#!/bin/bash
#$ -N T0_INPUT_F_2016_02_19_full_chipseq-16.02
#$ -q long-sl65
#$ -l virtual_free=40G
#$ -l h_rt=20:00:00
#$ -o /users/GR/mb/jquilez/pipelines/chipseq-16.02/job_out/T0_INPUT_F_2016_02_19_full_chipseq-16.02_$JOB_ID.out
#$ -e /users/GR/mb/jquilez/pipelines/chipseq-16.02/job_out/T0_INPUT_F_2016_02_19_full_chipseq-16.02_$JOB_ID.err
#$ -j y
#$ -M javier.quilez@crg.eu
#$ -m abe
#$ -pe smp 8

submitted_on=2016_02_19
pipeline_version=16.02
sample_id=T0_INPUT_F
data_type=chipseq
pipeline_run_mode=full
io_mode=standard
CUSTOM_IN=/users/GR/mb/jquilez/pipelines/chipseq-16.02/test
read1_fname=test_read1.fastq.gz
read2_fname=test_read2.fastq.gz
CUSTOM_OUT=/users/GR/mb/jquilez/misc
submit_to_cluster=yes
queue=long-sl65
memory=40G
max_time=20:00:00
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
species=homo_sapiens
version=hg19
peak_caller=macs2
use_control=no
control_bam=
macs2_qvalue=0.05
PIPELINE=/users/GR/mb/jquilez/pipelines/chipseq-16.02
config=pipelines/chipseq-16.02/chipseq.config
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
# alignment
BWA=$SAMPLE/alignments/bwa/$version
# tag directory
TAG_DIR=$SAMPLE/tag_directory/homer/$version
# reads per million (RPM) profiles
RPMS=$SAMPLE/rpms/$version
# peaks
PEAKS=$SAMPLE/peaks/$peak_caller/$version

# (2) Files

# input FASTQ
if [[ $io_mode == "custom" ]]; then
	ifq1=$CUSTOM_IN/$read1_fname
else
	ifq1=/users/GR/mb/jquilez/data/$data_type/raw/*/${sample_id}*read1.fastq.gz
fi

# tools
trimmomatic=`which trimmomatic`
# I use this newer BWA version instead of that pointed by `which bwa` (0.7.10-r789)
bwa=/software/mb/el6.3/bwa/bwa-0.7.12/bwa
genome_fasta=/users/GR/mb/jquilez/assemblies/$species/$version/ucsc/$version.fa
# + as of 2015-11-03, the `which samtools` version of samtools fails to sort the BAM of some samples
# Quique mentioned he had a similar problem and recommended me to use the '/software/mb/el6.3/samtools-1.2/samtools'
samtools=/software/mb/el6.3/samtools-1.2/samtools
bamToBed=`which bamToBed`
makeTagDirectory=`which makeTagDirectory`
bedToBigBed=`which bedToBigBed`
bedtools=`which bedtools`
perl=`which perl`
bam2wig=`which bam2wig.pl`
bedGraphToBigWig=`which bedGraphToBigWig`

# genome fasta and chromosome sizes
genome_fasta=/users/GR/mb/jquilez/assemblies/$species/$version/ucsc/$version.fa
genome_chrom_sizes=/users/GR/mb/jquilez/assemblies/$species/$version/ucsc/$version.chrom.sizes.autosomes.*



# =================================================================================================
# CODE EXECUTION
# =================================================================================================

main() {

	echo
	if [[ $pipeline_run_mode == 'full' ]]; then
		trim_reads_single_end
		align_single_end
		make_tag_directory
		make_bigbed
		calculate_rpms
		call_peaks
	elif [[ $pipeline_run_mode == 'full_no_call_peaks' ]]; then
		trim_reads_single_end
		align_single_end
		make_tag_directory
		make_bigbed
		calculate_rpms
	elif [[ $pipeline_run_mode == 'trim_reads_single_end' ]]; then trim_reads_single_end
	elif [[ $pipeline_run_mode == 'align_single_end' ]]; then align_single_end
	elif [[ $pipeline_run_mode == 'make_tag_directory' ]]; then make_tag_directory
	elif [[ $pipeline_run_mode == 'make_bigbed' ]]; then make_bigbed
	elif [[ $pipeline_run_mode == 'calculate_rpms' ]]; then calculate_rpms
	elif [[ $pipeline_run_mode == 'call_peaks' ]]; then call_peaks
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

	step="trim_reads_single_end"
	time0=$(date +"%s")

	# Check that FASTQ file exists
	if [ -f $ifq1 ]; then
		mkdir -p $SAMPLE
		mkdir -p $SE_READS
		mkdir -p $LOGS
		step_log=$SAMPLE/logs/${sample_id}_${step}.log
	else
		message_error $step "$ifq1 not found. Exiting..."
	fi

	# adapter trimming: the trimmomatic program directory contains a folder with the adapter sequences for
	# the Illumina sequencers in use. 'TruSeq3-PE.fa' is used, which contains the adapter sequences for the HiSeq 
	message_info $step "trimming adapter sequences for HiSeq, NextSeq or HiSeq"
	message_info $step "trimming low-quality reads ends using trimmomatic's recommended practices"
	seqs=$ADAPTERS/TruSeq3-$sequencing_type.fa
	trimmed1=$SE_READS/${sample_id}_read1.fastq.gz
	$trimmomatic $sequencing_type \
 					$ifq1 \
 					$trimmed1 \
 					ILLUMINACLIP:$seqs:$seedMismatches:$palindromeClipThreshold:$simpleClipThreshold:$minAdapterLength:$keepBothReads \
 					LEADING:$leading \
 					TRAILING:$trailing \
 					MAXINFO:$targetLength:$strictness \
 					MINLEN:$minLength >$step_log 2>&1
 	n_reads_trimmed=`grep Surviving $step_log | cut -f3 -d':' | cut -f1 -d'(' | sed "s/ //g"`
 	echo $n_reads_trimmed > tmp.txt
	message_info $step "reads after trimming = $n_reads_trimmed"

	message_time_step $step $time0

}


# =================================================================================================
# aling single-end reads
# =================================================================================================

align_single_end() {

	step="align_single_end"
	time0=$(date +"%s")

	# align single-end reads with BWA
	message_info $step "align single-end reads with BWA"
	message_info $step "alignments are converted to BAM, sorted by genomic coordinates and multi-mappings filtered out"
	trimmed1=$SE_READS/${sample_id}_read1.fastq.gz
	step_log=$LOGS/${sample_id}_${step}.log
	TMP_DIR=$BWA/my_tmp
	mkdir -p $TMP_DIR	
	tbam=$BWA/${sample_id}_sorted.bam
	obam=$BWA/${sample_id}_sorted_unique.bam
	read_group="@RG\tID:'$sample_id'\tLB:'$sample_id'\tPL:illumina\tPU:'$sample_id'\tSM:'$sample_id'"
	$bwa mem -t $slots -M $genome_fasta -R $read_group $trimmed1 -v 0 |$samtools sort -o $tbam -O bam -T $TMP_DIR/$sample_id - >$step_log
	$samtools view -bq 1 $tbam > $obam
	n_reads_aligned=`$samtools view $tbam | wc -l`
	n_reads_unique=`$samtools view $obam | wc -l`
	message_info $step "reads aligned = $n_reads_aligned"
	message_info $step "reads unique = $n_reads_unique"
	rm -fr $TMP_DIR
	rm -f $tbam

	message_time_step $step $time0

}


# =================================================================================================
# Make tag directory with HOMER
# =================================================================================================

make_tag_directory() {

	step="make_tag_directory"
	time0=$(date +"%s")

	# to facilitate the analysis of ChIP-Seq (or any other type of short read re-sequencing data)
	# it is useful to first transform the sequence alignment into platform independent data structure representing the experiment
	# analogous to loading the data into a database.

	# Convert BAM to BED --required for making tag directory
	message_info $step "converting BAM to BED --required for making tag directory"
	ibam=$BWA/${sample_id}_sorted_unique.bam
	obed=$BWA/${sample_id}_sorted_unique.bed
	$bamToBed -i $ibam | awk '{OFS="\t"; $4="."; print $0}' > $obed

	# Make tag directory with HOMER
	message_info $step "Make tag directory with HOMER"
	mkdir -p $TAG_DIR
	step_log=$LOGS/${sample_id}_${step}.log
	$makeTagDirectory $TAG_DIR $obed -format bed -genome $genome_fasta -checkGC 2>$step_log

	# parse step log to extract generated metadata
	message_info $step "parse step log to extract generated metadata"
	# parse
	estimated_genome_size=`grep "Estimated genome size" $step_log | cut -f2 -d'=' | sed "s/ //g"`
	estimated_average_read_density_per_bp=`grep "Estimated average read density" $step_log | cut -f2 -d'=' |sed "s/ //g" |sed "s/perbp//g"`
	total_tags=`grep "Total Tags" $step_log | cut -f2 -d'=' | sed "s/ //g"`
	total_positions=`grep "Total Positions" $step_log | cut -f2 -d'=' | sed "s/ //g"`
	avg_tag_length=`grep "Average tag length" $step_log | cut -f2 -d'=' | sed "s/ //g"`
	median_tags_per_position=`grep "Median tags per position" $step_log | cut -f2 -d'=' | cut -f1 -d'(' |sed "s/ //g"`
	avg_tags_per_position=`grep "Average tags per position" $step_log | cut -f2 -d'=' | sed "s/ //g"`
	fragment_length_estimate=`grep "Fragment Length Estimate" $step_log | cut -f2 -d':' | sed "s/ //g"`
	peak_width_estimate=`grep "Peak Width Estimate" $step_log | cut -f2 -d':' | sed "s/ //g"`
	autocorrelation_same_strand_fold_enrichment=`grep "Same strand fold enrichment" $step_log | cut -f2 -d':' | sed "s/ //g"`
	autocorrelation_diff_strand_fold_enrichment=`grep "Diff strand fold enrichment" $step_log | cut -f2 -d':' | sed "s/ //g"`
	autocorrelation_same_to_diff_strand_fold_enrichment=`grep "Same / Diff fold enrichment" $step_log | cut -f2 -d':' | sed "s/ //g"`
	avg_fragment_gc=`grep "Avg Fragment GC" $step_log | cut -f2 -d'=' | sed "s/ //g" |sed "s/%//g"`
	avg_expected_gc=`grep "Avg Expected GC" $step_log | cut -f2 -d'=' | sed "s/ //g" |sed "s/%//g"`
	# print
	message_info $step "estimated genome size = $estimated_genome_size"
	message_info $step "estimated average read density per bp = $estimated_average_read_density_per_bp"
	message_info $step "total tags = $total_tags"
	message_info $step "total_positions = $total_positions"
	message_info $step "avg. tag length = $avg_tag_length"
	message_info $step "median tags per position = $median_tags_per_position"
	message_info $step "avg. tags per position = $avg_tags_per_position"
	message_info $step "fragment length estimate = $fragment_length_estimate"
	message_info $step "peak width estimate = $peak_width_estimate"
	message_info $step "autocorrelation: same strand fold enrichment = $autocorrelation_same_strand_fold_enrichment"
	message_info $step "autocorrelation: diff strand fold enrichment = $autocorrelation_diff_strand_fold_enrichment"
	message_info $step "autocorrelation: same-to-diff strand fold enrichment = $autocorrelation_same_to_diff_strand_fold_enrichment"
	message_info $step "avg. fragment GC% = $avg_fragment_gc"
	message_info $step "avg. expected GC% = $avg_expected_gc"

	message_time_step $step $time0

}


# =================================================================================================
# Make BigBed file
# =================================================================================================

make_bigbed() {

	step="make_bigbed"
	time0=$(date +"%s")

	# Make BigBed file
	message_info $step "make BigBed file"
	step_log=$LOGS/${sample_id}_${step}.log
	ibed=$BWA/${sample_id}_sorted_unique.bed
	obb=$BWA/${sample_id}_sorted_unique.bb
	$bedToBigBed $ibed $genome_chrom_sizes $obb 2>$step_log
	
	message_time_step $step $time0

}


# =================================================================================================
# Calculate reads per million (RPM)
# =================================================================================================

calculate_rpms() {

	step="calculate_rpms"
	time0=$(date +"%s")

	# Get fragment length estimate as calculated in the 'make_tag_directory' step
	message_info $step "get fragment length estimate (l) as calculated in the 'make_tag_directory' step"
	fragment_length_estimate=`grep "Fragment Length Estimate" $LOGS/${sample_id}_make_tag_directory.log | cut -f2 -d':' | sed "s/ //g"`
	fragment_length_estimate_corrected=`cat $TAG_DIR/tagInfo.txt | grep fragmentLengthEstimate |cut -f 2 -d"=" | sed 's/[^0-9]*//g'`
	message_info $step "fragment length estimate = $fragment_length_estimate"
	message_info $step "fragment length correction = $fragment_length_estimate_corrected"
	message_info $step "the correction will be use if the estimate is not reliable"

	# Generate RPM fragment profile
	message_info $step "generate reads pe million profile (RPM) fragment profile"
	ibam=$BWA/${sample_id}_sorted_unique.bam
	mkdir -p $RPMS
	orpm=$RPMS/$sample_id.rpm
	step_log=$LOGS/${sample_id}_${step}.log
	$perl $bam2wig --bw --bwapp $bedGraphToBigWig --pos extend --ext $fragment_length_estimate_corrected --rpm --in $ibam --out $orpm > $step_log
	
	message_time_step $step $time0

}


# =================================================================================================
# Call peaks
# =================================================================================================

call_peaks() {

	step="call_peaks"
	time0=$(date +"%s")

	# Get fragment length estimate as calculated in the 'make_tag_directory' step
	message_info $step "Get fragment length estimate (l) as calculated in the 'make_tag_directory' step"
	fragment_length_estimate_corrected=`cat $TAG_DIR/tagInfo.txt | grep fragmentLengthEstimate |cut -f 2 -d"=" | sed 's/[^0-9]*//g'`
	message_info $step "Fragment length (l) is $fragment_length_estimate_corrected bp (note this is not used if peak caller is zerone)"

	mkdir -p $PEAKS
	ibam=$BWA/${sample_id}_sorted_unique.bam
	ibed=$BWA/${sample_id}_sorted_unique.bed	

	# Peak calling with MACS2
	if [[ $peak_caller == "macs2" ]]; then
		
		macs2=`which $peak_caller`
		message_info $step "Peak calling with MACS2"
		if [[ $species == "homo_sapiens" ]]; then
			genome_size="hs"
		fi
		message_info $step "genome size for $species will be used"
		message_info $step "q-value cutoff = $macs2_q (default is 0.01)"
		message_info $step "--nomodel (MACS2 will not try to model peak length but use l = $fragment_length_estimate_corrected instead"
		message_info $step "--call-summits = MACS reanalyzes the shape of signal profile to deconvolve subpeaks within each peak called"
		# Peak calling with input DNA as control
		if [[ $use_control == "yes" ]]; then
			if [ ! -f $control_bam ]; then
				message_error $step "$control_bam not found. Exiting..."
			fi
		 	OUTDIR=$PEAKS/with_control
		 	mkdir -p $OUTDIR
		 	step_log=$LOGS/${sample_id}_${step}_${peak_caller}_with_control.log
		 	message_info $step "peak calling with input DNA ($control_bam) as control"
		 	$macs2 callpeak -t $ibam \
		 					-c $control_bam \
		 					-q $macs2_qvalue \
		 					--nomodel \
		 					--extsize $fragment_length_estimate_corrected \
		 					--outdir $OUTDIR \
		 					--call-summits \
		 					--name $sample_id 2>$step_log
		# Peak calling with the sample alone
		elif [[ $use_control == "no" ]]; then
			OUTDIR=$PEAKS/sample_alone
			mkdir -p $OUTDIR
			step_log=$LOGS/${sample_id}_${step}_${peak_caller}_sample_alone.log
			message_info $step "peak calling with the sample alone (i.e. no input)"
			$macs2 callpeak -t $ibam \
							-q $macs2_qvalue \
							--nomodel \
							--extsize $fragment_length_estimate_corrected \
							--outdir $OUTDIR \
		 					--call-summits \
							--name $sample_id 2>$step_log
		fi

	# Peak calling with Zerone
	elif [[ $peak_caller == "zerone" ]]; then
		zerone=`which $peak_caller`
		message_info $step "peak calling with zerone"
		if [[ $use_control == "yes" ]]; then
			if [ ! -f $control_bam ]; then
				message_error $step "$control_bam not found. Exiting..."
			fi
			OUTDIR=$PEAKS/with_control
			mkdir -p $OUTDIR
			step_log=$LOGS/${sample_id}_${step}_${peak_caller}_with_control.log
			message_info $step "peak calling with input DNA ($control_bam) as control"
			# Zerone parameters
			# -l = zerone produces an alternative output in which 
			# (i) only enriched windows (i.e. peaks) are shown
			# (ii) contiguous windows are merged
			$zerone -l -0 $control_bam -1 $ibam > $OUTDIR/${sample_id}_zerone.txt
		# A control is mandatory for zerone, thus exiting otherwise
		else
			message_error $step "a control is mandatory for zerone. Exiting..."
		fi
	fi

}

main


