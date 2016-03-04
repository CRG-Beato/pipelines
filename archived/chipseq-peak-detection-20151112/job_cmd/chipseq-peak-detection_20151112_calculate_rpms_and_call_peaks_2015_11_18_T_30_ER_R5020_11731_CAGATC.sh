#!/bin/bash
#$ -N chipseq-peak-detection_20151112_calculate_rpms_and_call_peaks_2015_11_18_T_30_ER_R5020_11731_CAGATC
#$ -q short-sl65
#$ -l virtual_free=40G
#$ -l h_rt=6:00:00
#$ -o /users/GR/mb/jquilez/pipelines/chipseq-peak-detection-20151112/job_out/chipseq-peak-detection_20151112_calculate_rpms_and_call_peaks_2015_11_18_T_30_ER_R5020_11731_CAGATC_$JOB_ID.out
#$ -e /users/GR/mb/jquilez/pipelines/chipseq-peak-detection-20151112/job_out/chipseq-peak-detection_20151112_calculate_rpms_and_call_peaks_2015_11_18_T_30_ER_R5020_11731_CAGATC_$JOB_ID.err
#$ -j y
#$ -M javier.quilez@crg.eu
#$ -m abe
#$ -pe smp 8

submitted_on=2015_11_18
pipeline_version=20151112
sample_id=T_30_ER_R5020_11731_CAGATC
pipeline_run_mode=calculate_rpms_and_call_peaks
submit_to_cluster=yes
queue=short-sl65
memory=40G
max_time=6:00:00
slots=8
peak_caller=macs2
control_bam=$HOME/non_existing_control.bam
macs2_q=0.01
sicer_redundancy_threshold=1
sicer_window_size=100
sicer_fdr=0.001
sicer_effective_genome_fraction=0.8
PIPELINE=/users/GR/mb/jquilez/pipelines/chipseq-peak-detection-20151112
config=/users/GR/mb/jquilez/pipelines/chipseq-peak-detection-20151112/chipseq-peak-detection.config
time_start=$(date +"%s")




#==================================================================================================
# PATHS
#==================================================================================================

# Directories
DATA=$HOME/data
SCRIPTS=$PIPELINE/scripts
SAMPLE=$HOME/data/chipseq/samples/$sample_id
LOGS=$SAMPLE/logs 
BAMS=$SAMPLE/alignments/bowtie
TAG_DIR=$SAMPLE/tag_directory
RPMS=$SAMPLE/rpms
PEAKS=$SAMPLE/peaks/$peak_caller

# Files 
fq=$HOME/data/chipseq/raw/*/$sample_id.fastq.gz
pigz=`which pigz`
bowtie=`which bowtie`
# + as of 2015-11-03, the `which samtools` version of samtools fails to sort the BAM of some samples
# Quique mentioned he had a similar problem and recommended me to use the '/software/mb/el6.3/samtools-1.2/samtools'
# version
samtools=/software/mb/el6.3/samtools-1.2/samtools
index=/users/GR/mb/dsoronellas/tracks/human_genome_19/bowtie_index_t47d/hg19-noMY-noMMTV
# There are certain regions of the genome (mostly repeats) on which reads pile up
# As such regions are known from different papers, Dani combined and subtracted them from the genome
# so that we will used them to only retain reads mapping on them
regions_to_include=$HOME/assemblies/homo_sapiens/hg19/hg19_20150817_complement_to_highly_duplicated_regions.bed
regions_to_exclude=$HOME/assemblies/homo_sapiens/hg19/hg19_20150817_highly_duplicated_regions.bed
bamToBed=`which bamToBed`
makeTagDirectory=`which makeTagDirectory`
sortBed=`which sortBed`
bedToBigBed=`which bedToBigBed`
perl=`which perl`
bam2wig=`which bam2wig.pl`
bedGraphToBigWig=`which bedGraphToBigWig`
pyicotrocol_template=$SCRIPTS/pyicotrocol_template.ini



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
	python $SCRIPTS/update_and_get_samples_metadata.py $metadata update_metrics $sample_id $field_name $length
	echo -e "TIME \t`date +"%Y-%m-%d %T"` \t[pipeline] \ttotal time for completion (seconds) = $length"
	echo
}


# =================================================================================================
# Mapping
# =================================================================================================

mapping() {

	step="mapping"
	time0=$(date +"%s")

	# Check that FASTQ file exists
	if [ -f $fq ]; then
		mkdir -p $SAMPLE
		mkdir -p $BAMS
		mkdir -p $LOGS
		log_mapping=$LOGS/${sample_id}_${step}.log
	else
		message_error $step "$fq not found. Exiting..."
	fi

	# Mapping with bowtie
	# -m 1 = only report uniquely mapped reads
	# -k 1 = report up to 1 good alignments per read
	# --best --strata = hits in sub-optimal strata aren't reported (requires --best)
	# -n 2 = maximum 2 mismatches in seed
	# --chunkmbs = The number of megabytes of memory a given thread is given to store path descriptors in --best mode.
	# Default is 64Mb but I noted that it was insufficient and produced the following error:
	# "Warning: Exhausted best-first chunk memory for read XXX (patid YYY); skipping read"
	# -F 1804 = excludes invalid mappings 
	# from https://broadinstitute.github.io/picard/explain-flags.html: 1804 corresponds which are:
	# read unmapped, mate unmapped, read reverse or mate reverse strand
	message_info $step "mapping with bowtie; only unique and valid mappings are reported"
	obam=$BAMS/$sample_id.bam
	obam_sorted=$BAMS/${sample_id}_sorted.bam
	$pigz -dc -p $slots $fq \
			| $bowtie -p $slots -m 1 -k 1 --best --strata -n 2 --chunkmbs 250 -S $index - 2>$log_mapping \
			| $samtools view -F 1804 -L $regions_to_include -Sb - > $obam 
	
	# Sort and create index with samtools
	message_info $step "sort and index BAM"
	$samtools sort -@ $slots -m 1G -o $obam_sorted -T $BAMS/tmp $obam
	$samtools index $obam_sorted
	#rm -f $obam

	message_time_step $step $time0

}


# =================================================================================================
# Make tag directory with HOMER
# =================================================================================================

make_tag_directory() {

	step="make_tag_directory"
	time0=$(date +"%s")

	# Convert BAM to BED --required for making tag directory
	message_info $step "converting BAM to BED --required for making tag directory"
	ibam=$BAMS/${sample_id}_sorted.bam
	obed=$BAMS/${sample_id}_sorted.bed
	$bamToBed -i $ibam | awk '{OFS="\t"; $4="."; print $0}' > $obed

	# Make tag directory with HOMER
	message_info $step "Make tag directory with HOMER"
	mkdir -p $TAG_DIR
	log_make_tag_directory=$LOGS/${sample_id}_${step}.log
	$makeTagDirectory $TAG_DIR $obed -format bed -genome hg19 -checkGC 2>$log_make_tag_directory

	message_time_step $step $time0

}



# =================================================================================================
# Make BigBed file
# =================================================================================================

make_bigbed() {

	step="make_bigbed"
	time0=$(date +"%s")

	# Sort BED file
	#message_info $step "Sort BED file"
	#ibed=$BAMS/${sample_id}_sorted.bed
	#tbed=$BAMS/${sample_id}_tmp.bed
	#obed=$BAMS/${sample_id}_sorted.bed
	#$sortBed -i $ibed > $tbed
	#mv $tbed $obed

	# Make BigBed file
	message_info $step "Make BigBed file"
	log_make_bigbed=$LOGS/${sample_id}_${step}.log
	obb=$BAMS/${sample_id}_sorted.bb 
	$bedToBigBed $obed <(cut -f1,2 assemblies/homo_sapiens/hg19/t47d_genome.fa.fai) $obb 2>$log_make_bigbed
	
	message_time_step $step $time0

}


# =================================================================================================
# Calculate reads per million (RPM)
# =================================================================================================

calculate_rpms() {

	step="calculate_rpms"
	time0=$(date +"%s")

	# Get fragment length estimate as calculated in the 'make_tag_directory' step
	message_info $step "Get fragment length estimate (l) as calculated in the 'make_tag_directory' step"
	fragment_length=`cat $TAG_DIR/tagInfo.txt | grep fragmentLengthEstimate |cut -f 2 -d"=" | sed 's/[^0-9]*//g'`

	# Generate RPM fragment profile
	message_info $step "Generate RPM fragment profile"
	ibam=$BAMS/${sample_id}_sorted.bam
	mkdir -p $RPMS
	orpm=$RPMS/$sample_id.rpm
	log_calculate_rpms=$LOGS/${sample_id}_${step}.log
	$perl $bam2wig --bw --bwapp $bedGraphToBigWig --pos extend --ext $fragment_length --rpm --in $ibam --out $orpm > $log_calculate_rpms
	
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
	fragment_length=`cat $TAG_DIR/tagInfo.txt | grep fragmentLengthEstimate |cut -f 2 -d"=" | sed 's/[^0-9]*//g'`
	message_info $step "Fragment length (l) is $fragment_length bp"

	mkdir -p $PEAKS
	ibam=$BAMS/${sample_id}_sorted.bam
	ibed=$BAMS/${sample_id}_sorted.bed	

	# Peak calling with MACS2
	if [[ $peak_caller == "macs2" ]]; then
		macs2=`which $peak_caller`
		message_info $step "Peak calling with MACS2"
		message_info $step "q-value cutoff = $macs2_q (default is 0.01)"
		message_info $step "--nomodel (MACS2 will not try to model peak length but use l = $fragment_length instead"
		# Peak calling with input DNA as control
		if [ -f $control_bam ]; then
			OUTDIR=$PEAKS/with_control
			mkdir -p $OUTDIR
			log_call_peaks=$LOGS/${sample_id}_${step}_${peak_caller}_with_control.log
			message_info $step "Peak calling with input DNA ($control_bam) as control"
			$macs2 callpeak -t $ibam \
							-c $control_bam \
							-q $macs2_q \
							--nomodel \
							--extsize $fragment_length \
							--outdir $OUTDIR \
							--name $sample_id 2>$log_call_peaks
		# Peak calling with the sample alone
		else
			OUTDIR=$PEAKS/sample_alone
			mkdir -p $OUTDIR
			log_call_peaks=$LOGS/${sample_id}_${step}_${peak_caller}_sample_alone.log
			message_info $step "Peak calling with the sample alone (i.e. no input)"
			$macs2 callpeak -t $ibam \
							-q $macs2_q \
							--nomodel \
							--extsize $fragment_length \
							--outdir $OUTDIR \
							--name $sample_id 2>$log_call_peaks
		fi

	# Peak calling with SICER
	elif [[ $peak_caller == "sicer" ]]; then
		gap_size=`expr 3 \* $sicer_window_size`
		message_info $step "Peak calling with SICER"
		message_info $step "species = hg19"
		message_info $step "redundancy threshold = 1 (1 = remove duplicates)"
		message_info $step "window size (bp) = $sicer_window_size"
		message_info $step "fragment length (bp) = $fragment_length"
		message_info $step "gap size (bp) = $gap_size (typically, 3 x window size)"
		message_info $step "FDR = $sicer_fdr (default = 0.1%)"
		# Peak calling with input DNA as control
		if [ -f $control_bam ]; then
			OUTDIR=$PEAKS/with_control
			mkdir -p $OUTDIR
			log_call_peaks=$LOGS/${sample_id}_${step}_${peak_caller}_with_control.log
			message_info $step "Peak calling with input DNA ($control_bam) as control"
			sicer=/software/mb/el6.3/SICER_V1.1/SICER/SICER.sh
			control_bed=$BAMS/tmp_input_control.bed
			$bamToBed -i $control_bam | awk '{OFS="\t"; $4="."; print $0}' > $control_bed
			read_length=`awk '{print $3-$2}' $ibed | head -n 1`
			message_info $step "read length (bp) = $read_length"	
			$sicer $BAMS \
					`basename $ibed` \
					`basename $control_bed`\
					$OUTDIR \
					hg19 \
					$sicer_redundancy_threshold \
					$sicer_window_size \
					$fragment_length \
					$read_length \
					$gap_size \
					$sicer_fdr
			rm -f $control_bed
		# Peak calling with the sample alone
		else
			OUTDIR=$PEAKS/sample_alone
			mkdir -p $OUTDIR
			log_call_peaks=$LOGS/${sample_id}_${step}_${peak_caller}_sample_alone.log
			message_info $step "Peak calling with the sample alone (i.e. no input)"
			sicer=`which SICER-rb.sh`
			message_info $step "Effective genome fraction = $sicer_effective_genome_fraction"		
			$sicer $BAMS \
					`basename $ibed` \
					`basename $control_bed`\
					$OUTDIR \
					hg19 \
					$sicer_redundancy_threshold \
					$sicer_window_size \
					$fragment_length \
					$sicer_effective_genome_fraction \
					$gap_size \
					$sicer_fdr
		fi

	# Peak calling with Pyicos
	elif [[ $peak_caller == "pyicos" ]]; then
		pyicos=`which $peak_caller`
		message_info $step "Peak calling with Pycos"
		# Peak calling with input DNA as control
		if [ -f $control_bam ]; then
			OUTDIR=$PEAKS/with_control
			mkdir -p $OUTDIR
			log_call_peaks=$LOGS/${sample_id}_${step}_${peak_caller}_with_control.log
			message_info $step "Peak calling with input DNA ($control_bam) as control"
			control_bed=$BAMS/tmp_input_control.bed
			$bamToBed -i $control_bam | awk '{OFS="\t"; $4="."; print $0}' > $control_bed
			# Make pycotrocol file
			pyicotrocol_file=$OUTDIR/${sample_id}_pyicotrocol_file.ini
			ofile=$OUTDIR/${sample_id}.bedpk
			cat $pyicotrocol_template | sed 's/PATH_TO_EXPERIMENT_BED/$ibed/g' > $pyicotrocol_file
			#cat $pycotrocol_file | sed "s/PATH_TO_EXPERIMENT_BED/$ibed/g"
			#cat $pycotrocol_file | sed "s/PATH_TO_CONTROL_BED/$control_bed/g"
			#cat $pycotrocol_file | sed "s/PATH_TO_REGIONS_TO_EXCLUDE_BED/$regions_to_exclude/g"
			#cat $pycotrocol_file | sed "s/PATH_TO_SIGNIFICANT_PEAKS_BEDPK/$ofile/g"

		fi
	fi

	message_time_step $step $time0

}

# =================================================================================================
# MAIN CODE
# =================================================================================================

echo
#mapping
#make_tag_directory
#make_bigbed
calculate_rpms
call_peaks
echo

