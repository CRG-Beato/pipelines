# README
---------------------------------------------------------------------------------------------------

**Pipeline to quantify transcript abundance using RNA-seq data**


# Workflow

1. trim sequencing adapters and low-quality ends from the reads using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
2. map with [STAR](https://github.com/alexdobin/STAR)
3. make read per million (RPM) profiles
4a. quantification with [HTSeq](http://www-huber.embl.de/HTSeq/doc/overview.html)
4b. quantification with [Kallisto](http://pachterlab.github.io/kallisto/)
4c. quantification with [featureCounts](http://bioinf.wehi.edu.au/featureCounts/)

The pipeline can be configured to run all steps or just one (see below).


# Scripts

- `rnaseq.sh`: most of the code
- `rnaseq_submit`: wrapper script that both:
	- retrieves configuration variables and parameter values from the `rnaseq.config` file
	- (if applies) submits jobs (one per sample) to execute the pipeline in the CRG cluster
- `rnaseq.config`: configuration file (see below)


# Execute pipeline

```
/pipeline_location/rnaseq_submit.sh *.config
```

Users other than me have no writting permissions for the `rnaseq.config` file, so they need to provide their own file


# Configuration file

*Not all parameters are described*

- samples: list of samples for which the pipeline will be run, no quoting, one empty espace between samples

- pipeline_run_mode: one of the following:
	- full = runs all steps in the pipeline (see above)
	- trim_reads = trim sequencing adapters and low-quality ends from the reads using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
	- map_star = map with [STAR](https://github.com/alexdobin/STAR)
	- make_profiles = make read per million (RPM) profiles
	- quantification_htseq = quantification with [HTSeq](http://www-huber.embl.de/HTSeq/doc/overview.html)
	- quantification_kallisto = quantification with [Kallisto](http://pachterlab.github.io/kallisto/)

- io_mode: one of the following:
	- standard: input is searched/directed to `/users/GR/mb/jquilez/data/..` (see `rnaseq.sh` for more details)
- CUSTOM_IN	= specify directory where input FASTQ files are located
- read1_fname = file name within `CUSTOM_IN` for read1
- read2_fname = file name within `CUSTOM_IN` for read2
- CUSTOM_OUT = specify directory for pipeline output

- submit_to_cluster	= yes/no; the following are only applied if submit_to_cluster=yes
- queue	= CRG cluster queue (e.g. short-sl65, long-sl65)
- memory = requested memory (e.g. 40G)
- max_time = maximum running time allowed (e.g. 06:00:00)
- slots = number of slots in the CRG cluster

- read_length = read length (bp) as of 2015-12-23 STAR genome index (required) are only available for this read length


# Test dataset

It contains 2 FASTQ files (paired-end) with 1,000 reads each. This dataset was generate with:

```
# make directory with today's date
RAW=$HOME/data/rnaseq/raw/2015-12-18
mkdir -p $RAW
# use on of Dani's RNA-seq FASTQ file
INDIR=/users/GR/mb/dsoronellas/RNAseq/raw_data/polyA-ss-siC-siPADI/
ifq1=$INDIR/siC-T0-R1.read1.fastq.gz
ifq2=$INDIR/siC-T0-R1.read2.fastq.gz
zcat $ifq1 | head -n 4000 | gzip -c > $RAW/test_read1.fastq.gz
zcat $ifq2 | head -n 4000 | gzip -c > $RAW/test_read2.fastq.gz
```



