# README
---------------------------------------------------------------------------------------------------

**Pipeline to call binding site peaks using ChIP-seq data**


## NEW
- allows alignment to hg19 and hg38 assembly vesions including the MMTV luciferase construct as a pseudo-chromosome


## Modules
1. `trim_reads_single_end`: trim sequencing adapters and low-quality ends from the reads using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
2. `align_single_end`: align single-end reads with [BWA](http://bio-bwa.sourceforge.net/bwa.shtml)
3. `make_tag_directory`: generate quality reads of the ChIP-seq fragments with [HOMER](http://homer.salk.edu/homer/)
4. `make_bigbed`: convert alignments from BAM to BED format
5. `calculate_rpms`: generate reads per million (RPM) profiles
6. `call_peaks`: call binding site peaks with (a) [MACS2](https://github.com/taoliu/MACS) (either using control or not) or (b) [Zerone](https://github.com/gui11aume/zerone) (control is required)


## Scripts

- `chipseq.sh`: most of the code
- `chipseq_submit`: wrapper script that both:
	- retrieves configuration variables and parameter values from the `chipseq.config` file
	- (if applies) submits jobs (one per sample) to execute the pipeline in the CRG cluster
- `chipseq.config`: configuration file (see below)


## Execute pipeline

```
/pipeline_location/chipseq_submit.sh <*.config>
```

**Users other than me have no writting permissions for the `chipseq.config` file, so they need to provide their own file**


## Configuration file

- `pipeline_run_mode`:
	- by setting `trim_reads_single_end`, `align_single_end`, `make_tag_directory`, `make_bigbed`, `calculate_rpms`, or `call_peaks` the corresponding module (see above) is executed (**note modules are sequential so each cannot be run unles the preciding one is completed**) 
	- `full`: all the steps above in the sequential order
	- `full_no_call_peaks`: all the modules except `call_peaks` (this is useful when a sample will be used as control)


# Test dataset

It contains 2 FASTQ files (paired-end) with 1,000 reads each. This dataset was generate with:
```
INDIR=/users/GR/mb/jquilez/data/chipseq/raw/2016-02-04
ifq1=$INDIR/rf_001_01_01_chipseq_read1.fastq.gz
OUTDIR=/users/GR/mb/jquilez/pipelines/chipseq-20160217/test
mkdir -p $OUTDIR
zcat $ifq1 | head -n 4000 | gzip -c > $OUTDIR/test_read1.fastq.gz
```
