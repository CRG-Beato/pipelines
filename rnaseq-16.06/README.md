# README
---------------------------------------------------------------------------------------------------

**Pipeline to quantify gene/transcript abundance using RNA-seq data**


## New features
- read per million (RPM) profiles are generated directly with STAR in the `align_star` module
- the `make_profiles` is deprecated
- added the `clean_up` module to delete relatively large intermediate files which can be re-generated


## Modules

1. `trim_reads_trimmomatic`: trim sequencing adapters and low-quality ends from the reads using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
2. `align_star`: align reads to genome with [STAR](https://github.com/alexdobin/STAR)
3. `quality_alignments`: quality control of the mappings using [qualimap](http://qualimap.bioinfo.cipf.es/)
4. `quantification_featurecounts`: quantification of reads counts per gene, using alignments from `align_star`, with [featureCounts](http://bioinf.wehi.edu.au/featureCounts/)
5. `quantification_kallisto`: pseudo-quantification + quantification of read counts per transcript with [Kallisto](http://pachterlab.github.io/kallisto/)
6. `clean_up`: to delete relatively large intermediate files which can be re-generated

The pipeline can be configured to run all steps or just one (see below) considering that:
- `trim_reads_trimmomatic` needs to be run at least once, either running the full pipeline (`full` mode) or this step alone (`trim_reads_trimmomatic` mode), as it checks that the input FASTQ files do exist and, if so, creates the sample directory
- `quantification_kallisto` can be run alone provided `trim_reads_trimmomatic` has also been run
- `quantification_featurecounts` require that `align_star` has been run


## Scripts

- `rnaseq.sh`: most of the code
- `rnaseq_submit`: wrapper script that both:
	- retrieves configuration variables and parameter values from the `rnaseq.config` file
	- (if applies) submits jobs (one per sample) to execute the pipeline in the CRG cluster
- `rnaseq.config`: configuration file (see below)


## Execute pipeline

```
/pipeline_location/rnaseq_submit.sh <*.config>
```

**Users other than me have no writting permissions for the `rnaseq.config` file, so they need to provide their own file**


## Configuration file

- `data_type`: `rnaseq`
- `samples`: espace-separated list of samples
- `pipeline_run_mode`: any of the 5 steps described in the **Modules** section of `full` to run all of them sequentially
- `io_mode`:
	- `standard`: in/out directories defined in my home directory `/users/GR/mb/jquilez`
	- `custom`:	in/out directories defined by the user in `CUSTOM_IN` and `CUSTOM_OUT`, and the input FASTQs defined with `read1_fname` and, if paired-end, `read2_fname`
- `[cluster options]`: see [CRG cluster](http://www.linux.crg.es/index.php/Main_Page)
- `[trimmomatic]`: see [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- `species` and `version`: as of 2016-03-01, STAR genome files are only generated for `homo_sapiens`, `version` hg19 and hg38, and `read_length` 50, 75 and 100 bp
- `[kallisto]`: see [Kallisto](http://pachterlab.github.io/kallisto/)
- `[featureCounts]`: see [featureCounts](http://bioinf.wehi.edu.au/featureCounts/)


## Test dataset

```
IDIR=data/rnaseq/raw/2016-02-02
ODIR=pipelines/rnaseq-16.06/test
zcat $IDIR/fd_001_01_01_rnaseq_12661_ACAGTG_read1.fastq.gz |head -n 4000 | gzip -c > $ODIR/test1_read1.fastq.gz
zcat $IDIR/fd_001_01_01_rnaseq_12661_ACAGTG_read2.fastq.gz |head -n 4000 | gzip -c > $ODIR/test1_read2.fastq.gz
zcat $IDIR/fd_002_01_01_rnaseq_12662_GCCAAT_read1.fastq.gz |head -n 4000 | gzip -c > $ODIR/test2_read1.fastq.gz
zcat $IDIR/fd_002_01_01_rnaseq_12662_GCCAAT_read2.fastq.gz |head -n 4000 | gzip -c > $ODIR/test2_read2.fastq.gz
```