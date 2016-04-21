# README
---------------------------------------------------------------------------------------------------

**Pipeline to quantify gene/transcript abundance using RNA-seq data**


## To-do
- convert alignments to reads per million (RPM) profiles, use rpm directory as in chipseq-16.03 **use cpu option**
- integrate metadata


## NEW
- trimming low-quality reads ends using trimmomatic's recommended practices and the Maxinfo approach
- pipeline accepts both single-end and paired-end data
- allows alignment to hg19 and hg38 assembly vesions including the MMTV luciferase construct as a pseudo-chromosome
- quantification with HTSeq is deprecated


## Modules

1. `trim_reads_trimmomatic`: trim sequencing adapters and low-quality ends from the reads using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
2. `align_star`: align reads to genome with [STAR](https://github.com/alexdobin/STAR)
3. `quantification_featurecounts`: quantification of reads counts per gene, using alignments from `align_star`, with [featureCounts](http://bioinf.wehi.edu.au/featureCounts/)
4. `make_profiles`: make read per million (RPM) profiles using alignments from `align_star`
5. `quantification_kallisto`: pseudo-quantification + quantification of read counts per transcript with [Kallisto](http://pachterlab.github.io/kallisto/)

The pipeline can be configured to run all steps or just one (see below) considering that:
- **`trim_reads_trimmomatic` needs to be run at least once, either running the full pipeline (`full` mode) or this step alone (`trim_reads_trimmomatic` mode), as it checks that the input FASTQ files do exist and, if so, creates the sample directory**
- `quantification_kallisto` can be run alone provided `trim_reads_trimmomatic` has also been run
- `make_profiles` and `quantification_featurecounts` require that `align_star` has been run


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

It contains 2 FASTQ files (paired-end) with 1,000 reads each.
