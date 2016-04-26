# README
---------------------------------------------------------------------------------------------------

**Pipeline to call binding site peaks using ChIP-seq data**

## New features
- genome reference sequence includes chromosomes 1-22, X, Y, M and un-placed/un-localized scaffolds 
- use of qualimap to generate quality metrics of the alignments
- `make_bigbed` is deprecated as it was the resulting `*.bb` file was not used
- values for the fields `read_length` and `sequencing_type` are retrieved from the metadata database of `integrate_metadata=yes`


## Modules
1. `trim_reads_trimmomatic`: trim sequencing adapters and low-quality ends from the reads using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
2. `align_bwa`: align reads with [BWA](http://bio-bwa.sourceforge.net/bwa.shtml)
3. `quality_alignments`: quality control of the mappings using [qualimap](http://qualimap.bioinfo.cipf.es/)
4. `make_tag_directory`: generate quality reads of the ChIP-seq fragments with [HOMER](http://homer.salk.edu/homer/)
5. `make_profiles`: generate reads per million (RPM) profiles
6. `call_peaks`: call binding site peaks with (a) [MACS2](https://github.com/taoliu/MACS) (either using control or not) or (b) [Zerone](https://github.com/gui11aume/zerone) (control is required)

The pipeline can be configured to run all steps or just one (see below) considering that to run the `n` step the `n-1` step has to be run at leas once (except for `n=1`, obviously).


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
	- by setting `trim_reads_single_end`, `align_single_end`, `make_tag_directory`, `make_bigbed`, `calculate_rpms`, or `call_peaks` the corresponding module (see above) is executed (**note modules are sequential so each cannot be run unless the preceding one is completed**) 
	- `full`: all the steps above in the sequential order
	- `full_no_call_peaks`: all the modules except `call_peaks` (this is useful when a sample will be used as control)
- `data_type`: chipseq, atacseq (if `io_mode=standard` this will be used to search for the input/output directory)
- `samples`: espace-separated list of samples
- `pipeline_run_mode`: any of the 5 steps described in the **Modules** section of `full` to run all of them sequentially
- `io_mode`:
	- `standard`: in/out directories defined in my home directory `/users/GR/mb/jquilez`
	- `custom`:	in/out directories defined by the user in `CUSTOM_IN` and `CUSTOM_OUT`, and the input FASTQs defined in `sample_to_fastqs.txt`
- `[cluster options]`: see [CRG cluster](http://www.linux.crg.es/index.php/Main_Page)
- **if `integrate_metadata=yes` (see `*.config`), do not execute any two steps of the pipeline simulateneously as it will block the database**
- `[trimmomatic]`: see [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- `species` and `version`: as of 2016-03-01, genome sequence FASTA files are available for `homo_sapiens`, `version` hg19 and hg38 (either with or without the MMTV construct)
- `peak_caller`: macs2, zerone
- `use_control`: yes/no to using a BAM as control or input when calling peaks
- `control_bam`: BAM file to be used as control or input; mandatory if `use_control=yes`


## Test dataset

2 samples downloaded from the SRA data repository:
```
IDIR=data/atacseq/raw/2016-03-01
ODIR=pipelines/chipseq-16.03/test
# 1860 reads
cp -v $IDIR/SRR1779699_read1.fastq.gz $ODIR/test1_read1.fastq.gz
cp -v $IDIR/SRR1779699_read2.fastq.gz $ODIR/test1_read2.fastq.gz
# 1962 reads
cp -v $IDIR/SRR1779697_read1.fastq.gz $ODIR/test2_read1.fastq.gz
cp -v $IDIR/SRR1779697_read2.fastq.gz $ODIR/test2_read2.fastq.gz
```