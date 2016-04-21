# README
---------------------------------------------------------------------------------------------------

Pipelines for the analysis of next-generation sequencing data (e.g. RNA-seq, ChIP-seq, NET-seq). Old versions of the pipelines are in `archived`. 


## GitHub repository
---------------------------------------------------------------------------------------------------

This GitHub repository was created with:
- created repository named `pipelines` at [https://github.com/CRG-Beato](https://github.com/CRG-Beato)
- move to the `pipelines` directory in my computer
- in the terminal:
```
git init
git remote add origin https://github.com/CRG-Beato/pipelines.git
git add "*"
git commit -m "create repo"
git push -u origin master
```

## Log
---------------------------------------------------------------------------------------------------

**2016-03-04**
- `chipseq-16.03` includes:
	- both single- (`SE`) and paired-end (`PE`) data are accepted in `sequencing_type`
	- when `io_mode=custom`, multiple samples are accepted in `samples`
	- input metadata (e.g. `SAMPLE_ID`, `READ_LENGTH`) and output metrics from the pipeline (e.g. `N_UNIQUE_READS`, `N_PEAKS`) can be added to the `beato_lab_metadata.db` database if `integrate_metadata=yes`
- `chipseq-16.02` is moved to `archived`

**2016-04-21**
- `rnaseq-16.04` is released, with new features:
	- possibility to integrate metadata
	- genome reference sequence (used by STAR) includes chromosomes 1-22, X, Y, M and un-placed/un-localized scaffolds 
	- use of qualimap to generate quality metrics of the alignments
	- Kallisto's option to perform sequence based bias correction (such option failed in Kallisto versions prior to 0.42.5)
	- true read per million (RPM) profiles are generated, in rnaseq-16.02 reads piled up profiles were generated instead
- `rnaseq-16.02` is moved to `archived`