#!/bin/bash
#$ -N chipseq_call_peaks_20151110_mapping_T_30_PR_R5020_11729_GTGAAA
#$ -q short-sl65
#$ -l virtual_free=40G
#$ -l h_rt=06:00:00
#$ -o /users/GR/mb/jquilez/pipelines/chipseq_call_peaks/20151110/job_out/chipseq_call_peaks_20151110_mapping_T_30_PR_R5020_11729_GTGAAA_$JOB_ID.out
#$ -e /users/GR/mb/jquilez/pipelines/chipseq_call_peaks/20151110/job_out/chipseq_call_peaks_20151110_mapping_T_30_PR_R5020_11729_GTGAAA_$JOB_ID.err
#$ -M javier.quilez@crg.eu
#$ -m ae
#$ -pe smp 8
/software/mb/bin/pigz -dc -p 8 /users/GR/mb/jquilez/data/chipseq/raw/2015-11-03/T_30_PR_R5020_11729_GTGAAA.fastq.gz | /software/mb/bin/bowtie -p 8 -m 1 -S /users/GR/mb/dsoronellas/tracks/human_genome_19/bowtie_index_t47d/hg19-noMY-noMMTV - | /software/mb/bin/samtools view -@ 8 -Sb - > /users/GR/mb/jquilez/data/chipseq/samples/T_30_PR_R5020_11729_GTGAAA/alignments/bowtie/T_30_PR_R5020_11729_GTGAAA.bam
/software/mb/bin/samtools sort -@ 8 -m 4G /users/GR/mb/jquilez/data/chipseq/samples/T_30_PR_R5020_11729_GTGAAA/alignments/bowtie/T_30_PR_R5020_11729_GTGAAA.bam -f /users/GR/mb/jquilez/data/chipseq/samples/T_30_PR_R5020_11729_GTGAAA/alignments/bowtie/T_30_PR_R5020_11729_GTGAAA_sorted.bam
