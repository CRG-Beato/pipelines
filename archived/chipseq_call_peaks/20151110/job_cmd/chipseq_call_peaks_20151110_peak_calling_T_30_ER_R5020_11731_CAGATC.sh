#!/bin/bash
#$ -N chipseq_call_peaks_20151110_peak_calling_T_30_ER_R5020_11731_CAGATC
#$ -q short-sl65
#$ -l virtual_free=40G
#$ -l h_rt=06:00:00
#$ -o /users/GR/mb/jquilez/pipelines/chipseq_call_peaks/20151110/job_out/chipseq_call_peaks_20151110_peak_calling_T_30_ER_R5020_11731_CAGATC_$JOB_ID.out
#$ -e /users/GR/mb/jquilez/pipelines/chipseq_call_peaks/20151110/job_out/chipseq_call_peaks_20151110_peak_calling_T_30_ER_R5020_11731_CAGATC_$JOB_ID.err
#$ -M javier.quilez@crg.eu
#$ -m ae
#$ -pe smp 1
/software/mb/bin/macs2 callpeak -t /users/GR/mb/jquilez/data/chipseq/samples/T_30_ER_R5020_11731_CAGATC/alignments/bowtie/T_30_ER_R5020_11731_CAGATC_sorted.bam -n T_30_ER_R5020_11731_CAGATC --outdir /users/GR/mb/jquilez/data/chipseq/samples/T_30_ER_R5020_11731_CAGATC/peaks/macs2 -g hs -B --call-summits --verbose 3
