#!/bin/bash
#$ -N chipseq_call_peaks_20151110_peak_calling_T_0_ER_MCF_7_11730_CGATGT
#$ -q short-sl65
#$ -l virtual_free=40G
#$ -l h_rt=06:00:00
#$ -o /users/GR/mb/jquilez/pipelines/chipseq_call_peaks/20151110/job_out/chipseq_call_peaks_20151110_peak_calling_T_0_ER_MCF_7_11730_CGATGT_$JOB_ID.out
#$ -e /users/GR/mb/jquilez/pipelines/chipseq_call_peaks/20151110/job_out/chipseq_call_peaks_20151110_peak_calling_T_0_ER_MCF_7_11730_CGATGT_$JOB_ID.err
#$ -M javier.quilez@crg.eu
#$ -m ae
#$ -pe smp 1
/software/mb/bin/macs2 callpeak -t /users/GR/mb/jquilez/data/chipseq/samples/T_0_ER_MCF_7_11730_CGATGT/alignments/bowtie/T_0_ER_MCF_7_11730_CGATGT_sorted.bam -n T_0_ER_MCF_7_11730_CGATGT --outdir /users/GR/mb/jquilez/data/chipseq/samples/T_0_ER_MCF_7_11730_CGATGT/peaks/macs2 -g hs -B --call-summits --verbose 3
