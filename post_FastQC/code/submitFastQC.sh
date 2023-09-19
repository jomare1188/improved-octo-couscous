#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -pe smp 2
#$ -t 1-18
#$ -tc 10

INPUT_R1=`ls -1 ../../bbduk/results/seqs/*_R1.bbduk.fastq | head -n $SGE_TASK_ID | tail -n 1`
INPUT_R2=${INPUT_R1/_R1.bbduk.fastq/_R2.bbduk.fa}


module load FastQC/0.11.8

fastqc --threads $NSLOTS --outdir ../results/ --nogroup $INPUT_R1 $INPUT_R2
