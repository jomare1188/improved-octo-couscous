#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -pe smp 2
#$ -t 1-18

module load FastQC/0.11.8

INPUT_R1=`ls -1 *_R1_001.fastq.gz | head -n $SGE_TASK_ID|tail -n 1`
INPUT_R2=${INPUT_R1/_R1_001.fastq.gz/_R2_001.fastq.gz}

fastqc --threads $NSLOTS --nogroup $INPUT_R1 $INPUT_R2
