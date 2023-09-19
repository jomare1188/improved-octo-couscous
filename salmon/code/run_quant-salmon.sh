#!/bin/bash
#$ -q all.q
#$ -cwd
#$ -pe smp 5
#$ -t 1-18
#$ -tc 10

module load salmon/1.8.0

INPUT_R1=`ls -1 ../../bbduk/results/seqs/*_R1.bbduk.fastq | head -n $SGE_TASK_ID | tail -n 1`
INPUT_R2=${INPUT_R1/_R1.bbduk.fastq/_R2.bbduk.fastq}
PREFIX_R1=$(echo $INPUT_R1 | cut -f6 -d/ | cut -f1,2,3 -d_)
PREFIX_R2=$(echo $INPUT_R2 | cut -f6 -d/ | cut -f1,2,3 -d_)
PREFIX=$(echo $PREFIX_R1 | cut -f1,2 -d_)

# run salmon
salmon quant --libType A --threads 5 --index ../data/salmon_index/ --validateMappings --seqBias --posBias --softclip -1 ${INPUT_R1} -2 ${INPUT_R2} -o ../results/${PREFIX}
