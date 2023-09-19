#!/bin/bash
#$ -q all.q
#$ -cwd
#$ -pe smp 2
#$ -t 1-18
#$ -tc 10

INPUT_R1=`ls -1 ../../RAWDATA/*_R1_001.fastq.gz | head -n $SGE_TASK_ID | tail -n 1`
INPUT_R2=${INPUT_R1/_R1_001.fastq.gz/_R2_001.fastq.gz}
PREFIX_R1=$(echo $INPUT_R1 | cut -f4 -d/ | cut -f1,2,3 -d_)
PREFIX_R2=$(echo $INPUT_R2 | cut -f4 -d/ | cut -f1,2,3 -d_)
PREFIX=$(echo $PREFIX_R1 | cut -f1,2 -d_)

module load bbmap/38.91

mkdir -p ../results/stats/ ../results/seqs/
bbduk.sh -Xmx40g threads=${NSLOTS} in=${INPUT_R1} in2=${INPUT_R2} refstats=../results/stats/${PREFIX}.refstats stats=../results/stats/${PREFIX}.stats out=../results/seqs/${PREFIX_R1}.bbduk.fastq out2=../results/seqs/${PREFIX_R2}.bbduk.fastq ref=../data/adapters.fa,../data/rfam-5.8s-database-id98.fasta,../data/silva-euk-18s-id95.fasta,../data/silva-arc-23s-id98.fasta,../data/silva-euk-28s-id98.fasta,../data/silva-bac-16s-id90.fasta,../data/rfam-5s-database-id98.fasta,../data/silva-bac-23s-id98.fasta,../data/silva-arc-16s-id95.fasta forcetrimleft=11 forcetrimright2=3 minlength=80 qtrim=w trimq=20 tpe=t tbo=t 
