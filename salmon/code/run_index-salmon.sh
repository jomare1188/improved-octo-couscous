#!/bin/bash
#$ -q all.q
#$ -cwd
#$ -pe smp 25
module load salmon/1.8.0
salmon index --threads $NSLOTS -t ../data/47_genotypes_SC-pan.fasta -i ../data/salmon_index/
