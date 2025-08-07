#!/bin/bash
#$ -N exploratory_RNASEQ
#$ -cwd
#$ -V
#$ -pe smp 5
#$ -q all.q

module load miniconda3
conda activate r_RNASEQ

#Rscript preliminary_analysis_nice.r
Rscript RemoveUnwantedVariation.r

