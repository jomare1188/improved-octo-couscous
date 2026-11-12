#!/bin/bash
#$ -N dolores_gtcheck
#$ -q all.q
#$ -cwd
#$ -S /bin/bash
#$ -pe smp 2
#$ -j y
#$ -o logs/genotype_check.$JOB_ID.log
source /Storage/progs/miniconda3/etc/profile.d/conda.sh
conda activate r_RNASEQ
Rscript code/qc_genotype_check.R
