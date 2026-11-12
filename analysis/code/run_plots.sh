#!/bin/bash
#$ -N dolores_plots
#$ -q all.q
#$ -cwd
#$ -S /bin/bash
#$ -pe smp 2
#$ -j y
#$ -o logs/plots.$JOB_ID.log
# Figures: PCA, UpSet of DEGs, GO BP top-20 plots.   cd analysis && qsub code/run_plots.sh
source /Storage/progs/miniconda3/etc/profile.d/conda.sh
conda activate r_RNASEQ
Rscript code/06_plots.R
