#!/bin/bash
#$ -N dolores_story
#$ -q all.q
#$ -cwd
#$ -S /bin/bash
#$ -pe smp 2
#$ -j y
#$ -o logs/story.$JOB_ID.log
# Explanatory figures for results/report_story.md.   cd analysis && qsub code/run_story.sh
source /Storage/progs/miniconda3/etc/profile.d/conda.sh
conda activate r_RNASEQ
Rscript code/08_story_figures.R
