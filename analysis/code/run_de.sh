#!/bin/bash
#$ -N dolores_de
#$ -q all.q
#$ -cwd
#$ -S /bin/bash
#$ -pe smp 4
#$ -j y
#$ -o logs/de.$JOB_ID.log
# DESeq2 per genotype at one level.   cd analysis && qsub -v LEVEL=group code/run_de.sh
set -euo pipefail
source /Storage/progs/miniconda3/etc/profile.d/conda.sh
conda activate r_RNASEQ
Rscript code/03_de.R ${LEVEL:-group}
