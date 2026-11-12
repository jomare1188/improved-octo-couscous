#!/bin/bash
#$ -N dolores_qc2
#$ -q all.q
#$ -cwd
#$ -S /bin/bash
#$ -pe smp 4
#$ -j y
#$ -o logs/qc.$JOB_ID.log
set -euo pipefail
source /Storage/progs/miniconda3/etc/profile.d/conda.sh
conda activate r_RNASEQ
Rscript code/02_qc.R ${LEVEL:-gene}
