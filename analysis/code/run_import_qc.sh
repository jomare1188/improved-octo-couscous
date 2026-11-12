#!/bin/bash
#$ -N dolores_qc
#$ -q all.q
#$ -cwd
#$ -S /bin/bash
#$ -pe smp 4
#$ -j y
#$ -o logs/import_qc.$JOB_ID.log
# Import + QC at one level (gene or group).  cd analysis && qsub -v LEVEL=gene code/run_import_qc.sh
set -euo pipefail
source /Storage/progs/miniconda3/etc/profile.d/conda.sh
conda activate r_RNASEQ
LEVEL=${LEVEL:-gene}
Rscript code/01_import.R $LEVEL
Rscript code/02_qc.R $LEVEL
