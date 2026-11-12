#!/bin/bash
#$ -q all.q
#$ -cwd
#$ -S /bin/bash
#$ -pe smp 4
#$ -j y
#$ -o logs/$JOB_NAME.$JOB_ID.log
# Run one R script of the follow-up analyses as an SGE job, from analysis/:
#   qsub -N h2 code/run_r.sh code/10_genotype_cycle.R
#   qsub -N wgcna -v RENV=/Storage/data2/jorge.munoz/conda_envs/r_wgcna code/run_r.sh code/15_wgcna.R
# RENV = conda env (name or path); default r_RNASEQ.
set -euo pipefail
source /Storage/progs/miniconda3/etc/profile.d/conda.sh
set +u; conda activate "${RENV:-r_RNASEQ}"; set -u
echo "host $(hostname)  start $(date)  env ${RENV:-r_RNASEQ}  script $*"
Rscript "$@"
echo "end $(date)"
