#!/bin/bash
#$ -N nf_POJ2878
#$ -q all.q
#$ -cwd
#$ -S /bin/bash
#$ -pe smp 2
#$ -j y
#$ -o logs/nf_POJ2878.$JOB_ID.log
#$ -l h=!telura*

# nf-core/rnaseq 3.27.0, Salmon-only, POJ2878 reference, all samples.
# This job is only the Nextflow "head": every pipeline step is sent to SGE as its own job
# (conf_cluster/sge_cluster.config), at most 10 at a time.
#
# Submit from this directory:
#   cd /Storage/data1/jorge.munoz/DOLORES/nf-rnaseq/run_POJ2878
#   qsub run_nfcore_POJ2878.sh                       # 51 clean libraries (samplesheet_POJ2878_clean.csv)
#   qsub run_nfcore_POJ2878.sh samplesheet_pilot.csv # another samplesheet
# Resubmitting the same command resumes from the last finished step (-resume).

set -euo pipefail
source /etc/profile.d/modules.sh
module load miniconda3 singularity-ce/4.3.3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate nextflow-26.04.7

BASE=/Storage/data1/jorge.munoz/DOLORES/nf-rnaseq
INPUT=${1:-samplesheet_POJ2878_clean.csv}

export NXF_OPTS='-Xms1g -Xmx8g'
export NXF_SINGULARITY_CACHEDIR=${BASE}/singularity_cache

nextflow run nf-core/rnaseq -r 3.27.0 -profile singularity \
    -params-file ${BASE}/params_proposal/params_POJ2878_salmon.yaml \
    -c ${BASE}/conf_cluster/sge_cluster.config \
    -c ${BASE}/params_proposal/process_overrides.config \
    --input "$INPUT" \
    -work-dir ${BASE}/run_POJ2878/work \
    -resume
