#!/bin/bash
#$ -N ips6_POJ2878
#$ -q all.q
#$ -cwd
#$ -S /bin/bash
#$ -pe smp 2
#$ -j y
#$ -o logs/ips6.$JOB_ID.log
# InterProScan 6 head job. Each batch x member database is sent to SGE as its own job
# (sge_ips6.config), so the proteome is split across nodes.
#   cd /Storage/data1/jorge.munoz/DOLORES/annotation/interproscan
#   qsub run_ips6.sh                                  # full unique proteome -> out/
#   qsub run_ips6.sh <fasta> <outdir>                 # e.g. the test set
# Resubmitting the same command resumes (-resume).
set -euo pipefail
source /etc/profile.d/modules.sh
module load miniconda3 singularity-ce/4.3.3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate nextflow-26.04.7

BASE=/Storage/data1/jorge.munoz/DOLORES/annotation/interproscan
INPUT=${1:-/Storage/data1/jorge.munoz/DOLORES/annotation/proteins/POJ2878_v1.proteins.uniq.faa}
OUT=${2:-$BASE/out}
BATCH=${3:-20000}

export NXF_OPTS='-Xms1g -Xmx8g'
export NXF_SINGULARITY_CACHEDIR=$BASE/singularity_cache

nextflow run ebi-pf-team/interproscan6 -r 6.0.3.0 -profile singularity \
    -c $BASE/sge_ips6.config \
    --input "$INPUT" \
    --datadir $BASE/ips6_data \
    --interpro 110.0 \
    --goterms --pathways \
    --formats tsv,json,gff3 \
    --batchSize $BATCH --subBatchSize 2000 \
    --cpus 4 \
    --outdir "$OUT" \
    -work-dir $BASE/work \
    -resume
