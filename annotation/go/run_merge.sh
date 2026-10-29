#!/bin/bash
#$ -N merge_go
#$ -q all.q
#$ -cwd
#$ -S /bin/bash
#$ -pe smp 2
#$ -j y
#$ -o merge.$JOB_ID.log
# Merge InterProScan + eggNOG into GO / gene annotation tables (see merge_annotation.py).
#   cd /Storage/data1/jorge.munoz/DOLORES/annotation/go && qsub run_merge.sh
set -euo pipefail
D=/Storage/data1/jorge.munoz/DOLORES
python3 $D/annotation/go/merge_annotation.py \
    --ips $D/annotation/interproscan/out/part_{1,2,3,4}.faa.tsv \
    --eggnog $D/annotation/eggnog/out/POJ2878.emapper.annotations \
    --uniq-map $D/annotation/proteins/uniq_map.tsv \
    --tx2gene $D/nf-rnaseq/run_POJ2878/results/salmon/salmon.merged.tx2gene.tsv \
    --outdir $D/annotation/go
