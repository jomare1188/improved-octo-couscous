#!/bin/bash
#$ -q all.q
#$ -cwd
#$ -S /bin/bash
#$ -pe smp 24
#$ -j y
#$ -o logs/$JOB_NAME.$JOB_ID.log
# Genome-wide lncRNA search, step 1 of 3: HISAT2 index of the whole POJ2878 genome (10.4 Gb, all allele copies).
# Extends the chromosome-1 pilot (18_lnc_pilot.sh) to the whole genome and to the 48 libraries of the analysis.
# Work files: /Storage/data2/jorge.munoz/lnc_genome.   Run from analysis/:  qsub -N lncidx code/21_lnc_index.sh
set -euo pipefail
export PATH=/Storage/data2/jorge.munoz/conda_envs/lnc/bin:$PATH
DOL=/Storage/data1/jorge.munoz/DOLORES
GDIR=$DOL/nf-rnaseq/run_POJ2878/results_pilot/genome
W=/Storage/data2/jorge.munoz/lnc_genome
mkdir -p $W/bam $W/st; cd $W
echo "host $(hostname) start $(date)"
# the 48 libraries used in the analysis (samples.tsv minus excluded_samples.tsv)
awk -F'\t' 'NR==FNR{e[$1]=1;next} FNR>1 && !($1 in e){print $2}' $DOL/analysis/data/excluded_samples.tsv $DOL/analysis/data/samples.tsv > libs.txt
[ -s idx.1.ht2l ] || [ -s idx.1.ht2 ] || hisat2-build -p 24 $GDIR/POJ2878_1.0.genome.fa idx > build.log 2>&1
echo "index ready $(date)"
