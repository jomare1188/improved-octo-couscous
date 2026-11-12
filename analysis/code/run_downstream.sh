#!/bin/bash
#$ -N dolores_down
#$ -q all.q
#$ -cwd
#$ -S /bin/bash
#$ -pe smp 4
#$ -j y
#$ -o logs/downstream.$JOB_ID.log
# Everything after OrthoFinder, in order. Each step is skipped if its output already exists, so the job can be
# resubmitted after a failure.   cd analysis && qsub -hold_jid of_orth code/run_downstream.sh
set -euo pipefail
source /Storage/progs/miniconda3/etc/profile.d/conda.sh
D=/Storage/data1/jorge.munoz/DOLORES
step() { echo; echo "=== $(date '+%F %T') $*"; }

WD=/Storage/data2/jorge.munoz/DOLORES_orthofinder/run/Results_Oct29/WorkingDirectory
[ -n "$(find $WD -type d -name Orthologues_POJ2878 -newer $WD/OrthoFinder/Results_Oct29/Log.txt)" ] || [ -s $D/orthology/allele_groups/tx2group.tsv ] || { echo "OrthoFinder orthologues not found"; exit 1; }

if [ ! -s $D/orthology/allele_groups/tx2group.tsv ]; then
  step allele groups; conda activate ORTHOFINDER
  python $D/orthology/make_allele_groups.py; conda deactivate
fi
conda activate r_RNASEQ
[ -s data/txi_group.rds ] || { step import group; Rscript code/01_import.R group; }
[ -s results/qc/group/sample_correlation.png ] || { step QC group; Rscript code/02_qc.R group; }
[ -s results/de/group/deg_summary.tsv ] || { step DESeq2 group; Rscript code/03_de.R group; }
[ -s results/memory_genes/camera_literature_sets.tsv ] || { step memory genes; Rscript code/04_memory_genes.R; }
[ -s results/enrichment/go_enrichment.tsv ] || { step enrichment; Rscript code/05_enrichment.R; }
step all done
