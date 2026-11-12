#!/bin/bash
#$ -q all.q
#$ -cwd
#$ -S /bin/bash
#$ -pe smp 8
#$ -j y
#$ -o logs/$JOB_NAME.$JOB_ID.$TASK_ID.log
# Genome-wide lncRNA search, step 4 of 5: Salmon quantification of one library on the annotated transcripts
# plus the lncRNA candidates (index from 24_lnc_candidates.sh). Raw reads, library type ISR, as in the main run.
#   qsub -N lncquant -hold_jid lnccand -t 1-48 -tc 6 code/25_lnc_quant.sh        (from analysis/)
set -euo pipefail
SALMON=/home/jorge.munoz/.conda/envs/R2C/bin/salmon
DOL=/Storage/data1/jorge.munoz/DOLORES
W=/Storage/data2/jorge.munoz/lnc_genome
cd $W
L=$(sed -n "${SGE_TASK_ID}p" libs.txt)
echo "host $(hostname) library $L start $(date)"
[ -s quant/$L/quant.sf ] || $SALMON quant -i salmon_idx -l ISR -p 8 --seqBias --gcBias --validateMappings \
  -1 $DOL/RAWDATA/${L}_R1_001.fastq.gz -2 $DOL/RAWDATA/${L}_R2_001.fastq.gz -o quant/$L > quant/$L.log 2>&1
grep -h "Mapping rate" quant/$L/logs/salmon_quant.log | tail -1
echo "end $(date)"
