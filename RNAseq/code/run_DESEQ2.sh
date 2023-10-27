#!/bin/bash
#$ -q all.q
#$ -cwd
#$ -V
#$ -pe smp 6

echo $NSLOTS    "CORES"
module load R/4.0.0
Rscript make_DEA_by_groups.r
