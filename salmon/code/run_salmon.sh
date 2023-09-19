#!/bin/bash
#$ -q all.q
#$ -cwd
#$ -pe smp 1

qsub run_index-salmon.sh
qsub run_quant-salmon.sh
