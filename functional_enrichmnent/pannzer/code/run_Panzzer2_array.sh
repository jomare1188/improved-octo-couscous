#!/bin/bash
#SBATCH --export=ALL
#SBATCH --partition=long
#SBATCH --array=1-50
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8gb
#SBATCH --job-name=PANNZER2

GENOTYPE_FILE="../data/transdecoder_pep/parts/all.part-${SLURM_ARRAY_TASK_ID}.pep"

PANZZER=../SANSPANZ.3/runsanspanz.py
SPECIE="Saccharum officinarum x Saccharum spontaneum"

#/home/dmpachon/miniconda3/condabin/conda activate panzzer

source ~/.bashrc
conda activate panzzer
python ${PANZZER} -R -o ",../results/part_${SLURM_ARRAY_TASK_ID}_DE.out,../results/part_${SLURM_ARRAY_TASK_ID}_GO.out,../results/part_${SLURM_ARRAY_TASK_ID}_anno.out" -s "${SPECIE}" < ${GENOTYPE_FILE}



