#!/bin/bash
#SBATCH --export=ALL
#SBATCH --partition=long
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20gb
#SBATCH --job-name=PANNZER2

GENOTYPE="LongestsTranscript"
GENOTYPE_FILE="../data/reformated_fasta/no_empty_lines.fasta"
PANZZER="../SANSPANZ.3/runsanspanz.py"
SPECIE="Saccharum officinarum x Saccharum spontaneum"

#/home/dmpachon/miniconda3/condabin/conda activate panzzer

source ~/.bashrc
conda activate panzzer
python ${PANZZER} -R -o ",../results/${GENOTYPE}_DE.out,../results/${GENOTYPE}_GO.out,../results/${GENOTYPE}_anno.out" -s "${SPECIE}" < ${GENOTYPE_FILE}
