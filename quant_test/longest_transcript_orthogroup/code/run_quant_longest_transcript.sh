#!/bin/bash
#SBATCH --export=ALL
#SBATCH --partition=long
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=60
#SBATCH --mem=150gb
#SBATCH --job-name=quant_longest_transcript

source ~/.bashrc
conda activate salmon
## Make index
threads=60
transcriptome=../data/longest_trans.fa
index=../results/index/

## Make transcriptome index
salmon index -t ${transcriptome} -p ${threads} -i ${index}

## Quantification with salmon
rm -f salmon_commands
for id in $(cat ../../bbduk_seqs/ids_samples)
do
	echo salmon quant -i ${index} -l A -1 ../../bbduk_seqs/${id}_R1.bbduk.fastq -2 ../../bbduk_seqs/${id}_R2.bbduk.fastq -p 10 -o ../results/${id} >> salmon_commands
done
cat salmon_commands | parallel -j ${threads}

