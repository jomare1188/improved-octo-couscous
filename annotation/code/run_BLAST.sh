#!/bin/bash
#$ -q all.q
#$ -cwd
#$ -V
#$ -pe smp 30

module load miniconda3
conda activate DIAMOND
# yes, im using the blast used to run diamond.
sugarcane_pantranscriptome_coding=../../quant_test/protein_coding_transcripts/data/transcripts_in_orthogroups.ids.fasta
sorghum_db=/Storage/data1/jorge.munoz/DOLORES/refs/sorghum_db
maize_db=/Storage/data1/jorge.munoz/DOLORES/refs/maize_db
arabidopsis_db=/Storage/data1/jorge.munoz/DOLORES/refs/arabidopsis_db

mkdir -p ../results/


echo @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ running sorghum @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

blastn -db ${sorghum_db} -query ${sugarcane_pantranscriptome_coding} -out ../results/sorghum_blast.tbl -outfmt "6 qseqid sseqid evalue bitscore score" -max_target_seqs 1 -max_hsps 1 -num_threads $NSLOTS 

echo @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ running maize @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
blastn -db ${maize_db} -query ${sugarcane_pantranscriptome_coding} -out ../results/maize_blast.tbl -outfmt "6 qseqid  sseqid evalue bitscore score" -max_target_seqs 1 -max_hsps 1 -num_threads $NSLOTS

echo @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ running arabidopsis  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
blastn -db ${arabidopsis_db} -query ${sugarcane_pantranscriptome_coding} -out ../results/arabidopsis_blast.tbl -outfmt "6 qseqid sseqid evalue bitscore score" -max_target_seqs 1 -max_hsps 1 -num_threads $NSLOTS
# make one table relating all to sugarcane
# lets do it with R! :D
module load R/4.2.1
Rscript merge_tables.r
## 
