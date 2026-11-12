#!/bin/bash
#$ -q all.q
#$ -cwd
#$ -S /bin/bash
#$ -pe smp 8
#$ -j y
#$ -o logs/$JOB_NAME.$JOB_ID.$TASK_ID.log
# Genome-wide lncRNA search, step 2 of 3: align one library to the whole genome and assemble its transcripts.
# Libraries are poly-A, stranded (RF), 2 x 100 bp (results/qc/library_type/). The BAM is deleted after the
# assembly; expression is measured later with Salmon. Array job, one task per library:
#   qsub -N lncaln -hold_jid lncidx -t 1-48 -tc 8 code/22_lnc_align.sh        (from analysis/)
set -euo pipefail
export PATH=/Storage/data2/jorge.munoz/conda_envs/lnc/bin:$PATH
DOL=/Storage/data1/jorge.munoz/DOLORES
GTF=$DOL/nf-rnaseq/run_POJ2878/results_pilot/genome/Sh_POJ2878_v1.gtf
W=/Storage/data2/jorge.munoz/lnc_genome
cd $W
L=$(sed -n "${SGE_TASK_ID}p" libs.txt)
echo "host $(hostname) library $L start $(date)"
if [ ! -s st/$L.gtf ]; then
  hisat2 -p 8 --dta --rna-strandness RF --no-unal -x idx \
    -1 $DOL/RAWDATA/${L}_R1_001.fastq.gz -2 $DOL/RAWDATA/${L}_R2_001.fastq.gz 2> bam/$L.hisat2.log |
    samtools sort -@ 4 -m 2G -T bam/$L.tmp -o bam/$L.bam -
  stringtie bam/$L.bam -p 8 --rf -G $GTF -o st/$L.gtf.tmp -l $L
  mv st/$L.gtf.tmp st/$L.gtf
  rm -f bam/$L.bam
fi
tail -1 bam/$L.hisat2.log
echo "end $(date)"
