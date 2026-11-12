#!/bin/bash
#$ -q all.q
#$ -cwd
#$ -S /bin/bash
#$ -pe smp 4
#$ -j y
#$ -o logs/$JOB_NAME.$JOB_ID.log
# What kind of libraries are these? Poly-A selected or total RNA (rRNA-depleted)?
# Salmon already says they are stranded (ISR: read 1 antisense, dUTP-type) and paired 2 x 100 bp.
# Here: map 1M read pairs of six libraries to two chromosomes (Chr01g1, Chr02g1) with HISAT2 and count where
# read 1 lands: exon, intron (inside a gene, sense strand, no exon), antisense to a gene, or between genes.
# Poly-A libraries carry few intron reads (mature mRNA, typically < 5-10% of genic reads); total-RNA
# libraries carry many (pre-mRNA, typically 20-50%). This decides what a lncRNA annotation could find.
# Run from analysis/:  qsub -N libtype code/17_library_type.sh
set -euo pipefail
# tools by absolute path: the cluster modules do not load in batch shells, and Samtools/1.21 lacks libcrypto on some nodes
export PATH=/Storage/progs/hisat2-2.2.1:$PATH
samtools() { /home/jorge.munoz/.conda/envs/minimap2/bin/samtools "$@"; }   # samtools 1.9
bedtools() { /home/jorge.munoz/.conda/envs/quast/bin/bedtools "$@"; }      # bedtools 2.31.1
DOL=/Storage/data1/jorge.munoz/DOLORES
GENOME=$DOL/nf-rnaseq/run_POJ2878/results_pilot/genome/POJ2878_1.0.genome.fa
GFF=$DOL/refs/POJ2878/Sh_POJ2878_v1.gff3.gz
TMP=/Storage/data2/jorge.munoz/libtype_tmp
OUT=$DOL/analysis/results/qc/library_type
CHR="Chr01g1 Chr02g1"
NPAIRS=1000000
LIBS="168_S25 250_S2 Dolores_257_S5 79_S6 F28C1_S9 FIDH1800_S30"
mkdir -p $TMP $OUT
echo "host $(hostname) start $(date)"

# reference: two chromosomes, HISAT2 index, gene and exon BED (6 columns, strand kept)
if [ ! -s $TMP/ref.1.ht2 ]; then
  samtools faidx $GENOME $CHR > $TMP/ref.fa
  hisat2-build -p 4 $TMP/ref.fa $TMP/ref > $TMP/build.log 2>&1
fi
zcat $GFF | awk -F'\t' -v OFS='\t' 'BEGIN{split("'"$CHR"'",c," "); for(i in c) k[c[i]]=1}
  ($1 in k) && $3=="gene" {print $1,$4-1,$5,"g",".",$7 > "'$TMP'/genes.bed"}
  ($1 in k) && $3=="exon" {print $1,$4-1,$5,"e",".",$7 > "'$TMP'/exons.bed"}'

echo -e "library\tpairs_used\tread1_mapped_unique\tpct_mapped\tpct_spliced\texonic\tintronic\tantisense\tintergenic\tpct_intronic_of_genic\tpct_intergenic" > $OUT/read_locations.tsv
for L in $LIBS; do
  R1=$DOL/RAWDATA/${L}_R1_001.fastq.gz; R2=$DOL/RAWDATA/${L}_R2_001.fastq.gz
  zcat $R1 | head -n $((NPAIRS*4)) > $TMP/$L.1.fq || true
  zcat $R2 | head -n $((NPAIRS*4)) > $TMP/$L.2.fq || true
  hisat2 -p 4 --rna-strandness RF --no-unal -x $TMP/ref -1 $TMP/$L.1.fq -2 $TMP/$L.2.fq 2> $TMP/$L.hisat2.log |
    samtools view -b -f 64 -F 0x904 -q 60 - > $TMP/$L.r1.bam      # read 1, primary, unique (MAPQ 60)
  n=$(samtools view -c $TMP/$L.r1.bam)
  sp=$(samtools view $TMP/$L.r1.bam | awk '$6 ~ /N/' | wc -l)
  bedtools bamtobed -split -i $TMP/$L.r1.bam > $TMP/$L.r1.bed
  # read 1 is antisense to the RNA (RF), so "sense to a gene" = opposite strand: -S
  bedtools intersect -S -u -a $TMP/$L.r1.bed -b $TMP/exons.bed | cut -f4 | sort -u > $TMP/$L.exon
  bedtools intersect -S -u -a $TMP/$L.r1.bed -b $TMP/genes.bed | cut -f4 | sort -u > $TMP/$L.gene
  bedtools intersect -s -u -a $TMP/$L.r1.bed -b $TMP/genes.bed | cut -f4 | sort -u > $TMP/$L.anti
  cut -f4 $TMP/$L.r1.bed | sort -u > $TMP/$L.all
  ex=$(wc -l < $TMP/$L.exon)
  intr=$(comm -23 $TMP/$L.gene $TMP/$L.exon | wc -l)
  anti=$(comm -23 $TMP/$L.anti $TMP/$L.gene | wc -l)
  tot=$(wc -l < $TMP/$L.all)
  inter=$((tot - ex - intr - anti))
  awk -v L=$L -v N=$NPAIRS -v n=$n -v sp=$sp -v ex=$ex -v intr=$intr -v anti=$anti -v inter=$inter -v tot=$tot 'BEGIN{
    printf "%s\t%d\t%d\t%.1f\t%.1f\t%d\t%d\t%d\t%d\t%.1f\t%.1f\n", L, N, n, 100*n/N, 100*sp/n, ex, intr, anti, inter,
      100*intr/(ex+intr), 100*inter/tot}' >> $OUT/read_locations.tsv
  rm -f $TMP/$L.1.fq $TMP/$L.2.fq
  echo "$L done $(date)"
done
cp $TMP/*.hisat2.log $OUT/ 2>/dev/null || true
column -t $OUT/read_locations.tsv
echo "end $(date)"
