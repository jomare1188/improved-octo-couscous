#!/bin/bash
#$ -N name_screen
#$ -q all.q
#$ -cwd
#$ -S /bin/bash
#$ -pe smp 4
#$ -j y
#$ -o name_screen.$JOB_ID.log
# Read-only screen of every library in the samplesheet: reads, unique read names (R1, R2),
# whether R1/R2 names pair up in order, and the first-read barcode.
SHEET=/Storage/data1/jorge.munoz/DOLORES/nf-rnaseq/run_POJ2878/samplesheet_POJ2878.csv
OUT=name_screen.tsv
echo -e "sample\treads_R1\tunique_R1\treads_R2\tunique_R2\tpairs_in_order\tbarcode_R1" > $OUT
tail -n +2 $SHEET | while IFS=, read -r s f1 f2 st; do
  zcat "$f1" | awk 'NR%4==1{print $1}' > r1.txt
  zcat "$f2" | awk 'NR%4==1{print $1}' > r2.txt
  n1=$(wc -l < r1.txt); n2=$(wc -l < r2.txt)
  u1=$(sort -u -S 6G --parallel=4 r1.txt | wc -l); u2=$(sort -u -S 6G --parallel=4 r2.txt | wc -l)
  cmp -s r1.txt r2.txt && inorder=YES || inorder=NO
  bc=$(zcat "$f1" | head -1 | awk '{print $2}' | cut -d: -f4)
  echo -e "$s\t$n1\t$u1\t$n2\t$u2\t$inorder\t$bc" >> $OUT
done
rm -f r1.txt r2.txt
echo done
