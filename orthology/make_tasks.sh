#!/bin/bash
# Split the POJ2878 query (Species2) into chunks and write one DIAMOND task per line to tasks.txt.
# Each line: query_fasta db_index out_file
set -euo pipefail
WD=/Storage/data2/jorge.munoz/DOLORES_orthofinder/run/Results_Oct29/WorkingDirectory
NCHUNK=${NCHUNK:-20}
mkdir -p $WD/chunks
awk -v n=$NCHUNK -v d=$WD/chunks 'BEGIN{c=0} /^>/{f=sprintf("%s/Species2.c%02d.fa",d,(c++)%n)} {print > f}' $WD/Species2.fa
: > tasks.txt
for i in 0 1 2 3 4; do for j in 0 1 2 3 4; do
  if [ $i -eq 2 ]; then
    for k in $(seq -f %02g 0 $((NCHUNK-1))); do echo "$WD/chunks/Species2.c$k.fa $j $WD/chunks/Blast2_$j.c$k.txt" >> tasks.txt; done
  else
    echo "$WD/Species$i.fa $j $WD/Blast${i}_$j.txt" >> tasks.txt
  fi
done; done
wc -l tasks.txt
