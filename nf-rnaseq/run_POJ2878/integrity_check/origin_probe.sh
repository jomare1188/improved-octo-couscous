#!/bin/bash
#$ -N origin_probe
#$ -q all.q
#$ -cwd
#$ -S /bin/bash
#$ -pe smp 2
#$ -j y
#$ -o origin_probe.$JOB_ID.log
#$ -l h=!telura*
# Read-only. Phase 1: sample ~2000 reads from each suspect library.
# Phase 2: scan every library's R1 for those names and exact sequences, and record
# CRLF line endings, run/flowcell IDs and the top barcodes of each file.
set -euo pipefail
RAW=/Storage/data1/jorge.munoz/DOLORES/RAWDATA
SHEET=../samplesheet_POJ2878.csv
SUSPECT="159_S5 222_S6 255_S8 9_S15 13_S17 15_S18 D54_S5 D69_S47 79_S2 102_S3 101_S20 D158_S48"

# target list: library id <tab> R1 path (original 159_S5, not the repaired copy)
tail -n +2 $SHEET | tr -d '\r' | while IFS=, read -r s f1 f2 st; do
  lib=$(basename "$f1" | sed -E 's/(_repaired)?_R1_001.fastq.gz//'); echo -e "$lib\t$f1"; done \
  | sed "s#\t.*repair_159_S5/.*#\t$RAW/159_S5_R1_001.fastq.gz#" > targets.tsv

: > probes.tsv
for lib in $SUSPECT; do
  f=$(awk -F'\t' -v l=$lib '$1==l{print $2}' targets.tsv)
  n=$(zcat "$f" | awk 'END{print NR/4}')
  zcat "$f" | awk -v l=$lib -v step=$(( n/2000 )) '{sub(/\r$/,"")} NR%4==1{nm=$1} NR%4==2 && int((NR-2)/4)%step==0 && length($0)>=50 && $0!~/N/{print l"\t"nm"\t"$0}' >> probes.tsv
done
echo "probes: $(wc -l < probes.tsv)"

echo -e "target\tprobe_lib\tname_hits\tseq_hits" > hits.tsv
echo -e "target\treads\tcrlf_reads\truns\ttop_barcodes" > fileinfo.tsv
while IFS=$'\t' read -r t f; do
  zcat "$f" | awk -v t=$t -v P=probes.tsv -v H=hits.tsv -v I=fileinfo.tsv -F'\t' '
  BEGIN{while((getline l < P)>0){split(l,a,"\t"); nm[a[2]]=nm[a[2]]" "a[1]; sq[a[3]]=sq[a[3]]" "a[1]; libs[a[1]]}; FS=" "}
  {cr=sub(/\r$/,"")}
  NR%4==1{reads++; crl+=cr; split($1,h,":"); rc[h[1]":"h[2]":"h[3]":"h[4]]++; split($2,b,":"); bc[b[4]]++
          if($1 in nm){k=split(nm[$1],L," "); for(i=1;i<=k;i++) hn[L[i]]++}}
  NR%4==2{if($0 in sq){k=split(sq[$0],L," "); for(i=1;i<=k;i++) hs[L[i]]++}}
  END{for(p in libs) print t"\t"p"\t"(hn[p]+0)"\t"(hs[p]+0) >> H
      r=""; for(k in rc) r=r k"="rc[k]";"
      top=""; for(j=1;j<=3;j++){m=0; mk=""; for(k in bc) if(bc[k]>m){m=bc[k]; mk=k}; if(mk=="") break; top=top mk"="sprintf("%.3f",m/reads)";"; delete bc[mk]}
      print t"\t"reads"\t"crl"\t"r"\t"top >> I}'
  echo "$t done"
done < targets.tsv
echo finished
