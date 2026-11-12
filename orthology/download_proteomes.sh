#!/bin/bash
#$ -N dl_prot
#$ -q all.q
#$ -cwd
#$ -S /bin/bash
#$ -pe smp 1
#$ -j y
#$ -o download_proteomes.$JOB_ID.log
# Download reference proteomes (Ensembl Plants, current release) for OrthoFinder.
#   cd /Storage/data1/jorge.munoz/DOLORES/orthology && qsub download_proteomes.sh
set -euo pipefail
B=https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/fasta
cd proteomes
for sp in sorghum_bicolor zea_mays oryza_sativa arabidopsis_thaliana; do
  f=$(curl -s $B/$sp/pep/ | grep -o 'href="[^"]*pep.all.fa.gz"' | sed 's/href="//;s/"//' | head -1)
  for try in 1 2 3 4 5; do
    curl -sS --retry 5 --retry-all-errors -fLo $f $B/$sp/pep/$f && gzip -t $f && break
    sleep 60
  done
  echo "$sp $f $(md5sum < $f | cut -c1-32)"
done
curl -s https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/ | head -c 0
date > DOWNLOADED
