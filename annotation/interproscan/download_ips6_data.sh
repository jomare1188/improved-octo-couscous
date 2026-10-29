#!/bin/bash
#$ -N ips6_download
#$ -q all.q
#$ -cwd
#$ -S /bin/bash
#$ -pe smp 2
#$ -j y
#$ -o logs/download.$JOB_ID.log
#$ -l h=!telura*
# Download the InterProScan 6.0 member databases that the pipeline's 1 h download tasks could not
# finish (InterPro 110.0 versions, from ips6_data/interpro/110.0/databases.json).
# Same steps as the pipeline's DOWNLOAD process (curl, md5sum -c, tar into the datadir) but with
# no time limit, resumable curl, and no slow `chmod -R` over NFS. A partial target dir is removed
# before extracting. Already complete databases are skipped.
set -euo pipefail
DATA=/Storage/data1/jorge.munoz/DOLORES/annotation/interproscan/ips6_data
URL=https://ftp.ebi.ac.uk/pub/software/unix/iprscan/6/6.0
mkdir -p $DATA/.dl && cd $DATA/.dl
for nv in panther:19.0 cath:4.3.0 superfamily:1.75 ncbifam:20.0 prosite:2026_01 hamap:2026_01 pirsr:2025_05; do
  n=${nv%%:*}; v=${nv#*:}; a=$n-$v.tar.gz
  [ -f $DATA/$n/$v/.complete ] && { echo "$n $v already complete"; continue; }
  echo "$(date '+%F %T') $n $v: download"
  # resumable; EBI sometimes drops the connection (curl 56), so retry the whole call too
  ok=0
  for try in 1 2 3 4 5 6 7 8 9 10; do
    curl -fL --retry 10 --retry-all-errors --retry-delay 30 -C - -o $a $URL/$n/$a && { ok=1; break; }
    echo "$(date '+%F %T') $n $v: curl failed (try $try), resuming in 60 s"; sleep 60
  done
  [ $ok = 1 ] || { echo "ERROR: $a download failed"; exit 1; }
  curl -fL --retry 10 --retry-all-errors -o $a.md5 $URL/$n/$a.md5
  md5sum -c $a.md5
  echo "$(date '+%F %T') $n $v: extract"
  rm -rf $DATA/$n/$v
  tar -C $DATA -zxf $a
  [ -d $DATA/$n/$v ] || { echo "ERROR: $a did not create $n/$v"; exit 1; }
  touch $DATA/$n/$v/.complete
  rm -f $a $a.md5
  echo "$(date '+%F %T') $n $v: done ($(du -sh $DATA/$n/$v | cut -f1))"
done
echo "all downloads done"
