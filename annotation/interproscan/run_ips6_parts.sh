#!/bin/bash
#$ -N ips6_parts
#$ -q all.q
#$ -cwd
#$ -S /bin/bash
#$ -pe smp 2
#$ -j y
#$ -o logs/ips6_parts.$JOB_ID.log
# Full InterProScan 6 run on the 307,787 unique POJ2878 proteins, in 4 parts run one after another.
# Why parts: the work dir needs ~1.4 GB per 1,000 proteins (~430 GB for all), and /Storage/data1 is
# 99% full. Each part (~77k proteins, ~110 GB of work) gets its own work dir, which is deleted once
# that part's TSV is written. Within a part, each 3,000-protein batch x member database is its own
# SGE job (sge_ips6.config), so the work spreads over the nodes.
#   cd /Storage/data1/jorge.munoz/DOLORES/annotation/interproscan
#   qsub run_ips6_parts.sh          # resubmit to resume: finished parts are skipped,
#                                   # an unfinished part resumes from its work dir (-resume)
# Output: out/part_N.faa.tsv (+ .gff3). JSON is skipped (~120 GB; the TSV has GO and pathways).
set -euo pipefail
source /etc/profile.d/modules.sh
module load miniconda3 singularity-ce/4.3.3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate nextflow-26.04.7

BASE=/Storage/data1/jorge.munoz/DOLORES/annotation/interproscan
export NXF_OPTS='-Xms1g -Xmx8g'
export NXF_SINGULARITY_CACHEDIR=$BASE/singularity_cache
mkdir -p $BASE/out $BASE/logs
cd $BASE

for i in 1 2 3 4; do
  part=part_$i
  if [ -f out/$part.done ]; then echo "$(date '+%F %T') $part already done"; continue; fi
  echo "$(date '+%F %T') $part: start ($(grep -c '>' parts/$part.faa) proteins)"
  # Nextflow sometimes aborts with "KryoException: Buffer underflow": 4 member DBs split the same
  # batch file at once and race on its .chunks cache file. The file is complete afterwards, so a
  # -resume gets past it; retry the part up to 5 times.
  ok=0
  for try in 1 2 3 4 5; do
  echo "$(date '+%F %T') $part: nextflow try $try"
  if nextflow -log logs/nextflow.$part.try$try.log run ebi-pf-team/interproscan6 -r 6.0.3.0 -profile singularity \
      -c $BASE/sge_ips6.config \
      --input $BASE/parts/$part.faa \
      --datadir $BASE/ips6_data --interpro 110.0 \
      --goterms --pathways --formats tsv,gff3 \
      --batchSize 3000 --subBatchSize 3000 --cpus 4 \
      --outdir $BASE/out \
      -work-dir $BASE/work_$part \
      -with-trace $BASE/logs/trace.$part.try$try.$JOB_ID.txt \
      -resume
  then ok=1; break; fi
  echo "$(date '+%F %T') $part: nextflow failed (try $try)"; sleep 120
  done
  [ $ok = 1 ] || { echo "ERROR: $part failed 5 times"; exit 1; }
  [ -s out/$part.faa.tsv ] || { echo "ERROR: out/$part.faa.tsv missing"; exit 1; }
  echo "$(date '+%F %T') $part: $(cut -f1 out/$part.faa.tsv | sort -u | wc -l) proteins with matches; removing work_$part"
  rm -rf $BASE/work_$part
  touch out/$part.done
done
echo "$(date '+%F %T') all parts done"
