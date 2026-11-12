#!/bin/bash
#$ -N of_finish
#$ -q all.q
#$ -cwd
#$ -S /bin/bash
#$ -pe smp 16
#$ -j y
#$ -o logs/of_finish.$JOB_ID.log
# After the DIAMOND array: join the POJ2878 chunk results and run OrthoFinder from the precomputed hits.
#   qsub -hold_jid of_diamond run_orthofinder_finish.sh
set -euo pipefail
source /Storage/progs/miniconda3/etc/profile.d/conda.sh
conda activate ORTHOFINDER
WD=/Storage/data2/jorge.munoz/DOLORES_orthofinder/run/Results_Oct29/WorkingDirectory
n=0
while read q db out; do [ -e $out.done ] || { echo "missing: $out"; n=$((n+1)); }; done < tasks.txt
[ $n -eq 0 ] || { echo "$n DIAMOND tasks not done; resubmit run_diamond_array.sh (it skips finished tasks)"; exit 1; }
for j in 0 1 2 3 4; do
  cat $WD/chunks/Blast2_$j.c*.txt.gz > $WD/Blast2_$j.txt.gz
done
echo "$(date) hits joined; running OrthoFinder on $(hostname)"
orthofinder -b $WD -t 16 -a 16
echo "$(date) done"
ls -d $WD/OrthoFinder/Results_*
