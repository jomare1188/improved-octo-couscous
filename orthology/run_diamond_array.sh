#!/bin/bash
#$ -N of_diamond
#$ -q all.q
#$ -cwd
#$ -S /bin/bash
#$ -pe smp 8
#$ -j y
#$ -o logs/diamond.$JOB_ID.$TASK_ID.log
# One OrthoFinder DIAMOND search per task (line of tasks.txt), same settings as OrthoFinder 2.5.4
# but 8 threads. Skips tasks already done (resumable).
#   qsub -t 1-$(wc -l < tasks.txt) -tc 10 run_diamond_array.sh
set -euo pipefail
source /Storage/progs/miniconda3/etc/profile.d/conda.sh
conda activate ORTHOFINDER
WD=/Storage/data2/jorge.munoz/DOLORES_orthofinder/run/Results_Oct29/WorkingDirectory
read q db out < <(sed -n "${SGE_TASK_ID}p" tasks.txt)
[ -s $out.gz ] && [ -e $out.done ] && { echo "already done: $out"; exit 0; }
rm -f $out $out.gz
echo "$(date) $(hostname) $q vs Species$db"
diamond blastp -d $WD/diamondDBSpecies$db -q $q -o $out --more-sensitive -p 8 --quiet -e 0.001 --compress 1 \
    --tmpdir ${TMPDIR:-/tmp}
touch $out.done
echo "$(date) done"
