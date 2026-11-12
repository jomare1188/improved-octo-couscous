#!/bin/bash
#$ -N of_orth
#$ -q all.q
#$ -cwd
#$ -S /bin/bash
#$ -pe smp 16
#$ -j y
#$ -o logs/of_orth.$JOB_ID.log
# OrthoFinder 2.5.4 (conda env) built the orthogroups but crashed in the species-tree step (Python 3.11 rejects
# open(..., 'rU')). Resume from those orthogroups with the cluster's OrthoFinder 2.5.5, which has the fix.
#   qsub run_orthofinder_orthologues.sh
set -euo pipefail
source /etc/profile.d/modules.sh 2>/dev/null || true
module load OrthoFinder/2.5.5
WD=/Storage/data2/jorge.munoz/DOLORES_orthofinder/run/Results_Oct29/WorkingDirectory
orthofinder -fg $WD/OrthoFinder/Results_Oct29 -t 16 -a 16
echo "$(date) done"
find $WD -maxdepth 4 -type d -name Orthologues_POJ2878
