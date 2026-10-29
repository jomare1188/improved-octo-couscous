#!/bin/bash
#$ -N emapper_search
#$ -q all.q
#$ -cwd
#$ -S /bin/bash
#$ -pe smp 8
#$ -j y
#$ -o logs/search.$JOB_ID.$TASK_ID.log
#$ -l h=!telura*
# eggNOG-mapper step 1: DIAMOND search, one array task per protein chunk (no annotation).
# Submit:  N=$(ls ../proteins/chunks/chunk_*.faa | wc -l)
#          qsub -t 1-$N -tc 6 run_emapper_search.sh            (at most 6 chunks at once)
set -euo pipefail
source /etc/profile.d/modules.sh
module load miniconda3 eggnog-mapper/2.1.13
CHUNKDIR=${CHUNKDIR:-/Storage/data1/jorge.munoz/DOLORES/annotation/proteins/chunks}
c=$(printf "chunk_%03d" $SGE_TASK_ID)
[ -s out/$c.emapper.seed_orthologs ] && { echo "$c already done"; exit 0; }
emapper.py -m diamond --itype proteins --no_annot \
    -i $CHUNKDIR/$c.faa -o $c --output_dir out \
    --temp_dir tmp --cpu $NSLOTS --override
echo "$c done"
