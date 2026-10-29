#!/bin/bash
#$ -N emapper_annot
#$ -q all.q
#$ -cwd
#$ -S /bin/bash
#$ -pe smp 8
#$ -j y
#$ -o logs/annot.$JOB_ID.log
#$ -l h=!telura*,virtual_free=64G
# eggNOG-mapper step 2: annotate all search hits at once with eggnog.db loaded in RAM (~45 GB).
# Submit after the search array:  qsub -hold_jid emapper_search run_emapper_annot.sh
set -euo pipefail
source /etc/profile.d/modules.sh
module load miniconda3 eggnog-mapper/2.1.13
N=$(ls ${CHUNKDIR:-/Storage/data1/jorge.munoz/DOLORES/annotation/proteins/chunks}/chunk_*.faa | wc -l)
D=$(ls out/chunk_*.emapper.seed_orthologs | wc -l)
[ "$N" = "$D" ] || { echo "only $D of $N chunks searched; resubmit the missing array tasks"; exit 1; }
{ head -n 5 out/chunk_001.emapper.seed_orthologs | grep '^#' ; cat out/chunk_*.emapper.seed_orthologs | grep -v '^#'; } > out/POJ2878.seed_orthologs
emapper.py -m no_search --annotate_hits_table out/POJ2878.seed_orthologs \
    -o POJ2878 --output_dir out --dbmem --cpu $NSLOTS --override
echo done
