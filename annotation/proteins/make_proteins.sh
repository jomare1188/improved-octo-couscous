#!/bin/bash
#$ -N make_proteins
#$ -q all.q
#$ -cwd
#$ -S /bin/bash
#$ -pe smp 2
#$ -j y
#$ -o make_proteins.$JOB_ID.log
#$ -l h=!telura*
# Build the POJ2878 proteome from the GFF + genome, clean it, collapse identical sequences,
# and split the unique set into chunks of 20,000 for eggNOG.
set -euo pipefail
source /etc/profile.d/modules.sh; module load singularity-ce/4.3.3
D=/Storage/data1/jorge.munoz/DOLORES
GENOME=$D/nf-rnaseq/run_POJ2878/results_pilot/genome/POJ2878_1.0.genome.fa
GFF=$D/refs/POJ2878/Sh_POJ2878_v1.gff3.gz
IMG=$D/nf-rnaseq/singularity_cache/depot.galaxyproject.org-singularity-gffread-0.12.7--hdcf5f25_4.img
CHUNK=20000

zcat $GFF > tmp.gff3
singularity exec -B $D $IMG gffread -y raw.faa -g $GENOME tmp.gff3
rm -f tmp.gff3

# one line per sequence, strip terminal stop, mask internal stops/dots with X
awk '/^>/{if(id)print id"\t"s; id=substr($1,2); s=""; next}{s=s $0} END{print id"\t"s}' raw.faa \
 | awk -F'\t' '{sub(/[*.]+$/,"",$2); n=gsub(/[*.]/,"X",$2); if(n>0) bad++; print $1"\t"$2 > "all.tsv"} END{print "proteins with internal stops: " bad+0}'
awk -F'\t' '{print ">"$1"\n"$2}' all.tsv > POJ2878_v1.proteins.faa

# identical sequences -> first ID is the representative
awk -F'\t' '!($2 in rep){rep[$2]=$1; print ">"$1"\n"$2 > "POJ2878_v1.proteins.uniq.faa"} {print rep[$2]"\t"$1 > "uniq_map.tsv"}' all.tsv

rm -f chunks/chunk_*.faa
awk -v n=$CHUNK '/^>/{if(c%n==0){f=sprintf("chunks/chunk_%03d.faa", c/n+1)}; c++} {print > f}' POJ2878_v1.proteins.uniq.faa

echo "all proteins:    $(grep -c '>' POJ2878_v1.proteins.faa)"
echo "unique proteins: $(grep -c '>' POJ2878_v1.proteins.uniq.faa)"
echo "empty sequences: $(awk -F'\t' '$2==""' all.tsv | wc -l)"
echo "chunks:          $(ls chunks | wc -l)"
echo "stops left:      $(grep -v '>' POJ2878_v1.proteins.uniq.faa | grep -c '[*.]' || true)"
rm -f raw.faa all.tsv
md5sum POJ2878_v1.proteins.faa POJ2878_v1.proteins.uniq.faa uniq_map.tsv > md5sum.txt
echo done
