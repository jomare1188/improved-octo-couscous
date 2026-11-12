#!/bin/bash
#$ -q all.q
#$ -cwd
#$ -S /bin/bash
#$ -pe smp 16
#$ -j y
#$ -o logs/$JOB_NAME.$JOB_ID.log
# Genome-wide lncRNA search, step 3 of 5: from the 48 assemblies to lncRNA candidates and a Salmon index.
# Same filters as the chromosome-1 pilot (18_lnc_pilot.sh):
#  1. StringTie merge of the 48 assemblies (guided by the POJ2878 annotation); gffcompare against the annotation.
#  2. Candidates: new transcripts of class u (between genes), x (antisense to a gene exon), i (inside an intron),
#     >= 200 nt. Class u within 1 kb of a gene on the same strand is removed (possible UTR extension).
#  3. Coding potential: CPC2 "noncoding", longest ORF < 100 aa, no DIAMOND blastx hit (1e-5) to the POJ2878 proteome.
#  4. Allele copies: candidates clustered with MMseqs2 (>= 90% identity over >= 50% of the shorter one);
#     one cluster = one lncRNA, the analogue of an allele group.
#  5. Salmon index of the annotated transcripts plus the candidates, so that reads of coding genes do not
#     land on lncRNAs (Salmon 1.10.3; the main analysis used 2.7.0 inside nf-core, so numbers of coding
#     genes from this index are not mixed with the main analysis).
# Run from analysis/:  qsub -N lnccand -hold_jid lncaln code/24_lnc_candidates.sh
set -euo pipefail
export PATH=/Storage/data2/jorge.munoz/conda_envs/lnc/bin:$PATH
SALMON=/home/jorge.munoz/.conda/envs/R2C/bin/salmon
MMSEQS=/home/jorge.munoz/.conda/envs/ORTHOFINDER/bin/mmseqs
DOL=/Storage/data1/jorge.munoz/DOLORES
GDIR=$DOL/nf-rnaseq/run_POJ2878/results_pilot/genome
REF=$GDIR/POJ2878_1.0.genome.fa
GTF=$GDIR/Sh_POJ2878_v1.gtf
PROT=$DOL/annotation/proteins/POJ2878_v1.proteins.uniq.faa
W=/Storage/data2/jorge.munoz/lnc_genome
OUT=$DOL/analysis/results/lncrna_genome
T=16
mkdir -p $OUT; cd $W
echo "host $(hostname) start $(date)"
n=$(ls st/*.gtf | wc -l); [ "$n" -eq 48 ] || { echo "expected 48 assemblies, found $n"; exit 1; }

# ---- 1. merge and compare
ls st/*.gtf > gtf_list.txt
[ -s merged.gtf ] || stringtie --merge -p $T -G $GTF -m 200 -o merged.gtf gtf_list.txt
[ -s cmp.merged.gtf.tmap ] || gffcompare -r $GTF -o cmp merged.gtf
echo "merged $(date)"

# ---- 2. candidates by class
awk -F'\t' 'NR>1 && ($3=="u"||$3=="x"||$3=="i") && $10>=200 {print $5"\t"$3"\t"$6"\t"$10}' cmp.merged.gtf.tmap > cand_class.tsv
cut -f1 cand_class.tsv > cand_ids.txt
awk -F'\t' 'NR==FNR{k[$1]=1; next} {match($9,/transcript_id "[^"]+"/); id=substr($9,RSTART+15,RLENGTH-16); if (id in k) print}' cand_ids.txt merged.gtf > cand.gtf
gffread -w cand.fa -g $REF cand.gtf
awk -F'\t' -v OFS='\t' '$3=="transcript"{match($9,/transcript_id "[^"]+"/); print $1,$4-1,$5,substr($9,RSTART+15,RLENGTH-16),".",$7}' cand.gtf | sort -k1,1 -k2,2n > cand.bed
awk -F'\t' -v OFS='\t' '$3=="transcript"{match($9,/gene_id "[^"]+"/); print $1,$4-1,$5,substr($9,RSTART+9,RLENGTH-10),".",$7}' $GTF | sort -k1,1 -k2,2n > refgenes.bed
bedtools closest -s -d -t first -a cand.bed -b refgenes.bed | awk -F'\t' -v OFS='\t' '{print $4,$NF}' > cand_dist.tsv
echo "candidates $(wc -l < cand_class.tsv) $(date)"

# ---- 3. coding potential and protein homology
[ -s cpc2.txt ] || CPC2.py -i cand.fa -o cpc2 > cpc2.log 2>&1
[ -s prot.dmnd ] || diamond makedb --in $PROT -d prot -p $T > /dev/null 2>&1
[ -s blastx.tsv ] || diamond blastx -q cand.fa -d prot -p $T -e 1e-5 -k 1 --outfmt 6 qseqid sseqid pident length evalue -o blastx.tsv > /dev/null 2>&1
echo "coding filters $(date)"

python3 - "$W" <<'PY'
import sys, csv
W = sys.argv[1]
cls = {r[0]: (r[1], int(r[2]), int(r[3])) for r in csv.reader(open(f"{W}/cand_class.tsv"), delimiter="\t")}
dist = {r[0]: int(r[1]) for r in csv.reader(open(f"{W}/cand_dist.tsv"), delimiter="\t")}
cpc = {r["#ID"]: (r["label"], int(r["peptide_length"])) for r in csv.DictReader(open(f"{W}/cpc2.txt"), delimiter="\t")}
hit = {r[0] for r in csv.reader(open(f"{W}/blastx.tsv"), delimiter="\t")}
with open(f"{W}/candidates.tsv", "w") as fh, open(f"{W}/lnc_ids.txt", "w") as ok:
    fh.write("transcript\tclass_code\texons\tlength\tcpc2\torf_aa\tprotein_hit\tdist_same_strand_gene\tlncRNA\tfilters_failed\n")
    for t, (c, ne, ln) in cls.items():
        lab, pep = cpc.get(t, ("NA", -1)); d = dist.get(t, -1); fl = []
        if lab != "noncoding": fl.append("CPC2_coding")
        if pep >= 100: fl.append("ORF>=100aa")
        if t in hit: fl.append("protein_hit")
        if c == "u" and 0 <= d <= 1000: fl.append("near_gene_same_strand")
        fh.write(f"{t}\t{c}\t{ne}\t{ln}\t{lab}\t{pep}\t{int(t in hit)}\t{d}\t{int(not fl)}\t{';'.join(fl) or 'pass'}\n")
        if not fl: ok.write(t + "\n")
PY
awk 'NR==FNR{k[">"$1]=1; next} /^>/{p=($1 in k)} p' lnc_ids.txt cand.fa > lnc.fa
awk -F'\t' 'NR==FNR{k[$1]=1; next} {match($9,/transcript_id "[^"]+"/); id=substr($9,RSTART+15,RLENGTH-16); if (id in k) print}' lnc_ids.txt cand.gtf > lnc.gtf
echo "lncRNA candidates $(grep -c '>' lnc.fa) $(date)"

# ---- 4. allele copies -> clusters
if [ ! -s clu_cluster.tsv ]; then
  rm -rf mmtmp
  $MMSEQS easy-cluster lnc.fa clu mmtmp --min-seq-id 0.9 -c 0.5 --cov-mode 1 --cluster-mode 2 --threads $T > mmseqs.log 2>&1
  rm -rf mmtmp
fi
awk -F'\t' -v OFS='\t' 'BEGIN{print "transcript_id","group_id"} {print $2,"LNC:"$1}' clu_cluster.tsv > lnc_tx2group.tsv
echo "clusters $(cut -f1 clu_cluster.tsv | sort -u | wc -l) $(date)"

# ---- 5. Salmon index: annotated transcripts + lncRNA candidates
if [ ! -s salmon_idx/versionInfo.json ]; then
  cat $GDIR/genome.transcripts.fa lnc.fa > tx_plus_lnc.fa
  $SALMON index -t tx_plus_lnc.fa -i salmon_idx -k 31 -p $T > salmon_index.log 2>&1
fi
mkdir -p quant
cp candidates.tsv lnc_tx2group.tsv cmp.stats $OUT/
gzip -c lnc.gtf > $OUT/lncRNA_candidates.gtf.gz
for L in $(cat libs.txt); do echo -e "$L\t$(tail -1 bam/$L.hisat2.log)"; done > $OUT/alignment_rates.tsv
echo "end $(date)"
