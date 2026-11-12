#!/bin/bash
#$ -q all.q
#$ -cwd
#$ -S /bin/bash
#$ -pe smp 8
#$ -j y
#$ -o logs/$JOB_NAME.$JOB_ID.log
# lncRNA pilot: how many lncRNA candidates do our libraries give, and how reliable are they?
# Reference: all 13 copies of chromosome 1 (Chr01g1-g13, 1.49 Gb, 14% of the genome), so reads from each
# allele can map to their own copy. Eight libraries covering both cultivars, leaves and buds.
# Libraries are poly-A, stranded (RF), 2 x 100 bp (results/qc/library_type/).
#  1. HISAT2 index of Chr01; align (--dta, strand RF); sort.
#  2. StringTie per library, guided by the POJ2878 annotation; merge; gffcompare against the annotation.
#  3. Candidates: new transcripts of class u (between genes), x (antisense to a gene exon), i (inside an intron),
#     >= 200 nt. Class u within 1 kb of a gene on the same strand is flagged as a possible UTR extension
#     (about half of the annotated mRNAs have no UTR).
#  4. Coding potential: CPC2 "noncoding", longest ORF < 100 aa (CPC2 peptide length), and no DIAMOND blastx
#     hit (e-value 1e-5) to the POJ2878 proteome (an unannotated copy of a coding gene is not a lncRNA).
#  5. Expression: StringTie -e on the merged models, TPM per library.
# Re-running skips finished steps. Run from analysis/:  qsub -N lncpilot code/18_lnc_pilot.sh
set -euo pipefail
ENV=/Storage/data2/jorge.munoz/conda_envs/lnc
export PATH=$ENV/bin:$PATH
DOL=/Storage/data1/jorge.munoz/DOLORES
GDIR=$DOL/nf-rnaseq/run_POJ2878/results_pilot/genome
PROT=$DOL/annotation/proteins/POJ2878_v1.proteins.uniq.faa
W=/Storage/data2/jorge.munoz/lnc_pilot
OUT=$DOL/analysis/results/lncrna_pilot
LIBS="Dolores_257_S5 168_S25 230_S15 89_S30 FIDH1800_S30 250_S2 221_S19 21DH5_S21"
T=8
mkdir -p $W/bam $W/st $W/ste $OUT
cd $W
echo "host $(hostname) start $(date)"
for x in hisat2 samtools stringtie gffcompare gffread CPC2.py diamond bedtools; do command -v $x >/dev/null || { echo "missing $x"; exit 1; }; done

# ---- 0. reference
if [ ! -s ref.fa.fai ]; then
  samtools faidx $GDIR/POJ2878_1.0.genome.fa $(awk '$1 ~ /^Chr01g/ {print $1}' $GDIR/POJ2878_1.0.genome.fa.fai) > ref.fa
  samtools faidx ref.fa
fi
awk -F'\t' '$1 ~ /^Chr01g/' $GDIR/Sh_POJ2878_v1.gtf > ref.gtf
[ -s idx.1.ht2 ] || hisat2-build -p $T ref.fa idx > build.log 2>&1
echo "index ready $(date)"

# ---- 1-2. align and assemble per library
for L in $LIBS; do
  if [ ! -s bam/$L.bam.bai ]; then
    hisat2 -p $T --dta --rna-strandness RF --no-unal -x idx \
      -1 $DOL/RAWDATA/${L}_R1_001.fastq.gz -2 $DOL/RAWDATA/${L}_R2_001.fastq.gz 2> bam/$L.hisat2.log |
      samtools sort -@ 4 -m 1G -o bam/$L.bam -
    samtools index bam/$L.bam
  fi
  [ -s st/$L.gtf ] || stringtie bam/$L.bam -p $T --rf -G ref.gtf -o st/$L.gtf -l ${L}
  echo "$L aligned and assembled $(date)"
done
ls st/*.gtf > gtf_list.txt
stringtie --merge -p $T -G ref.gtf -m 200 -o merged.gtf gtf_list.txt
gffcompare -r ref.gtf -o cmp merged.gtf
echo "merged $(date)"

# ---- 3. candidate transcripts by class
# tmap columns: ref_gene_id ref_id class_code qry_gene_id qry_id num_exons FPKM TPM cov len major_iso_id ref_match_len
awk -F'\t' 'NR>1 && ($3=="u"||$3=="x"||$3=="i") && $10>=200 {print $5"\t"$3"\t"$6"\t"$10}' cmp.merged.gtf.tmap > cand_class.tsv
cut -f1 cand_class.tsv > cand_ids.txt
awk -F'\t' 'NR==FNR{k[$1]=1; next} {match($9,/transcript_id "[^"]+"/); id=substr($9,RSTART+15,RLENGTH-16); if (id in k) print}' cand_ids.txt merged.gtf > cand.gtf
gffread -w cand.fa -g ref.fa cand.gtf
# distance to the nearest annotated gene on the same strand (for class u)
awk -F'\t' -v OFS='\t' '$3=="transcript"{match($9,/transcript_id "[^"]+"/); print $1,$4-1,$5,substr($9,RSTART+15,RLENGTH-16),".",$7}' cand.gtf | sort -k1,1 -k2,2n > cand.bed
awk -F'\t' -v OFS='\t' '$3=="transcript"{match($9,/gene_id "[^"]+"/); print $1,$4-1,$5,substr($9,RSTART+9,RLENGTH-10),".",$7}' ref.gtf | sort -k1,1 -k2,2n > refgenes.bed
bedtools closest -s -d -t first -a cand.bed -b refgenes.bed | awk -F'\t' -v OFS='\t' '{print $4,$NF}' > cand_dist.tsv

# ---- 4. coding potential and protein homology
[ -s cpc2.txt ] || CPC2.py -i cand.fa -o cpc2 > cpc2.log 2>&1
[ -s prot.dmnd ] || diamond makedb --in $PROT -d prot -p $T > /dev/null 2>&1
[ -s blastx.tsv ] || diamond blastx -q cand.fa -d prot -p $T -e 1e-5 -k 1 --outfmt 6 qseqid sseqid pident length evalue -o blastx.tsv > /dev/null 2>&1

# ---- 5. expression on the merged models
for L in $LIBS; do
  [ -s ste/$L.gtf ] || stringtie bam/$L.bam -p $T --rf -e -G merged.gtf -o ste/$L.gtf
done

# ---- summary tables
python3 - "$W" "$OUT" $LIBS <<'EOF'
import sys, re, csv, statistics, collections
W, OUT, libs = sys.argv[1], sys.argv[2], sys.argv[3:]
cls = {}; nex = {}; ln = {}
for r in csv.reader(open(f"{W}/cand_class.tsv"), delimiter="\t"):
    cls[r[0]] = r[1]; nex[r[0]] = int(r[2]); ln[r[0]] = int(r[3])
dist = {r[0]: int(r[1]) for r in csv.reader(open(f"{W}/cand_dist.tsv"), delimiter="\t")}
cpc = {}
for r in csv.DictReader(open(f"{W}/cpc2.txt"), delimiter="\t"):
    cpc[r["#ID"]] = (r["label"], int(r["peptide_length"]))
hit = {r[0] for r in csv.reader(open(f"{W}/blastx.tsv"), delimiter="\t")}
tpm = collections.defaultdict(dict)        # transcript -> lib -> TPM
for L in libs:
    for line in open(f"{W}/ste/{L}.gtf"):
        f = line.split("\t")
        if len(f) > 8 and f[2] == "transcript":
            tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
            tpm[tid][L] = float(re.search(r'TPM "([^"]+)"', f[8]).group(1))
rows = []
for t in cls:
    lab, pep = cpc.get(t, ("NA", -1))
    d = dist.get(t, -1)
    flags = []
    if lab != "noncoding": flags.append("CPC2_coding")
    if pep >= 100: flags.append("ORF>=100aa")
    if t in hit: flags.append("protein_hit")
    near = cls[t] == "u" and 0 <= d <= 1000
    if near: flags.append("near_gene_same_strand")
    v = [tpm[t].get(L, 0.0) for L in libs]
    rows.append(dict(transcript=t, class_code=cls[t], exons=nex[t], length=ln[t], cpc2=lab, orf_aa=pep,
                     protein_hit=int(t in hit), dist_same_strand_gene=d,
                     lncRNA=int(not flags or flags == []), filters_failed=";".join(flags) or "pass",
                     libs_tpm_ge1=sum(x >= 1 for x in v), max_tpm=round(max(v), 2),
                     **{f"TPM_{L}": round(x, 2) for L, x in zip(libs, v)}))
with open(f"{OUT}/candidates.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t"); w.writeheader(); w.writerows(rows)
lnc = [r for r in rows if r["lncRNA"]]
def summ(sel, name):
    if not sel: return [name, 0] + [""] * 6
    multi = sum(r["exons"] > 1 for r in sel); exp = sum(r["libs_tpm_ge1"] >= 2 for r in sel)
    return [name, len(sel), multi, round(100 * multi / len(sel), 1), int(statistics.median(r["length"] for r in sel)),
            exp, round(100 * exp / len(sel), 1), round(statistics.median(r["max_tpm"] for r in sel), 2)]
tab = [["set", "transcripts", "multi_exon", "pct_multi_exon", "median_length", "tpm_ge1_in_2plus_libs", "pct", "median_max_tpm"]]
tab.append(summ(rows, "all u/x/i candidates >= 200 nt"))
for f in ["CPC2_coding", "ORF>=100aa", "protein_hit", "near_gene_same_strand"]:
    tab.append(summ([r for r in rows if f in r["filters_failed"]], f"removed: {f}"))
tab.append(summ(lnc, "lncRNA candidates (pass all)"))
for c, n in [("u", "lincRNA (between genes)"), ("x", "antisense"), ("i", "intronic")]:
    tab.append(summ([r for r in lnc if r["class_code"] == c], f"  {n}"))
tab.append(summ([r for r in lnc if r["exons"] > 1 and r["libs_tpm_ge1"] >= 2], "  high confidence: multi-exon, TPM >= 1 in 2+ libs"))
with open(f"{OUT}/summary.tsv", "w") as fh:
    for r in tab: fh.write("\t".join(map(str, r)) + "\n")
for r in tab: print("\t".join(map(str, r)))
EOF

# alignment rates and gffcompare stats
for L in $LIBS; do echo -e "$L\t$(tail -1 bam/$L.hisat2.log)"; done > $OUT/alignment_rates.tsv
cp cmp.stats $OUT/gffcompare.stats
cp merged.gtf $OUT/merged.gtf
awk -F'\t' 'NR==FNR{if ($9=="1") k[$1]=1; next} {match($9,/transcript_id "[^"]+"/); id=substr($9,RSTART+15,RLENGTH-16); if (id in k) print}' $OUT/candidates.tsv cand.gtf > $OUT/lncRNA_candidates.gtf
gzip -f $OUT/merged.gtf
echo "end $(date)"
