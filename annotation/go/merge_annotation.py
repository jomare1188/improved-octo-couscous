#!/usr/bin/env python3
"""Merge InterProScan 6 and eggNOG-mapper annotations of POJ2878 into GO and gene tables.

Both tools ran on the unique protein sequences (one representative per identical sequence).
Each representative's annotation is copied to every transcript with that sequence (uniq_map.tsv),
then transcripts are grouped into genes with the Salmon tx2gene table, so gene IDs match the counts.

Outputs (in the output directory):
  transcript2go.tsv        transcript_id, gene_id, go_id, source   (source: ips | eggnog | ips,eggnog)
  gene2go.tsv              gene_id, go_id, source                  (union over the gene's isoforms)
  gene2go_list.tsv         gene_id, GO ids comma-separated          (e.g. topGO readMappings)
  transcript_annotation.tsv / gene_annotation.tsv
                           name and description (eggNOG), InterPro, Pfam, KEGG KO / pathway, EC
  coverage.tsv             how many transcripts and genes have each kind of annotation
GO terms are taken as reported, without adding parent terms (enrichment tools do that).
InterProScan pathways are not used: they are mostly Reactome IDs for animal species.
"""
import argparse
import collections
import csv
import glob
import os
import sys

p = argparse.ArgumentParser()
p.add_argument('--ips', nargs='+', required=True, help='InterProScan TSV files')
p.add_argument('--eggnog', required=True, help='emapper .annotations file')
p.add_argument('--uniq-map', required=True, help='rep_id <tab> member_id')
p.add_argument('--tx2gene', required=True, help='transcript_id <tab> gene_id (with header)')
p.add_argument('--outdir', required=True)
a = p.parse_args()
os.makedirs(a.outdir, exist_ok=True)


def log(*x):
    print(*x, file=sys.stderr, flush=True)


def dash(v):
    return v if v not in ('', '-') else ''


# --- InterProScan (per representative protein) ---
ips_go = collections.defaultdict(set)
ips_ipr = collections.defaultdict(dict)    # rep -> {IPR: description}
ips_pfam = collections.defaultdict(set)
ips_any = set()
for f in a.ips:
    log('reading', f)
    with open(f) as fh:
        for line in fh:
            c = line.rstrip('\n').split('\t')
            if len(c) < 14:
                continue
            q = c[0]
            ips_any.add(q)
            if c[3] == 'Pfam':
                ips_pfam[q].add(c[4])
            if dash(c[11]):
                ips_ipr[q][c[11]] = c[12]
            if 'GO:' in c[13]:
                for g in c[13].split('|'):
                    ips_go[q].add(g.split('(')[0])
log('InterProScan: %d proteins with matches, %d with GO' % (len(ips_any), len(ips_go)))

# --- eggNOG (per representative protein) ---
egg = {}
egg_go = {}
with open(a.eggnog) as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        c = line.rstrip('\n').split('\t')
        q = c[0]
        egg[q] = dict(name=dash(c[8]), desc=dash(c[7]), ec=dash(c[10]), ko=dash(c[11]).replace('ko:', ''),
                      pathway=','.join(x for x in dash(c[12]).split(',') if x.startswith('map')),
                      cog=dash(c[6]), pfam=dash(c[20]))
        if dash(c[9]):
            egg_go[q] = set(c[9].split(','))
log('eggNOG: %d proteins annotated, %d with GO' % (len(egg), len(egg_go)))

# --- transcript -> representative, transcript -> gene ---
rep_of = {}
with open(a.uniq_map) as fh:
    for line in fh:
        r, m = line.rstrip('\n').split('\t')
        rep_of[m] = r
gene_of = {}
with open(a.tx2gene) as fh:
    rd = csv.reader(fh, delimiter='\t')
    next(rd)
    for row in rd:
        gene_of[row[0]] = row[1]
missing = [t for t in rep_of if t not in gene_of]
log('%d protein-coding transcripts, %d in tx2gene, %d not in tx2gene' % (len(rep_of), len(rep_of) - len(missing), len(missing)))
if missing:
    log('  e.g.', missing[:5])

transcripts = sorted(t for t in rep_of if t in gene_of)
genes_all = sorted(set(gene_of.values()))


def join(xs):
    return ','.join(sorted(xs))


# --- transcript tables ---
tx_go = {}
gene_go = collections.defaultdict(lambda: collections.defaultdict(set))
gene_tx = collections.defaultdict(list)
with open(os.path.join(a.outdir, 'transcript2go.tsv'), 'w') as go_out, \
        open(os.path.join(a.outdir, 'transcript_annotation.tsv'), 'w') as an_out:
    go_out.write('transcript_id\tgene_id\tgo_id\tsource\n')
    an_out.write('transcript_id\tgene_id\trepresentative\tpreferred_name\tdescription\tinterpro\tinterpro_desc\t'
                 'pfam\tkegg_ko\tkegg_pathway\tec\tcog_category\tn_go\n')
    for t in transcripts:
        g = gene_of[t]
        r = rep_of[t]
        gene_tx[g].append(t)
        src = collections.defaultdict(set)
        for go in ips_go.get(r, ()):
            src[go].add('ips')
        for go in egg_go.get(r, ()):
            src[go].add('eggnog')
        for go in sorted(src):
            s = ','.join(sorted(src[go], key=['ips', 'eggnog'].index))
            go_out.write('%s\t%s\t%s\t%s\n' % (t, g, go, s))
            gene_go[g][go] |= src[go]
        tx_go[t] = len(src)
        e = egg.get(r, {})
        ipr = ips_ipr.get(r, {})
        pfam = set(ips_pfam.get(r, ()))
        an_out.write('\t'.join([t, g, r, e.get('name', ''), e.get('desc', ''), join(ipr),
                                '; '.join(ipr[k] for k in sorted(ipr)), join(pfam), e.get('ko', ''),
                                e.get('pathway', ''), e.get('ec', ''), e.get('cog', ''), str(len(src))]) + '\n')

# --- gene tables ---
with open(os.path.join(a.outdir, 'gene2go.tsv'), 'w') as go_out, \
        open(os.path.join(a.outdir, 'gene2go_list.tsv'), 'w') as ls_out:
    go_out.write('gene_id\tgo_id\tsource\n')
    for g in sorted(gene_go):
        for go in sorted(gene_go[g]):
            go_out.write('%s\t%s\t%s\n' % (g, go, ','.join(sorted(gene_go[g][go], key=['ips', 'eggnog'].index))))
        ls_out.write('%s\t%s\n' % (g, ', '.join(sorted(gene_go[g]))))


def uniq_join(vals, sep=','):
    out = []
    for v in vals:
        for x in v.split(sep) if v else ():
            if x and x not in out:
                out.append(x)
    return sep.join(sorted(out))


with open(os.path.join(a.outdir, 'gene_annotation.tsv'), 'w') as out:
    out.write('gene_id\tn_transcripts\tpreferred_name\tdescription\tinterpro\tinterpro_desc\tpfam\t'
              'kegg_ko\tkegg_pathway\tec\tn_go\n')
    for g in genes_all:
        txs = gene_tx.get(g, [])
        reps = [rep_of[t] for t in txs]
        es = [egg[r] for r in reps if r in egg]
        # name/description: the most common one among isoforms (first transcript wins ties)
        names = collections.Counter(e['name'] for e in es if e['name'])
        descs = collections.Counter(e['desc'] for e in es if e['desc'])
        ipr = {}
        for r in reps:
            ipr.update(ips_ipr.get(r, {}))
        pfam = set()
        for r in reps:
            pfam |= ips_pfam.get(r, set())
        out.write('\t'.join([g, str(len(txs)),
                             names.most_common(1)[0][0] if names else '',
                             descs.most_common(1)[0][0] if descs else '',
                             join(ipr), '; '.join(ipr[k] for k in sorted(ipr)), join(pfam),
                             uniq_join([e['ko'] for e in es]), uniq_join([e['pathway'] for e in es]),
                             uniq_join([e['ec'] for e in es]), str(len(gene_go.get(g, ())))]) + '\n')

# --- coverage ---
def frac(n, d):
    return '%d\t%d\t%.1f' % (n, d, 100.0 * n / d if d else 0)


nt = len(transcripts)
coding_genes = sorted(gene_tx)
ng = len(coding_genes)
rows = []
tx_ips_go = sum(1 for t in transcripts if rep_of[t] in ips_go)
tx_egg_go = sum(1 for t in transcripts if rep_of[t] in egg_go)
tx_any_go = sum(1 for t in transcripts if tx_go[t])
rows += [('transcripts', 'GO (InterProScan)', frac(tx_ips_go, nt)),
         ('transcripts', 'GO (eggNOG)', frac(tx_egg_go, nt)),
         ('transcripts', 'GO (either)', frac(tx_any_go, nt)),
         ('transcripts', 'InterPro entry', frac(sum(1 for t in transcripts if ips_ipr.get(rep_of[t])), nt)),
         ('transcripts', 'eggNOG hit', frac(sum(1 for t in transcripts if rep_of[t] in egg), nt)),
         ('transcripts', 'KEGG KO', frac(sum(1 for t in transcripts if egg.get(rep_of[t], {}).get('ko')), nt))]


def gene_has(fn):
    return sum(1 for g in coding_genes if any(fn(rep_of[t]) for t in gene_tx[g]))


rows += [('genes', 'GO (InterProScan)', frac(gene_has(lambda r: r in ips_go), ng)),
         ('genes', 'GO (eggNOG)', frac(gene_has(lambda r: r in egg_go), ng)),
         ('genes', 'GO (either)', frac(sum(1 for g in coding_genes if gene_go.get(g)), ng)),
         ('genes', 'InterPro entry', frac(gene_has(lambda r: bool(ips_ipr.get(r))), ng)),
         ('genes', 'eggNOG hit', frac(gene_has(lambda r: r in egg), ng)),
         ('genes', 'KEGG KO', frac(gene_has(lambda r: bool(egg.get(r, {}).get('ko'))), ng)),
         ('genes', 'any annotation', frac(gene_has(lambda r: r in egg or r in ips_any), ng))]
n_pairs = sum(len(v) for v in gene_go.values())
n_terms = len({go for v in gene_go.values() for go in v})
with open(os.path.join(a.outdir, 'coverage.tsv'), 'w') as out:
    out.write('level\tannotation\tn\ttotal\tpercent\n')
    for r in rows:
        out.write('\t'.join(r) + '\n')
    out.write('genes\tgene-GO pairs\t%d\t\t\n' % n_pairs)
    out.write('genes\tdistinct GO terms\t%d\t\t\n' % n_terms)
log('genes in tx2gene: %d; protein-coding genes: %d' % (len(genes_all), ng))
with open(os.path.join(a.outdir, 'coverage.tsv')) as fh:
    sys.stderr.write(fh.read())
log('done')
