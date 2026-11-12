#!/usr/bin/env python3
"""Make OrthoFinder input: one protein per gene (longest isoform), FASTA headers = gene IDs.

POJ2878: longest isoform per gene from annotation/proteins/POJ2878_v1.proteins.faa (gene = ID without .tN).
Genes whose longest protein is identical are collapsed to one representative (POJ2878 is allele-resolved,
many alleles code the same protein); input/POJ2878_gene_rep.tsv maps rep_gene <tab> member_gene.
Ensembl proteomes: longest isoform per 'gene:' tag; gene symbols kept in input/<species>_symbols.tsv.
"""
import gzip, os, re, sys
D = '/Storage/data1/jorge.munoz/DOLORES/orthology'
out = os.path.join(D, 'input')


def read_fa(path):
    op = gzip.open if path.endswith('.gz') else open
    h, s = None, []
    with op(path, 'rt') as fh:
        for line in fh:
            if line.startswith('>'):
                if h:
                    yield h, ''.join(s)
                h, s = line[1:].rstrip('\n'), []
            else:
                s.append(line.strip())
    if h:
        yield h, ''.join(s)


def clean(seq):
    return seq.rstrip('*').replace('*', 'X').replace('.', 'X')


def write(name, best):
    with open(os.path.join(out, name + '.fa'), 'w') as fo:
        for g in sorted(best):
            fo.write('>%s\n%s\n' % (g, best[g]))
    print(name, len(best), file=sys.stderr)


# --- POJ2878 ---
best = {}
for h, s in read_fa('/Storage/data1/jorge.munoz/DOLORES/annotation/proteins/POJ2878_v1.proteins.faa'):
    g = re.sub(r'\.t\d+$', '', h.split()[0])
    s = clean(s)
    if g not in best or len(s) > len(best[g]):
        best[g] = s
rep_of_seq, rep = {}, {}
for g in sorted(best):
    r = rep_of_seq.setdefault(best[g], g)
    rep[g] = r
with open(os.path.join(out, 'POJ2878_gene_rep.tsv'), 'w') as fo:
    for g in sorted(rep):
        fo.write('%s\t%s\n' % (rep[g], g))
print('POJ2878 genes', len(best), 'unique proteins', len(rep_of_seq), file=sys.stderr)
write('POJ2878', {r: s for s, r in rep_of_seq.items()})

# --- Ensembl ---
for sp, f in [('Sbicolor', 'Sorghum_bicolor.Sorghum_bicolor_NCBIv3.pep.all.fa.gz'),
              ('Zmays', 'Zea_mays.Zm-B73-REFERENCE-NAM-5.0.pep.all.fa.gz'),
              ('Osativa', 'Oryza_sativa.IRGSP-1.0.pep.all.fa.gz'),
              ('Athaliana', 'Arabidopsis_thaliana.TAIR10.pep.all.fa.gz')]:
    best, sym = {}, {}
    for h, s in read_fa(os.path.join(D, 'proteomes', f)):
        m = re.search(r' gene:(\S+)', h)
        if not m:
            continue
        g = m.group(1)
        s = clean(s)
        sm = re.search(r' gene_symbol:(\S+)', h)
        ds = re.search(r' description:(.*?)( \[Source|$)', h)
        if sm or ds:
            sym[g] = (sm.group(1) if sm else '', ds.group(1) if ds else '')
        if g not in best or len(s) > len(best[g]):
            best[g] = s
    write(sp, best)
    with open(os.path.join(out, sp + '_symbols.tsv'), 'w') as fo:
        fo.write('gene_id\tsymbol\tdescription\n')
        for g in sorted(sym):
            fo.write('%s\t%s\t%s\n' % (g, *sym[g]))
