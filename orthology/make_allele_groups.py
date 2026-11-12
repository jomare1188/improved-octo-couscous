#!/usr/bin/env python3
"""Group POJ2878 alleles by their sorghum orthologs (OrthoFinder) for allele-group level DE.

Group rules, per POJ2878 gene (OrthoFinder ran on one representative per identical protein;
input/POJ2878_gene_rep.tsv gives the members):
  SB:<sorghum gene(s)>   POJ2878 genes in one row of POJ2878__v__Sbicolor.tsv (orthologues of the same sorghum gene(s))
  OG:<orthogroup>        genes in an orthogroup but with no sorghum orthologue (grouped per orthogroup)
  U:<gene>               genes in no orthogroup (kept alone)
Outputs (allele_groups/):
  gene2group.tsv   gene_id, group_id, group_type
  tx2group.tsv     transcript_id, group_id                 (for tximport)
  group_info.tsv   group_id, type, n_genes, Sbicolor, Zmays, Osativa, Athaliana orthologues, Athaliana symbols
  ortholog_pairs.tsv  POJ2878 gene, species, orthologue      (all pairs, for mapping literature genes)
  group2go_list.tsv   group_id, GO ids (union over member genes) for topGO
"""
import collections, csv, glob, os, sys
D = '/Storage/data1/jorge.munoz/DOLORES'
WD = '/Storage/data2/jorge.munoz/DOLORES_orthofinder/run/Results_Oct29/WorkingDirectory'
# newest OrthoFinder results folder that has orthologues (a -fg rerun writes a new Results_* folder)
cands = [os.path.dirname(os.path.dirname(d)) for d in glob.glob(WD + '/**/Orthologues/Orthologues_POJ2878', recursive=True)]
res = max(cands, key=os.path.getmtime)
out = D + '/orthology/allele_groups'
print('OrthoFinder results:', res, file=sys.stderr)


def split(x):
    return [g for g in x.replace(' ', '').split(',') if g]


members = collections.defaultdict(list)            # rep -> member genes
with open(D + '/orthology/input/POJ2878_gene_rep.tsv') as fh:
    for line in fh:
        r, m = line.split()
        members[r].append(m)

# orthogroup of each POJ rep
og_of = {}
ogf = res + '/Orthogroups/Orthogroups.tsv'
if not os.path.exists(ogf):
    ogf = WD + '/OrthoFinder/Results_Oct29/Orthogroups/Orthogroups.tsv'
with open(ogf) as fh:
    rd = csv.reader(fh, delimiter='\t'); hdr = next(rd); ip = hdr.index('POJ2878')
    for row in rd:
        if ip < len(row):
            for g in split(row[ip]):
                og_of[g] = row[0]

# orthologues of each POJ rep in each species
orth = {sp: collections.defaultdict(set) for sp in ['Sbicolor', 'Zmays', 'Osativa', 'Athaliana']}
sb_row = {}
for sp in orth:
    f = '%s/Orthologues/Orthologues_POJ2878/POJ2878__v__%s.tsv' % (res, sp)
    with open(f) as fh:
        rd = csv.reader(fh, delimiter='\t'); next(rd)
        for row in rd:
            pj, other = split(row[1]), split(row[2])
            for g in pj:
                orth[sp][g].update(other)
                if sp == 'Sbicolor':
                    key = 'SB:' + '|'.join(sorted(other))
                    if g in sb_row and sb_row[g] != key:   # gene in two rows: keep first, report
                        print('warning: %s in two sorghum rows' % g, file=sys.stderr)
                        continue
                    sb_row[g] = key

sym = {}
with open(D + '/orthology/input/Athaliana_symbols.tsv') as fh:
    next(fh)
    for line in fh:
        c = line.rstrip('\n').split('\t')
        if c[1]:
            sym[c[0]] = c[1]

gene_group, gtype = {}, {}
for r in members:
    if r in sb_row:
        grp, t = sb_row[r], 'sorghum_orthologue'
    elif r in og_of:
        grp, t = 'OG:' + og_of[r], 'orthogroup_no_sorghum'
    else:
        grp, t = 'U:' + r, 'unassigned'
    for m in members[r]:
        gene_group[m] = grp
        gtype[grp] = t

with open(out + '/gene2group.tsv', 'w') as fo:
    fo.write('gene_id\tgroup_id\tgroup_type\n')
    for g in sorted(gene_group):
        fo.write('%s\t%s\t%s\n' % (g, gene_group[g], gtype[gene_group[g]]))

n_tx = 0
with open(D + '/nf-rnaseq/run_POJ2878/results/salmon/salmon.merged.tx2gene.tsv') as fh, \
        open(out + '/tx2group.tsv', 'w') as fo:
    next(fh); fo.write('transcript_id\tgroup_id\n')
    for line in fh:
        t, g = line.split('\t')[:2]
        if g in gene_group:
            fo.write('%s\t%s\n' % (t, gene_group[g])); n_tx += 1

grp_genes = collections.defaultdict(list)
for g, grp in gene_group.items():
    grp_genes[grp].append(g)
rep_of = {m: r for r in members for m in members[r]}
with open(out + '/group_info.tsv', 'w') as fo, open(out + '/ortholog_pairs.tsv', 'w') as fp:
    fo.write('group_id\tgroup_type\tn_genes\tSbicolor\tZmays\tOsativa\tAthaliana\tAthaliana_symbols\n')
    fp.write('gene_id\tspecies\torthologue\n')
    for grp in sorted(grp_genes):
        reps = {rep_of[g] for g in grp_genes[grp]}
        cols = []
        for sp in ['Sbicolor', 'Zmays', 'Osativa', 'Athaliana']:
            s = set().union(*(orth[sp].get(r, set()) for r in reps))
            cols.append(','.join(sorted(s)))
        ath = cols[3].split(',') if cols[3] else []
        fo.write('\t'.join([grp, gtype[grp], str(len(grp_genes[grp]))] + cols +
                           [','.join(sorted({sym[a] for a in ath if a in sym}))]) + '\n')
    for g in sorted(gene_group):
        for sp in orth:
            for o in sorted(orth[sp].get(rep_of[g], ())):
                fp.write('%s\t%s\t%s\n' % (g, sp, o))

go = collections.defaultdict(set)
with open(D + '/annotation/go/gene2go.tsv') as fh:
    next(fh)
    for line in fh:
        g, t = line.split('\t')[:2]
        if g in gene_group:
            go[gene_group[g]].add(t)
with open(out + '/group2go_list.tsv', 'w') as fo:
    for grp in sorted(go):
        fo.write('%s\t%s\n' % (grp, ', '.join(sorted(go[grp]))))

c = collections.Counter(gtype[gene_group[g]] for g in gene_group)
print('genes:', len(gene_group), dict(c), file=sys.stderr)
print('groups:', len(grp_genes), dict(collections.Counter(gtype.values())), file=sys.stderr)
print('transcripts mapped:', n_tx, '; groups with GO:', len(go), file=sys.stderr)
sizes = collections.Counter(len(v) for k, v in grp_genes.items() if k.startswith('SB:'))
print('sorghum-anchored group sizes (n genes: n groups):', sorted(sizes.items())[:20], file=sys.stderr)
