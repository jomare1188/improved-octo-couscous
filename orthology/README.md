# Orthology: POJ2878 allele groups and orthologues (OrthoFinder)

POJ2878 is assembled allele by allele (for example Chr01g1…Chr01g13), so each of its 397,650 genes is one allele.
To work at the gene level, alleles are grouped by their sorghum orthologue. The same orthologues map the literature
memory genes from Arabidopsis and maize onto POJ2878.

## Inputs (`download_proteomes.sh`, `make_input.py`)

- Ensembl Plants proteomes, current release (downloaded 2026-10-29; MD5s are in `download_proteomes.*.log`):
  - *Sorghum bicolor* NCBIv3
  - *Zea mays* B73 RefGen v5 (Zm-B73-REFERENCE-NAM-5.0)
  - *Oryza sativa* IRGSP-1.0
  - *Arabidopsis thaliana* TAIR10
- One protein per gene: the longest isoform.
- POJ2878 genes with an identical longest protein are collapsed: 397,650 genes become 285,755 sequences.
  `input/POJ2878_gene_rep.tsv` maps each representative to its members.

## Run

- OrthoFinder 2.5.4 (`ORTHOFINDER` conda env), with DIAMOND `--more-sensitive -e 0.001`, as OrthoFinder sets it.
- Working directory, on data2 for space: `/Storage/data2/jorge.munoz/DOLORES_orthofinder/run/Results_Oct29/WorkingDirectory`
- Steps:
  1. `orthofinder -f fasta -op`: databases and the list of searches.
  2. `make_tasks.sh` → `tasks.txt`. That gives 120 DIAMOND searches: the POJ2878 query is split into 20 chunks.
  3. `run_diamond_array.sh`: an SGE array, 8 threads per task. It is resumable, since finished tasks are skipped.
  4. `run_orthofinder_finish.sh`: joins the chunk results, then runs `orthofinder -b` (16 threads).
  5. `run_orthofinder_orthologues.sh`: 2.5.4 built the orthogroups (`Results_Oct29`), then crashed in the species-tree
     step: the conda env's Python 3.11 rejects `open(..., 'rU')`. Orthologues were computed from those orthogroups with
     the cluster module OrthoFinder 2.5.5 (`orthofinder -fg`), giving `Results_Oct29_1`.
- Results are in `WorkingDirectory/OrthoFinder/Results_Oct29` (orthogroups) and `Results_Oct29_1` (orthologues).
- **Orthogroups:** 85.7% of POJ2878 proteins fall in an orthogroup, against 83–91% for the other species.
- **Allele groups:**
  - 397,650 genes become 78,613 groups.
  - 24,210 groups are sorghum-anchored, and they hold 218,424 genes.
  - 13,664 are orthogroup-only, holding 124,868 genes.
  - 40,739 are unassigned, holding 54,358 genes.
  - The most common sorghum-anchored group size is 12 alleles.

## Allele groups (`make_allele_groups.py` → `allele_groups/`)

- `SB:<sorghum gene>`: POJ2878 genes that are orthologues of the same sorghum gene(s).
- `OG:<orthogroup>`: genes in an orthogroup with no sorghum orthologue.
- `U:<gene>`: genes in no orthogroup, kept alone.
- Output files:
  - `gene2group.tsv`
  - `tx2group.tsv` (for tximport)
  - `group_info.tsv`, which lists the sorghum, maize, rice and Arabidopsis orthologues of each group
  - `ortholog_pairs.tsv`
  - `group2go_list.tsv`: the union of the member genes' GO terms, in topGO format.
