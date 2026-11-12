# lncRNA pilot (chromosome 1)

`../../code/18_lnc_pilot.sh` (conda env `/Storage/data2/jorge.munoz/conda_envs/lnc`; BAMs and work files in `/Storage/data2/jorge.munoz/lnc_pilot`).

**Setup.** All 13 copies of chromosome 1 (Chr01g1–g13, 1.49 Gb, 14% of the genome). Eight libraries: 8008 Mt_W (Dolores_257_S5), Mt_D1 (168_S25), Mt_D3 (230_S15), Ft_D1 (89_S30), T bud (FIDH1800_S30); 5000 Mt_W (250_S2), Mt_D3 (221_S19), M bud (21DH5_S21). HISAT2 (strand RF, `--dta`) → StringTie per library, guided by the POJ2878 annotation → merge → gffcompare.

**Filters.** New transcripts of class u (between genes), x (antisense) or i (intronic), ≥ 200 nt; CPC2 noncoding; longest ORF < 100 aa; no DIAMOND blastx hit (1e-5) to the POJ2878 proteome; class u not within 1 kb of a gene on the same strand (possible UTR extension).

**Results** (`summary.tsv`, `candidates.tsv`, `lncRNA_candidates.gtf`):

| Step | Transcripts |
|---|---|
| New u/x/i transcripts ≥ 200 nt | 24,577 |
| … with a hit to a POJ2878 protein (unannotated gene copies, transposons) | 16,453 |
| … coding by CPC2 / ORF ≥ 100 aa / next to a gene | 3,061 / 6,457 / 1,301 (overlapping) |
| **lncRNA candidates** (pass all filters) | **6,773** (6,237 between genes, 112 antisense, 424 intronic) |
| **High confidence:** ≥ 2 exons and TPM ≥ 1 in ≥ 2 libraries | **2,975** transcripts in 2,161 loci |
| … distinct after collapsing allele copies (≥ 90% identity over ≥ 50% of the shorter) | **~1,425** |

- High-confidence candidates: median 605 nt, 2 exons; found on all 13 copies (117–466 per copy).
- Expressed in leaf and bud: 2,385; leaf only: 342; bud only: 248. Both cultivars: 2,257; 8008 only: 510; 5000 only: 208 (5 vs 3 libraries, so not a fair comparison).
- 24–36% of read pairs aligned to chromosome 1; reads from other chromosomes can be pulled in, which the protein filter partly handles.
- Rough genome-wide extrapolation (×7): about 10,000 distinct lncRNAs.
- No drought test here: one library per group.
