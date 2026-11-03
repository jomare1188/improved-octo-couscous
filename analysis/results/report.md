# DOLORES: drought memory in sugarcane, results report

2026-10-29, updated 2026-11-03 with the follow-up analyses (H2 interaction test, bud batch check, rank-based H5 tests, allele-level view, WGCNA). The full narrative is [`report_story.md`](report_story.md). Results are at the allele-group level (78,613 groups) unless stated otherwise. 48 libraries are in the models.

## In short

- **Leaves of 8008 barely "learn".** The response to the 3rd drought repeats the response to the 1st:
  - log2FCs correlate at r = 0.90;
  - 5,673 groups respond in the same direction both times and only 1 reverses;
  - only 84 groups differ between D3 and D1.

  The few groups that respond more strongly at D3 are mostly desiccation-protection genes: the RAB18/XERO1 dehydrin ortholog, two LEA proteins, a seed-storage protein and an embryo-specific protein. RAB18 is the classic trainable memory gene in Arabidopsis, and its maize orthologue is a memory gene too.
- **5000 responds more strongly at D3, but its D1 group is a single plant.** The whole drought response grows, non-memory gene sets included, and translation genes are shut down. It does not look like memory specific to a few genes.
- **H2, that the tolerant 8008 has more memory, is not supported.** A joint genotype × cycle model plus downsampling show that 5000's larger D3-vs-D1 change can be produced by its single D1 plant: refitting 8008 with one D1 plant gives up to 1,171 D3-vs-D1 DEGs (5000: 1,229). 8008 responds more strongly at D1.
- **Buds:**
  - Tillering-stressed and maturation-stressed buds were sequenced on the same run, so T vs M is clean. It shows a timing-of-stress signature led by PAL (phenylpropanoid) genes, which go in opposite directions in the two genotypes.
  - The watered buds were sequenced on a different run and their libraries were prepared differently (fragment length 137 vs 182 bp, adapter 67% vs 23%). Stable genes shift across the run boundary and RUVg's first factor is the treatment itself (r = −0.97), so stressed vs never-stressed buds cannot be rescued.
- **Propagules (exploratory):** in both genotypes, propagules from stressed mothers repress the genes that drought represses less strongly than propagules from watered mothers do. Photosynthesis genes are the clearest case. This fits the better drought performance in Marcos et al. 2018, but each comparison rests on one reference plant.
- **H5:** no chain. Rank-based tests find a modest mother leaf → bud (T vs M) link in 8008 only, and no bud → propagule link. The propagule differences are a milder drought response, and the Ft_D1 vs Fw_D1 profiles of the two cultivars agree (rho 0.20).
- **Gene copies and modules:** alleles within a group move together (median I² = 0), and single-copy counts are too uncertain to go further (median InfRV 2.2). No WGCNA leaf module follows repeated drought; one translation/chloroplast module (M07) is less repressed in propagules of stressed mothers.
- **Next steps, in order:**
  1. Ask Dolores three things: whether the 3 libraries are mislabelled, the time of day each group was sampled, and how stressed each plant was.
  2. Re-sequence a few existing bud and leaf libraries together on one run, to measure the batch effect.
  3. Validate a short gene list by RT-qPCR.

## 1. Data and methods

| Item | Detail |
|---|---|
| Libraries | 51 original libraries. The 18 Drive files were copies and were dropped. 3 libraries carry the other genotype's alleles and are excluded (section 2), leaving **48 in the models**. |
| Quantification | nf-core/rnaseq 3.27.0, Salmon on the POJ2878 genome (397,650 allele-resolved genes) |
| Allele groups | OrthoFinder with sorghum, maize, rice and Arabidopsis. Alleles that share a sorghum orthologue are summed into one group, giving 78,613 groups. Allele-level results are kept as well (`de/gene/`). |
| Differential expression | DESeq2, each genotype separately. Buds: `~condition` (W, T, M). Leaves: mother and propagule leaves in one model, `~condition`, so groups with a single library borrow variance from the others. A DEG has padj < 0.05 and \|log2FC\| ≥ 1. |
| Contrast status | **conclusive**: both groups have ≥ 2 libraries. **exploratory**: one group has a single library. **descriptive**: maturation mother leaves. |
| Memory classes | As in Ding et al. 2013/2014. The first sign is D1 vs W, the second is D3 vs D1. `++`/`--` memory, `+-`/`-+` revised response, `+=`/`-=` no memory. We added `0+`/`0-` for late responders (no D1 response, but a response at D3). |
| Known memory genes | Ding 2013 (Arabidopsis) and Ding 2014 (maize) memory classes, plus 40 curated genes (ABA, LEA/dehydrins, osmolytes, aquaporins, chromatin and DNA-methylation machinery), all mapped to groups by orthology. Tested with limma `cameraPR` in every contrast. |
| Enrichment | topGO weight01 (BP and MF; terms kept at p < 0.01, **not corrected for multiple testing**) and KEGG with clusterProfiler (BH-adjusted). |

## 2. Limits on interpretation (read first)

1. **Three libraries carry the other genotype's alleles.** They are excluded until Dolores confirms (`../data/excluded_samples.tsv`). The check is in `qc/genotype_check/`.
   - `FIDH35000_S29`: labelled 5000 tillering bud, matches 8008.
   - `21DH8_S24`: labelled 8008 maturation bud, matches 5000. With `FIDH35000_S29` it looks like a swap.
   - `251_S3`: labelled 5000 Mt_W leaf, matches 8008.
2. **Batch, from the read headers.**

   | Comparison | Status |
   |---|---|
   | Bud T vs W, M vs W | **confounded**: every W bud was sequenced on run 44, every stressed bud on run 45 |
   | Bud T vs M | clean (same run) |
   | Leaf D1 vs W, D3 vs W | partly confounded: Mt_W on runs 61/63/73, drought leaves on 51/52 |
   | Leaf D3 vs D1 (the memory contrast) | clean |
   | Propagules | clean (all on runs 51/52) |
3. **Single-library groups.** The groups with a single library are:
   - 5000: Mt_D1, Mm_D1, Mm_D3, Fw_D1 and Ft_D1;
   - 8008: Fw_W and Fw_D1.

   Every propagule comparison against Fw_D1 therefore rests on one plant. In 8008 that plant (`79_S6`) is also the weakest library: about 76% rRNA (median 25%), 1.05 million reads on genes (median 6.2 million) and about 15,000 groups detected (median about 29,000).
4. **Age.** Control leaves were sampled on day 39, the same day as D3. D1 was sampled on day 9, so D1 vs W, and D3 vs D1 as well, include about 30 days of growth on top of the drought history.
5. **Possible time-of-day effect.** In 8008, D1 leaves differ from W leaves in clock genes:
   - morning gene CCA1/LHY up (+1.4);
   - PRR7, TOC1, GI and ELF3 down (−1.2 to −3.5);
   - the clock-controlled THI1 strongly down (−5.2).

   Drought does shift the clock, but this pattern also fits sampling at different hours. CCA1/LHY and THI1 also sit in the 8008 memory classes.
6. **Extreme fold changes in single-allele groups.** Many `U:` groups (one unassigned allele) switch on or off between plants, with |log2FC| of 6 to 27. Most of the 5000 `-+` class looks like this: about −8 at D1, then +8 at D3. It points to the one D1 plant rather than to memory.
7. **Reading the enrichment.** GO terms are uncorrected p < 0.01 and often rest on 2 to 5 genes. Some KEGG maps are animal pathways reached through shared gene families: "Estrogen signalling" is HSP70/HSP90, and "MicroRNAs in cancer" is a mix of transporters and kinases. Read them as gene families.

## 3. Hypotheses: answers

| # | Hypothesis | Answer | Strength |
|---|---|---|---|
| H1 | Leaves "learn": the 3rd drought response differs from the 1st | **8008: mostly no.** The response repeats; only a few dozen genes respond more strongly, mainly LEA proteins and dehydrins. **5000: the D3 response is larger, but not resolved.** | 8008 conclusive; 5000 exploratory |
| H2 | 8008 (tolerant) has more memory than 5000 (productive) | **Not supported.** 5000's larger change between cycles is within what one D1 plant can produce (downsampling of 8008). 8008 responds more strongly at D1. | joint model + downsampling |
| H3 | Drought leaves a lasting mark in buds | **Partly.** Long-ago and recent stress leave different bud states (T vs M, clean). Stressed vs never-stressed buds differ in library preparation as well as run, and cannot be rescued. | T vs M conclusive; T/M vs W confounded |
| H4 | Propagules from stressed mothers respond differently to drought | **Suggestive, the same in both genotypes:** drought-repressed genes, photosynthesis above all, are repressed less. | exploratory |
| H5 | Memory is the same genes from mother to bud to propagule | **No chain.** Mother leaf → bud link in 8008 only; no bud → propagule link. | rank-based tests |
| H6 | Known drought-memory genes behave as memory genes in sugarcane | **Partly.** The literature genes are strongly drought-responsive, but as sets they show no memory. A few individual orthologues do, led by the RAB18 dehydrin. | good (set tests) |

### H1. Do tillering leaves respond differently at the 3rd drought?

DEG counts for the tillering mother leaves (up/down):

| Genotype | D1 vs W | D3 vs W | D3 vs D1 | Status of D3 vs D1 |
|---|---|---|---|---|
| 8008 | 3,905 / 4,196 | 3,585 / 3,978 | **35 / 49** | conclusive (D3 n=2, D1 n=3) |
| 5000 | 2,091 / 1,750 | 4,646 / 4,657 | 500 / 729 | exploratory (D1 n=1) |

How similar the two responses are, on the groups that are DEGs in either contrast (`report_tables/h1_d1_vs_d3_concordance.tsv`):

| Genotype | r (log2FC D1 vs D3) | Shared DEGs, same direction | Opposite | Median \|log2FC\| D1 → D3 |
|---|---|---|---|---|
| 8008 | 0.90 | 5,673 | 1 | 1.97 → 2.11 |
| 5000 | 0.82 | 3,272 | 5 | 1.42 → 2.19 |

**Memory classes** (`de/group/memory_class_counts_*`):

| Genotype | `++` | `--` | `+-` | `-+` | `+=` | `-=` | `0+` | `0-` |
|---|---|---|---|---|---|---|---|---|
| 8008 | 22 | 10 | 7 | 12 | 3,876 | 4,174 | 657 | 1,232 |
| 5000 | 117 | 45 | 46 | 94 | 1,928 | 1,611 | 2,772 | 3,254 |

**8008: memory is limited to a small, coherent set of genes.**
- **`++`**, a stronger induction at D3: mostly desiccation and seed-maturation proteins.

  | Gene | D1 vs W | D3 vs D1 |
  |---|---|---|
  | RAB18/XERO1 dehydrin | +10.4 | +1.4 |
  | LEA3/SAG21 | +2.5 | +1.4 |
  | another LEA | +7.8 | +2.0 |
  | seed-storage (nutrient reservoir) protein | +7.1 | +4.9 |
  | embryo-specific protein 3 | +9.3 | +1.8 |

  Also in this class: a sugar transporter (+9.0, then +2.3), cytokinin oxidase CKX5, tyrosine/DOPA decarboxylase and CCA1/LHY.
- **`--`**, a stronger repression at D3: two THI1 groups (thiamine thiazole synthase: −5.2 then −2.6, and −6.2 then −7.8) and 2-isopropylmalate synthase (leucine synthesis). GO and KEGG for `--` and for D3 vs D1 down: thiamine biosynthesis.
- **Late responders:**
  - `0+`: flavonoid and phenylpropanoid biosynthesis, response to hypoxia and jasmonic acid.
  - `0-`: circadian rhythm, protein folding and cell cycle.
- **Overlap with the literature:** our `++` overlaps the Arabidopsis `++` class (3 groups: RAB18, LEA3 and PIT1; OR 25, FDR 0.003).

**5000: the D3 response is stronger overall.**
- **Not memory-specific:** in D3 vs D1, literature memory sets move in the memory direction, but so do non-memory sets (Arabidopsis and maize `+=` up, `-=` down; FDR 2e-5 to 6e-11).
- **Translation shut down:** the `0-` class and D3 vs D1 down are dominated by translation and ribosome genes (GO translation p = 4e-29; KEGG Ribosome, 102 groups, padj = 3e-41). Strong repression of translation usually means stronger stress.
- **Reading:** either D3 was a harsher stress, or the single D1 plant was mildly stressed. The data cannot tell these apart, so stress intensity has to be measured (section 5).
- **Overlap with the literature:**
  - `++` overlaps maize `++`: 8 groups, among them clade A PP2C HAI1-3 and GBF4; OR 24, FDR 4e-8.
  - `-+` overlaps maize `+-` (16 groups, FDR 3e-14), which is the inverted pattern. That is another sign the single D1 plant is atypical.

### H2. Does 8008 have more memory than 5000?

Not supported:
- **8008 is the stable genotype:** 84 groups change between D3 and D1, and the responses correlate at r = 0.90.
- **5000 changes more:** 1,229 groups, r = 0.82, but on one D1 plant.

Follow-up (`h2_genotype_cycle/`, `code/10_genotype_cycle.R`): in one model of both genotypes, the D3-vs-D1 interaction leans to "5000 stronger" (camera FDR 6e-5), but every test of it uses the one 5000 D1 plant. Refitting 8008 with one D1 and two W plants (9 fits) gives 17–1,171 D3-vs-D1 DEGs depending on the D1 plant; with plant 168_S25 it looks like 5000 (r(D1, D3) 0.81). At D1, 8008 induces drought genes more strongly than 5000 (FDR 1.5e-4).

### H3. Does drought leave a lasting mark in buds?

| Genotype | T vs W | M vs W | T vs M (clean) |
|---|---|---|---|
| 8008 | 3,748 / 2,014 | 2,773 / 1,873 | **261 / 146** |
| 5000 | 2,576 / 1,377 | 3,164 / 1,278 | **34 / 113** |

- **T vs W and M vs W:** thousands of DEGs, but watered and stressed buds were sequenced on different runs, so drought and run cannot be separated. In the PCA, watered buds separate from stressed buds along PC2, which is the same axis. Two signals stand out in both genotypes:
  - DNA-demethylation genes are up in stressed buds vs W: ROS1/DME set, FDR 2e-5 to 6e-11.
  - In 5000, RdDM genes are up as well.

  They would matter for transmission, but they carry the batch caveat.
- **T vs M (same run), tillering-stressed (135 days since stress) vs maturation-stressed (6 days since stress).**
  - **Both genotypes:** the top terms come from the **PAL (phenylalanine ammonia-lyase) family**. GO lists them as "salicylic acid catabolic process", "drought recovery", "lignin catabolic process" and "L-phenylalanine catabolic process", and KEGG as phenylalanine metabolism and phenylpropanoid biosynthesis. The direction is **opposite**:
    - 8008: PAL is higher in T buds (GO p = 2e-7).
    - 5000: PAL is higher in M buds (GO p = 2e-16; KEGG padj = 2e-6).
  - **8008 only:**
    - Lower in T buds: thiamine biosynthesis (THI1, p = 4e-8), circadian rhythm and response to chitin.
    - The T vs M DEGs are enriched for genes that respond to drought in mother leaves (140 groups, OR 2.4, FDR 6e-13). 90 of the 140 go in the same direction, so T buds look slightly more "drought-like" than M buds.

### H4. Do propagules from stressed mothers respond differently?

| Genotype | Fw_D1 vs Fw_W | Ft_D1 vs Fw_D1 | Fm_D1 vs Fw_D1 |
|---|---|---|---|
| 8008 | 1,622 / 1,072 (n 1 vs 1) | 136 / 93 (n 2 vs 1) | not sequenced |
| 5000 | 2,672 / 2,851 (n 1 vs 2) | 231 / 145 (n 1 vs 1) | 49 / 36 (n 3 vs 1) |

Drought-stressed propagules from stressed mothers differ little from those of watered mothers, but the difference has the same direction in both genotypes:
- **Gene-set test:** the Arabidopsis `--` set (genes drought represses) is **Up** in Ft_D1 vs Fw_D1:
  - 8008, FDR 1e-6;
  - 5000, FDR 7e-9;
  - and in 5000 Fm_D1 vs Fw_D1, FDR 1e-5.
- **Response size on independently chosen genes.** The genes are the DEGs of the mother-leaf D1 response, so they are not picked from the propagule data (`report_tables/h4_propagule_response.tsv`). Median propagule log2FC vs Fw_W:

  | Genotype | Propagule | Mother-induced genes | Mother-repressed genes |
  |---|---|---|---|
  | 8008 | Fw (watered mother) | +1.06 | −1.07 |
  | 8008 | Ft (tillering-stressed mother) | +1.27 | **−0.55** |
  | 5000 | Fw | +1.73 | −1.71 |
  | 5000 | Ft | +1.39 | **−1.21** |
  | 5000 | Fm (maturation-stressed mother) | +1.46 | **−1.26** |

  Repression is milder in propagules from stressed mothers in both genotypes. Induction is similar.
- **GO, higher in Ft than Fw_D1:**
  - 8008: photosystem I assembly (p = 4e-7), photosynthetic electron transport and photosynthetic acclimation;
  - 5000: photosynthesis (KEGG padj = 0.003), response to heat and to high light.
- **Literature:** this fits Marcos et al. 2018, where propagules from stressed mothers kept better photosynthesis under drought.
- **Caveat:** each test uses the same single Fw_D1 plant per genotype. If that plant was more stressed than the others, it would produce the same pattern.
  - **8008:** the reference is also the weakest library (`79_S6`, section 2). In a thin library many genes look repressed, which would inflate the Fw_D1 repression, so the 8008 result may be partly technical.
  - **5000:** the reference (`42_S39`) is a normal library, and Fm_D1 (3 plants) agrees with Ft_D1, so the 5000 result is the better supported one.

### H5. Is it the same genes from mother to bud to propagule?

This is a preliminary test: overlap of DEG lists (`report_tables/h5_overlaps.tsv`, `code/07_report_stats.R`).

| Overlap | 8008 | 5000 |
|---|---|---|
| Mother-leaf drought genes ∩ Ft_D1 vs Fw_D1 | 88 groups, OR 1.4, FDR 0.02; no direction bias | **152 groups, OR 2.8, FDR 2e-19; 118 of 152 opposite** (milder response) |
| Mother-leaf drought genes ∩ Fm_D1 vs Fw_D1 | — | 30 groups, OR 3.0, FDR 2e-5; 22 of 30 opposite |
| Bud T vs M ∩ Ft_D1 vs Fw_D1 | 5 groups, FDR 0.26 | 2 groups, FDR 0.69 |
| Bud T vs W ∩ Ft_D1 vs Fw_D1 (confounded) | 44 groups, OR 1.6, FDR 0.02 | 46 groups, OR 1.6, FDR 0.01 |
| Leaf D3 vs D1 ∩ bud T vs M | 4 groups, OR 4.9, FDR 0.02, all same direction | 6 groups, FDR 0.37 |

- **Reading:** the propagule effect sits on genes that respond to drought, and mostly makes that response milder. It does not reuse the genes that separate the clean bud groups.
- **Rank-based follow-up** (`h5_rank/`, `code/12_h5_rank.R`; camera on the top 300 genes in both directions, plus RRHO): 8008 mother-leaf drought genes and D3 > D1 genes are higher in T than M buds (FDR 0.01 and 0.001); 5000 shows no such link. Bud T vs M does not predict Ft vs Fw_D1 or Ft vs Fm_D1 (FDR ≥ 0.06). Mother drought genes respond more mildly in Ft/Fm than in Fw_D1 (5000 FDR down to 6e-14; 8008 repressed genes FDR 6e-4).
- **What it can't test:** whether propagules inherit a shifted resting state. That needs watered propagules from stressed mothers (T/W, M/W), which were never sequenced.

### H6. Do known drought-memory genes behave as memory genes?

- **As drought genes, yes.** In 8008 D1 vs W the maize `++` set is Up (FDR 5e-43), and so is the Arabidopsis `++` set (FDR 7e-17). Clade A PP2Cs are Up as well (FDR 2e-9). The orthologues of known drought genes respond in sugarcane.
- **As memory genes, mostly no.**
  - In 8008 D3 vs D1, no `++` set is induced further.
  - The Arabidopsis `--` set goes the other way, less repressed at D3 (FDR 0.005).
  - In 5000, the memory sets do move in the memory direction, but non-memory sets move as well (H1).
- **Individual orthologues that do behave as memory genes:**
  - 8008: RAB18/XERO1 (a `++` gene in both Arabidopsis and maize), LEA3/SAG21 and PIT1.
  - 5000: HAI1-3 and GBF4 (maize `++`).
  - The other canonical markers, RD29A and RD29B, have no grass orthologue.
- **Curated genes** (`memory_genes/curated_heatmap_*.png`). In 8008 D1 vs W (shrunk log2FC):
  - **up:** NCED3 (+1.8), ABF2 (+4.1), HAI1 (+7.4), P5CS1 (+1.8), GolS2 (+3.4) and PIP2;1 (+3.5);
  - **down:** TPS1 (−4.9);
  - **chromatin and methylation genes:** at least one group each of MET1, ROS1, ATX1 and DCL3 goes down by 1.2 to 3.2.

  Few change much between D3 and D1. RAB18 rises another +1.4 (padj 0.006). HSFA2 and AGO4 drop by about 1.5, but neither is significant.

### Maturation mother leaves (descriptive only)

- **5000 drought response:**
  - Mm_D1 vs Mm_W: 3,000 up / 2,563 down.
  - Mm_D3 vs Mm_W: 2,337 / 2,093.
  - Mm_D3 vs Mm_D1: 710 / 865.

  Every drought group is a single plant.
- **8008:** only watered maturation leaves exist.
- **Development stage (Mm_W vs Mt_W):** 915 groups differ in 5000 and 2,246 in 8008. In 5000, ABA-PP2C and aquaporin sets are lower in maturation leaves.

## 4. Figures and tables

| What | File |
|---|---|
| PCA, all libraries (excluded ones marked ×) | `plots/pca_all.png` |
| PCA per sample set and genotype | `plots/pca_panels.png`, `plots/pca_<set>_<genotype>.png` |
| UpSet of DEGs shared between contrasts, split up/down | `plots/upset_degs_5000.png`, `plots/upset_degs_8008.png` |
| Top 20 GO BP terms, up and down, per contrast (50 plots) | `plots/go/go_bp_<genotype>_<contrast>_<up\|down>.png` |
| Curated memory genes, log2FC heatmaps | `memory_genes/curated_heatmap_<genotype>.png` |
| QC: library sizes, PCA, correlation, genotype check | `qc/` |
| DE tables (all groups or alleles) | `de/group/de_<genotype>_<contrast>.tsv.gz`, `de/gene/…` |
| DEG counts | `de/group/deg_summary.tsv` |
| Memory classes per group | `de/group/memory_classes_<genotype>_tillering.tsv.gz` |
| Literature set tests, class overlaps, mapping | `memory_genes/camera_literature_sets.tsv`, `memory_class_overlap_literature.tsv`, `literature_to_groups.tsv` |
| GO and KEGG enrichment | `enrichment/go_enrichment.tsv`, `enrichment/kegg_enrichment.tsv` |
| Numbers added for this report (H1 similarity, H4 response size, H5 overlaps) | `report_tables/` (from `../code/07_report_stats.R`) |

Colours come from scico palettes and were checked for colour-vision deficiency with `../code/palette_check.R`.

## 5. Where to go now

The main uncertainties are about the samples, not the analysis, so the cheapest high-value steps come first.

### Questions for Dolores and the lab (no cost, blocks interpretation)

1. **The three mislabelled libraries** (`FIDH35000_S29`, `21DH8_S24`, `251_S3`): check tube and plate records. Relabel them if the records agree; otherwise drop them for good.
2. **Sampling time of day** for each leaf group and for the buds. This decides whether the clock and thiamine genes (CCA1/LHY, PRR7, TOC1, GI, THI1), which appear among the "memory" genes, are biology or time of sampling.
3. **Stress intensity for each sampled plant:** leaf water potential or RWC, soil moisture, gas exchange at D1, at D3 and in the propagules. This is what separates real memory from "this plant was simply drier":
   - the stronger 5000 response at D3;
   - the milder propagule repression.

   Marcos et al. 2018 may already hold the propagule data.
4. **Library preparation batches.** The data say the watered buds were prepared separately (shorter fragments, more adapter). Confirm from the lab notes, and check whether tissue or RNA is left to re-prepare them.
5. **Material left over.** Is there RNA or tissue for the groups with one library? The key ones are 5000 Mt_D1, Fw_D1 of both genotypes, 5000 Ft_D1 and 8008 Fw_W. Do the T/W and M/W propagules exist?

### Wet lab, if there is budget (highest return first)

1. **Re-sequence a few existing libraries together on one run.** For example, 2 W, 2 T and 2 M buds per genotype, plus 2 Mt_W and 2 Mt_D leaves. Shallow depth is enough. This measures the run effect directly and can rescue bud T vs W, the strongest test of a lasting mark.
2. **Add replicates where conclusions rest on one plant:**
   - 5000 Mt_D1, for H1 and H2;
   - Fw_D1 of both genotypes, for H4;
   - if possible, T/W propagules, so H5 can test an inherited resting state.
3. **RT-qPCR on a short list in every available plant**, sequenced or not:
   - RAB18/XERO1 dehydrin, LEA3/SAG21, the seed-storage protein and HAI1 (memory candidates);
   - THI1 and CCA1/LHY (memory or time of day);
   - PAL (bud signature);
   - two or three photosynthesis genes (propagule signature).
4. **Mechanism, only after the qPCR holds:**
   - H3K4me3 ChIP-qPCR at the RAB18 and LEA groups at D1 vs D3; this is the mark behind transcriptional memory in Arabidopsis.
   - DNA methylation in buds, if T vs W survives the batch check, since the demethylation genes point that way.

### Analysis with the current data: done

All five items are done; details and figures in [`report_story.md`](report_story.md).

1. **H2 interaction test** (`h2_genotype_cycle/`): no evidence that either genotype remembers more.
2. **Bud batch check** (`bud_batch/`): run 44 vs 45 is a library-prep difference (fragment length, adapter); stable genes shift (16–34% padj < 0.05 vs 0.3–0.4% within run); RUVg's first factor is the treatment (r = −0.97). Not rescuable. Side finding: 8008 contrasts against Mt_W carry a GC trend (run fingerprint).
3. **Rank-based H5 tests** (`h5_rank/`): see H5 above.
4. **Allele-level view** (`alleles/`): copies move together (median I² = 0; 4–11% of groups heterogeneous); single-copy InfRV median 2.2 vs 0.23 for groups. THI1 copy Chr02g1 in 8008 is the one possible single-allele case.
5. **WGCNA** (`wgcna/`): 13 leaf modules; 11 follow current drought, none follows repeated drought; M07 (translation/chloroplast) is less repressed in Ft/Fm (8008 FDR 0.04), on the single Fw_D1 plants.

### Possible story, as the data stand

- **Leaves:** sugarcane leaves mostly reset their drought response between cycles. The tolerant 8008 keeps a small, Arabidopsis-like memory module of desiccation-protection genes (RAB18 dehydrin, LEA).
- **Buds:** they carry a stress-timing signature built on phenylpropanoid genes, which differs between genotypes.
- **Propagules:** propagules from stressed mothers repress photosynthesis less under drought, which agrees with the physiology.

The first point is the most solid. The other two need the batch check and the extra replicates before they can be claimed.
