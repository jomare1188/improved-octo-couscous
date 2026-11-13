# Does sugarcane remember a drought? The DOLORES RNA-seq story, in short

*November 2026. This is the short version: the main insights, what was done with the problematic samples, and which comparisons can be trusted. The full documentation, with every analysis, figure, method and concept explained, is the project [`README.md`](../../README.md); section numbers in brackets point to it. Each section ends with a **Where** line: the files that hold its numbers, relative to `analysis/results/`, and the README section with the full account.*

---

## 1. The question

A plant that survives a drought is sometimes better prepared for the next one: **drought memory**. Sugarcane is planted from stalk pieces, so a mark left in the mother's **buds** could pass to the new plants, the **propagules**. DOLORES looks for that memory in gene expression in two cultivars, the productive **5000** (IACSP95-5000) and the drought-tolerant **8008** (IACCTC07-8008).

![Experiment 1: mother plants](../../figs/dry_cycle_1.png)

*Figure 1. Mother plants. T plants go through three drought cycles while young (days 0–45), M plants late in the season (days 135–180), W plants are always watered. Green arrows: leaf samples, at the peak of the 1st (D1) and 3rd (D3) drought. Red arrows: bud samples, all on day 180.*

Three hypotheses carry the story:

| # | Hypothesis | Answer |
|---|---|---|
| **H1** | Leaves respond differently to the 3rd drought than to the 1st. | **Yes, for a defined set of genes.** In 8008, among about 8,000 groups that respond the same way both times, 84 respond differently at the 3rd drought. They are led by drying-protection genes such as the RAB18 dehydrin. |
| **H3** | Drought leaves a lasting mark in the mother's buds. | **Yes, for the timing of the stress.** Buds stressed 135 days before sampling differ from buds stressed 6 days before in 407 groups (8008) and 147 (5000), led by PAL genes. |
| **H4** | Propagules of stressed mothers respond differently to drought. | **Yes, in both cultivars.** They respond more mildly and spare photosynthesis and growth genes. Exploratory: one reference plant per cultivar. |

The other hypotheses (H2, H5, H6) are too weak to conclude from with the present samples (section 7 below).

**Words used below.** A **group** is one gene: sugarcane has 8–12 copies of most genes, and the copies are summed (78,613 groups). A **DEG** is a group that changes at least two-fold with padj < 0.05. `A_vs_B` means A compared with B. `Mt_W`, `Mt_D1`, `Mt_D3` are tillering mother leaves: watered, 1st drought, 3rd drought. T and M buds come from tillering-stressed and maturation-stressed mothers. `Fw`, `Ft`, `Fm` are propagules of watered, T and M mothers; `_D1` means under their own drought.

---

## 2. The samples, and what we did with the problematic ones

78 samples were planned (26 groups of 3 plants). **48 libraries are in the analysis.**

| Problem | Libraries | What we did |
|---|---|---|
| Never available | 9 samples | Nothing to analyse. |
| Copies of other libraries | 18 files from the Google Drive | **Dropped.** Each one repeats the reads of one of the 51 original libraries, often under another sample's name. |
| Alleles of the other cultivar | `FIDH35000_S29` (labelled 5000 T bud), `21DH8_S24` (labelled 8008 M bud), `251_S3` (labelled 5000 watered leaf) | **Excluded, not relabelled.** Their expression suggests each keeps the treatment on its label, in the other cultivar. But every group was designed with 3 plants, relabelling would give groups of 4, and for two of the three the group cannot be assigned with confidence. |
| Failed poly-A selection | `79_S6` (76% rRNA, 1 million usable reads) | **Kept and flagged.** It is the only 8008 propagule of a watered mother under drought. |
| Groups with one plant | 5000: Mt_D1, Mm_D1, Mm_D3, Fw_D1, Ft_D1. 8008: Fw_W, Fw_D1 | **Kept.** Contrasts that use them are called **exploratory**: they can be computed, but one unusual plant would go unnoticed. |

The sample table from the lab is correct, so the cultivar mix-up happened when the tissue was sampled or the tubes were handled. **Where:** list of libraries `../data/samples.tsv`, excluded ones `../data/excluded_samples.tsv`; cultivar check `qc/genotype_check/`; which group each excluded library resembles `qc/relabel_check/`; poly-A check `qc/library_type/`. README sections 3.1 and 5.1.

---

## 3. Which contrasts can be biased by sequencing run or sampling time

The libraries were sequenced on several runs, and leaves of different ages were sampled on different days. Some contrasts compare groups that differ in run or day as well as in treatment.

| Contrast | Sequencing run | Sampling day | Can it be trusted? |
|---|---|---|---|
| Leaf D1 vs watered | drought leaves on runs 51/52, watered on 61/63/73 | day 9 vs day 39 | **Biased by run and by 30 days of age.** The drought response is large, so its main features are safe. Small differences are not. |
| Leaf D3 vs watered | same split of runs | both day 39 | **Biased by run.** Same reading. |
| **Leaf D3 vs D1** (H1) | both on runs 51/52 | day 39 vs day 9 | **Clean for run.** It carries 30 days of age, which no design can remove: the 3rd drought always comes later. |
| Bud T or M vs watered | stressed buds on run 45, watered on run 44, with a different library preparation | all day 180 | **Not usable.** Run and treatment are the same thing here, and no correction separates them. |
| **Bud T vs M** (H3) | both on run 45 | both day 180 | **Clean.** The stress happened 135 days (T) or 6 days (M) before, which is the comparison itself. |
| **Propagules** (H4) | all on runs 51/52; the reference plant is on the other run from the plants compared with it | same drought cycle | **Clean for run:** runs 51 and 52 behave as one batch, and the arrangement is reversed between the cultivars. Exploratory because of the single reference plant. |
| Maturation leaves | runs 51/52 | D1 day 144, D3 and watered day 174 | Descriptive only: one plant per drought group, 5000 only. |

**The hour of sampling is very important, and every contrast is confounded with it.** The lab confirmed that each sampling started early in the morning and in some cases lasted until past midday, and the hour of each plant was not recorded. Two groups in one contrast can therefore differ by several hours of the day, and this cannot be corrected. What it touches:
- **What the clock genes say.** We read the internal time of each library from six clock genes: the morning gene CCA1/LHY, which falls after dawn, and five day and evening genes (PRR73, PRR37, PRR59/95, GI, TOC1), which rise towards midday. A later hour lowers the first and raises the others.
  - Plants of one group share their internal time (99% of the variation is between groups), so the groups were sampled one after another. Hour is therefore tied to the group, and no model can remove it.
  - **1st and 3rd drought leaves** have the same level of the five day and evening genes, in both cultivars. A difference of several hours between them is unlikely. Only CCA1/LHY is higher at the 3rd drought (3 times in 8008, 2 times in 5000).
  - **T and M buds** and the **propagule groups** differ little in these genes.
  - **Drought against watered leaves** is where the clock genes differ most (in 8008 the watered leaves look later in the day), so genes of the daily clock are not interpreted in those contrasts.
- **The main results, less so.** The H1 memory genes are drying-protection genes, not clock genes, and an hour effect would not pick protection genes out of 8,000. T and M buds (H3) were collected on the same day, in an order we do not know. The propagule result (H4) is a shift of thousands of genes in the same direction in two cultivars. An hour effect on a single gene can never be excluded.

**Where:** runs per library `qc/sequencing_runs.tsv` and `plots/story/design_runs.png`; distance between libraries within and across runs `qc/run_pair_distances.tsv`; bud run effect `bud_batch/`; clock genes `clock/clock_by_contrast.tsv`, `clock/clock_scores.png`, made by `../code/23_clock_check.R`. README sections 5.2, 5.3 and 5.5.

---

## 4. Leaves: a set of genes remembers (H1), among conserved response.

### 4.1 What drought does to a leaf

| Cultivar | D1 vs W (up / down) | D3 vs W (up / down) | D3 vs D1 (up / down) |
|---|---|---|---|
| 8008 | 3,905 / 4,196 | 3,585 / 3,978 | **35 / 49** |
| 5000 | 2,091 / 1,750 | 4,646 / 4,657 | 500 / 729 (one D1 plant) |

*DEG counts in tillering mother leaves.*

What these genes do, from GO and KEGG enrichment:

| Contrast | Higher (up) | Lower (down) |
|---|---|---|
| **8008** D1 vs W | Response to water deprivation, heat and hydrogen peroxide (dehydrins, chaperones). Phenylpropanoid, flavonoid and stilbenoid synthesis. Breakdown of branched-chain amino acids and tryptophan; galactose metabolism; glyoxylate cycle. | Chloroplast organisation, chlorophyll synthesis, photosystem II assembly, antenna proteins, Calvin cycle. Thiamine synthesis. rRNA processing, ribosome biogenesis, tRNA charging. |
| **8008** D3 vs W | The same programme: heat, high light, water deprivation. Flavonoids and phenylpropanoids even stronger, lignin, cutin and wax. | Thiamine synthesis, circadian rhythm, chloroplast organisation, rRNA processing, tRNA charging, sulfur and nitrogen metabolism. |
| **8008** D3 vs D1 | Drying protection: RAB18 dehydrin, two LEA proteins, embryo-specific protein ATS3, a seed-storage protein. Also a sugar transporter, cytokinin oxidase, tyrosine/DOPA decarboxylase, glutamine synthetase and the clock gene CCA1/LHY. | Thiamine synthesis (THI1), papain-like cysteine proteases, branched-chain amino acid synthesis, rRNA maturation. |
| **5000** D1 vs W | Chlorophyll breakdown, breakdown of branched-chain amino acids, phenylpropanoids, tryptophan metabolism and auxin synthesis, response to water deprivation. | Photosynthesis, antenna proteins, Calvin cycle, chlorophyll synthesis, water transport (aquaporins), cytokinin signalling. |
| **5000** D3 vs W | Phenylpropanoids and flavonoids, amino acid breakdown, fatty acid β-oxidation, glutathione metabolism, chlorophyll breakdown. | Translation (482 groups), photosynthesis, chloroplast organisation, chlorophyll synthesis, tRNA charging, ribosome biogenesis. |
| **5000** D3 vs D1 | ABA response, PP2C phosphatases, MAPK and hormone signalling. | Translation and ribosome (102 ribosome groups), protein folding (HSP70/HSP90 chaperones). |

Both cultivars run the same programme: protection and secondary metabolism go up, photosynthesis and protein synthesis go down.

**Where:** DEG counts `de/group/deg_summary.tsv`; one table per contrast `de/group/de_<cultivar>_<contrast>.tsv.gz`; GO and KEGG `enrichment/go_enrichment.tsv`, `enrichment/kegg_enrichment.tsv`, plots in `plots/go/`. README section 6.1.

### 4.2 The 3rd drought replays the 1st, and some genes stand out

![D1 vs D3 response](plots/story/d1_vs_d3.png)

*Figure 2. Each dot is a group that responds to drought at D1 or at D3. x = response at the 1st drought (log2FC, D1 vs watered), y = response at the 3rd (D3 vs watered). On the dashed diagonal the two are equal; the further a dot is from it, the more the 3rd response differs from the 1st.*

Each group gets a **memory class** of two signs (Ding et al.). The first sign is the 1st drought against watered: induced `+`, repressed `-`, no response `0`. The second sign is the 3rd drought against the 1st: higher `+`, lower `-`, the same `=`. In Figure 2:

| Class | Meaning | Where in Figure 2 | 8008 | 5000 |
|---|---|---|---|---|
| `+=` | induced, and the same at the 3rd drought: no memory | grey, upper right, along the diagonal | 3,876 | 1,928 |
| `-=` | repressed, and the same at the 3rd drought: no memory | grey, lower left, along the diagonal | 4,174 | 1,611 |
| `++` | induced, and induced more at the 3rd drought: **memory** | red dots, upper right, above the diagonal | 22 | 117 |
| `--` | repressed, and repressed more at the 3rd drought: **memory** | blue dots, lower left, below the diagonal | 10 | 45 |
| `+-` | induced at the 1st drought, back towards the watered level at the 3rd: revised response | dark triangles, right half, below the diagonal | 7 | 46 |
| `-+` | repressed at the 1st drought, back towards the watered level at the 3rd: revised response | dark triangles, left half, above the diagonal | 12 | 94 |
| `0+` | no response at the 1st drought, induced at the 3rd: late responder | grey, near x = 0, upper half | 657 | 2,772 |
| `0-` | no response at the 1st drought, repressed at the 3rd: late responder | grey, near x = 0, lower half | 1,232 | 3,254 |

*Number of groups per class in tillering mother leaves. The 5000 numbers rest on one D1 plant.*

In 8008 the two responses correlate at **r = 0.90**, and of 5,674 groups that respond at both droughts, 5,673 respond in the same direction. This reproducibility is the sea against which the memory genes are visible:
- **84 groups** respond differently at the 3rd drought (35 higher, 49 lower).
- **32 are memory genes** in the sense of Ding et al.: 22 respond to drought and respond more strongly the 3rd time (`++`), 10 are repressed and more repressed the 3rd time (`--`).
- **1,889 more respond at the 3rd drought only** (`0+` and `0-`; flavonoid and phenylpropanoid synthesis up; circadian rhythm and protein folding down). They can be a late response or an effect of the 30 days of age between the two samplings, so they are not counted as memory.

The `++` genes are **drying-protection genes**:

| Gene | D1 vs W (log2FC) | D3 vs D1 (log2FC) |
|---|---|---|
| RAB18/XERO1 dehydrin | +10.4 | +1.4 |
| LEA3/SAG21 | +2.5 | +1.4 |
| another LEA protein | +7.8 | +2.0 |
| seed-storage protein | +7.1 | +4.9 |
| embryo-specific protein 3 | +9.3 | +1.8 |
| sugar transporter | +9.0 | +2.3 |

RAB18 is a known *Arabidopsis* gene that remembers drought, and its maize orthologue is a memory gene too. We compared our list with the published one:
- **"Shares" means the same gene is in both lists.** Ding et al. 2013 list the *Arabidopsis* `++` genes. Their sugarcane orthologues are 240 groups. Three of them (RAB18, LEA3 and PIT1) are also among our 22 `++` groups.
- **By chance we would expect about 0.14 shared groups, and we find 3.** Our 22 groups are 0.06% of the 37,767 groups tested, so 240 groups picked at random would hold 0.14 of them.
- **The odds ratio of 25 says the same thing as one number.** A sugarcane group whose *Arabidopsis* orthologue is a `++` memory gene has 25 times higher odds of being `++` here than any other group. An odds ratio of 1 would mean no link. The FDR of 0.003 is the probability of an overlap this large by chance, corrected for the number of class pairs tested.

![Example genes](plots/story/example_genes.png)

*Figure 3. Four genes plant by plant in tillering mother leaves (y = log2 expression). Left 8008, right 5000. RAB18 and LEA3 keep rising at D3 in 8008; HAI1 rises and stays; THI1 drops and drops further.*

On the repressed side, thiamine synthesis (THI1) is shut down by drought and shut down further at the 3rd drought. THI1 follows the daily clock, so we checked it against the clock genes (section 3):
- **THI1 is a drought-repressed gene, whatever the hour.** In 8008, 11 of the 13 THI1 groups fall 17 to 170 times at the 1st drought (the other two are barely expressed), and 5000 represses them too, although the clock genes of the two cultivars move in opposite directions in that contrast.
- **Its further fall at the 3rd drought goes together with the rise of CCA1/LHY,** in both cultivars, while the five day and evening clock genes do not move. Either the 3rd-drought leaves were sampled a little earlier in the morning, or the 3rd drought raises CCA1/LHY. The data cannot tell which. So THI1 and CCA1/LHY count as one open signal, and not as two memory genes.

**5000** changes more between the droughts: 1,229 groups, 162 of them in the memory classes, with a broad further shutdown of translation and stronger ABA signalling. This is a real difference in the data, but its 1st-drought group is one plant, and refitting 8008 with a single D1 plant can produce a similar picture. With one D1 plant it stays exploratory.

**H1 in one line:** sugarcane leaves repeat their drought response faithfully, and on top of that the tolerant cultivar strengthens a set of seed-like drying-protection genes.

**Where:** class of every group `de/group/memory_classes_<cultivar>_tillering.tsv.gz`, counts per class `de/group/memory_class_counts_<cultivar>_tillering.tsv`; agreement of the two droughts `report_tables/h1_d1_vs_d3_concordance.tsv`; overlap with the published classes `memory_genes/memory_class_overlap_literature.tsv`; figures `plots/story/`, made by `../code/08_story_figures.R`. README sections 6 and 8.

---

## 5. Buds: the timing of the stress leaves a mark (H3)

Buds become the next plants, so this is where an inherited memory would have to sit. The clean comparison is **T vs M**: buds of mothers stressed 135 days before sampling against buds of mothers stressed 6 days before, on the same run and the same day.

| Cultivar | T vs M (higher in T / higher in M) | Higher in T buds | Higher in M buds |
|---|---|---|---|
| 8008 | **261 / 146** | PAL (5 groups) and cinnamyl alcohol dehydrogenase (lignin), expansins, aquaporins, potassium transporters, B-box and NF-Y transcription factors | Raffinose synthase (6 groups), thiamine synthesis (THI1, 5 groups), circadian rhythm, phosphate-sensing SPX proteins |
| 5000 | **34 / 113** | Aquaporins, a chlorophyll a/b-binding protein, clock genes | PAL (8 groups), tubulins, thiamine-pyrophosphate enzymes, LOB and TCP transcription factors, malic enzyme, callose synthase |

![PAL in buds](plots/story/bud_pal.png)

*Figure 4. PAL expression in buds. Sugarcane has 56 groups annotated as phenylalanine ammonia-lyase (PAL). For each bud library we added up the normalised read counts of all 56, so each dot is the total PAL expression of one plant (y = log2 of that sum). The bar is the mean of the group. The grey band marks the watered buds, which were sequenced on another run and cannot be compared with T and M.*

- **PAL leads the list in both cultivars, in opposite directions.** PAL is the entry enzyme of the phenylpropanoid pathway, which makes lignin, flavonoids and a route to salicylic acid. It is higher in T buds in 8008 and higher in M buds in 5000. So how long ago the stress happened changes the phenylpropanoid state of the bud, and each cultivar does it its own way.
- **8008 T buds keep part of the mother's leaf drought programme.** Genes that 8008 mother leaves switch on under drought, and the genes that rise more at the 3rd drought, are slightly higher in T buds than in M buds (FDR 0.01 and 0.001), although T buds are 135 days past their last drought. That is what a long-term mark would look like. 5000 does not show it.
- **M buds look like a tissue just out of stress:** raffinose-family sugars and thiamine are typical of it, which fits 6 days of recovery.

**What is still hidden:** whether stressed buds differ from never-stressed buds. That comparison shows thousands of DEGs, but the watered buds' libraries were prepared differently (shorter fragments, three times more adapter) and sequenced on another run. One hint survived the correction attempt: after removing the run factor, the leaf drought genes come out higher in 8008 stressed buds than in watered buds. This is a loose hint and is not reported as a result (README section 5.3).

**H3 in one line:** the bud transcriptome records when the mother was stressed, through phenylpropanoid genes, differently in each cultivar.

**Where:** `de/group/de_<cultivar>_T_vs_M.tsv.gz`; GO and KEGG in `enrichment/`; leaf drought genes tested in buds `h5_rank/h5_rank_tests.tsv`; bud run effect and the correction attempt `bud_batch/`. README sections 9, 5.3 and 11.

---

## 6. Propagules: a milder drought response when the mother was stressed (H4)

Propagules of stressed mothers (Ft, Fm) and of watered mothers (Fw) went through one drought. Their own drought response is the expected one: lignin, wax, proline and trehalose synthesis up; photosynthesis, chloroplast and chlorophyll synthesis down (8008: 1,622 up / 1,072 down; 5000: 2,672 / 2,851).

Against the propagules of watered mothers under the same drought:

| Contrast | DEGs (up / down) | Higher | Lower |
|---|---|---|---|
| 8008 Ft_D1 vs Fw_D1 | 136 / 93 | Photosystem I assembly, chloroplast organisation, glucose metabolism; raffinose synthase, ClpB chaperones | Defence proteins (glucanases, chitinase, thaumatin), HSP90, asparagine and glutamine synthetases |
| 5000 Ft_D1 vs Fw_D1 | 231 / 145 | Small heat-shock proteins and other chaperones; chlorophyll a/b-binding and light-harvesting proteins; response to high light; translation | ABC transporters, PP2C phosphatases, wax synthesis, MYB transcription factors |
| 5000 Fm_D1 vs Fw_D1 | 49 / 36 | Light harvesting in photosystem I, ferredoxin-NADP reductase, chaperonins | Aminotransferases, wax synthesis |

To measure the size of the response fairly, we chose the genes in other plants, the mother leaves, and then followed them in the propagules:
- **Induced genes:** the groups that are DEGs **up** in mother leaves at the 1st drought (`Mt_D1_vs_Mt_W`): 3,905 groups in 8008 and 2,091 in 5000.
- **Repressed genes:** the groups that are DEGs **down** in the same contrast: 4,196 in 8008 and 1,750 in 5000.

For each propagule group under drought (Fw_D1, Ft_D1, Fm_D1) we then took the log2FC of every one of these genes against the same watered propagules (Fw_W).

![Propagule response](plots/story/propagule_response.png)

*Figure 5. Drought response of propagules on the mother-leaf drought genes. Left panels: the induced genes; right panels: the repressed genes. y = log2FC of each gene in the propagule group against watered propagules (Fw_W). Each box holds the middle half of the genes and the label is the median. A median of +1 means that the typical gene of the list doubles under drought, −1 that it halves, 0 that it does not respond.*

- **Induced genes are induced in every group, by a similar amount.** The typical induced gene rises 2.1 to 2.4 times in 8008 (median log2FC +1.06 in Fw, +1.27 in Ft) and 2.6 to 3.3 times in 5000 (+1.73 in Fw, +1.39 in Ft, +1.46 in Fm). The mother's history does not change this side in a consistent direction.
- **Repressed genes are repressed less** in propagules of stressed mothers. In 8008 the typical repressed gene falls to 48% of its watered level in Fw and only to 68% in Ft:

  | Cultivar | Fw (watered mother) | Ft (T mother) | Fm (M mother) |
  |---|---|---|---|
  | 8008 | −1.07 | −0.55 | not sequenced |
  | 5000 | −1.71 | −1.21 | −1.26 (3 plants) |

  *Median log2FC of the mother-repressed genes.*
- **The two cultivars agree.** The Ft vs Fw profiles of 8008 and 5000 correlate (rho 0.20; set tests FDR 10⁻¹¹ to 10⁻²⁵), and gene-set tests give the same direction in both (FDR 10⁻⁶ in 8008, 7 × 10⁻⁹ in 5000).
- **The spared genes are photosynthesis and growth genes.** "Spared" means genes that drought switches down in propagules of watered mothers and switches down less, or not at all, in propagules of stressed mothers. They are the "Higher" column of the table above: photosystem I assembly, photosynthetic electron transport and chloroplast organisation in 8008, light-harvesting proteins and translation in 5000 (KEGG "Photosynthesis", padj 0.003).
- **A co-expression analysis points the same way, with low power.** WGCNA groups genes that rise and fall together across the 32 leaf libraries into 13 modules, without using the labels of the samples. We then asked which modules follow the mother's history.
  - One does: **M07**, 258 groups of translation, ribosome and chloroplast organisation genes. Drought switches it down, and less so in propagules of stressed mothers.
  - It is significant in 8008 only (FDR 0.04). 5000 goes the same direction without significance (t = 1.3), and the two cultivars together give FDR 0.16.
  - In both cultivars M07 is rich in the genes that are higher in Ft_D1 than in Fw_D1 (odds ratio 5.1 in 8008, 7.1 in 5000).
  - The limits: 32 libraries are few for a network (the usual quality criterion, a scale-free fit of 0.8, was not reached; the best was 0.76), and M07 depends on the same single Fw_D1 plants. So it is the same result seen with another method, and not an independent confirmation.

This fits the physiology: Marcos et al. 2018 found that propagules of stressed mothers kept better photosynthesis under drought.

**What it rests on.** The reference, Fw_D1, is one plant per cultivar, and in 8008 it is the failed library `79_S6`. A thin library makes many genes look repressed, which would exaggerate the difference in 8008. 5000 is the stronger case: its reference is a normal library and its three Fm plants agree with Ft. A reference plant that was simply drier than the others would give the same picture, and gene expression alone cannot tell these two apart.

**H4 in one line:** in both cultivars, propagules of drought-stressed mothers hold on to photosynthesis and growth genes under drought.

**Where:** medians of Figure 5 `report_tables/h4_propagule_response.tsv`; contrasts `de/group/de_<cultivar>_Ft_D1_vs_Fw_D1.tsv.gz` and `…Fm_D1_vs_Fw_D1.tsv.gz`; set tests and agreement of the cultivars `h5_rank/h5_rank_tests.tsv`, `h5_rank/contrast_correlations.tsv`, `memory_genes/camera_literature_sets.tsv`; modules `wgcna/module_trait_tests.tsv`, `wgcna/module_list_overlaps.tsv`, `wgcna/history_modules.png`. README sections 10 and 13.

---

## 7. The other hypotheses: too weak for now

| # | Hypothesis | Where it stands | What it needs |
|---|---|---|---|
| H2 | The tolerant 8008 remembers more than 5000 | 5000 changes more between droughts, on one D1 plant. 8008 responds more strongly at the 1st drought. | More 5000 D1 plants: not possible with these samples |
| H5 | The same genes carry the mark from leaf to bud to propagule | A leaf → bud link in 8008 (section 5). No bud → propagule link detected. | Watered propagules of stressed mothers (T/W, M/W), never sequenced |
| H6 | Known memory genes are memory genes in sugarcane | They are drought genes here too; RAB18, LEA3 and PIT1 keep their memory behaviour (section 4). | A side question, answered as far as these data go |

Two checks support the main results: gene copies within a group move together, so no single copy hides a memory, and co-expression modules follow current drought.

**Where:** H2 `h2_genotype_cycle/`; H5 `h5_rank/`; H6 `memory_genes/`; gene copies `alleles/`; modules `wgcna/`. README sections 7, 8, 11, 12 and 13.

---

## 8. A side result: non-coding RNAs point to the same bud difference

Everything above uses the annotated protein-coding genes. As a side analysis we searched the whole genome for long non-coding RNAs (lncRNAs), transcripts that do not make a protein, and ran the same drought tests on them.

![lncRNA search](plots/story/lncrna.png)

*Figure 6. The lncRNA search. **A.** The filters, from all transcripts assembled from the reads to the candidates: a transcript must be new (outside the annotated genes), at least 200 nt long, and show no sign of coding for a protein. Copies of the same sequence are then joined in clusters, the unit we test. **B.** The one sequence family that stands out. Each dot is one bud library; y is the log2 of the summed expression (TPM) of the whole family; the bar is the mean. The grey band marks the watered buds, sequenced on another run.*

- **lncRNAs respond to drought like coding genes do.** In 8008 leaves, 1,327 clusters change at the 1st drought; in 5000, 1,206 at the 3rd. 459 change in the same direction in both cultivars and 22 in opposite directions.
- **They add no evidence for leaf memory.** Only two clusters differ between the 3rd and the 1st drought in 8008.
- **In buds they repeat H3.** 18 lncRNAs differ between T and M buds in both cultivars, in the same direction, where chance predicts 0.2. Of these, 13 are copies of one sequence family.
- **The family is higher in M buds, the ones stressed 6 days before.** It has 37 clusters and 815 copies in the genome. It is about 7 times lower in T than in M buds in 8008 and 4 times in 5000, and in M buds it reaches about 1% of all transcripts in 8008. In 8008 leaves drought represses it (about 14 times).
- **We could not name it.** It matches only unannotated sugarcane DNA at NCBI (R570, SP80-3280, *S. officinarum*), and no known protein, repeat or RNA family. It is also transcribed in SP80-3280. It is specific to sugarcane and, with so many copies, probably derived from a repeat.
- **Three of the other five are small protein genes missing from the annotation.** Two are copies of a glycine-rich cell wall protein, called cold and drought-regulated (CORA) in maize, higher in M buds. One is a wound-induced protein, higher in T buds.

**In one line:** with non-coding RNAs, on transcripts the main analysis never saw, T and M buds again differ in both cultivars. This supports H3, but the function of the family is unknown.

**Where:** `lncrna_genome/lnc_deg_summary.tsv` (counts per contrast), `lncrna_genome/bud_shared_lncRNAs.tsv` (the 18), `lncrna_genome/bud_family_clusters.tsv` and `family_tpm_by_library.tsv` (the family), `lncrna_genome/family_identity.md` (database searches). README section 13.4.

---

## 9. Where this leaves us

There is no budget for new sequencing or new experiments, and the lab has already answered what it could about the sampling. These 48 libraries are the data, so each result is stated with the confidence the design allows:

| Result | Confidence | What limits it, and cannot be changed |
|---|---|---|
| H1 in 8008: drying-protection genes rise more at the 3rd drought | The firmest result: 2–3 plants per group, on runs 51/52 that behave as one batch | 30 days of age between D1 and D3; hour of sampling |
| H3: T and M buds differ, led by PAL | Firm: same run, same sampling day | 2–3 plants per group; hour of sampling; says nothing about stressed vs never-stressed buds |
| H4: milder response in propagules of stressed mothers | Exploratory, and seen in both cultivars | One Fw_D1 plant per cultivar; no measure of how dry each plant was |
| H1 in 5000, H2 | Exploratory | One D1 plant |
| Stressed vs watered buds | Not answerable | Run and library preparation are the same thing as treatment |
| lncRNAs: one sugarcane-specific family differs between T and M buds (section 8) | Firm as a difference, in both cultivars, same run | Its function is unknown |
| THI1 and CCA1/LHY as memory | Open, and narrowed: THI1 is repressed by drought for certain; its further fall at the 3rd drought follows CCA1/LHY | A small difference in sampling hour and a real rise of CCA1/LHY look the same |

**What can still be done with the data in hand:** use the physiology that already exists. Marcos et al. 2018 describe the physiology of such propagules; if the lab kept water status or gas exchange per plant, it would say whether the Fw_D1 reference plants were simply drier.

---

*Everything behind this summary (the other figures, the checks, the methods and the file locations) is in the project [`README.md`](../../README.md).*
