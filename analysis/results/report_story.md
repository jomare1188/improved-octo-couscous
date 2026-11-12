# Does sugarcane remember a drought? The DOLORES RNA-seq story, in short

*November 2026. This is the short version: the main insights, what was done with the problematic samples, and which comparisons can be trusted. The full documentation, with every analysis, figure, method and concept explained, is the project [`README.md`](../../README.md); section numbers in brackets point to it.*

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

The sample table from the lab is correct, so the cultivar mix-up happened when the tissue was sampled or the tubes were handled. Details and the evidence for each library are in the README (section 5.1).

---

## 3. Which contrasts can be biased by sequencing run or sampling time

The libraries were sequenced on several runs and the leaves were sampled on different days. Some contrasts compare groups that differ in run or day as well as in treatment.

| Contrast | Sequencing run | Sampling day | Can it be trusted? |
|---|---|---|---|
| Leaf D1 vs watered | drought leaves on runs 51/52, watered on 61/63/73 | day 9 vs day 39 | **Biased by run and by 30 days of age.** The drought response is large, so its main features are safe. Small differences are not. |
| Leaf D3 vs watered | same split of runs | both day 39 | **Biased by run.** Same reading. |
| **Leaf D3 vs D1** (H1) | both on runs 51/52 | day 39 vs day 9 | **Clean for run.** It carries 30 days of age, which no design can remove: the 3rd drought always comes later. |
| Bud T or M vs watered | stressed buds on run 45, watered on run 44, with a different library preparation | all day 180 | **Not usable.** Run and treatment are the same thing here, and no correction separates them. |
| **Bud T vs M** (H3) | both on run 45 | both day 180 | **Clean.** The stress happened 135 days (T) or 6 days (M) before, which is the comparison itself. |
| **Propagules** (H4) | all on runs 51/52; the reference plant is on the other run from the plants compared with it | same drought cycle | **Clean for run:** runs 51 and 52 behave as one batch, and the arrangement is reversed between the cultivars. Exploratory because of the single reference plant. |
| Maturation leaves | runs 51/52 | D1 day 144, D3 and watered day 174 | Descriptive only: one plant per drought group, 5000 only. |

**The hour of sampling is not recorded for any group.** It matters because clock genes (CCA1/LHY, PRR7, TOC1, GI) and the clock-controlled THI1 differ between D1 and watered leaves in 8008, and two of them appear among the memory genes. Drought can shift the clock, but so can sampling at a different hour (README sections 5.2, 5.3 and 5.5).

---

## 4. Leaves: among a sea of equal responses, a set of genes remembers (H1)

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

### 4.2 The 3rd drought replays the 1st, and some genes stand out

![D1 vs D3 response](plots/story/d1_vs_d3.png)

*Figure 2. Each dot is a group that responds to drought. x = response at the 1st drought, y = response at the 3rd. On the dashed diagonal the two are equal. Coloured dots are the memory classes.*

In 8008 the two responses correlate at **r = 0.90**, and of 5,674 groups that respond at both droughts, 5,673 respond in the same direction. This reproducibility is the sea against which the memory genes are visible:
- **84 groups** respond differently at the 3rd drought (35 higher, 49 lower).
- **32 are memory genes** in the sense of Ding et al.: 22 respond to drought and respond more strongly the 3rd time (`++`), 10 are repressed and more repressed the 3rd time (`--`).
- **1,889 more respond at the 3rd drought only** (flavonoid and phenylpropanoid synthesis up; circadian rhythm and protein folding down).

The `++` genes are **drying-protection genes of the kind seeds use**:

| Gene | D1 vs W (log2FC) | D3 vs D1 (log2FC) |
|---|---|---|
| RAB18/XERO1 dehydrin | +10.4 | +1.4 |
| LEA3/SAG21 | +2.5 | +1.4 |
| another LEA protein | +7.8 | +2.0 |
| seed-storage protein | +7.1 | +4.9 |
| embryo-specific protein 3 | +9.3 | +1.8 |
| sugar transporter | +9.0 | +2.3 |

RAB18 is the textbook trainable gene of *Arabidopsis*, and its maize orthologue is a memory gene too. Our `++` class shares RAB18, LEA3 and PIT1 with the *Arabidopsis* `++` class, far more than chance (odds ratio 25, FDR 0.003).

![Example genes](plots/story/example_genes.png)

*Figure 3. Four genes plant by plant in tillering mother leaves (y = log2 expression). Left 8008, right 5000. RAB18 and LEA3 keep rising at D3 in 8008; HAI1 rises and stays; THI1 drops and drops further.*

On the repressed side, thiamine synthesis (THI1) is shut down further at the 3rd drought. THI1 follows the clock, so the hour of sampling has to be known before it is called memory (section 3).

**5000** changes more between the droughts: 1,229 groups, 162 of them in the memory classes, with a broad further shutdown of translation and stronger ABA signalling. This is a real difference in the data, but its 1st-drought group is one plant, and refitting 8008 with a single D1 plant can produce a similar picture. It becomes a result when more 5000 D1 plants are sequenced.

**H1 in one line:** sugarcane leaves repeat their drought response faithfully, and on top of that the tolerant cultivar strengthens a set of seed-like drying-protection genes, the same kind of memory described in *Arabidopsis* and maize.

---

## 5. Buds: the timing of the stress leaves a mark (H3)

Buds become the next plants, so this is where an inherited memory would have to sit. The clean comparison is **T vs M**: buds of mothers stressed 135 days before sampling against buds of mothers stressed 6 days before, on the same run and the same day.

| Cultivar | T vs M (higher in T / higher in M) | Higher in T buds | Higher in M buds |
|---|---|---|---|
| 8008 | **261 / 146** | PAL (5 groups) and cinnamyl alcohol dehydrogenase (lignin), expansins, aquaporins, potassium transporters, B-box and NF-Y transcription factors | Raffinose synthase (6 groups), thiamine synthesis (THI1, 5 groups), circadian rhythm, phosphate-sensing SPX proteins |
| 5000 | **34 / 113** | Aquaporins, a chlorophyll a/b-binding protein, clock genes | PAL (8 groups), tubulins, thiamine-pyrophosphate enzymes, LOB and TCP transcription factors, malic enzyme, callose synthase |

![PAL in buds](plots/story/bud_pal.png)

*Figure 4. Total expression of the 56 PAL groups in buds. The grey band marks the watered buds, sequenced on another run.*

- **PAL leads the list in both cultivars, in opposite directions.** PAL is the entry enzyme of the phenylpropanoid pathway, which makes lignin, flavonoids and a route to salicylic acid. It is higher in T buds in 8008 and higher in M buds in 5000. So how long ago the stress happened changes the phenylpropanoid state of the bud, and each cultivar does it its own way.
- **8008 T buds keep part of the mother's leaf drought programme.** Genes that 8008 mother leaves switch on under drought, and the genes that rise more at the 3rd drought, are slightly higher in T buds than in M buds (FDR 0.01 and 0.001), although T buds are 135 days past their last drought. That is what a long-term mark would look like. 5000 does not show it.
- **M buds look like a tissue just out of stress:** raffinose-family sugars and thiamine are typical of it, which fits 6 days of recovery.

**What is still hidden:** whether stressed buds differ from never-stressed buds. That comparison shows thousands of DEGs, but the watered buds' libraries were prepared differently (shorter fragments, three times more adapter) and sequenced on another run. One hint survived the correction attempt: after removing the run factor, the leaf drought genes come out higher in 8008 stressed buds than in watered buds. It needs a re-sequencing of a few buds together to become a result (README section 5.3).

**H3 in one line:** the bud transcriptome records when the mother was stressed, through phenylpropanoid genes, differently in each cultivar.

---

## 6. Propagules: a milder drought response when the mother was stressed (H4)

Propagules of stressed mothers (Ft, Fm) and of watered mothers (Fw) went through one drought. Their own drought response is the expected one: lignin, wax, proline and trehalose synthesis up; photosynthesis, chloroplast and chlorophyll synthesis down (8008: 1,622 up / 1,072 down; 5000: 2,672 / 2,851).

Against the propagules of watered mothers under the same drought:

| Contrast | DEGs (up / down) | Higher | Lower |
|---|---|---|---|
| 8008 Ft_D1 vs Fw_D1 | 136 / 93 | Photosystem I assembly, chloroplast organisation, glucose metabolism; raffinose synthase, ClpB chaperones | Defence proteins (glucanases, chitinase, thaumatin), HSP90, asparagine and glutamine synthetases |
| 5000 Ft_D1 vs Fw_D1 | 231 / 145 | Small heat-shock proteins and other chaperones; chlorophyll a/b-binding and light-harvesting proteins; response to high light; translation | ABC transporters, PP2C phosphatases, wax synthesis, MYB transcription factors |
| 5000 Fm_D1 vs Fw_D1 | 49 / 36 | Light harvesting in photosystem I, ferredoxin-NADP reductase, chaperonins | Aminotransferases, wax synthesis |

To measure the size of the response fairly, we took the genes that drought induces or represses in the *mother leaves* and followed them in each propagule group.

![Propagule response](plots/story/propagule_response.png)

*Figure 5. Drought response of propagules on mother-leaf drought genes. Box = middle half of the genes, label = median log2FC. Left: genes drought induces; right: genes drought represses.*

- **Induced genes** respond about the same in every group.
- **Repressed genes are repressed less** in propagules of stressed mothers:

  | Cultivar | Fw (watered mother) | Ft (T mother) | Fm (M mother) |
  |---|---|---|---|
  | 8008 | −1.07 | −0.55 | not sequenced |
  | 5000 | −1.71 | −1.21 | −1.26 (3 plants) |

  *Median log2FC of the mother-repressed genes.*
- **The two cultivars agree.** The Ft vs Fw profiles of 8008 and 5000 correlate (rho 0.20; set tests FDR 10⁻¹¹ to 10⁻²⁵), and gene-set tests give the same direction in both (FDR 10⁻⁶ in 8008, 7 × 10⁻⁹ in 5000).
- **The genes spared are photosynthesis and growth genes.** A co-expression analysis finds the same on its own: the one module that follows the mother's history is a translation and chloroplast module, less repressed in propagules of stressed mothers.

This fits the physiology: Marcos et al. 2018 found that propagules of stressed mothers kept better photosynthesis under drought.

**What it rests on.** The reference, Fw_D1, is one plant per cultivar, and in 8008 it is the failed library `79_S6`. A thin library makes many genes look repressed, which would exaggerate the difference in 8008. 5000 is the stronger case: its reference is a normal library and its three Fm plants agree with Ft. A reference plant that was simply drier than the others would give the same picture, so a measure of plant water status is the missing piece.

**H4 in one line:** in both cultivars, propagules of drought-stressed mothers hold on to photosynthesis and growth genes under drought.

---

## 7. The other hypotheses: too weak for now

| # | Hypothesis | Where it stands | What it needs |
|---|---|---|---|
| H2 | The tolerant 8008 remembers more than 5000 | 5000 changes more between droughts, on one D1 plant. 8008 responds more strongly at the 1st drought. | More 5000 D1 plants |
| H5 | The same genes carry the mark from leaf to bud to propagule | A leaf → bud link in 8008 (section 5). No bud → propagule link detected. | Watered propagules of stressed mothers (T/W, M/W), never sequenced |
| H6 | Known memory genes are memory genes in sugarcane | They are drought genes here too; RAB18, LEA3 and PIT1 keep their memory behaviour (section 4). | A side question, answered as far as these data go |

Two checks support the main results: gene copies within a group move together, so no single copy hides a memory, and co-expression modules follow current drought (README sections 12 and 13).

---

## 8. What would make these results firm

**Ask the lab (free):**
1. **Hour of sampling** for every group: decides whether THI1 and the clock genes are memory or time of day.
2. **How stressed each sampled plant was** (leaf water potential or relative water content, soil moisture, gas exchange): separates memory from "this plant was drier", for 5000 leaves and for the propagules.
3. **How the watered buds' libraries were prepared,** and whether tissue or RNA is left.

**In the lab, if there is budget:**
1. **Re-prepare and re-sequence a few buds together** (two W, two T, two M per cultivar): opens the stressed vs never-stressed comparison of H3.
2. **Add plants where one carries everything:** Fw_D1 in both cultivars (H4), 5000 Mt_D1 (H1, H2).
3. **RT-qPCR in every available plant:** RAB18, LEA3 and the seed-storage protein (H1); PAL (H3); photosynthesis genes (H4); THI1 and CCA1/LHY (memory or clock).

---

*Everything behind this summary (the other figures, the checks, the methods and the file locations) is in the project [`README.md`](../../README.md).*
