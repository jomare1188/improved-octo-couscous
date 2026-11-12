# Drought-memory gene lists: sources

Downloaded 2026-10-29. The supplementary files came from Europe PMC, `literature/PMC*/`.

## Memory classes from repeated dehydration (Avramova lab)

The class is the response at the 1st stress compared with water, then the change at the 3rd stress compared with the 1st:
- `++` / `--`: memory genes; the 3rd response is stronger.
- `+-` / `-+`: revised-response memory genes.
- `+=` / `-=`: non-memory genes; same response each time.

| Species | Paper | File used | Genes per class (unique) |
|---|---|---|---|
| *Arabidopsis thaliana* | Ding Y, Liu N, Virlouvet L, Riethoven JJ, Fromm M, Avramova Z (2013). Four distinct types of dehydration stress memory genes in *Arabidopsis thaliana*. BMC Plant Biol 13:229. doi:10.1186/1471-2229-13-229 | Additional file 4 (`1471-2229-13-229-S4.xls`) | ++ 362, -- 310, +- 857, -+ 434, += 2177, -= 2439 |
| *Zea mays* (B73 v3 IDs) | Ding Y, Virlouvet L, Liu N, Riethoven JJ, Fromm M, Avramova Z (2014). Dehydration stress memory genes of *Zea mays*; comparison with *Arabidopsis thaliana*. BMC Plant Biol 14:141. doi:10.1186/1471-2229-14-141 | Additional file 4 (`1471-2229-14-141-S4.xlsx`) | ++ 159, -- 70, +- 531, -+ 48, += 933, -= 298 |

- `code/extract_ding_lists.R` writes `ding_memory_classes.tsv`, which keeps FPKM and log2FC from the papers.
- Maize v3 IDs are converted to the B73 v5 IDs used in our OrthoFinder run with MaizeGDB's `B73v3_to_B73v5.tsv` (`literature/maizegdb/`).

## Curated candidates (`curated_candidates.tsv`)

40 Arabidopsis genes, each ID checked against the TAIR10 symbols in the Ensembl proteome. They cover:
- the classic trainable markers: RD29B and RAB18 (`+/+`), with RD29A and COR15A as non-memory controls (Ding et al. 2012, Nat Commun 3:740);
- ABA synthesis and signalling: NCED3, the clade A PP2Cs, OST1/SnRK2.6 and the ABFs;
- LEA proteins and dehydrins;
- osmolytes: proline P5CS1, galactinol synthase and trehalose. Leaf proline was higher in propagules from drought-exposed mothers in Marcos et al. 2018, PLoS One 13:e0206716;
- aquaporin PIP2;1;
- stress-memory regulators: MYC2, HSFA2, BRM, FGT1 and ATX1;
- the DNA methylation and RdDM machinery: MET1, CMT3, DRM2, DDM1, ROS1, DME, AGO4, NRPD1, DCL3 and RDR2. These matter for the question of whether memory passes to the propagules.

## Sugarcane context

- Marcos FCC et al. (2018), PLoS One 13:e0206716, is the physiological basis of this experiment. It found that propagules from mother plants put through 3 water-deficit cycles performed better under drought.
- No published sugarcane drought-memory gene list was found, so sugarcane genes are reached through orthology. The Arabidopsis and maize genes are mapped to POJ2878 allele groups with OrthoFinder (`../orthology/`).
