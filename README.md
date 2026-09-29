# Exposure of Sugarcane Plants to Drought Cycles at Different Phenological Stages

## Experimental Design

Three-month-old sugarcane seedlings of genotypes **IACSP95-5000** (productive material) and **IACCTC07-8008** (drought-tolerant material) developed by the Sugarcane Breeding Program of the Agronomic Institute (ProCana, IAC, Brazil) were transplanted into 20 L pots containing red-yellow latosol and fertilized. Cultivation was carried out in a greenhouse with daily irrigation.

## Treatment Groups

At five months of age (from transplant date), plants were divided into three experimental groups per genotype (total of six groups).

### a) W (Well-watered)
- **Treatment**: Maintained under daily irrigation
- **Purpose**: Control group

### b) T (Tillering)  
- **Treatment**: Subjected to three drought cycles during the tillering phase
- **Purpose**: Assess drought stress during vegetative growth

### c) M (Maturation)
- **Treatment**: Subjected to three drought cycles during the maturation phase  
- **Purpose**: Assess drought stress during reproductive development

## Drought Cycle Protocol

**Drought induction**: Irrigation suspension for **9 days**

**Recovery period**: Plants re-irrigated for **6 days** after each drought interval

**Total cycles**: 3 drought/recovery cycles per treatment group

This experimental design allows for the evaluation of genotype-specific responses to drought stress at two different developmental stages in sugarcane.

![Drought cycle 1](figs/dry_cycle_1.png)

*Figure 1. Schematic showing drought cycles in two sugarcane genotypes (IACCTC08-8008 and IACSP95-5000) and in two phenological phases (tillering and maturation). W (control), plants maintained with irrigation throughout the experiment; T, plants subjected to three drought cycles (D1, D2, D3) during the tillering phase; M, plants subjected to three drought cycles (D1, D2, D3) during the maturation phase. Green arrows indicate leaf sample collection (34 of 36 samples in `RAWDATA/`) and red arrows indicate bud collection (18 of 18 samples in `RAWDATA/`).*

# Propagule Experiment

At the end of the experiment, stem cuttings were collected and new seedlings were vegetatively propagated from the above experimental groups to carry out the propagule experiment where we evaluated the "heritability" of drought memory.

Three-month-old sugarcane seedlings of genotypes IACSP95-5000 and IACCTC07-8008 obtained through vegetative propagation from groups W, T and M were transplanted into 12 L pots containing red-yellow latosol and fertilized. Cultivation was carried out in a greenhouse with daily irrigation.

At two months of age (from transplant), plants were divided and drought treatment was initiated in both genotypes (total of 12 groups):

## Treatment Groups

### a) W/D
- **Origin**: Propagules obtained from non-stressed plants
- **Treatment**: Subjected to one drought cycle

### b) W/W  
- **Origin**: Propagules obtained from non-stressed plants
- **Treatment**: Well-watered

### c) T/D
- **Origin**: Propagules obtained from plants stressed during tillering
- **Treatment**: Subjected to one drought cycle

### d) T/W
- **Origin**: Propagules obtained from plants stressed during tillering  
- **Treatment**: Well-watered

### e) M/D
- **Origin**: Propagules obtained from plants stressed during maturation
- **Treatment**: Subjected to one drought cycle

### f) M/W
- **Origin**: Propagules obtained from plants stressed during maturation
- **Treatment**: Well-watered

> **Note:** only groups a) W/D, b) W/W, c) T/D and e) M/D were sequenced (8 groups, 24 samples). There are no RNA-seq samples for T/W and M/W.

## Drought Cycle Protocol

**Drought induction**: Irrigation suspension for **9 days**

**Recovery period**: Plants re-irrigated after the cycle completion (**6 days** recovery)

![Drought cycle 2](figs/dry_cycle_2.png)

*Figure 2. Experimental scheme of sugarcane propagules from two genotypes (IACCTC08-8008 and IACSP95-5000) subjected to one drought cycle. Green arrows indicate leaf sample collection for transcriptomic analysis*

# So, what do we have?

The original metadata ([experimental design](Esquema_experimental_RNAseq.pdf)) and the sample spreadsheet `arquivos_metaboloma.xlsx` describe **26 groups and 78 RNA-seq samples** (3 biological replicates per group) in three sample sets. The leaf samples labelled `Mt`/`Mm` belong to the **mother plants of experiment 1** (green arrows in Figure 1), not to the propagule experiment. Only the propagule leaves (`Fw`/`Ft`/`Fm`) come from experiment 2, and only 4 of its 6 groups per genotype were sequenced.

We have **69 of the 78 samples** in `RAWDATA/`: 51 we already had, plus 18 downloaded from Dolores's Google Drive on 2026-09-29 (MD5s checked, see `RAWDATA/drive_provenance.tsv`). The other 9 were never available. Missing samples and open doubts are summarized in [`INVENTORY.md`](INVENTORY.md).

| Sample set | Groups | Planned samples | In `RAWDATA/` | Missing (not available) |
|---|---|---|---|---|
| Mother plants, buds (end of experiment 1) | 6 | 18 | **18** | 0 |
| Mother plants, leaves (experiment 1, during the drought cycles) | 12 | 36 | **34** | 2 |
| Propagules, leaves (experiment 2) | 8 | 24 | **17** | 7 |
| **Total** | **26** | **78** | **69** | **9** |

## Replicates per group in `RAWDATA/`

Genotype 5000 = IACSP95-5000, 8008 = IACCTC07-8008. W = well-watered (C in the spreadsheet), D = drought, D1/D3 = sampled at the 1st/3rd drought cycle. Every group has at least 2 replicates.

**Mother buds**

| Genotype | W | Tillering D | Maturation D |
|---|---|---|---|
| 5000 | 3/3 | 3/3 | 3/3 |
| 8008 | 3/3 | 3/3 | 3/3 |

**Mother leaves**

| Genotype | Tillering W (Mt_W) | Tillering D1 (Mt_D1) | Tillering D3 (Mt_D3) | Maturation W (Mm_W) | Maturation D1 (Mm_D1) | Maturation D3 (Mm_D3) |
|---|---|---|---|---|---|---|
| 5000 | 3/3 | 3/3 | 3/3 | 3/3 | 3/3 | 3/3 |
| 8008 | 3/3 | 3/3 | 3/3 | 3/3 | 2/3 | 2/3 |

**Propagule leaves (mother group / propagule treatment)**

| Genotype | W/W (Fw_W) | W/D (Fw_D1) | T/D (Ft_D1) | M/D (Fm_D1) |
|---|---|---|---|---|
| 5000 | 2/3 | 2/3 | 2/3 | 3/3 |
| 8008 | 2/3 | 2/3 | 2/3 | 2/3 |

The T/W and M/W propagule groups were not sequenced, so the propagule set can compare drought-stressed propagules from stressed and non-stressed mothers, but it has no watered controls for propagules from stressed mothers.

**Before analysis:** the `159_S5` pair is out of sync (R2 has 100,000 fewer reads than R1) and must be re-paired. See [`INVENTORY.md`](INVENTORY.md) for this and the other open doubts.

## Sample table

Full details, with read counts and notes, are in [`sample_inventory.csv`](sample_inventory.csv). Group names follow the nf-core sample sheet: Mt = mother, tillering; Mm = mother, maturation; Fw / Ft / Fm = propagule from a well-watered / tillering-stressed / maturation-stressed mother.

| Sample | Set | Genotype | Group | Replicate | In `RAWDATA/` | Source |
|---|---|---|---|---|---|---|
| F25C1_S6 | Mother bud | 5000 | 5000_Maturation_W | 1 | YES | original |
| F25C2_S7 | Mother bud | 5000 | 5000_Maturation_W | 2 | YES | original |
| F25C3_S8 | Mother bud | 5000 | 5000_Maturation_W | 3 | YES | original |
| FIDH15000_S27 | Mother bud | 5000 | 5000_Tillering_D | 1 | YES | original |
| FIDH25000_S28 | Mother bud | 5000 | 5000_Tillering_D | 2 | YES | original |
| FIDH35000_S29 | Mother bud | 5000 | 5000_Tillering_D | 3 | YES | original |
| 21DH5_S21 | Mother bud | 5000 | 5000_Maturation_D | 1 | YES | original |
| 22DH5_S22 | Mother bud | 5000 | 5000_Maturation_D | 2 | YES | original |
| 23DH5_S23 | Mother bud | 5000 | 5000_Maturation_D | 3 | YES | original |
| F28C1_S9 | Mother bud | 8008 | 8008_Maturation_W | 1 | YES | original |
| F28C2_S10 | Mother bud | 8008 | 8008_Maturation_W | 2 | YES | original |
| F28C3_S11 | Mother bud | 8008 | 8008_Maturation_W | 3 | YES | original |
| FIDH1800_S30 | Mother bud | 8008 | 8008_Tillering_D | 1 | YES | original |
| FIDH28000_S31 | Mother bud | 8008 | 8008_Tillering_D | 2 | YES | original |
| FIDH38000_S32 | Mother bud | 8008 | 8008_Tillering_D | 3 | YES | original |
| 21DH8_S24 | Mother bud | 8008 | 8008_Maturation_D | 1 | YES | original |
| 22DH8_S25 | Mother bud | 8008 | 8008_Maturation_D | 2 | YES | original |
| 23DH8_S26 | Mother bud | 8008 | 8008_Maturation_D | 3 | YES | original |
| Dolores_249_S3 | Mother leaf | 5000 | 5000_Mt_W | 1 | YES | original |
| 250_S2 | Mother leaf | 5000 | 5000_Mt_W | 2 | YES | original |
| 251_S3 | Mother leaf | 5000 | 5000_Mt_W | 3 | YES | original |
| 157_S12 | Mother leaf | 5000 | 5000_Mt_D1 | 2 | YES | original |
| D158_S48 | Mother leaf | 5000 | 5000_Mt_D1 | 2 | YES | Drive (transcriptome_3) |
| 159_S5 | Mother leaf | 5000 | 5000_Mt_D1 | 3 | YES | Drive (transcriptome_3) |
| 221_S19 | Mother leaf | 5000 | 5000_Mt_D3 | 2 | YES | original |
| 222_S6 | Mother leaf | 5000 | 5000_Mt_D3 | 2 | YES | Drive (transcriptome_3) |
| 223_S42 | Mother leaf | 5000 | 5000_Mt_D3 | 3 | YES | original |
| 156_S22 | Mother leaf | 5000 | 5000_Mm_W | 1 | YES | original |
| 154_S28 | Mother leaf | 5000 | 5000_Mm_W | 2 | YES | original |
| 155_S11 | Mother leaf | 5000 | 5000_Mm_W | 3 | YES | original |
| 253_S41 | Mother leaf | 5000 | 5000_Mm_D1 | 1 | YES | original |
| 254_S7 | Mother leaf | 5000 | 5000_Mm_D1 | 2 | YES | Drive (transcriptome_2) |
| 255_S8 | Mother leaf | 5000 | 5000_Mm_D1 | 3 | YES | Drive (transcriptome_3) |
| 9_S15 | Mother leaf | 5000 | 5000_Mm_D3 | 1 | YES | Drive (transcriptome_1) |
| 10_S16 | Mother leaf | 5000 | 5000_Mm_D3 | 2 | YES | Drive (transcriptome_2) |
| 11_S1 | Mother leaf | 5000 | 5000_Mm_D3 | 3 | YES | original |
| Dolores_257_S5 | Mother leaf | 8008 | 8008_Mt_W | 1 | YES | original |
| 258_S5 | Mother leaf | 8008 | 8008_Mt_W | 2 | YES | original |
| Dolores_259_S6 | Mother leaf | 8008 | 8008_Mt_W | 3 | YES | original |
| 165_a_S40 | Mother leaf | 8008 | 8008_Mt_D1 | 2 | YES | original |
| 166_S23 | Mother leaf | 8008 | 8008_Mt_D1 | 2 | YES | original |
| 168_S25 | Mother leaf | 8008 | 8008_Mt_D1 | 3 | YES | original |
| 230_S15 | Mother leaf | 8008 | 8008_Mt_D3 | 2 | YES | original |
| D231_D46 | Mother leaf | 8008 | 8008_Mt_D3 | 2 | YES | Drive (transcriptome_2) |
| 232_S8 | Mother leaf | 8008 | 8008_Mt_D3 | 3 | YES | original |
| 161_S13 | Mother leaf | 8008 | 8008_Mm_W | 1 | YES | original |
| 162_S20 | Mother leaf | 8008 | 8008_Mm_W | 2 | YES | original |
| 163_S7 | Mother leaf | 8008 | 8008_Mm_W | 3 | YES | original |
| 261_S9 | Mother leaf | 8008 | 8008_Mm_D1 | 1 | YES | Drive (transcriptome_2) |
| 262_S10 | Mother leaf | 8008 | 8008_Mm_D1 | 2 | YES | Drive (transcriptome_2) |
| 263_S11 | Mother leaf | 8008 | 8008_Mm_D1 | 3 | **NO** | not available |
| 13_S17 | Mother leaf | 8008 | 8008_Mm_D3 | 1 | YES | Drive (transcriptome_2) |
| 14 | Mother leaf | 8008 | 8008_Mm_D3 | 2 | **NO** | not available |
| 15_S18 | Mother leaf | 8008 | 8008_Mm_D3 | 3 | YES | Drive (transcriptome_2) |
| 33_S2 | Propagule leaf | 5000 | 5000_Fw_W | 1 | YES | original |
| 34_S19 | Propagule leaf | 5000 | 5000_Fw_W | 2 | **NO** | not available |
| 35_S9 | Propagule leaf | 5000 | 5000_Fw_W | 3 | YES | original |
| D33_S40 | Propagule leaf | 5000 | 5000_Fw_D1 | 1 | YES | Drive (transcriptome_2) |
| 42_S39 | Propagule leaf | 5000 | 5000_Fw_D1 | 2 | YES | original |
| D43_S52 | Propagule leaf | 5000 | 5000_Fw_D1 | 3 | **NO** | not available |
| 53_S10 | Propagule leaf | 5000 | 5000_Ft_D1 | 1 | YES | original |
| D54_S5 | Propagule leaf | 5000 | 5000_Ft_D1 | 2 | YES | Drive (transcriptome_2) |
| D55_S5 | Propagule leaf | 5000 | 5000_Ft_D1 | 3 | **NO** | not available |
| 65_S3 | Propagule leaf | 5000 | 5000_Fm_D1 | 1 | YES | original |
| 66_S4 | Propagule leaf | 5000 | 5000_Fm_D1 | 2 | YES | original |
| 67_S5 | Propagule leaf | 5000 | 5000_Fm_D1 | 3 | YES | original |
| D69_S47 | Propagule leaf | 8008 | 8008_Fw_W | 1 | YES | Drive (transcriptome_2) |
| 70_S34 | Propagule leaf | 8008 | 8008_Fw_W | 2 | YES | original |
| D71_S54 | Propagule leaf | 8008 | 8008_Fw_W | 3 | **NO** | not available |
| 77 | Propagule leaf | 8008 | 8008_Fw_D1 | 1 | **NO** | not available |
| 79_S2 | Propagule leaf | 8008 | 8008_Fw_D1 | 2 | YES | Drive (transcriptome_1) |
| 79_S6 | Propagule leaf | 8008 | 8008_Fw_D1 | 3 | YES | original |
| 89_S30 | Propagule leaf | 8008 | 8008_Ft_D1 | 1 | YES | original |
| D66_S41 | Propagule leaf | 8008 | 8008_Ft_D1 | 2 | **NO** | not available |
| 91_S37 | Propagule leaf | 8008 | 8008_Ft_D1 | 3 | YES | original |
| 101_S20 | Propagule leaf | 8008 | 8008_Fm_D1 | 1 | YES | Drive (transcriptome_1) |
| 102_S3 | Propagule leaf | 8008 | 8008_Fm_D1 | 2 | YES | Drive (transcriptome_1) |
| 103_S4 | Propagule leaf | 8008 | 8008_Fm_D1 | 3 | **NO** | not available |

# Results

We can divide our analysis in three steps: data cleaning, preliminary analysis and analysis.

## Data cleaning

We assessed the quality of the raw sequencing data using **FastQC**, and summarized the results across all samples with **MultiQC**.  

Reads were then cleaned with **BBDuk**, using the SILVA and Rfam databases to remove rRNA contamination. The following parameters were applied:  

- `forcetrimleft=11` → removes the first 11 bases from the 5′ end.  
- `forcetrimright2=3` → removes 3 bases from the 3′ end.  
- `minlength=80` → discards reads shorter than 80 bp after trimming.  
- `qtrim=w` → trims both ends of the read based on quality.  
- `trimq=20` → sets the quality threshold for trimming (Phred score 20).  
- `tpe=t` → trims both reads in a pair equally.  
- `tbo=t` → trims adapters based on pair overlap.  






