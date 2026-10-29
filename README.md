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

**We can use 51 of the 78 samples.** `RAWDATA/` holds 69 libraries: the 51 we already had, plus 18 downloaded from Dolores's Google Drive on 2026-09-29 (MD5s checked, see `RAWDATA/drive_provenance.tsv`). The other 9 were never available.

> **The 18 Drive libraries are not independent samples (checked 2026-10-22).** Every one is a copy of one of the 51 original libraries, often relabelled as another sample, so all 18 were dropped. The evidence:
> - Their read names are found in the original libraries.
> - They share the originals' barcode error spectra.
> - Some have repeated read names, R1/R2 out of step, or Windows line endings.
>
> Which original each one copies is in the sample table below and in [`nf-rnaseq/run_POJ2878/integrity_check/FINDINGS.md`](nf-rnaseq/run_POJ2878/integrity_check/FINDINGS.md).

| Sample set | Groups | Planned samples | In `RAWDATA/` | Copies (dropped) | **Used** | Missing (not available) |
|---|---|---|---|---|---|---|
| Mother plants, buds (end of experiment 1) | 6 | 18 | 18 | 0 | **18** | 0 |
| Mother plants, leaves (experiment 1, during the drought cycles) | 12 | 36 | 34 | 12 | **22** | 2 |
| Propagules, leaves (experiment 2) | 8 | 24 | 17 | 6 | **11** | 7 |
| **Total** | **26** | **78** | **69** | **18** | **51** | **9** |

## Usable replicates per group

Genotype 5000 = IACSP95-5000, 8008 = IACCTC07-8008. W = well-watered (C in the spreadsheet), D = drought, D1/D3 = sampled at the 1st/3rd drought cycle. Counts are the 51 libraries used. **Three groups have no usable sample** (8008 Mm_D1, 8008 Mm_D3, 8008 Fm_D1), and 7 more have only one.

**Mother buds**

| Genotype | W | Tillering D | Maturation D |
|---|---|---|---|
| 5000 | 3/3 | 3/3 | 3/3 |
| 8008 | 3/3 | 3/3 | 3/3 |

**Mother leaves**

| Genotype | Tillering W (Mt_W) | Tillering D1 (Mt_D1) | Tillering D3 (Mt_D3) | Maturation W (Mm_W) | Maturation D1 (Mm_D1) | Maturation D3 (Mm_D3) |
|---|---|---|---|---|---|---|
| 5000 | 3/3 | 1/3 | 2/3 | 3/3 | 1/3 | 1/3 |
| 8008 | 3/3 | 3/3 | 2/3 | 3/3 | **0/3** | **0/3** |

**Propagule leaves (mother group / propagule treatment)**

| Genotype | W/W (Fw_W) | W/D (Fw_D1) | T/D (Ft_D1) | M/D (Fm_D1) |
|---|---|---|---|---|
| 5000 | 2/3 | 1/3 | 1/3 | 3/3 |
| 8008 | 1/3 | 1/3 | 2/3 | **0/3** |

The T/W and M/W propagule groups were not sequenced, so the propagule set can compare drought-stressed propagules from stressed and non-stressed mothers, but it has no watered controls for propagules from stressed mothers.

Missing samples and other doubts are listed in [`INVENTORY.md`](INVENTORY.md). The out-of-sync `159_S5` pair noted there is one of the dropped copies.

## Sample table

Full details, with read counts and notes, are in [`sample_inventory.csv`](sample_inventory.csv). Group names follow the nf-core sample sheet: Mt = mother, tillering; Mm = mother, maturation; Fw / Ft / Fm = propagule from a well-watered / tillering-stressed / maturation-stressed mother.

| Sample | Set | Genotype | Group | Replicate | In `RAWDATA/` | Source | Used |
|---|---|---|---|---|---|---|---|
| F25C1_S6 | Mother bud | 5000 | 5000_Maturation_W | 1 | YES | original | yes |
| F25C2_S7 | Mother bud | 5000 | 5000_Maturation_W | 2 | YES | original | yes |
| F25C3_S8 | Mother bud | 5000 | 5000_Maturation_W | 3 | YES | original | yes |
| FIDH15000_S27 | Mother bud | 5000 | 5000_Tillering_D | 1 | YES | original | yes |
| FIDH25000_S28 | Mother bud | 5000 | 5000_Tillering_D | 2 | YES | original | yes |
| FIDH35000_S29 | Mother bud | 5000 | 5000_Tillering_D | 3 | YES | original | yes |
| 21DH5_S21 | Mother bud | 5000 | 5000_Maturation_D | 1 | YES | original | yes |
| 22DH5_S22 | Mother bud | 5000 | 5000_Maturation_D | 2 | YES | original | yes |
| 23DH5_S23 | Mother bud | 5000 | 5000_Maturation_D | 3 | YES | original | yes |
| F28C1_S9 | Mother bud | 8008 | 8008_Maturation_W | 1 | YES | original | yes |
| F28C2_S10 | Mother bud | 8008 | 8008_Maturation_W | 2 | YES | original | yes |
| F28C3_S11 | Mother bud | 8008 | 8008_Maturation_W | 3 | YES | original | yes |
| FIDH1800_S30 | Mother bud | 8008 | 8008_Tillering_D | 1 | YES | original | yes |
| FIDH28000_S31 | Mother bud | 8008 | 8008_Tillering_D | 2 | YES | original | yes |
| FIDH38000_S32 | Mother bud | 8008 | 8008_Tillering_D | 3 | YES | original | yes |
| 21DH8_S24 | Mother bud | 8008 | 8008_Maturation_D | 1 | YES | original | yes |
| 22DH8_S25 | Mother bud | 8008 | 8008_Maturation_D | 2 | YES | original | yes |
| 23DH8_S26 | Mother bud | 8008 | 8008_Maturation_D | 3 | YES | original | yes |
| Dolores_249_S3 | Mother leaf | 5000 | 5000_Mt_W | 1 | YES | original | yes |
| 250_S2 | Mother leaf | 5000 | 5000_Mt_W | 2 | YES | original | yes |
| 251_S3 | Mother leaf | 5000 | 5000_Mt_W | 3 | YES | original | yes |
| 157_S12 | Mother leaf | 5000 | 5000_Mt_D1 | 2 | YES | original | yes |
| D158_S48 | Mother leaf | 5000 | 5000_Mt_D1 | 2 | YES | Drive (transcriptome_3) | **no**, copy of `223_S42` |
| 159_S5 | Mother leaf | 5000 | 5000_Mt_D1 | 3 | YES | Drive (transcriptome_3) | **no**, copy of `157_S12` |
| 221_S19 | Mother leaf | 5000 | 5000_Mt_D3 | 2 | YES | original | yes |
| 222_S6 | Mother leaf | 5000 | 5000_Mt_D3 | 2 | YES | Drive (transcriptome_3) | **no**, copy of `223_S42` |
| 223_S42 | Mother leaf | 5000 | 5000_Mt_D3 | 3 | YES | original | yes |
| 156_S22 | Mother leaf | 5000 | 5000_Mm_W | 1 | YES | original | yes |
| 154_S28 | Mother leaf | 5000 | 5000_Mm_W | 2 | YES | original | yes |
| 155_S11 | Mother leaf | 5000 | 5000_Mm_W | 3 | YES | original | yes |
| 253_S41 | Mother leaf | 5000 | 5000_Mm_D1 | 1 | YES | original | yes |
| 254_S7 | Mother leaf | 5000 | 5000_Mm_D1 | 2 | YES | Drive (transcriptome_2) | **no**, copy of `253_S41` |
| 255_S8 | Mother leaf | 5000 | 5000_Mm_D1 | 3 | YES | Drive (transcriptome_3) | **no**, copy of `157_S12` |
| 9_S15 | Mother leaf | 5000 | 5000_Mm_D3 | 1 | YES | Drive (transcriptome_1) | **no**, copy of `11_S1` |
| 10_S16 | Mother leaf | 5000 | 5000_Mm_D3 | 2 | YES | Drive (transcriptome_2) | **no**, copy of `251_S3` |
| 11_S1 | Mother leaf | 5000 | 5000_Mm_D3 | 3 | YES | original | yes |
| Dolores_257_S5 | Mother leaf | 8008 | 8008_Mt_W | 1 | YES | original | yes |
| 258_S5 | Mother leaf | 8008 | 8008_Mt_W | 2 | YES | original | yes |
| Dolores_259_S6 | Mother leaf | 8008 | 8008_Mt_W | 3 | YES | original | yes |
| 165_a_S40 | Mother leaf | 8008 | 8008_Mt_D1 | 2 | YES | original | yes |
| 166_S23 | Mother leaf | 8008 | 8008_Mt_D1 | 2 | YES | original | yes |
| 168_S25 | Mother leaf | 8008 | 8008_Mt_D1 | 3 | YES | original | yes |
| 230_S15 | Mother leaf | 8008 | 8008_Mt_D3 | 2 | YES | original | yes |
| D231_D46 | Mother leaf | 8008 | 8008_Mt_D3 | 2 | YES | Drive (transcriptome_2) | **no**, copy of `230_S15` |
| 232_S8 | Mother leaf | 8008 | 8008_Mt_D3 | 3 | YES | original | yes |
| 161_S13 | Mother leaf | 8008 | 8008_Mm_W | 1 | YES | original | yes |
| 162_S20 | Mother leaf | 8008 | 8008_Mm_W | 2 | YES | original | yes |
| 163_S7 | Mother leaf | 8008 | 8008_Mm_W | 3 | YES | original | yes |
| 261_S9 | Mother leaf | 8008 | 8008_Mm_D1 | 1 | YES | Drive (transcriptome_2) | **no**, copy of `258_S5` |
| 262_S10 | Mother leaf | 8008 | 8008_Mm_D1 | 2 | YES | Drive (transcriptome_2) | **no**, copy of `157_S12` |
| 263_S11 | Mother leaf | 8008 | 8008_Mm_D1 | 3 | **NO** | not available | no |
| 13_S17 | Mother leaf | 8008 | 8008_Mm_D3 | 1 | YES | Drive (transcriptome_2) | **no**, copy of `230_S15` |
| 14 | Mother leaf | 8008 | 8008_Mm_D3 | 2 | **NO** | not available | no |
| 15_S18 | Mother leaf | 8008 | 8008_Mm_D3 | 3 | YES | Drive (transcriptome_2) | **no**, copy of `230_S15` |
| 33_S2 | Propagule leaf | 5000 | 5000_Fw_W | 1 | YES | original | yes |
| 34_S19 | Propagule leaf | 5000 | 5000_Fw_W | 2 | **NO** | not available | no |
| 35_S9 | Propagule leaf | 5000 | 5000_Fw_W | 3 | YES | original | yes |
| D33_S40 | Propagule leaf | 5000 | 5000_Fw_D1 | 1 | YES | Drive (transcriptome_2) | **no**, copy of `65_S3` |
| 42_S39 | Propagule leaf | 5000 | 5000_Fw_D1 | 2 | YES | original | yes |
| D43_S52 | Propagule leaf | 5000 | 5000_Fw_D1 | 3 | **NO** | not available | no |
| 53_S10 | Propagule leaf | 5000 | 5000_Ft_D1 | 1 | YES | original | yes |
| D54_S5 | Propagule leaf | 5000 | 5000_Ft_D1 | 2 | YES | Drive (transcriptome_2) | **no**, copy of `65_S3` |
| D55_S5 | Propagule leaf | 5000 | 5000_Ft_D1 | 3 | **NO** | not available | no |
| 65_S3 | Propagule leaf | 5000 | 5000_Fm_D1 | 1 | YES | original | yes |
| 66_S4 | Propagule leaf | 5000 | 5000_Fm_D1 | 2 | YES | original | yes |
| 67_S5 | Propagule leaf | 5000 | 5000_Fm_D1 | 3 | YES | original | yes |
| D69_S47 | Propagule leaf | 8008 | 8008_Fw_W | 1 | YES | Drive (transcriptome_2) | **no**, copy of `70_S34` |
| 70_S34 | Propagule leaf | 8008 | 8008_Fw_W | 2 | YES | original | yes |
| D71_S54 | Propagule leaf | 8008 | 8008_Fw_W | 3 | **NO** | not available | no |
| 77 | Propagule leaf | 8008 | 8008_Fw_D1 | 1 | **NO** | not available | no |
| 79_S2 | Propagule leaf | 8008 | 8008_Fw_D1 | 2 | YES | Drive (transcriptome_1) | **no**, copy of `91_S37` |
| 79_S6 | Propagule leaf | 8008 | 8008_Fw_D1 | 3 | YES | original | yes |
| 89_S30 | Propagule leaf | 8008 | 8008_Ft_D1 | 1 | YES | original | yes |
| D66_S41 | Propagule leaf | 8008 | 8008_Ft_D1 | 2 | **NO** | not available | no |
| 91_S37 | Propagule leaf | 8008 | 8008_Ft_D1 | 3 | YES | original | yes |
| 101_S20 | Propagule leaf | 8008 | 8008_Fm_D1 | 1 | YES | Drive (transcriptome_1) | **no**, copy of `91_S37` |
| 102_S3 | Propagule leaf | 8008 | 8008_Fm_D1 | 2 | YES | Drive (transcriptome_1) | **no**, copy of `91_S37` |
| 103_S4 | Propagule leaf | 8008 | 8008_Fm_D1 | 3 | **NO** | not available | no |

# Results

We can divide our analysis in three steps: data cleaning, preliminary analysis and analysis.

## Data cleaning (earlier analysis)

This section describes the earlier analysis, done before the POJ2878 reference was available (`bbduk/`, `RNAseq/`). The current quantification is described below.

We assessed the quality of the raw sequencing data using **FastQC**, and summarized the results across all samples with **MultiQC**.  

Reads were then cleaned with **BBDuk**, using the SILVA and Rfam databases to remove rRNA contamination. The following parameters were applied:  

- `forcetrimleft=11` → removes the first 11 bases from the 5′ end.  
- `forcetrimright2=3` → removes 3 bases from the 3′ end.  
- `minlength=80` → discards reads shorter than 80 bp after trimming.  
- `qtrim=w` → trims both ends of the read based on quality.  
- `trimq=20` → sets the quality threshold for trimming (Phred score 20).  
- `tpe=t` → trims both reads in a pair equally.  
- `tbo=t` → trims adapters based on pair overlap.

## Reference genome: POJ2878

Reads are quantified against the **POJ2878** genome. POJ2878 is a *S. officinarum* × *S. spontaneum* hybrid and a founder of modern sugarcane cultivars (Nature 2026, doi:10.1038/s41586-026-10576-7).
- **Genome:** NCBI GCA_060231685.1. The headers are renamed from GenBank accessions to the sequence names used in the GFF (`refs/POJ2878/rename_headers.sh`).
- **Gene models:** `Sh_POJ2878_v1.gff3`, from Figshare: 397,650 genes and 424,890 transcripts.
- **Details:** sources, MD5s and checks are in [`refs/POJ2878/SOURCE.txt`](refs/POJ2878/SOURCE.txt). The genome files themselves are not in git.

## Quantification: nf-core/rnaseq 3.27.0 (Salmon)

The 51 usable libraries were run through **nf-core/rnaseq 3.27.0** on the SGE cluster, in Salmon-only mode (no genome alignment), on 2026-10-23.
- **Run script:** [`nf-rnaseq/run_POJ2878/run_nfcore_POJ2878.sh`](nf-rnaseq/run_POJ2878/run_nfcore_POJ2878.sh)
- **Parameters:** [`nf-rnaseq/params_proposal/params_POJ2878_salmon.yaml`](nf-rnaseq/params_proposal/params_POJ2878_salmon.yaml)
- **Cluster config:** [`nf-rnaseq/conf_cluster/sge_cluster.config`](nf-rnaseq/conf_cluster/sge_cluster.config) and [`process_overrides.config`](nf-rnaseq/params_proposal/process_overrides.config)
- **Sample sheets:** `samplesheet_POJ2878_clean.csv` and `sample_metadata_POJ2878_clean.csv`

Main settings:
- Trim Galore 2.3.0 (`--2colour 20 --length 36`).
- rRNA removed with Bowtie2 2.5.4 against the pipeline's default rRNA database (SortMeRNA `smr_v4.3_fast_db`, built from SILVA and Rfam).
- Salmon 2.7.0 with a decoy-aware index (genome as decoy), `--seqBias --gcBias --numBootstraps 30`. The 30 bootstraps are for allele-level analysis, e.g. with Swish.
- Salmon keeps one copy of each set of identical transcript sequences (its default): 65,316 transcripts were collapsed this way, which is common in this polyploid genome.

Results (QC report: [`nf-rnaseq/run_POJ2878/results/multiqc/multiqc_report.html`](nf-rnaseq/run_POJ2878/results/multiqc/multiqc_report.html)):
- **Strandedness:** all 51 libraries are stranded, library type ISR.
- **Salmon mapping rate:** median 69.3% (range 54.6–76.4%).
- **rRNA:** median 14.4% of reads.
- **Weakest libraries:**
  - `79_S6`: 1.05 M mapped fragments, 77% rRNA.
  - `66_S4`, `67_S5`, `65_S3`: 1.9–2.8 M mapped fragments, 52–58% rRNA.
- **Count tables:** the merged tables (`salmon.merged.gene_counts.tsv`, TPM, transcript counts and `tx2gene`) are on the server in `nf-rnaseq/run_POJ2878/results/salmon/`. They are not in git: each is 20–250 MB.

## Functional annotation (GO, InterPro, KEGG)

The POJ2878 gene models come with no functional annotation, so we annotated the predicted proteins ourselves (2026-10-24 to 2026-10-26). Steps, versions and commands are in [`annotation/SOURCE.md`](annotation/SOURCE.md).

1. **Proteins:** gffread on the genome and GFF gives 424,890 proteins. Identical sequences were collapsed to 307,787 unique proteins for annotation (`annotation/proteins/make_proteins.sh`).
2. **eggNOG-mapper 2.1.13** (eggNOG 5.0.2) as an SGE array over 16 chunks (`annotation/eggnog/`). It gives a name, description, GO terms, KEGG KO and pathways, and EC numbers.
3. **InterProScan 6.0.3.0** with InterPro 110.0 and all free member databases, with InterPro2GO (`annotation/interproscan/`). It ran in 4 parts on the SGE cluster.
4. **Merge** (`annotation/go/merge_annotation.py`): copies each unique protein's annotation back to every transcript with that sequence, then groups transcripts into genes using Salmon's `tx2gene`. The gene IDs therefore match the count tables.

Coverage of the 397,650 genes ([`annotation/go/coverage.tsv`](annotation/go/coverage.tsv)):

| Annotation | Genes | % |
|---|---|---|
| Any annotation | 396,632 | 99.7 |
| eggNOG hit | 392,089 | 98.6 |
| InterPro entry | 369,349 | 92.9 |
| GO, InterProScan | 246,830 | 62.1 |
| GO, eggNOG | 137,404 | 34.6 |
| **GO, combined** | **272,954** | **68.6** |
| KEGG KO | 160,589 | 40.4 |

Output tables, in `annotation/go/` on the server (not in git, 90–450 MB each):
- `gene2go_list.tsv`: one line per gene, ready for topGO `readMappings`.
- `gene2go.tsv` and `transcript2go.tsv`: one row per pair, with a `source` column (ips, eggnog or both).
- `gene_annotation.tsv` and `transcript_annotation.tsv`: name, description, InterPro, Pfam, KEGG KO and pathway, EC.

Notes for enrichment:
- eggNOG GO lists already include parent terms; InterPro2GO lists only the specific terms. topGO and clusterProfiler add the parents anyway.
- For pathways, use the eggNOG KEGG pathways. InterProScan's pathway column is mostly Reactome IDs for animal species, so it was left out.

## What is in git and what is not

Only scripts, configs, sample sheets, documentation and small summary tables are in this repository. The following stay on the server, under `/Storage/data1/jorge.munoz/DOLORES`, and are excluded in `.gitignore`:
- raw reads (`RAWDATA/*.fastq.gz`)
- the reference genome (`refs/POJ2878/*.gz`)
- the Salmon results, apart from the MultiQC report
- the protein FASTA files
- the eggNOG and InterProScan outputs and databases
- the GO tables
- Singularity images, Nextflow work directories and logs
