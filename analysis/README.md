# Drought-memory analysis: plan, steps and scripts

This file lists the comparisons and the scripts. The results are documented in full in the project [`README.md`](../README.md), and in short in [`results/report_story.md`](results/report_story.md).

Plan agreed on 2026-10-29. Each genotype is analysed separately, because the two cultivars are expected to respond
differently. Groups with one replicate stay in the models, where they help estimate variance, but no conclusions rest
on them. Maturation mother leaves are descriptive only. Alleles are grouped through OrthoFinder orthologues to sorghum.

## Hypotheses and comparisons

| # | Question | Contrasts (per genotype) | Status |
|---|---|---|---|
| H1 | Do leaves "learn": is the 3rd drought response different from the 1st? | `Mt_D1_vs_Mt_W`, `Mt_D3_vs_Mt_W`, `Mt_D3_vs_Mt_D1`, plus memory classes | 8008 conclusive; 5000 exploratory (Mt_D1 n=1) |
| — | Maturation mother leaves | `Mm_D1_vs_Mm_W`, `Mm_D3_vs_Mm_W`, `Mm_D3_vs_Mm_D1`, `Mm_W_vs_Mt_W` | descriptive (5000 n=1; 8008 has W only) |
| H3 | Does drought leave a lasting mark in buds? | `T_vs_W`, `M_vs_W`, `T_vs_M` | conclusive (but see batch note) |
| H4 | Do propagules from stressed mothers respond differently? | `Fw_D1_vs_Fw_W`, `Ft_D1_vs_Fw_D1`, `Fm_D1_vs_Fw_D1`, `Ft_D1_vs_Fw_W`, `Fm_D1_vs_Fw_W` | mostly exploratory (n=1 groups) |
| H5 | Is it the same genes from mother to bud to propagule? | overlap of H1 / H3 / H4 results | preliminary, `results/report_tables/h5_overlaps.tsv` |
| H6 | Do known memory genes behave as memory genes here? | gene-set tests of literature lists in every contrast | `results/memory_genes/` |

**Memory classes** follow Ding et al. 2013 and 2014: the first sign is D1 vs W, the second is D3 vs D1.
- `++` / `--`: memory genes, with a stronger response at D3.
- `+-` / `-+`: revised-response genes.
- `+=` / `-=`: non-memory genes.
- `0+` / `0-`: late responders, with no D1 response but a D3 response.

A DEG is padj < 0.05 and |log2FC| ≥ 1.

**To read the documents on a laptop:** `results/DOLORES_report.zip` holds both as HTML, plus every file they name (rebuild with `code/09_report_package.py`). [`results/report.md`](results/report.md) is the older technical summary of 2026-11-03; it is kept for reference and is no longer updated.

## Steps and scripts (`code/`)

| Step | Script | SGE wrapper | Output |
|---|---|---|---|
| Import Salmon → genes or allele groups | `01_import.R <gene\|group>` | `run_import_qc.sh` | `data/txi_*.rds`, `data/samples.tsv` |
| QC: library sizes, PCA, correlation | `02_qc.R <level>` | `run_qc.sh` | `results/qc/<level>/` |
| Genotype identity check | `qc_genotype_check.R` | `run_genotype_check.sh` | `results/qc/genotype_check/` |
| DESeq2 per genotype, memory classes | `03_de.R <level>` | `run_de.sh` | `results/de/<level>/` |
| Known memory genes | `04_memory_genes.R` | in `run_downstream.sh` | `results/memory_genes/` |
| GO / KEGG enrichment | `05_enrichment.R` | in `run_downstream.sh` | `results/enrichment/` |
| Plots: PCA, UpSet, top GO terms (scico palettes) | `06_plots.R` | `run_plots.sh` | `results/plots/` |
| Extra numbers for the report (H1, H4, H5) | `07_report_stats.R` | runs in seconds, no job needed | `results/report_tables/` |
| Explanatory figures for the narrative report | `08_story_figures.R` | `run_story.sh` | `results/plots/story/`, `results/qc/sequencing_runs.tsv` |
| HTML pages (short story + full documentation) and a zip with every file they reference, to read off the server | `09_report_package.py` | runs in seconds in `r_RNASEQ` (needs pandoc), no job needed | `results/DOLORES_report.zip` |
| Shared settings and helpers for steps 10–16 | `common.R` | `run_r.sh` (generic: `qsub -N name [-v RENV=env] code/run_r.sh code/X.R`) | |
| H2: genotype × cycle interaction, downsampling check | `10_genotype_cycle.R` | `run_r.sh` | `results/h2_genotype_cycle/` |
| Bud batch check: run QC metrics, GC/length bias, control genes, RUVg | `11_bud_batch.R` | `run_r.sh` | `results/bud_batch/` |
| H5 rank-based tests (camera, RRHO) | `12_h5_rank.R` | `run_r.sh` | `results/h5_rank/` |
| Allele-level look at the candidates (I², bootstrap InfRV) | `13_alleles.R` | `run_r.sh` | `results/alleles/` |
| WGCNA input (leaf VST, genotype removed; traits) | `14_wgcna_prep.R` | `run_r.sh` | `results/wgcna/input_*` |
| WGCNA network and modules | `15_wgcna.R` | `run_r.sh` with `RENV=/Storage/data2/jorge.munoz/conda_envs/r_wgcna` | `results/wgcna/modules.tsv`, `eigengenes.tsv` |
| Module–trait tests, GO, overlaps, buds | `16_wgcna_traits.R` | `run_r.sh` | `results/wgcna/` |
| Which group the three excluded libraries resemble (nearest centroid, leave-one-out calibration) | `19_relabel_check.R` | `run_r.sh` | `results/qc/relabel_check/` |
| Distance between libraries of a group, within and across sequencing runs | `20_run_check.R` | runs in seconds, no job needed | `results/qc/run_pair_distances*.tsv` |
| Clock genes: internal time of each library, and THI1 (sampling hour) | `23_clock_check.R` | runs in 2 minutes, no job needed | `results/clock/` |
| Library type: poly-A vs total RNA (intron reads, HISAT2 on 2 chromosomes) | `17_library_type.sh` | itself (`qsub -N libtype code/17_library_type.sh`) | `results/qc/library_type/` |
| lncRNA pilot on chromosome 1 (HISAT2, StringTie, CPC2, DIAMOND) | `18_lnc_pilot.sh` | itself (`qsub -N lncpilot code/18_lnc_pilot.sh`; env `/Storage/data2/jorge.munoz/conda_envs/lnc`) | `results/lncrna_pilot/` |
| lncRNA search on the whole genome, 48 libraries (side analysis; story section 8, README 13.4): index, align and assemble, candidates, Salmon, drought test, bud family and figure | `21_lnc_index.sh`, `22_lnc_align.sh`, `24_lnc_candidates.sh`, `25_lnc_quant.sh`, `26_lnc_de.R`, `27_lnc_family.R` | chain of jobs, see the header of each script (work files in `/Storage/data2/jorge.munoz/lnc_genome`) | `results/lncrna_genome/` |

- **`run_downstream.sh`** runs everything after OrthoFinder: allele groups → import → QC → DE → memory genes → enrichment.
- **DESeq2 runs serially:** with forked BiocParallel workers the group-level run hung.
- **Resubmitting:** steps whose output already exists are skipped, so the job can simply be resubmitted after a failure.
- **Allele groups** are built in `../orthology/`, see `../orthology/README.md`.
- **Literature lists** are in `../memory_genes/`, see `SOURCES.md`.

## QC findings (2026-10-29, allele level)

1. **Three libraries carry the other genotype's alleles.** The check uses on/off expression of 31,000 genes, which reflects which alleles each cultivar carries. Every other library matches its label clearly. None of the three duplicates another library: each looks like a real sample of the other genotype from the same tissue and treatment type, which suggests label or tube mix-ups.
   - `FIDH35000_S29`: labelled 5000 tillering bud, matches 8008.
   - `21DH8_S24`: labelled 8008 maturation bud, matches 5000. This and `FIDH35000_S29` look like a swap.
   - `251_S3`: labelled 5000 Mt_W leaf, matches 8008.

   **They are excluded from every model** (`data/excluded_samples.tsv`; decided 2026-11-04). Dolores's sample table is correct, so the mix-up happened at sampling. An expression check (`results/qc/relabel_check/`) suggests each keeps the treatment on its label, in the other cultivar, but every group was designed with 3 plants and relabelling would give 4. To relabel anyway, correct the cultivar and group in `data/samples.tsv`, remove the library from `data/excluded_samples.tsv` and rerun from `03_de.R`.
2. **Batch confounding, from the read headers (instrument VH00942):**
   - **Buds:** all 6 watered buds were sequenced on run 44 and all 12 stressed buds on run 45. Bud `T_vs_W` and `M_vs_W` are therefore confounded with sequencing run. `T_vs_M` is clean, since both are on run 45. In the PCA, watered buds separate from stressed buds on PC2 in both genotypes, while T and M buds overlap.
   - **Mother leaves:** the tillering controls (Mt_W) were sequenced on later runs (61, 63, 73), and every drought leaf on runs 51/52. `Mt_D1_vs_Mt_W` and `Mt_D3_vs_Mt_W` are therefore partly confounded with run. The memory contrast `Mt_D3_vs_Mt_D1` is within runs 51/52, so it is clean. The maturation controls (Mm_W) are on runs 51/52.
   - **Propagules:** all on runs 51/52. Each Fw_D1 reference plant is on the other run from the plants compared with it, in opposite arrangements in the two cultivars. Runs 51 and 52 behave as one batch (`results/qc/run_pair_distances_summary.tsv`).
3. **Sample structure:** in leaves, drought dominates PC1. In buds, genotype is PC1.

## Sampling times (from the experiment figures; to confirm with Dolores)

- Leaves were sampled at the peak of each drought (day 9 of the cycle).
- Controls: Mt_W on day 39, Mm_W on day 174.
- Buds on day 180. Tillering-stressed buds were then 135 days past their last drought; maturation-stressed buds were 6 days into recovery.
- Propagules were sampled at the peak of their one drought cycle.
