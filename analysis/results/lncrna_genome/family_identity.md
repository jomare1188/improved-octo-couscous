# What is the bud lncRNA family? Database searches (2026-11-13)

Query: the 20 lncRNAs of interest (18 shared T-vs-M bud lncRNAs and the 2 leaf D3-vs-D1 ones); for NCBI, only the
family representative `MSTRG.182280.1` (2,623 nt). Search files: `/Storage/data2/jorge.munoz/lnc_genome/family/`.

| Database | Search | Result for the family (13 core sequences) |
|---|---|---|
| NCBI core_nt (online, 2026-11-12 release) | megablast | Only *Saccharum* genomic clones: R570 BACs (KF184738.1, 90% identity; KF184809.1; KF184962.1), SP80-3280 genomic (MF737072.1), *S. officinarum* BAC (MH182569.1). No annotated feature. Outside *Saccharum* only 43-base matches to *Eragrostis tef*. Full list: `family_ncbi_core_nt.txt` |
| NCBI nr proteins (local copy of 2018/2019) | DIAMOND blastx, very sensitive, e < 1e-3 | No hit |
| RepBase 23.10 (57,798 repeats), RepeatMasker library, RepeatPeps | blastn, blastx | No hit (one 68-base match to a rice Gypsy LTR) |
| Rfam | cmscan, gathering cut-off | No hit |
| Sugarcane pan-transcriptome, SP80-3280 transcripts | blastn | Up to 98% identity to `SP80-3280_Mattielo_SP803280_c118623_g1_i1` (7,964 nt), in its non-coding stretch (positions about 1,800 to 5,500); the 302 aa protein predicted there lies outside our sequences |
| Sorghum transcripts (Phytozome v5.1) | blastn | No hit |

Reading: a sugarcane-specific, multi-copy (815 genomic copies), non-coding transcript that no database names.
"Repeat-derived" is our inference from the copy number, not a database result.

Three of the five shared lncRNAs outside the family are small protein genes missing from the annotation (protein
shorter than the 100 aa filter), found by the nr search:

| Transcript | Best protein match | T vs M |
|---|---|---|
| `MSTRG.84729.1`, `MSTRG.95658.1` | glycine-rich cell wall structural protein, sorghum XP_002457317.1 (maize: cold and drought-regulated protein CORA) | higher in M, both genotypes |
| `MSTRG.136860.1` | wound-induced basic protein, sorghum XP_002449445.1 | higher in T, both genotypes |
