# POJ2878 functional annotation: sources and steps

The reference (Sh_POJ2878_v1.gff3, Nature 2026, Figshare) ships gene structures only, with no GO,
InterPro or KEGG. The annotation here was made from its predicted proteins, on the SGE cluster,
on 2026-10-24 to 2026-10-26.

## 1. Proteins: `proteins/`
- `make_proteins.sh`: gffread 0.12.7 (`-y`) on the genome `nf-rnaseq/run_POJ2878/results_pilot/genome/POJ2878_1.0.genome.fa`
  and Sh_POJ2878_v1.gff3.
- Terminal stops are removed; 8 proteins with internal stops have those positions replaced by X.
- 424,890 proteins (397,650 genes); identical sequences collapsed to 307,787 unique proteins.
  - `POJ2878_v1.proteins.uniq.faa`
  - `uniq_map.tsv`: representative <tab> member transcript.
- `chunks/`: 16 chunks of 20,000 proteins, used for eggNOG.

## 2. eggNOG-mapper: `eggnog/`
- eggNOG-mapper 2.1.13 (cluster module), database eggNOG 5.0.2 (`/Storage/databases/eggnog-mapper-2.1.13`).
- `run_emapper_search.sh`: SGE array, DIAMOND search per chunk (`-m diamond --no_annot`).
- `run_emapper_annot.sh`: `-m no_search --annotate_hits_table` on the merged seed orthologs, `--dbmem`.
- Output: `out/POJ2878.emapper.annotations`. 302,790 of 307,787 unique proteins annotated
  (GO 34.0%, KEGG KO 39.6%).

## 3. InterProScan: `interproscan/`
- InterProScan 6.0.3.0 (ebi-pf-team/interproscan6, Nextflow 26.04.7, Singularity) with InterPro 110.0.
- Member databases are in `ips6_data/`; the large ones were fetched by `download_ips6_data.sh`.
- Licensed tools were not run: SignalP, Phobius, DeepTMHMM; `--runMl` was off.
- `sge_ips6.config`: SGE executor, up to 16 jobs at once.
- `run_ips6_parts.sh`: 4 parts of about 77k proteins, run one after another, each part's work dir
  deleted at the end (disk space).
- Settings: `--goterms --pathways --formats tsv,gff3 --batchSize 3000 --cpus 4`.
- Output: `out/part_{1..4}.faa.tsv` (+ `.gff3`). 304,945 proteins with a match, 283,335 with an
  InterPro entry, 188,995 with GO.
- `test/`: 1,000-protein test.

## 4. GO and gene tables: `go/`
- `merge_annotation.py` (`run_merge.sh`): copies each representative's annotation to every
  transcript with that sequence, and groups transcripts into genes with
  `nf-rnaseq/run_POJ2878/results/salmon/salmon.merged.tx2gene.tsv`.
- `transcript2go.tsv` and `gene2go.tsv` have a `source` column: ips, eggnog, or ips,eggnog.
  Gene level = union over isoforms.
- `gene2go_list.tsv`: one line per gene, for topGO (`readMappings`).
- `gene_annotation.tsv` and `transcript_annotation.tsv`: eggNOG name and description, InterPro, Pfam
  (from InterProScan), KEGG KO, pathway, EC.
- `coverage.tsv`: coverage summary.
- GO terms are used as reported. eggNOG lists already include parent terms; InterPro2GO does not.
  Enrichment tools such as topGO add parents themselves.
- InterProScan pathways are not used: they are mostly Reactome IDs for animal species. KEGG pathways
  come from eggNOG.
