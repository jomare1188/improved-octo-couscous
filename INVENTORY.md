# RNA-seq inventory (2026-09-29)

> **Update 2026-10-22: the 18 Drive libraries are copies.** Every library downloaded from Dolores's Drive is a copy of one of the 51 original libraries, often relabelled. All 18 were dropped, and **51 libraries are used**.
> - Evidence and the copy → original table: [`nf-rnaseq/run_POJ2878/integrity_check/FINDINGS.md`](nf-rnaseq/run_POJ2878/integrity_check/FINDINGS.md).
> - Usable replicates per group: README, section "Usable replicates per group".
>
> The counts below (69 samples, at least 2 replicates per group) describe the files, not the usable samples.

Per-sample details: [`sample_inventory.csv`](sample_inventory.csv).

## What we have: 69 / 78 samples, all in `RAWDATA/`

| Set | Have | Groups with 3/3 | Groups with 2/3 |
|---|---|---|---|
| Mother buds | 18/18 | 6 | 0 |
| Mother leaves | 34/36 | 10 | `8008_Mm_D1`, `8008_Mm_D3` |
| Propagule leaves | 17/24 | `5000_Fm_D1` | the other 7 |

No group has fewer than 2 replicates.

- 51 samples were already in `RAWDATA/`.
- 18 were downloaded from Dolores's Drive on 2026-09-29. Their MD5s match Drive and their read counts match her spreadsheet. Where each file came from is recorded in `RAWDATA/drive_provenance.tsv`.

## Missing: 9 samples (they do not exist, red in Dolores's spreadsheet)

`263_S11`, `14` (mother leaves); `34_S19`, `D43_S52`, `D55_S5`, `D71_S54`, `77`, `D66_S41`, `103_S4` (propagule leaves).

The propagule groups T/W and M/W were never sequenced.

## Doubts

1. **`159_S5` is out of sync.** R2 has 100,000 fewer reads than R1. The last read matches, so the reads are missing from the middle. The file on Drive is the same, so this is not a download problem. It needs re-pairing (e.g. BBTools `repair.sh`) or the original from Dolores.
2. **Duplicate replicate labels.** In `Mt_D1` and `Mt_D3` of both genotypes, two samples are labelled replicate 02 and none 01: `157_S12`/`D158_S48`, `221_S19`/`222_S6`, `165_a_S40`/`166_S23`, `230_S15`/`D231_D46`.
3. **Sample names.** `D54_S5` was called `D54_S53` in the old metadata. `79_S2` and `79_S6` are two replicates of `8008_Fw_D1` that share plant number 79.
4. **`261_S9` and `258_S5`** have exactly the same number of reads (57,208,988) but are different files with different barcodes. It is probably a coincidence; it is noted here just in case.
