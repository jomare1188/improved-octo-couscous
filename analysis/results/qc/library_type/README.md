# Library type

`../../../code/17_library_type.sh`: 1M read pairs of six libraries mapped with HISAT2 (strand RF) to Chr01g1 and Chr02g1; read 1, unique (MAPQ 60), classified against the POJ2878 GFF. Table: `read_locations.tsv`.

- **Stranded, dUTP-type:** Salmon library type ISR in all 51 libraries (≤ 2.6% reads on the wrong strand).
- **Paired-end 2 × 100 bp**, about 11 M pairs per library.
- **Poly-A selected:** only 3.4–7.8% of genic reads fall in introns in five of six libraries (leaves on runs 52, 63, 73; buds on runs 44 and 45). Total-RNA libraries usually give 20–50%.
- **Exception, `79_S6` (8008 Fw_D1):** 48% of genic reads are intronic and only 3% are spliced, even after removing an rRNA-like pile-up (Chr01g1 24–25 Mb, 56% of its mapped reads). Genic reads are sense to genes (85k sense vs 1.5k antisense), so it is RNA, not DNA: the poly-A selection of this library seems to have failed (it also has 76% rRNA).
- Reads outside annotated genes (29–38% in good libraries) are mostly unannotated UTRs (about half of the mRNAs have no UTR in the GFF) and reads from homeologous copies that are missing from this two-chromosome reference.
