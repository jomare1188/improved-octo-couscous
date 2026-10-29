# Copied libraries from Dolores's Google Drive (checked 2026-10-22)

Every library in the sample sheet that came from Dolores's Drive (transcriptome_1/2/3) is a copy of
a library already in the original RAWDATA, often relabelled as another sample. 18 of 18. None are
independent samples; all are dropped. The 51 original RAWDATA libraries show no problems.

Evidence (name_screen.tsv, fileinfo.tsv, hits.tsv; jobs 28337 and 28341):
- Repeated read names inside the file, and/or R1/R2 out of step.
- ~2000 sampled read names per suspect: nearly all found in the original listed below.
- Same barcode error spectrum as the original; where the main barcode was rewritten,
  the original's barcode variants remain. Some headers were rewritten to other runs.
- 79_S2, 101_S20, 102_S3 have Windows (CRLF) line endings.

| original (RAWDATA) | label        | copies from Drive (label)                                   |
|--------------------|--------------|-------------------------------------------------------------|
| 91_S37             | 8008_Ft_D1   | 79_S2 (8008_Fw_D1), 101_S20, 102_S3 (8008_Fm_D1)            |
| 230_S15            | 8008_Mt_D3   | 13_S17, 15_S18 (8008_Mm_D3), D231_D46 (8008_Mt_D3)          |
| 157_S12            | 5000_Mt_D1   | 159_S5 (5000_Mt_D1), 255_S8 (5000_Mm_D1), 262_S10 (8008_Mm_D1) |
| 11_S1              | 5000_Mm_D3   | 9_S15 (5000_Mm_D3)                                          |
| 223_S42            | 5000_Mt_D3   | D158_S48 (5000_Mt_D1), 222_S6 (5000_Mt_D3)                  |
| 65_S3              | 5000_Fm_D1   | D54_S5 (5000_Ft_D1), D33_S40 (5000_Fw_D1)                   |
| 70_S34             | 8008_Fw_W    | D69_S47 (8008_Fw_W)                                         |
| 258_S5             | 8008_Mt_W    | 261_S9 (8008_Mm_D1)                                         |
| 251_S3             | 5000_Mt_W    | 10_S16 (5000_Mm_D3)                                         |
| 253_S41            | 5000_Mm_D1   | 254_S7 (5000_Mm_D1)                                         |

Clean sheets: ../samplesheet_POJ2878_clean.csv, ../sample_metadata_POJ2878_clean.csv (51 libraries).
