# Data dictionary

| Field | Meaning |
|---|---|
| requested target | User input string |
| official_symbol | Hugo gene symbol after alias + provider confirm |
| entrez_gene_id | NCBI Entrez from cBioPortal `/api/genes` |
| molecularProfileId | cBioPortal profile actually returned |
| value (DepMap 25Q3 RNA) | `log2(TPM+1)` as published — no extra log |
| value (CCLE 2019 RNA) | RPKM in tables; `log2(RPKM+1)` in plots |
| value (CCLE protein) | Normalized whole-cell MS relative log2 (Nusinow); negatives valid |
| not_measured | Gene absent from profile or zero returned rows — ≠ antigen absent |
| measured 0 | True numeric zero from source |
| delta_mean_log2p1 | mean log2(expr+1) high-quartile − low-quartile of reference TAA |
| data_fidelity | Pass/fail of ID/profile/value transport vs source |
| biology_concordance | Optional direction check vs literature surface assays — independent of fidelity |
