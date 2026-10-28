# TAA Immune Profiler — change log (vs prior Immuno Target Explorer / v13 handoff)

## Fixed
- Alias failures for Mesothelin, TROP-2, Nectin-4, PVRL4, CD20, BCMA, CD269, CEA, CD66e, CLDN18.2, Claudin18.2
- HTTP 429 on CD19 protein, BCMA protein, DLL3 RNA 2019 — Retry-After + bounded backoff + concurrency ≤3 + TTL cache; revalidated successfully
- GD2 no longer mapped to a gene; CLDN18.2 isoform warning; FAP stromal warning
- Failure states are not rendered as empty/zero expression charts

## Added
- cBioPortal-backed cancer bulk RNA/RPPA tab
- Cell-line assays (25Q3 RNA, 2019 RNA, 2019 protein) with filters, pagination, compare, export
- TAA Top 20 co-elevation ranking with explicit metric definition
- Single-cell dataset catalog + analysis controls
- Validation tab + CLI with fidelity vs biology columns
- Private auth gate (token / allowlist); public mode off by default
- Google Drive project folder for artifacts

## Preserved
- Existing Open Targets / HPA / GTEx / CPTAC / character modules remain in the repo for reuse
