# TAA Immune Profiler — Backtest Report

- Run id folder: `validation_results/backtest_20261028T040445Z/`
- Machine UTC stamp in folder name; Asia/Seoul = UTC+9
- Public URL at time of UI backtest: https://peninsula-fairfield-lucas-cleared.trycloudflare.com

## 1. Data-fidelity suites (source transport)

| Suite | pass | mismatch | fetch_failed | data_insufficient | held | fidelity_pass | fidelity_fail |
|---|---:|---:|---:|---:|---:|---:|---:|
| Aliases (25) | 25 | 0 | 0 | 0 | 0 | 25 | 0 |
| HER2 regression (3 assays) | 3 | 0 | 0 | 0 | 0 | 3 | 0 |
| 17×3 cell-line profiles | 47 | 0 | 0 | 1 | 1 | 48 | 0 |
| Unit tests `test_gene_resolution` | 4/4 | — | — | — | — | — | — |

Notes:
- GPRC5D protein: `data_insufficient` (n_measured=0) — not treated as antigen absent.
- GD2: `held` (glycolipid, no gene mapping).
- Prior v13 HTTP 429 targets (CD19 protein, BCMA protein, DLL3 RNA 2019) completed with fidelity pass.

## 2. Biology concordance (separate from fidelity)

From `literature_direction.json`:
- biology_pass_direction: **26**
- biology_mismatch_direction: **1** — HER2 CCLE 2019 protein ZR-75-1 (0.8146) < MDA-MB-231 (1.2413) vs surface literature ZR>MDA. Values **not** swapped.
- biology_data_insufficient: 11 (missing comparison lines / unmeasured in that assay)
- held: 1 (GD2)

## 3. Public URL user-flow (Playwright)

7/7 pass:
1. home_load
2. her2_resolve (ERBB2 / 2064)
3. her2_cell_lines (n_measured table/chart)
4. trop2_alias (TROP-2 → TACSTD2)
5. cldn18_2_warning
6. gd2_glycolipid
7. validation_section

Screenshots: `/opt/cursor/artifacts/screenshots/bt_*.png`

## 4. Not claimed

- Full TISCH2/UCSC matrix marker computation without local H5AD was not marked pass.
- Statistical significance of Top 20 was not claimed.
- Approval-status of therapies was not re-checked against regulators this run.

## 5. Single-cell analysis (added after backtest)

- UCSC `head-neck` live matrix analysis: status ok, barcode join, composition + marker stats + coexpression.
- Playwright public UI: Run single-cell analysis → n_total / lineage / markers / coexpression pass.
- Analyzable datasets: head-neck, quake-gbm, gbm, colorectal-cancer, human-pancreas (first run caches matrix).
