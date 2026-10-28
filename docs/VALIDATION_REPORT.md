# TAA Immune Profiler — Validation Report

- Run timezone notes below use Asia/Seoul (UTC+9) when converted from UTC stamps.
- Data fidelity and biological concordance are reported as separate columns.
- Source values were never altered to match expected biology.

## Suites executed (this build)

### 1. Alias / official-symbol mapping
- File: `validation_results/aliases_latest.json`
- Result: **25/25 pass**, 0 mismatch
- Includes all v13 failed inputs: Mesothelin, TROP-2, Nectin-4, PVRL4, CD20, BCMA, CD269, CEA, CD66e, CLDN18.2, Claudin18.2
- GD2 correctly unresolved as glycolipid
- CLDN18.2 maps to CLDN18 with isoform_limited warning

### 2. HER2 regression (2026-10-20 snapshot)
- File: `validation_results/her2_latest.json`
- DepMap 25Q3 RNA: ZR751_BREAST **6.7664**, MDAMB231_BREAST **4.5967** — fidelity **pass**; biology direction ZR>MDA **pass**
- CCLE 2019 RNA RPKM: 53.7175 vs 12.8311 — fidelity **pass**
- CCLE 2019 protein: 0.8146 vs 1.2413 — fidelity **pass**; biology vs surface literature **mismatch_direction** (documented; do not swap values)
- Entrez ERBB2 = 2064 confirmed

### 3. 17-target × 3 assay cell-line profiles
- File: `validation_results/profiles_latest.json`
- Summary: pass=47, mismatch=0, fetch_failed=**0**, data_insufficient=1, held=1
- Fidelity pass=48, fidelity fail=0
- Biology pass_direction=26, mismatch_direction=1 (HER2 protein vs surface flow)

Previously failed HTTP 429 queries (v13 handoff) that now complete:
| Query | Assay | n_measured | status |
|---|---|---:|---|
| CD19 | protein | 54 | pass |
| BCMA | protein | 18 | pass |
| DLL3 | RNA 2019 | 1156 | pass |

Notable availability (not absence):
| Target | Assay | n_measured | interpretation |
|---|---|---:|---|
| PSMA | protein | 53 | measured subset only |
| GPRC5D | protein | 0 | **not_measured** / data_insufficient — not “antigen absent” |
| GD2 | all gene assays | — | held (glycolipid) |

## Not claimed as done
- Full browser E2E click-through was not automated as PASS.
- Full TISCH2/UCSC matrix marker stats require local H5AD or authorized matrix fetch; catalog + analysis controls are live, fabricated expression values are not.
- Regulatory approval status of each TAA was not re-verified against FDA/EMA this run; categories are indicative labels in the panel table.

## Re-run

```bash
python3 scripts/run_validation.py --suite aliases --out validation_results/aliases_latest.json
python3 scripts/run_validation.py --suite her2 --out validation_results/her2_latest.json
python3 scripts/run_validation.py --suite profiles --out validation_results/profiles_latest.json
python3 scripts/run_validation.py --suite full --out validation_results/full_latest.json
```
