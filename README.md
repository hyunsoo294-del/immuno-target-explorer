# TAA Immune Profiler

Private research web app for exploring public TAA (tumor-associated antigen) expression across TCGA bulk tumors, human cancer cell lines (DepMap/CCLE), and curated single-cell catalogs.

**Not for clinical decisions.** Expression values are never converted into treatment-response or affinity scores.

Google Drive workspace: [TAA Immune Profiler](https://drive.google.com/drive/folders/1nU-dk6ShTsskcue7Mo9SIL9hsNZfXVRu)

## Features

1. **Target expression** — TCGA PanCancer Atlas 2018 RNA (+ RPPA when present), subtype notes, relative TAA Top 20 (quartile co-elevation).
2. **Cancer cell lines** — DepMap 25Q3 RNA `log2(TPM+1)`, CCLE 2019 RNA RPKM, Nusinow 2020 whole-cell protein; filter, paginate, compare ≤12 lines, CSV/JSON export.
3. **Single-cell lab** — UCSC/TISCH2 curated catalog, lineage/marker controls, co-expression math helpers; matrix analysis when a local H5AD is provided.
4. **Evidence & limits** — assay separation, isoform/glycolipid/stromal warnings.
5. **Validation** — live check log with separate **data fidelity** vs **biology concordance** columns.

## Private access (default)

```bash
cp .env.example .env
# set TAA_PROFILER_ACCESS_TOKEN and/or TAA_PROFILER_ALLOWLIST
# TAA_PROFILER_PUBLIC=0  (do not enable without owner approval)
# local only: TAA_PROFILER_DEV_OPEN=1
```

All routes are gated in `profiler/auth.py`. Obscure URLs are not treated as privacy.

## Run locally

```bash
python3 -m pip install -r requirements.txt
export TAA_PROFILER_DEV_OPEN=1   # or set ACCESS_TOKEN
streamlit run app.py --server.port 8502
```

Open http://127.0.0.1:8502

Legacy Open Targets / HPA helpers remain in `data_sources.py`, `analysis.py`, etc., and are still importable.

## Validation

```bash
python3 scripts/run_validation.py --suite aliases
python3 scripts/run_validation.py --suite her2
python3 scripts/run_validation.py --suite profiles
python3 -m pytest tests/test_gene_resolution.py -v
```

See `docs/VALIDATION_REPORT.md` and `validation_results/`.

## Data sources

- cBioPortal API — https://www.cbioportal.org/api
- TCGA PanCancer Atlas 2018 studies
- `ccle_broad_2025` / `ccle_broad_2019` (RNA + protein)
- UCSC Cell Browser / TISCH2 (catalog + optional local H5AD)

## Architecture

```
profiler/
  gene_resolution.py   # aliases, GD2/CLDN18.2 semantics
  http_client.py       # 429 Retry-After, concurrency, cache
  cbioportal.py        # studies/profiles/molecular-data
  bulk_expression.py   # cancer/subtype RNA+RPPA
  cell_lines.py        # DepMap/CCLE + literature direction
  taa_ranking.py       # Top 20 co-elevation
  single_cell.py       # catalog + marker/coexpression helpers
  validation.py        # 17-target suite
  auth.py              # private allowlist/token gate
  exports.py           # CSV/JSON with provenance
app.py                 # Streamlit UI (5 tabs)
```
