"""DepMap bulk RNA and model metadata.

Primary expression is log2(TPM+1) from OmicsExpressionProteinCodingGenesTPMLogp1.
This is a prioritization heuristic input, not a surface-protein measurement.
"""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

import pandas as pd

from immunoscore.src.config_loader import DATA_DIR, load_configs

RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"


def manifest() -> dict:
    path = RAW_DIR / "MANIFEST.json"
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def _require_processed() -> None:
    expr = PROCESSED_DIR / "expression_log2tpm.parquet"
    models = PROCESSED_DIR / "models.parquet"
    if expr.exists() and models.exists():
        return
    raise FileNotFoundError(
        "DepMap parquet cache is missing. Run scripts/fetch_depmap.py --release 24Q4"
    )


@lru_cache(maxsize=1)
def load_models() -> pd.DataFrame:
    _require_processed()
    models = pd.read_parquet(PROCESSED_DIR / "models.parquet")
    params = load_configs()["params"]
    excluded = set(params.get("exclude_lineages") or [])
    lineage = models["OncotreeLineage"].fillna("")
    keep = ~lineage.isin(excluded) & lineage.ne("")
    disease = models["OncotreePrimaryDisease"].fillna("").astype(str).str.strip()
    models = models.loc[keep].copy()
    models["cancer"] = disease.loc[models.index]
    missing = models["cancer"].eq("")
    models.loc[missing, "cancer"] = models.loc[missing, "OncotreeLineage"].astype(str)
    models["ModelID"] = models["ModelID"].astype(str)
    return models.reset_index(drop=True)


@lru_cache(maxsize=1)
def load_expression() -> pd.DataFrame:
    """Expression indexed by ModelID. Columns are gene symbols. Values are log2(TPM+1)."""
    _require_processed()
    frame = pd.read_parquet(PROCESSED_DIR / "expression_log2tpm.parquet")
    frame.index = frame.index.astype(str)
    frame.index.name = "ModelID"
    return frame


@lru_cache(maxsize=1)
def load_gene_ids() -> pd.Series:
    path = DATA_DIR / "reference" / "gene_ids.csv"
    if not path.exists():
        return pd.Series(dtype=str)
    table = pd.read_csv(path)
    return table.drop_duplicates("gene_symbol").set_index("gene_symbol")["gene_id"]


def cohort_expression() -> pd.DataFrame:
    """Expression restricted to the cancer-model cohort."""
    models = load_models()
    expr = load_expression()
    ids = [model_id for model_id in models["ModelID"] if model_id in expr.index]
    return expr.loc[ids]


def gene_log2(gene: str) -> pd.Series:
    expr = cohort_expression()
    if gene not in expr.columns:
        return pd.Series(index=expr.index, dtype="float64")
    return expr[gene].astype("float64")
