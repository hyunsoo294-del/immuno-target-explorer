"""Frozen full-cohort anchors for expression logistics.

x50 and k are computed once from the cancer-cell cohort and stored in
data/reference/expression_reference.parquet. Scoring reads that file.
It does not recompute percentiles from the lines the user selected.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from immunoscore.src.config_loader import DATA_DIR

REFERENCE_PATH = DATA_DIR / "reference" / "expression_reference.parquet"


def build_reference_frame(expression: pd.DataFrame, span: float) -> pd.DataFrame:
    """Per-gene median, p10, p90, and k = span / (p90 - p10).

    expression values are log2(TPM+1). A gene with no spread gets k NaN and
    drops out of the weighted logistic.
    """
    quantiles = expression.quantile([0.1, 0.5, 0.9], interpolation="linear")
    frame = pd.DataFrame(
        {
            "gene": quantiles.columns.astype(str),
            "p10": quantiles.loc[0.1].to_numpy(),
            "median": quantiles.loc[0.5].to_numpy(),
            "p90": quantiles.loc[0.9].to_numpy(),
            "n": expression.notna().sum().to_numpy(),
        }
    )
    width = frame["p90"] - frame["p10"]
    frame["k"] = np.where(width > 0, float(span) / width, np.nan)
    return frame


def write_reference(frame: pd.DataFrame, metadata: dict) -> None:
    REFERENCE_PATH.parent.mkdir(parents=True, exist_ok=True)
    table = pa.Table.from_pandas(frame, preserve_index=False)
    existing = dict(table.schema.metadata or {})
    for key, value in metadata.items():
        existing[f"immunoscore.{key}".encode()] = str(value).encode()
    table = table.replace_schema_metadata(existing)
    pq.write_table(table, REFERENCE_PATH)


def reference_metadata() -> dict:
    if not REFERENCE_PATH.exists():
        return {}
    schema = pq.read_schema(REFERENCE_PATH)
    raw = schema.metadata or {}
    out = {}
    prefix = b"immunoscore."
    for key, value in raw.items():
        if key.startswith(prefix):
            out[key[len(prefix):].decode()] = value.decode()
    return out


@lru_cache(maxsize=1)
def load_expression_reference() -> pd.DataFrame:
    if not REFERENCE_PATH.exists():
        raise FileNotFoundError(
            "Frozen expression reference is missing. Run immunoscore/scripts/build_expression_reference.py"
        )
    frame = pd.read_parquet(REFERENCE_PATH)
    frame["gene"] = frame["gene"].astype(str)
    return frame.drop_duplicates("gene").set_index("gene")
