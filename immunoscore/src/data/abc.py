"""Surface ABC estimate.

Order: measured ABC, else the RNA fallback curve. CCLE proteomics is not loaded
in this build, so the protein step is skipped and the tier stays rna_only.
Whole-cell lysate MS would not be surface density even if it were loaded.
This is a prioritization heuristic, not a calibrated ABC.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
import pandas as pd

from immunoscore.src.config_loader import DATA_DIR, load_configs


@lru_cache(maxsize=1)
def measured_abc_table() -> pd.DataFrame:
    path = DATA_DIR / "curated" / "measured_abc.csv"
    if not path.exists():
        return pd.DataFrame(columns=["gene_symbol", "ModelID", "abc"])
    table = pd.read_csv(path)
    if table.empty:
        return table
    table["ModelID"] = table["ModelID"].astype(str)
    table["gene_symbol"] = table["gene_symbol"].astype(str)
    return table


def tpm_from_log2(log2_tpm_plus_1: pd.Series) -> pd.Series:
    return np.exp2(log2_tpm_plus_1) - 1.0


def estimate_surface_abc(gene: str, log2_values: pd.Series, surface_fraction: float) -> pd.DataFrame:
    """Return ModelID, abc, evidence_tier, source_detail."""
    params = load_configs()["params"]["abc_rna_fallback"]
    scale = float(params["molecules_per_tpm"])
    measured = measured_abc_table()
    rows = []
    gene_measured = measured[measured["gene_symbol"] == gene] if not measured.empty else measured
    measured_map = {}
    if gene_measured is not None and not gene_measured.empty:
        for record in gene_measured.itertuples(index=False):
            tier = getattr(record, "evidence_tier", "protein_measured")
            detail = getattr(record, "source_detail", "measured_abc.csv")
            if tier is None or (isinstance(tier, float) and pd.isna(tier)) or str(tier) == "nan":
                tier = "protein_measured"
            if detail is None or (isinstance(detail, float) and pd.isna(detail)) or str(detail) == "nan":
                detail = "measured_abc.csv"
            measured_map[str(record.ModelID)] = (float(record.abc), str(tier), str(detail))
    for model_id, log_value in log2_values.items():
        if model_id in measured_map and pd.notna(measured_map[model_id][0]):
            abc, tier, detail = measured_map[model_id]
            rows.append(
                {
                    "ModelID": model_id,
                    "abc": abc,
                    "evidence_tier": tier,
                    "source_detail": detail,
                }
            )
            continue
        if pd.isna(log_value):
            rows.append(
                {
                    "ModelID": model_id,
                    "abc": np.nan,
                    "evidence_tier": "rna_only",
                    "source_detail": "expression missing",
                }
            )
            continue
        tpm = max(float(np.exp2(log_value) - 1.0), 0.0)
        abc = scale * tpm * float(surface_fraction)
        rows.append(
            {
                "ModelID": model_id,
                "abc": abc,
                "evidence_tier": "rna_only",
                "source_detail": f"RNA fallback ABC={scale:.0f}*TPM*{surface_fraction:.2f}; unvalidated",
            }
        )
    return pd.DataFrame(rows)
