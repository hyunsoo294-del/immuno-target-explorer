"""Shared numeric helpers. Cutoffs are arguments, not constants baked into callers."""

from __future__ import annotations

import numpy as np
import pandas as pd


def empty_factor(cell_lines: pd.DataFrame, tier: str, detail: str) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "ModelID": cell_lines["ModelID"].astype(str).tolist(),
            "sub_score": np.nan,
            "evidence_tier": tier,
            "source_detail": detail,
        }
    )


def cohort_z(values: pd.Series, ddof: int) -> pd.Series:
    numeric = values.astype("float64")
    mean = numeric.mean(skipna=True)
    std = numeric.std(skipna=True, ddof=ddof)
    if pd.isna(std) or std == 0:
        return pd.Series(np.nan, index=numeric.index)
    return (numeric - mean) / std


def weighted_z(expr: pd.DataFrame, weights: dict, ddof: int) -> pd.Series:
    """Renormalize gene weights over genes that exist. Missing genes are not zeroes."""
    present = {gene: weight for gene, weight in weights.items() if gene in expr.columns and weight}
    if not present:
        return pd.Series(np.nan, index=expr.index)
    total = float(sum(present.values()))
    acc = pd.Series(0.0, index=expr.index)
    coverage = pd.Series(0.0, index=expr.index)
    for gene, weight in present.items():
        z = cohort_z(expr[gene], ddof)
        share = float(weight) / total
        valid = z.notna()
        acc = acc.add(z.fillna(0.0) * share, fill_value=0.0)
        coverage = coverage.add(valid.astype(float) * share, fill_value=0.0)
    acc = acc.where(coverage > 0, np.nan)
    # Put the mass back on the genes that were actually observed for that row.
    return acc / coverage.where(coverage > 0, np.nan)


def sigmoid_scaled(values: pd.Series, center: float, width: float) -> pd.Series:
    width = width if width else 1.0
    return 1.0 / (1.0 + np.exp(-(values - center) / width))


def clip_score(values: pd.Series) -> pd.Series:
    return values.clip(lower=0.0, upper=100.0)
