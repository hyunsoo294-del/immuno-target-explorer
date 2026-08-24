# -*- coding: utf-8 -*-
"""Validate and repair RNA fold, CPTAC log-ratio scale, and IHC ordinal stats."""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import pandas as pd

RNA_PSEUDOCOUNT = 0.25
CDAP_LOG_RATIO_ABS_MAX = 5.0
CDAP_DETECT_ABS_MAX = 8.0
SANITY_ATOL = 1e-6
SANITY_RTOL = 1e-9


def rna_log2_fold(tumor: float | None, normal: float | None, pseudocount: float = RNA_PSEUDOCOUNT) -> dict[str, Any]:
    """Fold = (Tumor + 0.25) / (Normal + 0.25); return linear fold and log2 fold."""
    t = _as_float(tumor)
    n = _as_float(normal)
    if t is None or n is None:
        return {"fold": None, "log2_fold": None}
    fold = (t + pseudocount) / (n + pseudocount)
    if fold <= 0 or math.isnan(fold) or math.isinf(fold):
        return {"fold": None, "log2_fold": None}
    return {"fold": float(fold), "log2_fold": float(math.log2(fold))}


def tcga_log2_plus1(raw_median: float | None) -> float | None:
    value = _as_float(raw_median)
    if value is None or value < 0:
        return None
    return float(math.log2(value + 1.0))


def invert_log2_plus1(log2_value: float | None) -> float | None:
    value = _as_float(log2_value)
    if value is None:
        return None
    return float((2.0 ** value) - 1.0)


def sanity_check_tcga_inverse(raw_median: float | None, log2_value: float | None | None = None) -> dict[str, Any]:
    """Confirm 2^(log2(raw+1)) - 1 recovers the raw tumor median."""
    raw = _as_float(raw_median)
    stored_log2 = _as_float(log2_value) if log2_value is not None else None
    expected_log2 = tcga_log2_plus1(raw)
    reconstructed = invert_log2_plus1(expected_log2 if stored_log2 is None else stored_log2)
    ok = False
    if raw is not None and reconstructed is not None:
        ok = math.isclose(reconstructed, raw, rel_tol=SANITY_RTOL, abs_tol=max(SANITY_ATOL, abs(raw) * 1e-9))
    if stored_log2 is not None and expected_log2 is not None:
        ok = ok and math.isclose(stored_log2, expected_log2, rel_tol=1e-9, abs_tol=1e-6)
    return {
        "raw_tumor_median": raw,
        "log2_rsem_plus1": expected_log2,
        "reconstructed_from_2pow_minus1": reconstructed,
        "sanity_ok": bool(ok),
        "looks_like_already_logged": bool(raw is not None and 0 < raw < 30),
    }


def remap_rna_arrays(df: pd.DataFrame) -> pd.DataFrame:
    """Map HPA nTPM -> Normal and TCGA median -> Tumor, then recompute log2 fold."""
    out = df.copy()
    tumor = pd.to_numeric(out.get("cancer_rna_tcga_median"), errors="coerce")
    if tumor is None:
        tumor = pd.to_numeric(out.get("cancer_rna_abs"), errors="coerce")
    hpa = pd.to_numeric(out.get("normal_rna_hpa_ntpm"), errors="coerce")
    gtex = pd.to_numeric(out.get("normal_rna_gtex_ntpm"), errors="coerce")
    normal = hpa.where(hpa.notna(), gtex)

    folds: list[float] = []
    log2s: list[float] = []
    log2_tumor: list[float] = []
    reconstructed: list[float] = []
    sanity: list[bool] = []
    sources: list[str] = []

    for i in out.index:
        t = None if pd.isna(tumor.loc[i]) else float(tumor.loc[i])
        n = None if pd.isna(normal.loc[i]) else float(normal.loc[i])
        calc = rna_log2_fold(t, n)
        check = sanity_check_tcga_inverse(t)
        folds.append(calc["fold"] if calc["fold"] is not None else np.nan)
        log2s.append(calc["log2_fold"] if calc["log2_fold"] is not None else np.nan)
        log2_tumor.append(check["log2_rsem_plus1"] if check["log2_rsem_plus1"] is not None else np.nan)
        reconstructed.append(
            check["reconstructed_from_2pow_minus1"] if check["reconstructed_from_2pow_minus1"] is not None else np.nan
        )
        sanity.append(bool(check["sanity_ok"]) if t is not None else False)
        if t is None or n is None:
            sources.append("RNA fold unavailable")
        else:
            src = "HPA nTPM" if pd.notna(hpa.loc[i]) else "GTEx nTPM"
            sources.append(f"TCGA tumor median / {src}; Fold=(T+0.25)/(N+0.25)")

    out["normal_rna_abs"] = normal
    out["cancer_rna_abs"] = tumor
    out["rna_fold_induction"] = folds
    out["rna_log2_fold"] = log2s
    out["rna_fold_normalized"] = folds
    out["log2_fold_normalized"] = log2s
    out["rna_log2_cancer"] = log2_tumor
    out["rna_log2_normal"] = (normal + RNA_PSEUDOCOUNT).map(lambda x: math.log2(x) if pd.notna(x) and x > 0 else np.nan)
    out["tcga_inverse_2pow_minus1"] = reconstructed
    out["tcga_inverse_sanity_ok"] = sanity
    out["rna_fold_source"] = sources
    return out


def to_cdap_log_ratio(values: pd.Series) -> tuple[pd.Series, str]:
    """Keep CDAP log-ratio if already in ~[-5, 5]; otherwise convert raw intensity."""
    series = pd.to_numeric(values, errors="coerce")
    valid = series.dropna()
    if valid.empty:
        return series, "empty"
    median = float(valid.median())
    p05 = float(valid.quantile(0.05))
    p95 = float(valid.quantile(0.95))
    abs_median = float(valid.abs().median())

    already_ratio = abs_median <= CDAP_LOG_RATIO_ABS_MAX and p95 <= CDAP_DETECT_ABS_MAX and p05 >= -CDAP_DETECT_ABS_MAX
    if already_ratio:
        return series, "cdap_log_ratio"

    if median > CDAP_DETECT_ABS_MAX:
        # Log2 TMT intensity stored as if it were a ratio (e.g. ERBB2 PAAD ~24).
        return series - median, "median_centered_log2_intensity"

    if median > 50:
        centered = np.log2(series / median)
        return pd.Series(centered, index=series.index), "log2_ratio_vs_cohort_median"

    return series, "unclassified_kept"


def normalize_cptac_sample_table(cptac: pd.DataFrame) -> pd.DataFrame:
    if cptac is None or cptac.empty or "value" not in cptac.columns:
        return cptac if isinstance(cptac, pd.DataFrame) else pd.DataFrame()
    frames = []
    group_col = "studyId" if "studyId" in cptac.columns else "cancer_type"
    for _, subset in cptac.groupby(group_col, dropna=False):
        work = subset.copy()
        converted, scale = to_cdap_log_ratio(work["value"])
        work["value_raw"] = pd.to_numeric(work["value"], errors="coerce")
        work["value"] = converted
        work["cptac_scale"] = scale
        frames.append(work)
    return pd.concat(frames, ignore_index=True)


def expand_ihc_counts(high: Any, medium: Any, low: Any, not_detected: Any) -> np.ndarray:
    scores = (
        [3] * int(high or 0)
        + [2] * int(medium or 0)
        + [1] * int(low or 0)
        + [0] * int(not_detected or 0)
    )
    return np.asarray(scores, dtype=float)


def mannwhitney_ihc(tumor_scores: np.ndarray, normal_scores: np.ndarray) -> dict[str, Any]:
    """Two-sided Mann-Whitney U on ordinal IHC 0-3 scores. No pseudocount fold."""
    tumor = np.asarray(tumor_scores, dtype=float)
    normal = np.asarray(normal_scores, dtype=float)
    tumor = tumor[np.isfinite(tumor)]
    normal = normal[np.isfinite(normal)]
    empty = {
        "ihc_mwu_u": None,
        "ihc_mwu_pvalue": None,
        "ihc_mwu_effect": None,
        "ihc_mwu_n_tumor": int(tumor.size),
        "ihc_mwu_n_normal": int(normal.size),
        "ihc_mwu_significant": None,
    }
    if tumor.size == 0 or normal.size == 0:
        return empty
    try:
        from scipy.stats import mannwhitneyu

        result = mannwhitneyu(tumor, normal, alternative="two-sided")
        u_stat = float(result.statistic)
        p_value = float(result.pvalue)
    except Exception:
        u_stat, p_value = _mannwhitney_numpy(tumor, normal)
    n1 = float(tumor.size)
    n2 = float(normal.size)
    effect = (2.0 * u_stat) / (n1 * n2) - 1.0
    return {
        "ihc_mwu_u": u_stat,
        "ihc_mwu_pvalue": p_value,
        "ihc_mwu_effect": effect,
        "ihc_mwu_n_tumor": int(n1),
        "ihc_mwu_n_normal": int(n2),
        "ihc_mwu_significant": bool(p_value < 0.05) if p_value is not None and math.isfinite(p_value) else None,
    }


def _mannwhitney_numpy(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    from math import erfc

    n1 = x.size
    n2 = y.size
    ranks = pd.Series(np.concatenate([x, y])).rank(method="average").to_numpy()
    u = float(ranks[:n1].sum() - n1 * (n1 + 1) / 2.0)
    mu = n1 * n2 / 2.0
    sigma = math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12.0)
    if sigma == 0:
        return u, 1.0
    z = abs(u - mu) / sigma
    p = float(erfc(z / math.sqrt(2.0)))
    return u, min(max(p, 0.0), 1.0)


def drop_ihc_pseudocount_fold(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for col in ("protein_fold_induction", "protein_log2_fold", "ihc_fold_vs_negative"):
        if col in out.columns:
            out[col] = np.nan
    return out


def _as_float(value: Any) -> float | None:
    try:
        if value is None or value == "":
            return None
        out = float(value)
        if math.isnan(out) or math.isinf(out):
            return None
        return out
    except (TypeError, ValueError):
        return None
