# -*- coding: utf-8 -*-
"""Normalize expression and compute cancer-level IHC / fold-induction metrics."""

from __future__ import annotations

import io
import math
import re
import time
import zipfile
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import requests

from csv_io import read_csv_safe
from data_qc import expand_ihc_counts, mannwhitney_ihc

HPA_CANCER_IHC_ZIP = "https://www.proteinatlas.org/download/tsv/cancer_data.tsv.zip"
CACHE_DIR = Path(__file__).resolve().parent / "cache"
CANCER_IHC_CACHE = CACHE_DIR / "hpa_cancer_data.tsv"
TIMEOUT = 120

IHC_LEVELS = ("Not detected", "Low", "Medium", "High")
IHC_SCORE_MAP = {"Not detected": 0, "Low": 1, "Medium": 2, "High": 3}


def _normalize_label(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", str(text or "").lower()).strip()


def _safe_float(value: Any) -> float | None:
    try:
        if value is None or value == "":
            return None
        out = float(value)
        if math.isnan(out) or math.isinf(out):
            return None
        return out
    except (TypeError, ValueError):
        return None


def ensure_cancer_ihc_cache(session: requests.Session | None = None) -> Path:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    if CANCER_IHC_CACHE.exists() and CANCER_IHC_CACHE.stat().st_size > 1000:
        return CANCER_IHC_CACHE

    sess = session or requests.Session()
    sess.headers.setdefault("User-Agent", "Mozilla/5.0 (compatible; ImmunoTargetExplorer/1.0; research use)")
    sess.headers.setdefault("Connection", "close")
    last = None
    for attempt in range(4):
        try:
            response = sess.get(HPA_CANCER_IHC_ZIP, timeout=(12, 120))
            response.raise_for_status()
            break
        except (requests.ConnectionError, requests.Timeout, requests.exceptions.ChunkedEncodingError) as exc:
            last = exc
            time.sleep(0.7 * (2 ** attempt))
    else:
        raise ConnectionError("Could not download HPA IHC cache. Search again in a moment.") from last
    with zipfile.ZipFile(io.BytesIO(response.content)) as zf:
        inner = zf.namelist()[0]
        CANCER_IHC_CACHE.write_bytes(zf.read(inner))
    return CANCER_IHC_CACHE


def fetch_cancer_ihc_rows(ensembl_id: str, session: requests.Session | None = None) -> pd.DataFrame:
    path = ensure_cancer_ihc_cache(session)
    df = read_csv_safe(path, sep="\t", dtype=str)
    subset = df[df["Gene"].astype(str).str.upper() == ensembl_id.upper()].copy()
    if subset.empty:
        return pd.DataFrame(
            columns=[
                "cancer_type",
                "ihc_high",
                "ihc_medium",
                "ihc_low",
                "ihc_not_detected",
                "patients_total",
            ]
        )

    for col in ["High", "Medium", "Low", "Not detected"]:
        subset[col] = pd.to_numeric(subset[col], errors="coerce").fillna(0).astype(int)

    subset["patients_total"] = subset[["High", "Medium", "Low", "Not detected"]].sum(axis=1)
    subset = subset.rename(
        columns={
            "Cancer": "cancer_type",
            "High": "ihc_high",
            "Medium": "ihc_medium",
            "Low": "ihc_low",
            "Not detected": "ihc_not_detected",
        }
    )
    return subset[
        [
            "cancer_type",
            "ihc_high",
            "ihc_medium",
            "ihc_low",
            "ihc_not_detected",
            "patients_total",
        ]
    ].reset_index(drop=True)


def _dominant_ihc_score(row: pd.Series) -> int:
    counts = {
        0: int(row.get("ihc_not_detected") or 0),
        1: int(row.get("ihc_low") or 0),
        2: int(row.get("ihc_medium") or 0),
        3: int(row.get("ihc_high") or 0),
    }
    return max(counts, key=counts.get)


def _weighted_ihc_score(row: pd.Series) -> float:
    total = float(row.get("patients_total") or 0)
    if total <= 0:
        return 0.0
    weighted = (
        int(row.get("ihc_high") or 0) * 3
        + int(row.get("ihc_medium") or 0) * 2
        + int(row.get("ihc_low") or 0) * 1
    )
    return weighted / total


def _ihc_rank_stats(row: pd.Series, normal_scores) -> dict:
    tumor = expand_ihc_counts(
        row.get("ihc_high"),
        row.get("ihc_medium"),
        row.get("ihc_low"),
        row.get("ihc_not_detected"),
    )
    return mannwhitney_ihc(tumor, normal_scores)


def _extract_hpa_cancer_rna(hpa_summary: dict[str, Any]) -> dict[str, float]:
    mapping = hpa_summary.get("RNA cell line specific nTPM") or {}
    if not isinstance(mapping, dict):
        return {}
    out: dict[str, float] = {}
    for cancer, value in mapping.items():
        num = _safe_float(value)
        if num is not None:
            out[str(cancer)] = num
    return out


def _extract_normal_baseline(expr_df: pd.DataFrame) -> float:
    if expr_df is None or expr_df.empty:
        return 1.0
    work = expr_df.copy()
    work["value_num"] = pd.to_numeric(work["expression_value"], errors="coerce")
    work = work.dropna(subset=["value_num"])
    work = work[work["condition"].astype(str).str.contains("normal", case=False, na=False)]
    if work.empty:
        return 1.0
    baseline = float(work["value_num"].quantile(0.10))
    return max(baseline, 1.0)


def _match_cancer_rna(cancer_type: str, rna_map: dict[str, float]) -> float | None:
    target = _normalize_label(cancer_type)
    if not rna_map:
        return None
    for key, value in rna_map.items():
        key_norm = _normalize_label(key)
        if target == key_norm or target in key_norm or key_norm in target:
            return value
    return None


def _match_association_score(cancer_type: str, assoc_df: pd.DataFrame) -> float | None:
    if assoc_df is None or assoc_df.empty:
        return None
    target = _normalize_label(cancer_type)
    best_score = None
    for rec in assoc_df.to_dict("records"):
        disease = _normalize_label(rec.get("disease"))
        if not disease:
            continue
        if target in disease or disease in target:
            score = _safe_float(rec.get("association_score"))
            if score is not None and (best_score is None or score > best_score):
                best_score = score
    return best_score


def _zscore(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    std = float(values.std(ddof=0))
    if std == 0 or math.isnan(std):
        return pd.Series([0.0] * len(values), index=series.index)
    mean = float(values.mean())
    return (values - mean) / std


def build_cancer_analytics_table(
    *,
    symbol: str,
    ensembl_id: str,
    expr_df: pd.DataFrame,
    assoc_df: pd.DataFrame,
    hpa_summary: dict[str, Any],
    session: requests.Session | None = None,
) -> pd.DataFrame:
    ihc_df = fetch_cancer_ihc_rows(ensembl_id, session=session)
    if ihc_df.empty:
        return pd.DataFrame()

    from expression_compare import CANCER_CONTEXT, _ensure_hpa_file, _filter_gene_tsv, _normal_ihc_score_array

    rna_map = _extract_hpa_cancer_rna(hpa_summary)
    normal_baseline = _extract_normal_baseline(expr_df)
    sess = session or requests.Session()
    try:
        normal_ihc = _filter_gene_tsv(_ensure_hpa_file("normal_ihc", sess), ensembl_id)
    except Exception:
        normal_ihc = pd.DataFrame()

    rows: list[dict[str, Any]] = []
    for rec in ihc_df.to_dict("records"):
        row = pd.Series(rec)
        cancer_type = rec["cancer_type"]
        organ = (CANCER_CONTEXT.get(str(cancer_type).lower()) or {}).get("organ")
        normal_scores = _normal_ihc_score_array(normal_ihc, organ) if organ else np.asarray([], dtype=float)
        mwu = _ihc_rank_stats(row, normal_scores)
        rna_ntpm = _match_cancer_rna(cancer_type, rna_map)
        rna_log2 = math.log2(rna_ntpm + 1.0) if rna_ntpm is not None else None

        rows.append(
            {
                "gene_symbol": symbol,
                "ensembl_id": ensembl_id,
                "cancer_type": cancer_type,
                "ihc_score_0_3": round(_dominant_ihc_score(row), 2),
                "ihc_weighted_0_3": round(_weighted_ihc_score(row), 3),
                "ihc_high": int(rec["ihc_high"]),
                "ihc_medium": int(rec["ihc_medium"]),
                "ihc_low": int(rec["ihc_low"]),
                "ihc_not_detected": int(rec["ihc_not_detected"]),
                "patients_total": int(rec["patients_total"]),
                "ihc_mwu_u": mwu["ihc_mwu_u"],
                "ihc_mwu_pvalue": mwu["ihc_mwu_pvalue"],
                "ihc_mwu_effect": mwu["ihc_mwu_effect"],
                "ihc_mwu_significant": mwu["ihc_mwu_significant"],
                "rna_ntpm": round(rna_ntpm, 3) if rna_ntpm is not None else None,
                "rna_log2_tpm1": round(rna_log2, 3) if rna_log2 is not None else None,
                "normal_baseline_tpm": round(normal_baseline, 3),
                "ot_association_score": _match_association_score(cancer_type, assoc_df),
            }
        )

    out = pd.DataFrame(rows)
    out["rna_zscore"] = _zscore(out["rna_ntpm"]).round(3)
    out["ihc_zscore"] = _zscore(out["ihc_weighted_0_3"]).round(3)
    neglogp = (-np.log10(pd.to_numeric(out["ihc_mwu_pvalue"], errors="coerce").clip(lower=1e-300))).fillna(0)
    out["combined_priority"] = (
        out["ihc_weighted_0_3"].fillna(0) * 0.45
        + neglogp * 0.20
        + out["rna_zscore"].fillna(0).clip(lower=0) * 0.25
        + out["ot_association_score"].fillna(0) * 0.10
    ).round(3)
    out = out.sort_values(["combined_priority", "ihc_weighted_0_3"], ascending=False).reset_index(drop=True)
    return out
