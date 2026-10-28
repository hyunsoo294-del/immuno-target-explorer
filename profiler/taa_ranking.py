# -*- coding: utf-8 -*-
"""Relative co-elevation Top 20 TAA ranking within a cohort."""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Any

import numpy as np
import pandas as pd

from profiler.cbioportal import fetch_molecular_data, resolve_rna_profile, sample_list_id_for_study
from profiler.gene_resolution import TAA_PANEL_FOR_RANKING, fetch_genes_batch
from profiler.http_client import request_json

# Minimum samples per quartile group
MIN_GROUP_N = 10


def _log2p1(values: np.ndarray) -> np.ndarray:
    return np.log2(np.asarray(values, dtype=float) + 1.0)


def rank_coelevated_taas(
    study_id: str,
    reference_symbol: str,
    reference_entrez: int,
    *,
    panel: list[str] | None = None,
) -> dict[str, Any]:
    """
    Define high/low quartiles of the reference TAA in the same RNA profile,
    then rank other panel TAAs by mean(log2(expr+1)) difference (high − low).

    This is exploratory co-elevation, NOT treatment-induced increase, NOT
    statistical significance testing, and NOT fold-change of raw means.
    """
    retrieved = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    profile = resolve_rna_profile(study_id)
    if not profile:
        return {
            "status": "profile_missing",
            "error": f"No RNA profile for {study_id}",
            "rows": [],
            "retrieved_at": retrieved,
        }

    ref = fetch_molecular_data(
        profile,
        reference_entrez,
        sample_list_id_for_study(study_id),
        study_id=study_id,
        symbol=reference_symbol,
        assay_key="tcga_rna",
        unit="RNA-seq native; analysis uses log2(value+1)",
        transform="log2_plus1_for_ranking",
    )
    if ref.status == "fetch_failed":
        return {"status": "fetch_failed", "error": ref.error, "rows": [], "retrieved_at": retrieved}
    if ref.n_measured < 2 * MIN_GROUP_N:
        return {
            "status": "held",
            "error": f"Insufficient samples for quartile split (n={ref.n_measured}, need ≥{2*MIN_GROUP_N})",
            "rows": [],
            "retrieved_at": retrieved,
        }

    ref_df = ref.to_dataframe()[["sampleId", "value"]].dropna()
    ref_df = ref_df.drop_duplicates(subset=["sampleId"], keep="first")
    q1 = float(ref_df["value"].quantile(0.25))
    q3 = float(ref_df["value"].quantile(0.75))
    if q1 == q3:
        return {
            "status": "held",
            "error": "Upper and lower quartiles are not separable (tied boundaries)",
            "rows": [],
            "retrieved_at": retrieved,
            "q1": q1,
            "q3": q3,
        }

    # Inclusive boundaries with ties assigned consistently; exclude middle.
    low = ref_df[ref_df["value"] <= q1]
    high = ref_df[ref_df["value"] >= q3]
    # Remove overlap if any sample satisfies both (flat distribution)
    overlap = set(low["sampleId"]) & set(high["sampleId"])
    if overlap:
        low = low[~low["sampleId"].isin(overlap)]
        high = high[~high["sampleId"].isin(overlap)]
    if len(low) < MIN_GROUP_N or len(high) < MIN_GROUP_N:
        return {
            "status": "held",
            "error": f"Group sizes below minimum after tie handling (low={len(low)}, high={len(high)})",
            "rows": [],
            "retrieved_at": retrieved,
            "q1": q1,
            "q3": q3,
        }

    low_ids = set(low["sampleId"])
    high_ids = set(high["sampleId"])
    candidates = [s for s in (panel or TAA_PANEL_FOR_RANKING) if s.upper() != reference_symbol.upper()]

    # Resolve entrez IDs in batch
    genes = fetch_genes_batch(candidates)
    symbol_to_entrez = {
        g["hugoGeneSymbol"]: int(g["entrezGeneId"])
        for g in genes
        if g.get("hugoGeneSymbol") and g.get("entrezGeneId") is not None
    }

    rows: list[dict[str, Any]] = []
    for symbol in candidates:
        entrez = symbol_to_entrez.get(symbol)
        if entrez is None:
            continue
        other = fetch_molecular_data(
            profile,
            entrez,
            sample_list_id_for_study(study_id),
            study_id=study_id,
            symbol=symbol,
            assay_key="tcga_rna",
            unit="log2(expression+1)",
            transform="log2_plus1",
        )
        if other.status == "fetch_failed" or other.n_measured == 0:
            continue
        odf = other.to_dataframe()[["sampleId", "value"]].drop_duplicates("sampleId")
        omap = dict(zip(odf["sampleId"], odf["value"]))
        high_vals = [omap[s] for s in high_ids if s in omap]
        low_vals = [omap[s] for s in low_ids if s in omap]
        if len(high_vals) < MIN_GROUP_N or len(low_vals) < MIN_GROUP_N:
            continue
        high_log = _log2p1(np.array(high_vals))
        low_log = _log2p1(np.array(low_vals))
        delta = float(np.mean(high_log) - np.mean(low_log))
        if delta <= 0:
            continue
        rows.append(
            {
                "rank": 0,
                "candidate": symbol,
                "official_symbol": symbol,
                "entrez_gene_id": entrez,
                "n_high": len(high_vals),
                "n_low": len(low_vals),
                "mean_log2p1_high": float(np.mean(high_log)),
                "mean_log2p1_low": float(np.mean(low_log)),
                "delta_mean_log2p1": delta,
                "unit": "mean log2(expression+1) difference (high − low quartile of reference)",
                "cohort": study_id,
                "profile": profile,
                "reference_symbol": reference_symbol,
                "retrieved_at": retrieved,
            }
        )

    rows.sort(key=lambda r: r["delta_mean_log2p1"], reverse=True)
    for i, r in enumerate(rows[:20], start=1):
        r["rank"] = i
    top = rows[:20]

    return {
        "status": "ok",
        "rows": top,
        "n_candidates_tested": len(candidates),
        "n_positive_delta": len(rows),
        "reference_symbol": reference_symbol,
        "reference_entrez": reference_entrez,
        "study_id": study_id,
        "profile": profile,
        "q1": q1,
        "q3": q3,
        "n_high_ref": len(high_ids),
        "n_low_ref": len(low_ids),
        "retrieved_at": retrieved,
        "interpretation_limits": [
            "Shows genes relatively higher in samples where the reference TAA is high (same untreated cohort).",
            "Not a claim that the reference decreased and this target increased.",
            "Not post-Enhertu / treatment-induced dynamics.",
            "No statistical significance testing was performed.",
            "Metric is difference of mean log2(expression+1), not a fold-change of raw means.",
        ],
    }
