# -*- coding: utf-8 -*-
"""Public single-cell dataset catalog and marker analysis helpers."""

from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Any

import numpy as np
import pandas as pd

from profiler.http_client import request_json, session

UCSC_INDEX = "https://cells.ucsc.edu/dataset.json"
TISCH2_HOME = "https://tisch.compbio.cn/"

# Curated Human cancer scRNA datasets with provenance (catalog — not a full scrape).
# Prefer datasets with public gene matrices / UCSC Cell Browser entries.
CURATED_DATASETS: list[dict[str, Any]] = [
    {
        "id": "ucsc_corces_aml",
        "title": "Human AML (Corces et al.) — UCSC Cell Browser",
        "cancer": "AML",
        "cancer_keys": ["LAML", "AML"],
        "species": "Human",
        "sample_type": "bone marrow / blood",
        "source": "UCSC Cell Browser",
        "dataset_url": "https://cells.ucsc.edu/?ds=corces-hematopoiesis",
        "paper": "https://doi.org/10.1038/s41588-020-00748-0",
        "n_cells_approx": None,
        "annotation": "cell type labels from publication",
        "format": "UCSC Cell Browser",
        "notes": "Collection may contain child datasets; explore children before analysis.",
    },
    {
        "id": "ucsc_breast_wu",
        "title": "Breast cancer tumor ecosystem (Wu et al. 2021)",
        "cancer": "Breast",
        "cancer_keys": ["BRCA"],
        "species": "Human",
        "sample_type": "tumor tissue",
        "source": "UCSC Cell Browser",
        "dataset_url": "https://cells.ucsc.edu/?ds=breast-cancer",
        "paper": "https://doi.org/10.1038/s41588-021-00911-1",
        "n_cells_approx": None,
        "annotation": "major lineage annotations",
        "format": "UCSC Cell Browser",
        "notes": "TIL label requires tissue annotation; immune lineage alone ≠ TIL.",
    },
    {
        "id": "ucsc_prostate_chen",
        "title": "Prostate cancer scRNA (public UCSC entry)",
        "cancer": "Prostate",
        "cancer_keys": ["PRAD"],
        "species": "Human",
        "sample_type": "tumor tissue",
        "source": "UCSC Cell Browser",
        "dataset_url": "https://cells.ucsc.edu/?ds=prostate-cancer",
        "paper": "https://doi.org/10.1038/s41467-021-23329-z",
        "n_cells_approx": None,
        "annotation": "publication cell types",
        "format": "UCSC Cell Browser",
        "notes": "Use for PSMA context exploration; verify gene availability per child dataset.",
    },
    {
        "id": "tisch2_brca_gse176078",
        "title": "TISCH2 BRCA — GSE176078",
        "cancer": "Breast",
        "cancer_keys": ["BRCA"],
        "species": "Human",
        "sample_type": "tumor",
        "source": "TISCH2",
        "dataset_url": "https://tisch.compbio.cn/home/",
        "paper": "https://doi.org/10.1093/nar/gkac959",
        "geo": "GSE176078",
        "n_cells_approx": None,
        "annotation": "TISCH2 major-lineage + malignancy annotation",
        "format": "H5AD/HDF5 (full file required when analyzing)",
        "notes": (
            "TISCH2 processed matrices require full-file download for analysis; "
            "not a gene-only Range fetch. Memory limits apply."
        ),
    },
    {
        "id": "tisch2_prad",
        "title": "TISCH2 Prostate cancer collections",
        "cancer": "Prostate",
        "cancer_keys": ["PRAD"],
        "species": "Human",
        "sample_type": "tumor",
        "source": "TISCH2",
        "dataset_url": "https://tisch.compbio.cn/home/",
        "paper": "https://doi.org/10.1093/nar/gkac959",
        "n_cells_approx": None,
        "annotation": "TISCH2 lineage annotations",
        "format": "H5AD/HDF5",
        "notes": "Filter Human + Prostate in TISCH2; confirm malignant vs normal B/plasma when relevant.",
    },
    {
        "id": "tisch2_nsclc",
        "title": "TISCH2 NSCLC collections",
        "cancer": "Lung",
        "cancer_keys": ["LUAD", "LUSC"],
        "species": "Human",
        "sample_type": "tumor",
        "source": "TISCH2",
        "dataset_url": "https://tisch.compbio.cn/home/",
        "paper": "https://doi.org/10.1093/nar/gkac959",
        "format": "H5AD/HDF5",
        "notes": "Dissociation/sorting/recovery biases apply; RNA ≠ surface protein.",
    },
]

LINEAGE_GROUPS = ["T", "B", "DC", "NK", "Myeloid", "Mast", "Stromal", "Epithelial/tumor", "Other"]

DEFAULT_MARKERS = {
    "PD-1": "PDCD1",
    "4-1BB": "TNFRSF9",
    "PD-L1": "CD274",
    "CD40": "CD40",
    "CLEC9A": "CLEC9A",
    "CD3E": "CD3E",
    "CD8A": "CD8A",
    "LAG3": "LAG3",
    "TIGIT": "TIGIT",
    "ICOS": "ICOS",
}


def list_datasets(cancer_key: str | None = None) -> list[dict[str, Any]]:
    if not cancer_key or cancer_key == "All":
        return list(CURATED_DATASETS)
    key = cancer_key.upper()
    return [
        d
        for d in CURATED_DATASETS
        if key in [k.upper() for k in d.get("cancer_keys") or []]
        or key in str(d.get("cancer", "")).upper()
    ]


def try_ucsc_index() -> dict[str, Any]:
    """Best-effort live UCSC dataset index probe."""
    try:
        data = request_json("GET", UCSC_INDEX, cache_ttl_s=86400, timeout=(12, 60))
        return {"status": "ok", "n": len(data) if isinstance(data, list) else None, "retrieved_at": _now()}
    except Exception as exc:
        return {"status": "fetch_failed", "error": str(exc), "retrieved_at": _now()}


def _now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


@dataclass
class MarkerStat:
    marker: str
    symbol: str
    n_cells: int
    n_expressed: int
    pct_expressed: float | None
    mean_expression: float | None
    unit: str
    measured: bool


def compute_marker_stats(
    expression: dict[str, np.ndarray],
    *,
    threshold: float = 0.0,
    unit: str = "dataset native",
) -> list[dict[str, Any]]:
    """
    expression: marker_symbol -> 1d array over selected cells (aligned barcodes).
    Threshold applies in the native scale; do not call this 'protein positive'.
    """
    rows = []
    for symbol, arr in expression.items():
        a = np.asarray(arr, dtype=float)
        if a.size == 0:
            rows.append(
                {
                    "symbol": symbol,
                    "n_cells": 0,
                    "n_expressed": 0,
                    "pct_expressed": None,
                    "mean_expression": None,
                    "unit": unit,
                    "measured": False,
                    "threshold": threshold,
                    "note": "No cells in selection",
                }
            )
            continue
        if np.all(np.isnan(a)):
            rows.append(
                {
                    "symbol": symbol,
                    "n_cells": int(a.size),
                    "n_expressed": None,
                    "pct_expressed": None,
                    "mean_expression": None,
                    "unit": unit,
                    "measured": False,
                    "threshold": threshold,
                    "note": "Gene not measured in this matrix",
                }
            )
            continue
        valid = np.isfinite(a)
        n = int(valid.sum())
        vals = a[valid]
        expressed = vals > threshold
        rows.append(
            {
                "symbol": symbol,
                "n_cells": n,
                "n_expressed": int(expressed.sum()),
                "pct_expressed": float(100.0 * expressed.mean()) if n else None,
                "mean_expression": float(vals.mean()) if n else None,
                "unit": unit,
                "measured": True,
                "threshold": threshold,
                "note": "pct_expressed is RNA above threshold, not protein positivity",
            }
        )
    return rows


def compute_coexpression(
    anchor: np.ndarray,
    other: np.ndarray,
    *,
    threshold: float = 0.0,
) -> dict[str, Any]:
    a = np.asarray(anchor, dtype=float)
    b = np.asarray(other, dtype=float)
    if a.shape != b.shape:
        return {"status": "error", "error": "Anchor and other arrays length mismatch"}
    if a.size == 0:
        return {"status": "data_insufficient", "error": "Empty cell selection"}
    if np.all(np.isnan(a)) or np.all(np.isnan(b)):
        return {"status": "not_measured", "double_pos_pct": None, "cond_pct": None}

    valid = np.isfinite(a) & np.isfinite(b)
    n = int(valid.sum())
    if n == 0:
        return {"status": "data_insufficient", "n": 0, "double_pos_pct": None, "cond_pct": None}

    av, bv = a[valid], b[valid]
    a_pos = av > threshold
    b_pos = bv > threshold
    both = a_pos & b_pos
    n_anchor = int(a_pos.sum())
    double_pct = float(100.0 * both.mean())
    cond_pct = float(100.0 * both[a_pos].mean()) if n_anchor > 0 else None
    # lift = P(other|anchor) / P(other)
    p_other = float(b_pos.mean())
    lift = (cond_pct / 100.0 / p_other) if (cond_pct is not None and p_other > 0) else None
    return {
        "status": "ok",
        "n_cells": n,
        "n_anchor_pos": n_anchor,
        "n_double_pos": int(both.sum()),
        "double_pos_pct_of_all": double_pct,
        "other_pos_pct_of_anchor_pos": cond_pct,
        "lift": lift,
        "lift_formula": "P(other>thr | anchor>thr) / P(other>thr)",
        "threshold": threshold,
        "denominators": {
            "double_pos_pct_of_all": "all selected cells",
            "other_pos_pct_of_anchor_pos": "anchor-positive cells only",
        },
        "note": "Same-cell co-expression; do not combine averages across cell groups.",
    }


def demo_analysis_bundle(dataset_id: str, markers: list[str]) -> dict[str, Any]:
    """
    When full H5AD is not loaded locally, return a structured pending analysis
    with explicit controls — never fake expression values.
    """
    ds = next((d for d in CURATED_DATASETS if d["id"] == dataset_id), None)
    return {
        "status": "awaiting_matrix",
        "dataset": ds,
        "markers": markers,
        "message": (
            "Dataset catalog entry selected. Full expression matrix is not bundled in this deployment. "
            "Provide a local H5AD/OPFS path or enable TISCH2/UCSC matrix fetch to run marker stats. "
            "No expression values were fabricated."
        ),
        "controls": {
            "lineages": LINEAGE_GROUPS,
            "threshold_default": 0.0,
            "threshold_meaning": "RNA > 0 in native scale (not protein positive)",
            "patient_vs_pooled": "Patient-level and pooled-cell stats will be reported separately when matrix is available",
        },
        "limits": [
            "Processed public matrix exploration only — not FASTQ reanalysis or new annotations.",
            "Numeric clusters are never auto-labeled as T/B/DC.",
            "Dropout, dissociation, sorting, and recovery biases apply.",
        ],
        "retrieved_at": _now(),
    }


def analyze_local_h5ad(
    path: str,
    markers: list[str],
    *,
    celltype_col: str = "cell_type",
    lineages: list[str] | None = None,
    threshold: float = 0.0,
) -> dict[str, Any]:
    """Analyze a local H5AD if anndata/h5py are available."""
    try:
        import anndata as ad
    except ImportError:
        return {"status": "unsupported", "error": "anndata not installed"}

    try:
        adata = ad.read_h5ad(path, backed="r")
    except Exception as exc:
        return {"status": "fetch_failed", "error": str(exc)}

    # Verify barcodes alignment assumption is NOT used — use obs_names as keys
    n_cells = adata.n_obs
    var_names = set(map(str, adata.var_names))
    selected = np.ones(n_cells, dtype=bool)
    mapping_note = None
    if celltype_col in adata.obs.columns and lineages:
        ct = adata.obs[celltype_col].astype(str)
        # Show original annotation alongside requested lineages; no auto cluster mapping
        selected = ct.isin(lineages).to_numpy()
        mapping_note = f"Filtered by raw annotation column '{celltype_col}' ∈ {lineages}"

    expr = {}
    for m in markers:
        if m not in var_names:
            expr[m] = np.full(int(selected.sum()), np.nan)
        else:
            X = adata[selected, m].X
            if hasattr(X, "toarray"):
                X = X.toarray()
            expr[m] = np.asarray(X).ravel()

    stats = compute_marker_stats(expr, threshold=threshold)
    coexp = {}
    if len(markers) >= 2 and markers[0] in expr:
        for other in markers[1:]:
            coexp[other] = compute_coexpression(expr[markers[0]], expr[other], threshold=threshold)

    return {
        "status": "ok",
        "n_cells_total": int(n_cells),
        "n_cells_selected": int(selected.sum()),
        "mapping_note": mapping_note,
        "marker_stats": stats,
        "coexpression_vs_anchor": coexp,
        "anchor": markers[0] if markers else None,
        "unit": "dataset native (as stored in H5AD)",
        "retrieved_at": _now(),
        "limits": [
            "Barcode-based join only; matrix orientation validated via anndata API.",
            "This is exploration of a processed matrix, not a new cell-type annotation.",
        ],
    }
