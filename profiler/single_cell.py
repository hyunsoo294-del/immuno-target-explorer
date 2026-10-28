# -*- coding: utf-8 -*-
"""Public single-cell catalog + real UCSC Cell Browser matrix analysis."""

from __future__ import annotations

import gzip
import json
import re
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd

from profiler.http_client import request_json, session

UCSC_INDEX = "https://cells.ucsc.edu/dataset.json"
UCSC_BASE = "https://cells.ucsc.edu"
CACHE_DIR = Path(__file__).resolve().parent.parent / "cache" / "singlecell"
TISCH2_HOME = "https://tisch.compbio.cn/"

# Verified leaf datasets with public exprMatrix.tsv.gz + meta.tsv on UCSC.
CURATED_DATASETS: list[dict[str, Any]] = [
    {
        "id": "ucsc_head_neck",
        "ucsc_path": "head-neck",
        "title": "Head and Neck Cancer (Puram et al.) — UCSC",
        "cancer": "Head and neck",
        "cancer_keys": ["HNSC"],
        "species": "Human",
        "sample_type": "tumor tissue",
        "source": "UCSC Cell Browser",
        "dataset_url": "https://cells.ucsc.edu/?ds=head-neck",
        "paper": "https://doi.org/10.1016/j.cell.2017.10.044",
        "n_cells_approx": 3608,
        "annotation_field": "cell type",
        "lineage_field": "cell_category",
        "format": "UCSC exprMatrix.tsv.gz + meta.tsv",
        "matrix_mb_approx": 85,
        "analyzable": True,
        "notes": "Patient tumors; immune lineage ≠ automatically TIL without tissue context.",
    },
    {
        "id": "ucsc_quake_gbm",
        "ucsc_path": "quake-gbm",
        "title": "Glioblastoma infiltrating vs tumor core — UCSC",
        "cancer": "Glioblastoma",
        "cancer_keys": ["GBM"],
        "species": "Human",
        "sample_type": "tumor tissue",
        "source": "UCSC Cell Browser",
        "dataset_url": "https://cells.ucsc.edu/?ds=quake-gbm",
        "paper": "https://doi.org/10.1016/j.celrep.2019.04.063",
        "n_cells_approx": 3588,
        "annotation_field": "Cell Type",
        "format": "UCSC exprMatrix.tsv.gz + meta.tsv",
        "matrix_mb_approx": 19,
        "analyzable": True,
        "notes": "Processed public matrix exploration only.",
    },
    {
        "id": "ucsc_gbm",
        "ucsc_path": "gbm",
        "title": "Glioblastoma (Neftel et al.) — UCSC",
        "cancer": "Glioblastoma",
        "cancer_keys": ["GBM"],
        "species": "Human",
        "sample_type": "tumor tissue",
        "source": "UCSC Cell Browser",
        "dataset_url": "https://cells.ucsc.edu/?ds=gbm",
        "paper": "https://doi.org/10.1016/j.cell.2019.06.024",
        "n_cells_approx": 32877,
        "annotation_field": "Cell Type Assignment",
        "format": "UCSC exprMatrix.tsv.gz + meta.tsv",
        "matrix_mb_approx": 595,
        "analyzable": True,
        "notes": "Large matrix (~600MB). First run downloads once to local cache.",
    },
    {
        "id": "ucsc_colorectal",
        "ucsc_path": "colorectal-cancer",
        "title": "Colorectal Cancer — UCSC",
        "cancer": "Colorectal",
        "cancer_keys": ["COADREAD"],
        "species": "Human",
        "sample_type": "tumor tissue",
        "source": "UCSC Cell Browser",
        "dataset_url": "https://cells.ucsc.edu/?ds=colorectal-cancer",
        "paper": "https://cells.ucsc.edu/?ds=colorectal-cancer",
        "n_cells_approx": 63689,
        "annotation_field": "Cell_subtype",
        "format": "UCSC exprMatrix.tsv.gz + meta.tsv",
        "matrix_mb_approx": 1026,
        "analyzable": True,
        "notes": "Very large (~1GB). Prefer smaller datasets unless needed.",
    },
    {
        "id": "ucsc_human_pancreas",
        "ucsc_path": "human-pancreas",
        "title": "Human Pancreas — UCSC",
        "cancer": "Pancreas",
        "cancer_keys": ["PAAD"],
        "species": "Human",
        "sample_type": "pancreas (mixed / reference)",
        "source": "UCSC Cell Browser",
        "dataset_url": "https://cells.ucsc.edu/?ds=human-pancreas",
        "n_cells_approx": 7714,
        "annotation_field": "CellType",
        "format": "UCSC exprMatrix.tsv.gz + meta.tsv",
        "matrix_mb_approx": 68,
        "analyzable": True,
        "notes": "Confirm disease annotation before treating as tumor-only.",
    },
    {
        "id": "tisch2_catalog_brca",
        "ucsc_path": None,
        "title": "TISCH2 BRCA collections (catalog only until H5AD provided)",
        "cancer": "Breast",
        "cancer_keys": ["BRCA"],
        "species": "Human",
        "sample_type": "tumor",
        "source": "TISCH2",
        "dataset_url": TISCH2_HOME,
        "paper": "https://doi.org/10.1093/nar/gkac959",
        "format": "H5AD/HDF5",
        "analyzable": False,
        "notes": "Requires full processed H5AD download; use local path analysis when file is available.",
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

# Keyword → lineage. Numeric clusters are NEVER auto-mapped.
_LINEAGE_RULES: list[tuple[str, list[str]]] = [
    ("T", [" t cell", "tcell", "cd4", "cd8", "treg", "t helper", "cytotoxic t", "nk t", "nkt", " t ", "thymocyte"]),
    ("B", [" b cell", "bcell", "plasma", "plasmablast", "b_cell", "b-cell"]),
    ("DC", ["dendritic", " cd c", "cdc1", "cdc2", "pdc", "clec9"]),
    ("NK", ["nk cell", "natural killer", "nkcell"]),
    ("Myeloid", ["monocyte", "macrophage", "myeloid", "neutrophil", "mdsc", "kupffer", "microglia"]),
    ("Mast", ["mast"]),
    ("Stromal", ["fibroblast", "caf", "endotheli", "pericyte", "stromal", "smooth muscle", "mesenchymal"]),
    ("Epithelial/tumor", ["epithelial", "tumor", "malignant", "carcinoma", "cancer cell", "neoplastic", "ductal", "acinar", "hepatocyte", "keratinocyte"]),
]


def _now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


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


def get_dataset(dataset_id: str) -> dict[str, Any] | None:
    return next((d for d in CURATED_DATASETS if d["id"] == dataset_id), None)


def try_ucsc_index() -> dict[str, Any]:
    try:
        data = request_json("GET", UCSC_INDEX, cache_ttl_s=86400, timeout=(12, 60))
        n = len(data.get("datasets") or []) if isinstance(data, dict) else (len(data) if isinstance(data, list) else None)
        return {"status": "ok", "n": n, "retrieved_at": _now()}
    except Exception as exc:
        return {"status": "fetch_failed", "error": str(exc), "retrieved_at": _now()}


def map_lineage(label: str) -> str:
    """Map free-text annotation to coarse lineage; unknown → Other. No numeric-cluster guessing."""
    raw = str(label or "").strip()
    if not raw:
        return "Other"
    if re.fullmatch(r"\d+", raw):
        return "Other"  # never auto-label numeric clusters
    text = f" {raw.lower().replace('_', ' ').replace('-', ' ')} "
    for lineage, keys in _LINEAGE_RULES:
        if any(k in text for k in keys):
            return lineage
    return "Other"


def compute_marker_stats(
    expression: dict[str, np.ndarray],
    *,
    threshold: float = 0.0,
    unit: str = "dataset native",
) -> list[dict[str, Any]]:
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


def _download(url: str, dest: Path, progress: Callable | None = None) -> Path:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        return dest
    tmp = dest.with_suffix(dest.suffix + ".part")
    sess = session()
    with sess.get(url, stream=True, timeout=(30, 600)) as resp:
        resp.raise_for_status()
        total = int(resp.headers.get("content-length") or 0)
        done = 0
        with open(tmp, "wb") as f:
            for chunk in resp.iter_content(1024 * 1024):
                if not chunk:
                    continue
                f.write(chunk)
                done += len(chunk)
                if progress and total:
                    progress(f"Downloading {dest.name}: {100.0 * done / total:.1f}%")
    tmp.replace(dest)
    return dest


def ensure_ucsc_files(ucsc_path: str, progress: Callable | None = None) -> dict[str, Path]:
    meta = CACHE_DIR / f"{ucsc_path.replace('/', '_')}_meta.tsv"
    matrix = CACHE_DIR / f"{ucsc_path.replace('/', '_')}_exprMatrix.tsv.gz"
    dsjson = CACHE_DIR / f"{ucsc_path.replace('/', '_')}_dataset.json"
    _download(f"{UCSC_BASE}/{ucsc_path}/meta.tsv", meta, progress=progress)
    _download(f"{UCSC_BASE}/{ucsc_path}/exprMatrix.tsv.gz", matrix, progress=progress)
    try:
        _download(f"{UCSC_BASE}/{ucsc_path}/dataset.json", dsjson, progress=progress)
    except Exception:
        pass
    return {"meta": meta, "matrix": matrix, "dataset_json": dsjson}


def load_meta(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    # First column is usually cell barcode / index
    if df.columns[0].lower() in ("index", "cell", "cellid", "barcode", "cell_id"):
        df = df.rename(columns={df.columns[0]: "barcode"})
    else:
        df = df.rename(columns={df.columns[0]: "barcode"})
    return df


def _normalize_gene_token(token: str) -> str:
    t = str(token).strip()
    # Common UCSC forms: SYMBOL, SYMBOL|ENSG, ENSG|SYMBOL
    parts = re.split(r"[|_\s]+", t)
    for p in parts:
        if re.fullmatch(r"ENSG\d+", p, flags=re.I):
            continue
        if p:
            return p.upper()
    return t.upper()


def extract_gene_rows(
    matrix_path: Path,
    genes: list[str],
    progress: Callable | None = None,
) -> tuple[list[str], dict[str, np.ndarray], dict[str, Any]]:
    """Single-pass gzip scan; extract only requested gene rows. Genes are columns after barcode header."""
    wanted = {g.upper(): g for g in genes}
    found: dict[str, np.ndarray] = {}
    barcodes: list[str] = []
    t0 = time.time()
    with gzip.open(matrix_path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        if not header:
            raise RuntimeError("Empty expression matrix")
        # header[0] is gene label; remaining are cell barcodes
        barcodes = header[1:]
        n_cells = len(barcodes)
        for i, line in enumerate(f, start=1):
            if not line:
                continue
            gene_tok, _, rest = line.partition("\t")
            sym = _normalize_gene_token(gene_tok)
            if sym not in wanted:
                continue
            vals = np.fromstring(rest, sep="\t", dtype=float)
            if vals.size != n_cells:
                # pad/truncate defensively but record
                fixed = np.full(n_cells, np.nan, dtype=float)
                n = min(n_cells, vals.size)
                fixed[:n] = vals[:n]
                vals = fixed
            # Keep first occurrence; do not silently average duplicates
            key = wanted[sym]
            if key not in found:
                found[key] = vals
            if progress and i % 5000 == 0:
                progress(f"Scanning matrix row {i}; found {len(found)}/{len(wanted)} genes")
            if len(found) == len(wanted):
                break
    meta = {
        "n_cells": len(barcodes),
        "n_genes_requested": len(genes),
        "n_genes_found": len(found),
        "genes_missing": [g for g in genes if g not in found],
        "seconds": round(time.time() - t0, 2),
        "matrix_path": str(matrix_path),
    }
    # Fill missing genes with NaN vectors so downstream can mark not_measured
    for g in genes:
        if g not in found:
            found[g] = np.full(len(barcodes), np.nan, dtype=float)
    return barcodes, found, meta


def pick_annotation_column(meta: pd.DataFrame, preferred: str | None = None) -> str:
    if preferred and preferred in meta.columns:
        return preferred
    candidates = [
        "cell type", "Cell Type", "CellType", "cell_type", "Cell_subtype",
        "cell_category", "Cell Type Assignment", "cluster", "annotation",
        "manual_celltype_annotation_high", "author_cell_type",
    ]
    for c in candidates:
        if c in meta.columns:
            return c
    # fallback: first non-barcode categorical-like column
    for c in meta.columns:
        if c != "barcode":
            return c
    raise RuntimeError("No annotation column in meta")


def composition_tables(meta: pd.DataFrame, ann_col: str) -> dict[str, Any]:
    work = meta.copy()
    work["lineage"] = work[ann_col].map(map_lineage)
    n_total = len(work)
    by_ann = (
        work[ann_col]
        .fillna("Unknown")
        .value_counts(dropna=False)
        .rename_axis("annotation")
        .reset_index(name="n_cells")
    )
    by_ann["pct_of_all_cells"] = 100.0 * by_ann["n_cells"] / n_total if n_total else None
    by_ann["lineage"] = by_ann["annotation"].map(map_lineage)

    by_lin = (
        work["lineage"].value_counts().rename_axis("lineage").reset_index(name="n_cells")
    )
    by_lin["pct_of_all_cells"] = 100.0 * by_lin["n_cells"] / n_total if n_total else None
    immune = {"T", "B", "DC", "NK", "Myeloid", "Mast"}
    n_immune = int(work["lineage"].isin(immune).sum())
    immune_share = (
        work.loc[work["lineage"].isin(immune), "lineage"]
        .value_counts()
        .rename_axis("lineage")
        .reset_index(name="n_cells")
    )
    if n_immune:
        immune_share["pct_of_immune_cells"] = 100.0 * immune_share["n_cells"] / n_immune
    else:
        immune_share["pct_of_immune_cells"] = None

    return {
        "n_cells_total": n_total,
        "n_immune_cells": n_immune,
        "pct_immune_of_all": (100.0 * n_immune / n_total) if n_total else None,
        "by_annotation": by_ann.to_dict(orient="records"),
        "by_lineage": by_lin.to_dict(orient="records"),
        "immune_lineage_share": immune_share.to_dict(orient="records"),
        "denominator_notes": [
            "pct_of_all_cells uses all recovered cells in the processed matrix as denominator.",
            "pct_of_immune_cells uses cells mapped to T/B/DC/NK/Myeloid/Mast only.",
            "Lineage mapping is keyword-based on author annotations; numeric clusters are not auto-labeled.",
            "Immune lineage ≠ TIL unless tissue/compartment annotation says so.",
        ],
    }


def analyze_ucsc_dataset(
    dataset_id: str,
    markers: list[str],
    *,
    lineages: list[str] | None = None,
    threshold: float = 0.0,
    anchor: str | None = None,
    progress: Callable | None = None,
) -> dict[str, Any]:
    ds = get_dataset(dataset_id)
    if not ds:
        return {"status": "error", "error": f"Unknown dataset id: {dataset_id}"}
    if not ds.get("analyzable") or not ds.get("ucsc_path"):
        return {
            "status": "awaiting_matrix",
            "dataset": ds,
            "message": ds.get("notes") or "Dataset is catalog-only; provide local H5AD or choose an analyzable UCSC dataset.",
            "retrieved_at": _now(),
        }

    if progress:
        progress(f"Fetching UCSC files for {ds['ucsc_path']}…")
    files = ensure_ucsc_files(ds["ucsc_path"], progress=progress)
    meta = load_meta(files["meta"])
    # Align barcodes: matrix header cells must match meta barcodes
    # UCSC matrices use same cell ids as meta index
    ann_col = pick_annotation_column(meta, ds.get("annotation_field") or ds.get("lineage_field"))
    comp = composition_tables(meta, ann_col)

    if progress:
        progress(f"Extracting {len(markers)} gene rows from matrix…")
    barcodes, expr, gene_meta = extract_gene_rows(files["matrix"], markers, progress=progress)

    # Join by barcode — do NOT assume order
    meta_idx = meta.set_index("barcode", drop=False)
    order = []
    missing_bc = 0
    for bc in barcodes:
        if bc in meta_idx.index:
            order.append(bc)
        else:
            missing_bc += 1
            order.append(None)

    if missing_bc:
        # Fall back to positional align only if lengths match and missing is total (rare)
        if len(barcodes) == len(meta) and missing_bc == len(barcodes):
            order_flag = "positional_fallback_same_length"
            aligned_meta = meta.reset_index(drop=True)
            select_mask = np.ones(len(meta), dtype=bool)
            expr_aligned = expr
        else:
            return {
                "status": "error",
                "error": f"Barcode join failed: {missing_bc}/{len(barcodes)} matrix barcodes missing from meta",
                "dataset": ds,
            }
    else:
        order_flag = "barcode_join"
        aligned_meta = meta_idx.loc[order].reset_index(drop=True)
        select_mask = np.ones(len(aligned_meta), dtype=bool)
        expr_aligned = expr

    # Lineage filter using mapped lineage from annotation
    aligned_meta = aligned_meta.copy()
    aligned_meta["lineage"] = aligned_meta[ann_col].map(map_lineage)
    mapping_note = f"Annotation column '{ann_col}' mapped to lineages via keywords; originals retained."
    if lineages:
        select_mask = aligned_meta["lineage"].isin(lineages).to_numpy()
        mapping_note += f" Filter lineages={lineages}."

    expr_sel = {g: np.asarray(expr_aligned[g])[select_mask] for g in markers if g in expr_aligned}
    stats = compute_marker_stats(expr_sel, threshold=threshold, unit="UCSC matrix native")

    # Per-lineage mean table (pooled cells)
    lineage_rows = []
    if select_mask.any():
        lin_labels = aligned_meta.loc[select_mask, "lineage"].to_numpy()
        for lin in sorted(set(lin_labels.tolist())):
            m = lin_labels == lin
            for g, arr in expr_sel.items():
                a = np.asarray(arr)[m]
                if np.all(np.isnan(a)):
                    lineage_rows.append({"lineage": lin, "symbol": g, "n": int(m.sum()), "mean": None, "pct_gt_thr": None, "measured": False})
                else:
                    vv = a[np.isfinite(a)]
                    lineage_rows.append(
                        {
                            "lineage": lin,
                            "symbol": g,
                            "n": int(vv.size),
                            "mean": float(vv.mean()) if vv.size else None,
                            "pct_gt_thr": float(100.0 * (vv > threshold).mean()) if vv.size else None,
                            "measured": True,
                        }
                    )

    anchor_gene = anchor or (markers[0] if markers else None)
    coexp = {}
    if anchor_gene and anchor_gene in expr_sel:
        for other in markers:
            if other == anchor_gene:
                continue
            coexp[other] = compute_coexpression(expr_sel[anchor_gene], expr_sel[other], threshold=threshold)

    # Patient/sample fields if present
    patient_cols = [c for c in aligned_meta.columns if c.lower() in ("patient", "sample", "donor", "patient_id", "sample_id")]
    patient_note = (
        f"Patient/sample fields available: {patient_cols}. Pooled-cell stats weight cells equally; "
        "patients with more cells contribute more to pooled means."
        if patient_cols
        else "No clear patient ID column; do not assume barcodes encode patients."
    )

    return {
        "status": "ok",
        "dataset": ds,
        "annotation_column": ann_col,
        "mapping_note": mapping_note,
        "barcode_align": order_flag,
        "composition": comp,
        "marker_stats": stats,
        "marker_by_lineage": lineage_rows,
        "coexpression_vs_anchor": coexp,
        "anchor": anchor_gene,
        "threshold": threshold,
        "threshold_meaning": "RNA > threshold in native matrix units (not protein-positive)",
        "n_cells_total": int(len(aligned_meta)),
        "n_cells_selected": int(select_mask.sum()),
        "gene_extract": gene_meta,
        "patient_note": patient_note,
        "limits": [
            "Processed public UCSC matrix exploration — not FASTQ reanalysis or new cell-type annotation.",
            "Dropout, dissociation, sorting, and recovery biases apply.",
            "RNA threshold positivity ≠ surface protein positivity.",
            "Author annotations shown with lineage mapping; numeric clusters are never auto-named T/B/DC.",
        ],
        "retrieved_at": _now(),
        "source_files": {k: str(v) for k, v in files.items()},
    }


def analyze_local_h5ad(
    path: str,
    markers: list[str],
    *,
    celltype_col: str = "cell_type",
    lineages: list[str] | None = None,
    threshold: float = 0.0,
) -> dict[str, Any]:
    try:
        import anndata as ad
    except ImportError:
        return {"status": "unsupported", "error": "anndata not installed"}

    try:
        adata = ad.read_h5ad(path, backed="r")
    except Exception as exc:
        return {"status": "fetch_failed", "error": str(exc)}

    n_cells = adata.n_obs
    var_names = {str(v).upper(): str(v) for v in adata.var_names}
    selected = np.ones(n_cells, dtype=bool)
    mapping_note = None
    if celltype_col in adata.obs.columns and lineages:
        ct = adata.obs[celltype_col].astype(str)
        selected = ct.isin(lineages).to_numpy()
        mapping_note = f"Filtered by raw annotation column '{celltype_col}' ∈ {lineages}"

    expr = {}
    for m in markers:
        key = var_names.get(m.upper())
        if not key:
            expr[m] = np.full(int(selected.sum()), np.nan)
        else:
            X = adata[selected, key].X
            if hasattr(X, "toarray"):
                X = X.toarray()
            expr[m] = np.asarray(X).ravel()

    stats = compute_marker_stats(expr, threshold=threshold)
    coexp = {}
    if markers and markers[0] in expr:
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
            "Barcode-based join only via anndata API.",
            "Processed matrix exploration only.",
        ],
    }


# Back-compat name used by older UI
def demo_analysis_bundle(dataset_id: str, markers: list[str]) -> dict[str, Any]:
    return analyze_ucsc_dataset(dataset_id, markers)
