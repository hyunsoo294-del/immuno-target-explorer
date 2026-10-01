# -*- coding: utf-8 -*-
"""Molecule panel, detection rules, and healthy-blood reference labeling."""

from __future__ import annotations

import io
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests

MOLECULE_PANEL = (
    {"molecule": "4-1BB / CD137", "gene_symbol": "TNFRSF9"},
    {"molecule": "PD-1", "gene_symbol": "PDCD1"},
    {"molecule": "PD-L1", "gene_symbol": "CD274"},
    {"molecule": "CTLA-4", "gene_symbol": "CTLA4"},
    {"molecule": "TIGIT", "gene_symbol": "TIGIT"},
    {"molecule": "LAG-3", "gene_symbol": "LAG3"},
    {"molecule": "TIM-3", "gene_symbol": "HAVCR2"},
    {"molecule": "OX40", "gene_symbol": "TNFRSF4"},
    {"molecule": "ICOS", "gene_symbol": "ICOS"},
    {"molecule": "CD27", "gene_symbol": "CD27"},
    {"molecule": "CD40", "gene_symbol": "CD40"},
    {"molecule": "CD80", "gene_symbol": "CD80"},
    {"molecule": "CD86", "gene_symbol": "CD86"},
    {"molecule": "HLA-DR alpha transcript", "gene_symbol": "HLA-DRA"},
    {"molecule": "MHC-I-related", "gene_symbol": "HLA-A"},
    {"molecule": "MHC-I-related", "gene_symbol": "HLA-B"},
    {"molecule": "MHC-I-related", "gene_symbol": "HLA-C"},
    {"molecule": "MHC-I-related", "gene_symbol": "B2M"},
)

ENRICHMENT_SUMMARY_FIELDS = {
    "rna blood cell specific ntp m",
    "rna blood cell specificity",
    "rna blood cell specificity score",
    "rna tissue specific ntp m",
    "rna cancer specific ptp m",
    "rna single cell type specific ncpm",
}

INELIGIBLE_TME_TOKENS = (
    "heart",
    "cerebellum",
    "pbmc",
    "peripheral blood mononuclear",
    "normal blood",
    "healthy blood",
)


class DetectionError(ValueError):
    pass


def is_enrichment_summary_field(field_name: str) -> bool:
    compact = " ".join(str(field_name or "").lower().replace("_", " ").split())
    return any(token in compact for token in (
        "specificity score",
        "specific ntp",
        "specific ptpm",
        "specific ncpm",
        "enriched",
    )) or compact in ENRICHMENT_SUMMARY_FIELDS


def tme_catalog_eligible(tissue: str, specimen: str, species: str = "Homo sapiens") -> bool:
    if species.lower() not in {"homo sapiens", "human"}:
        return False
    blob = f"{tissue} {specimen}".lower()
    if any(token in blob for token in INELIGIBLE_TME_TOKENS):
        return False
    if "tumor" not in specimen.lower() and "tumour" not in specimen.lower():
        return False
    return True


def rna_detection_rate(values: pd.Series, layer: str) -> float:
    """Fraction of cells with raw count > 0. Scaled or imputed layers are refused."""
    name = layer.lower()
    if any(token in name for token in ("z-score", "zscore", "scaled", "imputed")):
        raise DetectionError(f"detection rate cannot be computed from layer {layer}")
    if "count" not in name:
        raise DetectionError(f"detection rate needs a raw count layer, got {layer}")
    numeric = pd.to_numeric(values, errors="coerce").dropna()
    if numeric.empty:
        raise DetectionError("no count values")
    if (numeric < 0).any():
        raise DetectionError("raw counts contain negatives")
    return float((numeric > 0).mean())


def summarize_cells_by_patient(frame: pd.DataFrame, value_column: str) -> dict:
    """Patient-level mean first. A donor with no selected cells is unavailable, not zero."""
    if frame.empty or "patient_id" not in frame.columns:
        return {
            "patient_mean": None,
            "pooled_mean": None,
            "n_patients_with_cells": 0,
            "status": "Unavailable",
            "missing_reason": "unavailable",
        }
    rows = []
    for patient_id, group in frame.groupby("patient_id"):
        if group.empty:
            rows.append({"patient_id": patient_id, "mean": None, "status": "Unavailable"})
            continue
        values = pd.to_numeric(group[value_column], errors="coerce").dropna()
        if values.empty:
            rows.append({"patient_id": patient_id, "mean": None, "status": "Unavailable"})
        else:
            rows.append({"patient_id": patient_id, "mean": float(values.mean()), "status": "Observed"})
    means = [row["mean"] for row in rows if row["mean"] is not None]
    pooled = pd.to_numeric(frame[value_column], errors="coerce").dropna()
    return {
        "patient_mean": float(sum(means) / len(means)) if means else None,
        "pooled_mean": float(pooled.mean()) if not pooled.empty else None,
        "n_patients_with_cells": len(means),
        "n_patients_without_cells": sum(row["status"] == "Unavailable" for row in rows),
        "status": "Observed" if means else "Unavailable",
        "missing_reason": "" if means else "unavailable",
        "per_patient": rows,
    }


def load_hpa_immune_cell_table(url: str, cache_dir: Path, genes: list[str] | None = None) -> pd.DataFrame:
    """Cell-type summary nTPM from healthy blood. This is reference only, never tumor TIL."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    dest = cache_dir / "rna_immune_cell.tsv"
    if not dest.exists() or dest.stat().st_size < 40:
        response = requests.get(url, timeout=(12, 180))
        response.raise_for_status()
        with zipfile.ZipFile(io.BytesIO(response.content)) as archive:
            name = archive.namelist()[0]
            dest.write_bytes(archive.read(name))
    table = pd.read_csv(dest, sep="\t")
    table["Gene name"] = table["Gene name"].astype(str)
    if genes:
        table = table[table["Gene name"].str.upper().isin({gene.upper() for gene in genes})].copy()
    retrieved = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    table["evidence_status"] = "Reference only"
    table["specimen_type"] = "immune cells isolated from healthy blood"
    table["tissue_origin"] = "healthy blood"
    table["cancer_type"] = pd.NA
    table["subtype_label"] = pd.NA
    table["unit"] = "nTPM"
    table["source_name"] = "Human Protein Atlas immune cell"
    table["source_url"] = url
    table["dataset_id"] = "rna_immune_cell.tsv"
    table["retrieved_at"] = retrieved
    table["limitation"] = (
        "One summary value per immune-cell category. Not a patient-level distribution, "
        "not a tumor TIL measurement, and not an all-cell tumor fraction."
    )
    table["protein_evidence"] = "Protein unavailable"
    return table.reset_index(drop=True)
