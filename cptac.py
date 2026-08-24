# -*- coding: utf-8 -*-
"""CPTAC mass-spectrometry protein abundance via cBioPortal (NCI/PDC-linked)."""

from __future__ import annotations

import math
from typing import Any

import pandas as pd
import requests

from expression_compare import fetch_entrez_id
from data_qc import normalize_cptac_sample_table

TIMEOUT = 90
PSEUDOCOUNT = 1.0

# Prospective CPTAC studies on cBioPortal with gene-level protein quantification.
# TCIA/IDC host imaging; they do not provide per-gene IHC scores, so quantified
# protein here comes from CPTAC mass spectrometry.
CPTAC_STUDIES = {
    "brca_cptac_2020": {
        "cancer_type": "breast cancer",
        "profile": "brca_cptac_2020_protein_quantification",
        "citation": "CPTAC Cell 2020",
    },
    "luad_cptac_2020": {
        "cancer_type": "lung adenocarcinoma",
        "profile": "luad_cptac_2020_protein_quantification",
        "citation": "CPTAC Cell 2020",
    },
    "lusc_cptac_2021": {
        "cancer_type": "lung squamous cell carcinoma",
        "profile": "lusc_cptac_2021_protein_quantification",
        "citation": "CPTAC Cell 2021",
    },
    "coad_cptac_2019": {
        "cancer_type": "colon cancer",
        "profile": "coad_cptac_2019_protein_quantification",
        "citation": "CPTAC Cell 2019",
    },
    "ucec_cptac_2020": {
        "cancer_type": "endometrial carcinoma",
        "profile": "ucec_cptac_2020_protein_quantification",
        "citation": "CPTAC Cell 2020",
    },
    "paad_cptac_2021": {
        "cancer_type": "pancreatic ductal adenocarcinoma",
        "profile": "paad_cptac_2021_protein_quantification",
        "citation": "CPTAC Cell 2021",
    },
    "gbm_cptac_2021": {
        "cancer_type": "glioblastoma",
        "profile": "gbm_cptac_2021_protein_quantification",
        "citation": "CPTAC Cell 2021",
    },
    "brain_cptac_2020": {
        "cancer_type": "pediatric brain cancer",
        "profile": "brain_cptac_2020_protein_quantification",
        "citation": "CPTAC/CHOP Cell 2020",
    },
}


def fetch_cptac_protein(session: requests.Session, entrez_id: int) -> pd.DataFrame:
    profiles = [meta["profile"] for meta in CPTAC_STUDIES.values()]
    try:
        response = session.post(
            "https://www.cbioportal.org/api/molecular-data/fetch?projection=SUMMARY",
            json={"entrezGeneIds": [entrez_id], "molecularProfileIds": profiles},
            timeout=TIMEOUT,
            headers={"Content-Type": "application/json"},
        )
        response.raise_for_status()
        rows = response.json() or []
    except Exception:
        return pd.DataFrame()

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    df = df.dropna(subset=["value"])
    study_to_meta = CPTAC_STUDIES
    df["cancer_type"] = df["studyId"].map(lambda sid: study_to_meta.get(sid, {}).get("cancer_type"))
    df["source"] = df["studyId"].map(lambda sid: study_to_meta.get(sid, {}).get("citation", "CPTAC"))
    df["sample_class"] = "tumor"
    df = df.dropna(subset=["cancer_type"]).reset_index(drop=True)
    return normalize_cptac_sample_table(df)


def build_cptac_summary(cptac: pd.DataFrame) -> pd.DataFrame:
    if cptac is None or cptac.empty:
        return pd.DataFrame()

    rows: list[dict[str, Any]] = []
    for cancer_type, subset in cptac.groupby("cancer_type"):
        values = pd.to_numeric(subset["value"], errors="coerce").dropna()
        if values.empty:
            continue
        median = float(values.median())
        mean = float(values.mean())
        q1 = float(values.quantile(0.25))
        q3 = float(values.quantile(0.75))
        scale = subset["cptac_scale"].iloc[0] if "cptac_scale" in subset.columns else "cdap_log_ratio"
        rows.append(
            {
                "cancer_type": cancer_type,
                "source": subset["source"].iloc[0] if "source" in subset.columns else "CPTAC",
                "n_tumor": int(values.count()),
                "cptac_scale": scale,
                "median_log2_ratio": round(median, 3),
                "mean_log2_ratio": round(mean, 3),
                "q1_log2_ratio": round(q1, 3),
                "q3_log2_ratio": round(q3, 3),
                "fold_vs_reference": round(float(2 ** median), 3) if abs(median) <= 8 else None,
                "pct_above_reference": round(float((values > 0).mean()) * 100.0, 1),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values("median_log2_ratio", ascending=False, na_position="last").reset_index(drop=True)


def collect_cptac_for_gene(
    session: requests.Session,
    symbol: str,
    ensembl_id: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, str | None]:
    entrez = fetch_entrez_id(session, symbol, ensembl_id)
    if not entrez:
        return pd.DataFrame(), pd.DataFrame(), f"Could not resolve gene ID for {symbol}."
    samples = fetch_cptac_protein(session, entrez)
    if samples.empty:
        return samples, pd.DataFrame(), "No CPTAC protein quantification found for this gene."
    return samples, build_cptac_summary(samples), None
