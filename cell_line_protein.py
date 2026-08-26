# -*- coding: utf-8 -*-
"""Cancer cell-line protein abundance from ProCan-DepMapSanger (Cell Model Passports)."""

from __future__ import annotations

import json
from typing import Any

import pandas as pd
import requests

CMP_API = "https://api.cellmodelpassports.sanger.ac.uk"
CMP_PASSPORT = "https://cellmodelpassports.sanger.ac.uk/passports/{model_id}"
HPA_CELLINE = "https://www.proteinatlas.org/{ensembl}-{symbol}/celline"
DEPMAP_GENE = "https://depmap.org/portal/gene/{symbol}?tab=characterization"
PROCAN_PAPER = "https://doi.org/10.1016/j.ccell.2022.06.010"


def _lineage_from_identifiers(identifiers: list[str]) -> str:
    for ident in identifiers:
        text = str(ident or "")
        if "_" in text and not text.startswith(("ACH-", "CVCL_", "SID")):
            lineage = text.split("_")[-1].replace("-", " ").strip()
            if lineage and not lineage.isdigit() and lineage.upper() != text.upper():
                return lineage.title()
    return ""


def _model_name(model_attrs: dict[str, Any], identifiers: list[str], model_id: str) -> str:
    names = model_attrs.get("names") or []
    if names:
        return str(names[0])
    for ident in identifiers:
        text = str(ident or "")
        if text and not text.startswith(("ACH-", "CVCL_", "SID")) and not text.isdigit():
            return text.split("_")[0]
    return model_id


def _resolve_uniprot(session: requests.Session, symbol: str, ensembl_id: str, accession: str = "") -> str:
    if accession:
        return accession
    from taa_character import fetch_uniprot

    rec = fetch_uniprot(session, symbol, ensembl_id)
    return str(rec.get("accession") or "")


def _fetch_proteomics_page(
    session: requests.Session,
    uniprot_id: str,
    page: int,
    page_size: int = 500,
) -> dict[str, Any]:
    from data_sources import _request

    filt = json.dumps([{"name": "uniprot_id", "op": "eq", "val": uniprot_id}])
    response = _request(
        session,
        "GET",
        f"{CMP_API}/datasets/proteomics",
        params={
            "filter": filt,
            "page[size]": page_size,
            "page[number]": page,
            "include": "model",
        },
        headers={"Accept": "application/vnd.api+json, application/json"},
    )
    return response.json() or {}


def collect_cell_line_protein(
    session: requests.Session,
    symbol: str,
    ensembl_id: str = "",
    uniprot_accession: str = "",
) -> tuple[pd.DataFrame, dict[str, Any]]:
    accession = _resolve_uniprot(session, symbol, ensembl_id, uniprot_accession)
    meta: dict[str, Any] = {
        "symbol": symbol,
        "ensembl_id": ensembl_id,
        "uniprot_id": accession,
        "source": "ProCan-DepMapSanger",
        "hpa_celline_url": HPA_CELLINE.format(ensembl=ensembl_id, symbol=symbol) if ensembl_id else "",
        "depmap_url": DEPMAP_GENE.format(symbol=symbol),
        "procan_url": PROCAN_PAPER,
        "passports_url": "https://cellmodelpassports.sanger.ac.uk/",
        "error": None,
        "n_measured": 0,
        "n_positive": 0,
    }
    if not accession:
        meta["error"] = "No UniProt accession, so cell-line proteomics could not be queried."
        return pd.DataFrame(), meta

    rows: list[dict[str, Any]] = []
    page = 1
    try:
        while page <= 8:
            payload = _fetch_proteomics_page(session, accession, page)
            data = payload.get("data") or []
            included = payload.get("included") or []
            models = {item.get("id"): item for item in included if item.get("type") == "model"}
            idents_by_id = {item.get("id"): item for item in included if item.get("type") == "model_identifier"}
            for rec in data:
                attrs = rec.get("attributes") or {}
                zscore = attrs.get("zscore")
                intensity = attrs.get("protein_intensity")
                if zscore is None or intensity is None:
                    continue
                try:
                    z_val = float(zscore)
                    int_val = float(intensity)
                except (TypeError, ValueError):
                    continue
                model_id = str(attrs.get("model_id") or "")
                model = models.get(model_id) or {}
                model_attrs = model.get("attributes") or {}
                ident_ids = [
                    rel.get("id")
                    for rel in (((model.get("relationships") or {}).get("identifiers") or {}).get("data") or [])
                    if rel.get("id")
                ]
                identifiers = [
                    str((idents_by_id.get(i) or {}).get("attributes", {}).get("identifier") or "")
                    for i in ident_ids
                ]
                rows.append(
                    {
                        "cell_line": _model_name(model_attrs, identifiers, model_id),
                        "lineage": _lineage_from_identifiers(identifiers),
                        "model_id": model_id,
                        "uniprot_id": accession,
                        "protein_intensity_log2": int_val,
                        "zscore": z_val,
                        "peptides": attrs.get("median_no_peptides"),
                        "passport_url": CMP_PASSPORT.format(model_id=model_id) if model_id else "",
                    }
                )
            if not (payload.get("links") or {}).get("next") or not data:
                break
            page += 1
    except Exception as exc:
        meta["error"] = str(exc)
        return pd.DataFrame(), meta

    if not rows:
        meta["error"] = f"No ProCan cell-line protein values for {symbol} ({accession})."
        return pd.DataFrame(), meta

    work = pd.DataFrame(rows)
    meta["n_measured"] = int(work["cell_line"].nunique())
    positive = work.loc[work["zscore"] > 0].copy()
    if positive.empty:
        meta["n_positive"] = 0
        return pd.DataFrame(), meta

    positive = (
        positive.sort_values("zscore", ascending=False)
        .drop_duplicates(subset=["cell_line"], keep="first")
        .reset_index(drop=True)
    )
    z = pd.to_numeric(positive["zscore"], errors="coerce")
    z_min = float(z.min())
    z_max = float(z.max())
    span = z_max - z_min
    if span <= 0:
        positive["protein_norm_0_1"] = 1.0
    else:
        positive["protein_norm_0_1"] = ((z - z_min) / span).round(3)
    n = len(positive)
    positive["rank"] = range(1, n + 1)
    positive["percentile"] = (((n - positive["rank"] + 1) / n) * 100).round(1)
    meta["n_positive"] = n
    return positive, meta
