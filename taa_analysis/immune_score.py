# -*- coding: utf-8 -*-
"""Immune-activation index for one TAA on human cancer cell lines.

Three measured factors:
- TAA expression: HPA nTPM of the searched gene
- Cell-to-cell contact: HPA nTPM of HLA-A, HLA-B, B2M, ICAM1, CD58 on the same lines
- Binding hindrance: UniProt extracellular length and N-glycans inside that region

A missing factor stays missing. It is not filled with zero.
"""

from __future__ import annotations

import math

import pandas as pd
import requests

from taa_analysis.cell_lines import (
    UNKNOWN_CANCER,
    ensure_hpa_celline_zip,
    extract_hpa_symbols,
    load_cell_line_result,
)
from taa_analysis.settings import load_settings

CONTACT_GENES = ("HLA-A", "HLA-B", "B2M", "ICAM1", "CD58")
UNSCORED_LABELS = {UNKNOWN_CANCER, "Non-Cancerous"}
LEVEL_SATURATION_NTPM = 1024.0
WEIGHT_EXPRESSION = 0.50
WEIGHT_CONTACT = 0.30
WEIGHT_ACCESS = 0.20
ECD_BULK_AA = 700.0
GLYCAN_CAP = 8.0


def formula_text() -> str:
    return (
        f"면역 활성 점수 = 100 × ({WEIGHT_EXPRESSION:.2f} × TAA 발현 + {WEIGHT_CONTACT:.2f} × 세포 간 접촉 + "
        f"{WEIGHT_ACCESS:.2f} × (1 − 결합 방해)). "
        f"발현과 접촉은 log2(중앙값 nTPM + 1) / log2({LEVEL_SATURATION_NTPM:.0f} + 1) 이고 1을 넘기면 1. "
        "접촉은 같은 세포주에서 HLA-A, HLA-B, B2M, ICAM1, CD58 중앙값의 중앙값이다. "
        "결합 방해는 세포외 도메인 길이와 그 안의 N-glycan 수를 0–1로 둔 평균이다. "
        "없는 항목은 0으로 넣지 않고 남은 항목만 다시 맞춘다."
    )


def _median(values) -> float | None:
    series = pd.to_numeric(pd.Series(values), errors="coerce")
    series = series[series >= 0]
    if series.empty:
        return None
    return float(series.median())


def level_unit(median_ntpm: float | None) -> float | None:
    if median_ntpm is None or median_ntpm < 0:
        return None
    span = math.log2(LEVEL_SATURATION_NTPM + 1.0)
    if span <= 0:
        return None
    return min(math.log2(median_ntpm + 1.0) / span, 1.0)


def _combine(expression: float | None, contact: float | None, hindrance: float | None) -> float | None:
    access = None if hindrance is None else 1.0 - hindrance
    parts = [
        (WEIGHT_EXPRESSION, expression),
        (WEIGHT_CONTACT, contact),
        (WEIGHT_ACCESS, access),
    ]
    present = [(weight, value) for weight, value in parts if value is not None]
    if not present:
        return None
    total = sum(weight for weight, _value in present)
    return 100.0 * sum(weight * value for weight, value in present) / total


def binding_hindrance(features: list[dict], location_class: str = "") -> dict:
    """Extracellular bulk and N-glycans. Intracellular proteins are not given a fake surface score."""
    spans: list[tuple[int, int]] = []
    glycans: list[int] = []
    for feature in features or []:
        kind = str(feature.get("type") or "")
        text = str(feature.get("description") or "")
        start = feature.get("start")
        end = feature.get("end")
        try:
            start_i = int(start) if start is not None else None
            end_i = int(end) if end is not None else None
        except (TypeError, ValueError):
            start_i = end_i = None
        low = f"{kind} {text}".lower()
        if "extracellular" in low and start_i and end_i and end_i >= start_i:
            spans.append((start_i, end_i))
        if kind.lower() == "glycosylation" and "n-linked" in low and start_i:
            glycans.append(start_i)
    ecd_aa = sum(end - start + 1 for start, end in spans) if spans else None
    if spans:
        extracellular_glycans = sum(1 for site in glycans if any(start <= site <= end for start, end in spans))
    else:
        extracellular_glycans = None
    surface = location_class in {"Cell surface", "Membrane-associated"}
    parts = []
    if ecd_aa is not None:
        parts.append(min(ecd_aa / ECD_BULK_AA, 1.0))
    if extracellular_glycans is not None:
        parts.append(min(extracellular_glycans / GLYCAN_CAP, 1.0))
    hindrance = float(sum(parts) / len(parts)) if parts and surface else None
    reason = ""
    if not surface:
        reason = "plasma-membrane annotation is missing, so binding hindrance is not scored"
    elif hindrance is None:
        reason = "extracellular length and N-glycan sites are missing"
    return {
        "ecd_aa": ecd_aa,
        "n_glycan_extracellular": extracellular_glycans,
        "location_class": location_class,
        "hindrance_0_1": hindrance,
        "reason": reason,
    }


def _cancer_lines(lines: pd.DataFrame) -> pd.DataFrame:
    work = lines.copy()
    work["group_label"] = work["group_label"].astype(str)
    return work[~work["group_label"].isin(UNSCORED_LABELS)].copy()


def score_lines(lines: pd.DataFrame, contact: pd.DataFrame, hindrance: float | None) -> dict:
    """Score the supplied cancer-cell-line rows. Caller already applied the user's selection."""
    if lines is None or lines.empty:
        return {"score": None, "per_cancer": pd.DataFrame(), "contact_genes": pd.DataFrame(), "factors": {}}
    selected = set(lines["cell_line"].astype(str))
    expression = _median(lines["nTPM"])
    contact_rows = []
    if contact is not None and not contact.empty:
        subset = contact[contact["cell_line"].astype(str).isin(selected)]
        for gene, group in subset.groupby("gene_symbol"):
            contact_rows.append({"gene_symbol": gene, "median_nTPM": _median(group["nTPM"]), "n_cell_lines": int(group["cell_line"].nunique())})
    contact_table = pd.DataFrame(contact_rows)
    contact_median = _median(contact_table["median_nTPM"]) if not contact_table.empty else None
    expression_unit = level_unit(expression)
    contact_unit = level_unit(contact_median)
    per_cancer = []
    for label, group in lines.groupby("group_label", dropna=False):
        cancer_expression = _median(group["nTPM"])
        names = set(group["cell_line"].astype(str))
        cancer_contact = None
        if not contact_table.empty and contact is not None:
            gene_medians = []
            scoped = contact[contact["cell_line"].astype(str).isin(names)]
            for _gene, gene_group in scoped.groupby("gene_symbol"):
                value = _median(gene_group["nTPM"])
                if value is not None:
                    gene_medians.append(value)
            cancer_contact = _median(gene_medians)
        per_cancer.append(
            {
                "cancer": str(label),
                "n_cell_lines": int(group["cell_line"].nunique()),
                "taa_median_nTPM": cancer_expression,
                "contact_median_nTPM": cancer_contact,
                "hindrance_0_1": hindrance,
                "score": _combine(level_unit(cancer_expression), level_unit(cancer_contact), hindrance),
            }
        )
    per_frame = pd.DataFrame(per_cancer)
    if not per_frame.empty:
        per_frame = per_frame.sort_values("score", ascending=False, na_position="last").reset_index(drop=True)
    return {
        "score": _combine(expression_unit, contact_unit, hindrance),
        "factors": {
            "taa_median_nTPM": expression,
            "taa_expression_0_1": expression_unit,
            "contact_median_nTPM": contact_median,
            "contact_0_1": contact_unit,
            "hindrance_0_1": hindrance,
            "n_cell_lines": int(lines["cell_line"].nunique()),
        },
        "per_cancer": per_frame,
        "contact_genes": contact_table,
    }


def select_lines(lines: pd.DataFrame, cancers: list[str] | None = None, cell_lines: list[str] | None = None) -> pd.DataFrame:
    """Cancer-wide selection, or an explicit cell-line list. Unknown and non-cancerous stay out."""
    work = _cancer_lines(lines)
    if cell_lines:
        chosen = {str(name) for name in cell_lines}
        return work[work["cell_line"].astype(str).isin(chosen)].copy()
    if cancers:
        chosen_cancers = {str(name) for name in cancers}
        return work[work["group_label"].astype(str).isin(chosen_cancers)].copy()
    return work


def load_contact_panel(force: bool = False) -> pd.DataFrame:
    from pathlib import Path

    cache = Path(load_settings()["TAA_CACHE_DIR"])
    cache.mkdir(parents=True, exist_ok=True)
    dest = cache / "hpa_contact_panel.csv"
    if dest.exists() and not force and dest.stat().st_size > 100:
        return pd.read_csv(dest)
    found = extract_hpa_symbols(ensure_hpa_celline_zip(cache), list(CONTACT_GENES), cache / "hpa_celline_genes.csv")
    frames = [frame for frame in found.values() if frame is not None and not frame.empty]
    panel = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=["gene_symbol", "cell_line", "nTPM"])
    panel.to_csv(dest, index=False)
    return panel


def fetch_binding_hindrance(symbol: str) -> dict:
    from taa_character import fetch_uniprot

    session = requests.Session()
    try:
        record = fetch_uniprot(session, symbol, "")
    finally:
        session.close()
    if not record:
        return binding_hindrance([], "")
    features = record.get("features") or []
    # The character fetch omits glycosylation. Read it from the same accession when present.
    accession = record.get("accession") or ""
    if accession:
        response = requests.get(
            f"https://rest.uniprot.org/uniprotkb/{accession}.json",
            timeout=40,
        )
        if response.ok:
            features = (response.json() or {}).get("features") or features
    parsed = []
    for feature in features:
        location = feature.get("location") or {}
        start = (location.get("start") or {}).get("value") if isinstance(location, dict) else feature.get("start")
        end = (location.get("end") or {}).get("value") if isinstance(location, dict) else feature.get("end")
        parsed.append(
            {
                "type": feature.get("type"),
                "description": feature.get("description"),
                "start": start if start is not None else feature.get("start"),
                "end": end if end is not None else feature.get("end"),
            }
        )
    result = binding_hindrance(parsed, str(record.get("location_class") or ""))
    result["accession"] = accession
    result["uniprot_url"] = record.get("uniprot_url") or ""
    return result


def prepare_gene(gene_query: str, force: bool = False) -> dict:
    loaded = load_cell_line_result(gene_query, force=force)
    contact = load_contact_panel(force=force)
    symbol = loaded["identity"].get("gene_symbol") or gene_query
    try:
        hindrance = fetch_binding_hindrance(symbol)
    except Exception as exc:
        hindrance = {"hindrance_0_1": None, "reason": str(exc), "ecd_aa": None, "n_glycan_extracellular": None}
    return {
        "identity": loaded["identity"],
        "meta": loaded["meta"],
        "lines": _cancer_lines(loaded["lines"]),
        "contact": contact,
        "hindrance": hindrance,
        "formula": formula_text(),
    }
