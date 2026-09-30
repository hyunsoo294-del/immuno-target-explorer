"""Excel export. Three sheets plus the request JSON on Provenance."""

from __future__ import annotations

import json
from io import BytesIO

import pandas as pd

from immunoscore.src.config_loader import load_configs

DISPLAY_COLUMNS = [
    ("cell_line", "세포주"),
    ("cancer", "암종"),
    ("F1_expression", "발현(게이트)"),
    ("F2_heterogeneity", "균일성"),
    ("F3_internalization", "내재화"),
    ("F4_checkpoint", "체크포인트"),
    ("F5_adhesion", "접촉"),
    ("accessibility_label", "접근성"),
    ("F6_glycocalyx", "당질층"),
    ("F7_epitope_proximity", "에피톱"),
    ("score", "점수"),
    ("emax_text", "예측 Emax"),
    ("confidence", "신뢰도"),
    ("flag_text", "플래그"),
]


def _access_label(row: pd.Series) -> str:
    labels = load_configs()["params"].get("accessibility_source_labels") or {}
    value = row.get("accessibility")
    source = labels.get(row.get("f5_accessibility_source"), row.get("f5_accessibility_source") or "")
    if pd.isna(value):
        return ""
    return f"{float(value):.2f} {source}".strip()


def display_frame(result: pd.DataFrame) -> pd.DataFrame:
    labels = load_configs()["params"]["confidence_labels"]
    frame = pd.DataFrame()
    for source, label in DISPLAY_COLUMNS:
        if source == "confidence":
            frame[label] = result[source].map(lambda value: labels.get(value, value))
        elif source in ("cell_line", "cancer", "flag_text", "emax_text"):
            frame[label] = result[source] if source in result.columns else ""
        elif source == "accessibility_label":
            frame[label] = result.apply(_access_label, axis=1)
        elif source == "F4_checkpoint":
            shown = []
            for _idx, row in result.iterrows():
                if row.get("F4_checkpoint_tier") == "not_applicable" or pd.isna(row.get(source)):
                    shown.append("—")
                else:
                    shown.append(round(float(row[source]), 1))
            frame[label] = shown
        else:
            frame[label] = pd.to_numeric(result[source], errors="coerce").round(1)
    return frame


def to_excel(result: pd.DataFrame, request: dict) -> bytes:
    ranked = display_frame(result)
    factor_cols = ["ModelID", "cell_line", "cancer", "score", "data_completeness", "confidence", "f2_source"]
    for factor, _label in DISPLAY_COLUMNS:
        if factor.startswith("F"):
            factor_cols.extend([factor, factor + "_tier", factor + "_detail"])
    breakdown = result[[column for column in factor_cols if column in result.columns]].copy()
    provenance_rows = [
        {"key": "request", "value": json.dumps(request, ensure_ascii=False)},
        {"key": "manifest", "value": json.dumps(request.get("manifest") or {}, ensure_ascii=False)},
        {
            "key": "effective_weights",
            "value": json.dumps(request.get("effective_weights") or {}, ensure_ascii=False),
        },
        {
            "key": "nominal_weights",
            "value": json.dumps(request.get("nominal_weights") or {}, ensure_ascii=False),
        },
        {
            "key": "calibration",
            "value": json.dumps(request.get("calibration") or {}, ensure_ascii=False),
        },
        {
            "key": "disclaimer",
            "value": "Preliminary panel-selection heuristic. Not a validated predictor of assay Emax. "
            "The Emax column is a monotone calibration of the score and does not re-rank lines.",
        },
    ]
    provenance = pd.DataFrame(provenance_rows)
    buffer = BytesIO()
    with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
        ranked.to_excel(writer, sheet_name="Ranked", index=False)
        breakdown.to_excel(writer, sheet_name="Factor_Breakdown", index=False)
        provenance.to_excel(writer, sheet_name="Provenance", index=False)
    return buffer.getvalue()
