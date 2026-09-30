"""Excel export. Three sheets plus the request JSON on Provenance."""

from __future__ import annotations

import json
from io import BytesIO

import pandas as pd

from immunoscore.src.config_loader import load_configs

DISPLAY_COLUMNS = [
    ("cell_line", "세포주"),
    ("cancer", "암종"),
    ("F1_expression", "발현"),
    ("F2_heterogeneity", "균일성"),
    ("F3_internalization", "내재화"),
    ("F4_checkpoint", "체크포인트"),
    ("F5_adhesion", "접촉"),
    ("F6_glycocalyx", "당질층"),
    ("F7_epitope_proximity", "에피톱"),
    ("score", "점수"),
    ("confidence", "신뢰도"),
    ("flag_text", "플래그"),
]


def display_frame(result: pd.DataFrame) -> pd.DataFrame:
    labels = load_configs()["params"]["confidence_labels"]
    frame = pd.DataFrame()
    for source, label in DISPLAY_COLUMNS:
        if source == "confidence":
            frame[label] = result[source].map(lambda value: labels.get(value, value))
        elif source in ("cell_line", "cancer", "flag_text"):
            frame[label] = result[source]
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
    provenance = pd.DataFrame(
        [
            {"key": "request", "value": json.dumps(request, ensure_ascii=False)},
            {"key": "manifest", "value": json.dumps(request.get("manifest") or {}, ensure_ascii=False)},
            {
                "key": "disclaimer",
                "value": "Preliminary panel-selection heuristic. Not a validated predictor of assay Emax.",
            },
        ]
    )
    buffer = BytesIO()
    with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
        ranked.to_excel(writer, sheet_name="Ranked", index=False)
        breakdown.to_excel(writer, sheet_name="Factor_Breakdown", index=False)
        provenance.to_excel(writer, sheet_name="Provenance", index=False)
    return buffer.getvalue()
