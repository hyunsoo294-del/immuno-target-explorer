# -*- coding: utf-8 -*-
"""Score one TAA on human cancer cell lines.

Each cancer row uses that gene's HPA nTPM only. Unknown and non-cancerous
lines stay out of the score. A missing factor stays missing.
"""

from __future__ import annotations

import math

import pandas as pd

from taa_analysis.cell_lines import UNKNOWN_CANCER, load_cell_line_result

DETECTED_NTPM = 1.0
LEVEL_SATURATION_NTPM = 1024.0
SELECTIVITY_CLIP = 2.0
WEIGHT_LEVEL = 0.50
WEIGHT_PREVALENCE = 0.30
WEIGHT_SELECTIVITY = 0.20
UNSCORED_LABELS = {UNKNOWN_CANCER, "Non-Cancerous"}


def formula_text() -> str:
    return (
        f"점수 = 100 × ({WEIGHT_LEVEL:.2f} × 발현 + {WEIGHT_PREVALENCE:.2f} × 검출비율 + "
        f"{WEIGHT_SELECTIVITY:.2f} × 선택성). "
        f"발현 = min(log2(중앙값 nTPM + 1) / log2({LEVEL_SATURATION_NTPM:.0f} + 1), 1). "
        f"검출비율 = 측정된 세포주 가운데 nTPM ≥ {DETECTED_NTPM:g} 인 비율. "
        f"선택성 = 다른 암종 중앙값 대비 log2 비를 ±{SELECTIVITY_CLIP:g}로 자른 뒤 0–1. "
        "선택성을 계산할 다른 암종이 없으면 발현과 검출비율만 다시 맞춰 합칩니다. "
        "Unknown과 Non-Cancerous는 점수에 넣지 않습니다."
    )


def _finite(value) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if math.isnan(number) or math.isinf(number):
        return None
    return number


def _level(median_ntpm: float | None) -> float | None:
    if median_ntpm is None or median_ntpm < 0:
        return None
    span = math.log2(LEVEL_SATURATION_NTPM + 1.0)
    if span <= 0:
        return None
    return min(math.log2(median_ntpm + 1.0) / span, 1.0)


def _selectivity_unit(log2_ratio: float | None) -> float | None:
    if log2_ratio is None:
        return None
    clipped = min(SELECTIVITY_CLIP, max(-SELECTIVITY_CLIP, log2_ratio))
    return (clipped + SELECTIVITY_CLIP) / (2.0 * SELECTIVITY_CLIP)


def _combine(level: float | None, prevalence: float | None, selectivity: float | None) -> float | None:
    parts = [
        (WEIGHT_LEVEL, level),
        (WEIGHT_PREVALENCE, prevalence),
        (WEIGHT_SELECTIVITY, selectivity),
    ]
    present = [(weight, value) for weight, value in parts if value is not None]
    if not present or level is None or prevalence is None:
        return None
    total = sum(weight for weight, _value in present)
    if total <= 0:
        return None
    return 100.0 * sum(weight * value for weight, value in present) / total


def score_cell_line_groups(lines: pd.DataFrame) -> pd.DataFrame:
    """One row per cancer type, with measured factors and a 0–100 score."""
    if lines is None or lines.empty or "group_label" not in lines.columns:
        return pd.DataFrame()
    work = lines.copy()
    work["nTPM"] = pd.to_numeric(work["nTPM"], errors="coerce")
    work.loc[work["nTPM"] < 0, "nTPM"] = pd.NA
    grouped: list[dict] = []
    for label, group in work.groupby("group_label", dropna=False):
        name = str(label)
        values = pd.to_numeric(group["nTPM"], errors="coerce").dropna()
        grouped.append(
            {
                "cancer": name,
                "n_cell_lines": int(group["cell_line"].nunique()) if "cell_line" in group.columns else int(len(group)),
                "n_measured": int(len(values)),
                "median_nTPM": float(values.median()) if not values.empty else None,
                "q1_nTPM": float(values.quantile(0.25)) if not values.empty else None,
                "q3_nTPM": float(values.quantile(0.75)) if not values.empty else None,
                "detected_n": int((values >= DETECTED_NTPM).sum()) if not values.empty else 0,
                "scored": name not in UNSCORED_LABELS,
            }
        )
    frame = pd.DataFrame(grouped)
    if frame.empty:
        return frame
    scored = frame[frame["scored"] & frame["median_nTPM"].notna()]
    other_pool = {
        row.cancer: row.median_nTPM
        for row in scored.itertuples(index=False)
    }
    rows = []
    for record in frame.itertuples(index=False):
        median = _finite(record.median_nTPM)
        measured = int(record.n_measured)
        prevalence = (int(record.detected_n) / measured) if measured else None
        others = [value for cancer, value in other_pool.items() if cancer != record.cancer]
        other_median = float(pd.Series(others).median()) if others else None
        selectivity_log2 = None
        if record.scored and median is not None and other_median is not None and other_median >= 0:
            selectivity_log2 = math.log2((median + 1.0) / (other_median + 1.0))
        level = _level(median) if record.scored else None
        selectivity_unit = _selectivity_unit(selectivity_log2) if record.scored else None
        score = _combine(level, prevalence if record.scored else None, selectivity_unit)
        rows.append(
            {
                "cancer": record.cancer,
                "n_cell_lines": int(record.n_cell_lines),
                "n_measured": measured,
                "median_nTPM": median,
                "q1_nTPM": _finite(record.q1_nTPM),
                "q3_nTPM": _finite(record.q3_nTPM),
                "detected_n": int(record.detected_n),
                "detected_fraction": prevalence,
                "other_median_nTPM": other_median if record.scored else None,
                "selectivity_log2": selectivity_log2,
                "level_0_1": level,
                "prevalence_0_1": prevalence if record.scored else None,
                "selectivity_0_1": selectivity_unit,
                "score": score,
                "unit": "nTPM",
                "scored": bool(record.scored),
            }
        )
    out = pd.DataFrame(rows)
    out = out.sort_values(
        ["scored", "score", "median_nTPM"],
        ascending=[False, False, False],
        na_position="last",
    ).reset_index(drop=True)
    return out


def score_gene(gene_query: str, force: bool = False) -> dict:
    loaded = load_cell_line_result(gene_query, force=force)
    scored = score_cell_line_groups(loaded["lines"])
    cancer_lines = scored[scored["scored"]] if not scored.empty else scored
    return {
        "identity": loaded["identity"],
        "meta": loaded["meta"],
        "lines": loaded["lines"],
        "scored": scored,
        "n_cancer_groups": int(cancer_lines["cancer"].nunique()) if not cancer_lines.empty else 0,
        "n_scored_lines": int(cancer_lines["n_cell_lines"].sum()) if not cancer_lines.empty else 0,
        "formula": formula_text(),
    }
