# -*- coding: utf-8 -*-
"""Immune fractions. Medians are not renormalized. Scores are not percentages."""

from __future__ import annotations

import pandas as pd

RELATIVE_ONLY_ALGORITHMS = {"CIBERSORT_LM22_relative", "CIBERSORT_absolute_score"}
SCORE_NOT_PERCENT = {"xCell", "ESTIMATE"}
ALL_CELL_ALGORITHMS = {"quanTIseq", "EPIC"}

BROAD_IMMUNE_GROUPS = (
    "T",
    "B",
    "Plasma",
    "NK",
    "Monocyte/Macrophage",
    "DC",
    "Granulocyte",
    "Mast",
    "Other/Unassigned immune",
)

SUM_TOLERANCE = 1e-6


class CompositionError(ValueError):
    pass


def complete_rows(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    """Drop incomplete composition vectors. Missing components are not replaced with zero."""
    work = frame.copy()
    for column in columns:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    return work.dropna(subset=columns)


def vector_sums(frame: pd.DataFrame, columns: list[str]) -> pd.Series:
    return frame[columns].sum(axis=1)


def assert_relative_vectors_sum_to_one(frame: pd.DataFrame, columns: list[str], tolerance: float = SUM_TOLERANCE) -> None:
    if frame.empty:
        return
    sums = vector_sums(frame, columns)
    if ((sums - 1.0).abs() > tolerance).any():
        raise CompositionError("complete relative vector does not sum to 1")


def mean_complete_composition(frame: pd.DataFrame, columns: list[str]) -> pd.Series:
    """Mean of complete per-patient vectors. The mean of relative vectors also sums to 1."""
    complete = complete_rows(frame, columns)
    assert_relative_vectors_sum_to_one(complete, columns)
    if complete.empty:
        return pd.Series({column: float("nan") for column in columns})
    means = complete[columns].mean(axis=0)
    total = float(means.sum())
    if abs(total - 1.0) > SUM_TOLERANCE:
        raise CompositionError(f"mean composition sums to {total}, not 1")
    return means


def component_medians(frame: pd.DataFrame, columns: list[str]) -> pd.Series:
    """Per-component medians. These are not scaled to 100% and have no tumor residual."""
    complete = complete_rows(frame, columns)
    if complete.empty:
        return pd.Series({column: float("nan") for column in columns})
    return complete[columns].median(axis=0)


def all_cell_percentage_status(algorithm: str) -> dict[str, str]:
    if algorithm in RELATIVE_ONLY_ALGORITHMS or algorithm in SCORE_NOT_PERCENT:
        return {
            "status": "Unavailable",
            "missing_reason": "unsupported",
            "detail": (
                f"{algorithm} cannot supply an all-cell immune percentage. "
                "Relative fractions, absolute scores, xCell, and ESTIMATE are not relabeled as whole-tumor percents."
            ),
        }
    if algorithm in ALL_CELL_ALGORITHMS:
        return {
            "status": "Estimated",
            "missing_reason": "",
            "detail": "Estimated fraction of modeled immune populations, not a direct cell count.",
        }
    return {
        "status": "Unavailable",
        "missing_reason": "unsupported",
        "detail": f"No all-cell percentage is defined for {algorithm}.",
    }


def same_sample_all_cell_fraction(immune_fraction: float, within_immune_fraction: float, *, same_sample: bool, same_definition: bool) -> float:
    """A = I * R only for the same sample, denominator, and cell definition."""
    if not same_sample or not same_definition:
        raise CompositionError(
            "refusing to multiply fractions from different samples, studies, or cell definitions"
        )
    return float(immune_fraction) * float(within_immune_fraction)


def whole_tumor_fraction_allowed(enrichment_protocol: str) -> bool:
    protocol = (enrichment_protocol or "").strip().lower()
    blocked = ("cd45", "t-cell sorted", "t cell sorted", "immune-only", "immune only")
    return not any(token in protocol for token in blocked)
