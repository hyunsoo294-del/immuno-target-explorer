"""Confidence labels. RNA-only antigen density cannot be HIGH or MEDIUM."""

from __future__ import annotations

import pandas as pd


def cap_confidence(label: str, accessibility_source: str, params: dict) -> str:
    """Unlisted morphology cannot support HIGH. Curated lines keep their label."""
    if accessibility_source == "curated":
        return label
    cap = str((params.get("confidence") or {}).get("non_curated_accessibility_cap") or "MEDIUM")
    rank = {"LOW": 0, "MEDIUM": 1, "HIGH": 2}
    if rank.get(label, 0) > rank.get(cap, 1):
        return cap
    return label


def row_confidence(completeness: float, f1_tier: str, params: dict) -> str:
    spec = params["confidence"]
    if f1_tier == "rna_only":
        return str(spec["rna_only_f1_caps_at"])
    if pd.isna(completeness):
        return "LOW"
    if completeness >= float(spec["high_min_completeness"]):
        return "HIGH"
    if completeness >= float(spec["medium_min_completeness"]):
        return "MEDIUM"
    return "LOW"


def completeness(weights: dict, tiers: dict, scores: dict, params: dict) -> float:
    """Fraction of intended weight that is backed by data. NaN factors contribute 0."""
    credits = params["tier_completeness"]
    intended = float(sum(weights.values()))
    if intended <= 0:
        return 0.0
    total = 0.0
    for factor, weight in weights.items():
        if factor not in scores or pd.isna(scores[factor]):
            continue
        tier = tiers.get(factor) or "default"
        total += float(weight) * float(credits.get(tier, 0.0))
    return total / intended
