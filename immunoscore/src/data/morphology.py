"""Hand-entered culture accessibility. Not inferred from RNA or GrowthPattern alone."""

from __future__ import annotations


def lookup_accessibility(cell_line: str, stripped: str, growth: str, morphology: dict) -> tuple[float, str]:
    """Return (accessibility 0-1, source).

    source is curated, growthpattern_default, or default.
    """
    table = morphology.get("cell_lines") or {}
    folded = {str(key).upper(): value for key, value in table.items()}
    for name in (cell_line, stripped):
        entry = table.get(name) or folded.get(str(name).upper())
        if entry and entry.get("accessibility") is not None:
            return float(entry["accessibility"]), "curated"
    defaults = morphology.get("default") or {}
    growth_name = str(growth or "").strip()
    if growth_name in defaults:
        return float(defaults[growth_name]), "growthpattern_default"
    if "Unknown" in defaults:
        return float(defaults["Unknown"]), "default"
    return 0.75, "default"
