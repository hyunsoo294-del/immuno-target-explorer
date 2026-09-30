"""Hand-entered culture accessibility. Not inferred from RNA or GrowthPattern alone."""

from __future__ import annotations


def _entry(cell_line: str, stripped: str, morphology: dict) -> dict | None:
    table = morphology.get("cell_lines") or {}
    folded = {str(key).upper(): value for key, value in table.items()}
    for name in (cell_line, stripped):
        entry = table.get(name) or folded.get(str(name).upper())
        if entry:
            return entry
    return None


def lookup_accessibility(cell_line: str, stripped: str, growth: str, morphology: dict) -> tuple[float, str]:
    """Return (accessibility 0-1, source).

    source is curated, growthpattern_default, or default.
    """
    entry = _entry(cell_line, stripped, morphology)
    if entry and entry.get("accessibility") is not None:
        return float(entry["accessibility"]), "curated"
    defaults = morphology.get("default") or {}
    growth_name = str(growth or "").strip()
    if growth_name in defaults:
        return float(defaults[growth_name]), "growthpattern_default"
    if "Unknown" in defaults:
        return float(defaults["Unknown"]), "default"
    return 0.75, "default"


def effective_accessibility(
    cell_line: str,
    stripped: str,
    growth: str,
    morphology: dict,
    params: dict | None,
) -> tuple[float, str]:
    """Accessibility, times (diameter / reference)^2 only when cell size is enabled.

    use_cell_size false, a null diameter, or a null reference leaves the value unchanged.
    """
    value, source = lookup_accessibility(cell_line, stripped, growth, morphology)
    spec = (params or {}).get("cell_size") or {}
    if not spec.get("use_cell_size"):
        return value, source
    entry = _entry(cell_line, stripped, morphology) or {}
    diameter = entry.get("mean_diameter_um")
    reference = spec.get("d_ref_um")
    if diameter is None or reference is None:
        return value, source
    reference = float(reference)
    if reference <= 0:
        return value, source
    return value * (float(diameter) / reference) ** 2, source
