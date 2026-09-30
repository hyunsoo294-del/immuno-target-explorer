"""F7 epitope membrane proximity.

Biological claim: a membrane-proximal epitope crosslinks more readily than a
distal one. Closer scores higher.

Data behind it: epitope_membrane_distance_nm in taa_properties.yaml when the user
or a curated template supplied it. Unknown epitopes use the default score in
scoring_params.yaml and confidence unknown. UniProt topology is not queried here
because an unknown epitope midpoint cannot be placed on the chain.

Limitation: a default of 50 does not rank TAAs against each other. The nanometer
bins are a coarse map, not a structure. Heuristic for panel selection, not a predictor.
"""

from __future__ import annotations

import pandas as pd

from immunoscore.src.factors.common import empty_factor


def _bin_score(distance_nm: float, bins: list) -> float:
    for item in bins:
        cap = item.get("max_nm")
        if cap is None or distance_nm < float(cap):
            return float(item["score"])
    return float(bins[-1]["score"])


def score(cell_lines: pd.DataFrame, taa: str, arm: str, effector: str, cfg: dict) -> pd.DataFrame:
    del taa, arm, effector
    entry = cfg.get("taa_entry") or {}
    distance = entry.get("epitope_membrane_distance_nm")
    params = cfg["params"]
    if distance in (None, ""):
        sub = float(params["f7_unknown_score"])
        frame = empty_factor(cell_lines, "default", "epitope unknown; default score; confidence unknown")
        frame["sub_score"] = sub
        return frame
    sub = _bin_score(float(distance), params["f7_bins"])
    frame = empty_factor(
        cell_lines,
        "default",
        f"curated distance {float(distance):.1f} nm; confidence {entry.get('confidence', 'unknown')}",
    )
    frame["sub_score"] = sub
    return frame
