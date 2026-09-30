"""F3 TAA internalization.

Biological claim: antibody that leaves the surface during a 24 h assay cannot
keep crosslinking. Slow internalizers score higher.

Data behind it: a curated per-TAA rate and surface half-life in taa_properties.yaml.
There is no per-cell-line antibody internalization database. A weak modifier from
endocytic-gene RNA (CLTC, AP2M1, DNM2, RAB5A, RAB7A, CAV1) is an unvalidated proxy
and is capped.

Limitation: the modifier has not been checked against measured internalization.
Missing curated entries stay NaN. Heuristic for panel selection, not a predictor.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from immunoscore.src.factors.common import clip_score, cohort_z, empty_factor


def score(cell_lines: pd.DataFrame, taa: str, arm: str, effector: str, cfg: dict) -> pd.DataFrame:
    del arm, effector, taa
    entry = cfg.get("taa_entry") or {}
    rate_name = entry.get("internalization_rate")
    half_life = entry.get("half_life_surface_h")
    rate_map = cfg["params"]["internalization_rate_score"]
    if rate_name not in rate_map or half_life in (None, ""):
        return empty_factor(cell_lines, "default", "no curated internalization entry")
    params = cfg["params"]
    duration = float(cfg.get("assay_duration_h") or params["assay_duration_h_default"])
    rate_score = float(rate_map[rate_name])
    retention = 0.5 ** (duration / float(half_life))
    alpha = float(params["retention_blend_alpha"])
    base = alpha * rate_score + (1.0 - alpha) * (100.0 * retention)
    ids = cell_lines["ModelID"].astype(str)
    expr = cfg["expression"]
    genes = [gene for gene in params["endocytic_genes"] if gene in expr.columns]
    modifier = pd.Series(0.0, index=ids)
    if genes:
        z = pd.concat({gene: cohort_z(expr[gene], int(params["zscore_ddof"])) for gene in genes}, axis=1)
        mean_z = z.reindex(ids).mean(axis=1, skipna=True)
        cap = float(params["internalization_modifier_cap"])
        z_at_cap = float(params["internalization_modifier_z_at_cap"])
        modifier = (mean_z / z_at_cap * cap).clip(-cap, cap).fillna(0.0)
    # High endocytic RNA lowers the score. The modifier is an unvalidated proxy.
    sub = clip_score(pd.Series(base, index=ids) - modifier)
    detail = (
        f"rate={rate_name} ({rate_score:.0f}); retention={retention:.3f}; "
        f"duration_h={duration:.0f}; endocytic modifier is an unvalidated proxy"
    )
    return pd.DataFrame(
        {
            "ModelID": ids.tolist(),
            "sub_score": sub.to_numpy(),
            "evidence_tier": "default",
            "source_detail": detail,
        }
    )
