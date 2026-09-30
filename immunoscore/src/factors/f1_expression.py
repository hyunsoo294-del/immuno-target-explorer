"""F1 TAA expression.

Biological claim: at saturating antibody, Emax rises with surface antigen density
and then saturates. A Hill function is the claim, not a straight line through TPM.

Data behind it: DepMap log2(TPM+1) passed through an unvalidated global RNA-to-ABC
curve, times a curated surface fraction when one exists. A measured ABC overrides
the curve when data/curated/measured_abc.csv has a row.

Limitation: RNA-to-surface correlation is modest. This tier is rna_only and must
not be read as a molecule count. Heuristic for panel selection, not a predictor.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from immunoscore.src.data.abc import estimate_surface_abc
from immunoscore.src.factors.common import empty_factor


def score(cell_lines: pd.DataFrame, taa: str, arm: str, effector: str, cfg: dict) -> pd.DataFrame:
    del arm
    params = cfg["params"]
    effectors = cfg["effectors"]["effectors"][effector]
    expr = cfg["expression"]
    if taa not in expr.columns:
        return empty_factor(cell_lines, "rna_only", "gene absent from expression matrix")
    ids = cell_lines["ModelID"].astype(str)
    log_values = expr[taa].reindex(ids)
    surface_fraction = float(
        (cfg.get("taa_entry") or {}).get("surface_fraction")
        or params["abc_rna_fallback"]["default_surface_fraction"]
    )
    estimated = estimate_surface_abc(taa, log_values, surface_fraction).set_index("ModelID")
    abc50_base = float(cfg.get("abc50_base") or params["ABC50_base"])
    abc50 = abc50_base * float(effectors["abc50_multiplier"])
    hill = cfg.get("hill_h_override")
    hill_h = float(effectors["hill_h"] if hill in (None, "") else hill)
    abc = estimated["abc"]
    # A real zero is a zero score. A missing value stays missing.
    power_abc = np.power(abc.clip(lower=0), hill_h)
    power_50 = abc50 ** hill_h
    sub = 100.0 * power_abc / (power_50 + power_abc)
    sub = sub.where(abc.notna(), np.nan)
    sub = sub.where(~(abc.notna() & (abc <= 0)), 0.0)
    detail = estimated["source_detail"].astype(str) + f"; abc50={abc50:.0f}; h={hill_h:.2f}"
    return pd.DataFrame(
        {
            "ModelID": ids.tolist(),
            "sub_score": sub.reindex(ids).to_numpy(),
            "evidence_tier": estimated["evidence_tier"].reindex(ids).fillna("rna_only").tolist(),
            "source_detail": detail.reindex(ids).fillna("").tolist(),
            "abc": estimated["abc"].reindex(ids).to_numpy(),
        }
    )
