"""F1 TAA expression, used as a multiplicative gate.

Biological claim: in a fixed effector-to-target reporter assay the signal ceiling
is set by how many effector cells form a productive conjugate. Antigen density is
a threshold. Below it there is no signal. Above it, density neither raises nor
lowers the score. This is deliberate. In the 2026-09 HER2 panel, density and
contact morphology run in opposite directions, so a descending or bell-shaped F1
would just memorize that panel. Antigen depletion is not the explanation: even
the highest expressor presents about 0.4 nM of sites at typical plating density,
two orders below the top of the dose-response.

Data behind it: approximate ABC for five HER2 lines in measured_abc.csv, otherwise
DepMap log2(TPM+1) through an unvalidated RNA-to-ABC curve. ABC50_base is 30000,
then multiplied by the effector's abc50_multiplier.

Limitation: RNA ABC is not a molecule count. The five panel ABCs are approximate,
not a QIFIKIT fit. Heuristic for panel selection, not a validated predictor.
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
    f1_params = params.get("f1_expression") or {}
    abc50_base = float(cfg.get("abc50_base") or f1_params.get("abc50_base"))
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
