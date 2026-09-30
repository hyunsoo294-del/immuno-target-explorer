"""F5 adhesion and synapse capacity.

Biological claim: ICAM1 and CD58 on the tumor cell support an LFA-1/CD2 synapse.
Higher adhesion RNA scores higher. The effector must actually carry the counter-receptor.

Data behind it: DepMap log2(TPM+1) z-scores, gene weights in gene_sets.yaml, then a
cohort percentile rank. synapse_competence is applied in scoring.py, not here.

Limitation: RNA is not surface ICAM1 density, and HEK293 has almost no counter-receptor.
The competence scaler is a judgment, not a measured binding constant.
Heuristic for panel selection, not a predictor.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from immunoscore.src.factors.common import weighted_z


def score(cell_lines: pd.DataFrame, taa: str, arm: str, effector: str, cfg: dict) -> pd.DataFrame:
    del taa, arm, effector
    params = cfg["params"]
    expr = cfg["expression"]
    ids = cell_lines["ModelID"].astype(str)
    index = weighted_z(expr, cfg["genes"]["adhesion"], int(params["zscore_ddof"]))
    # Rank against the full cohort so the score does not depend on which lines were picked.
    ranked = index.rank(pct=True, method="average") * 100.0
    sub = ranked.reindex(ids)
    return pd.DataFrame(
        {
            "ModelID": ids.tolist(),
            "sub_score": sub.to_numpy(),
            "evidence_tier": np.where(sub.notna(), "rna_only", "rna_only"),
            "source_detail": "adhesion RNA percentile of the weighted z-score; competence applied at assembly",
            "raw_score": sub.to_numpy(),
        }
    )
