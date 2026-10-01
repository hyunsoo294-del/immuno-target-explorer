"""F5 adhesion and synapse capacity.

Biological claim: ICAM1 and CD58 on the tumor cell support an LFA-1/CD2 synapse.
Higher adhesion RNA scores higher. The effector must actually carry the counter-receptor.

Data behind it: DepMap log2(TPM+1) mapped with a logistic whose median and slope
are frozen on the full cancer cohort, then combined with the gene weights in
gene_sets.yaml. synapse_competence and culture accessibility are applied in
scoring.py, not here.

Limitation: RNA is not surface ICAM1 density, and HEK293 has almost no counter-receptor.
The competence scaler is a judgment, not a measured binding constant. Culture accessibility
is a separate hand-entered multiplier (morphology.yaml). CDH1 RNA cannot tell a flat
monolayer from a dense clump, so that multiplier is not computed from expression.
Heuristic for panel selection, not a predictor.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from immunoscore.src.factors.common import logistic_context, weighted_logistic


def score(cell_lines: pd.DataFrame, taa: str, arm: str, effector: str, cfg: dict) -> pd.DataFrame:
    del taa, arm, effector
    expr = cfg["expression"]
    ids = cell_lines["ModelID"].astype(str)
    reference, scale = logistic_context(cfg)
    combined = weighted_logistic(expr, cfg["genes"]["adhesion"], reference, scale)
    sub = combined.reindex(ids)
    return pd.DataFrame(
        {
            "ModelID": ids.tolist(),
            "sub_score": sub.to_numpy(),
            "evidence_tier": np.where(sub.notna(), "rna_only", "rna_only"),
            "source_detail": "adhesion RNA logistic on the frozen cohort reference; competence applied at assembly",
            "raw_score": sub.to_numpy(),
        }
    )
