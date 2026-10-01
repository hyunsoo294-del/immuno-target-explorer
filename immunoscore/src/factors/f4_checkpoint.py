"""F4 checkpoint ligand load.

Biological claim: ligands for inhibitory receptors reduce a PBMC readout. Lower
ligand RNA scores higher.

Data behind it: DepMap log2(TPM+1) logistics for CD274, PDCD1LG2, CD276, VTCN1,
LGALS9, PVR, NECTIN2, and HLA-E, anchored to the frozen cohort reference.
Weights are in gene_sets.yaml. The factor score is 100 minus that mix.
JAK mutation status is applied later as a flag only when a mutation table is present.

Limitation: effector gating (parental Jurkat, HEK, Ramos) is not done here.
scoring.py sets this factor to not_applicable when the effector has no inhibitory
receptors. Baseline RNA is not IFN-gamma-induced protein. Heuristic, not a predictor.
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
    burden = weighted_logistic(expr, cfg["genes"]["checkpoint"], reference, scale).reindex(ids)
    sub = scale - burden
    sub = sub.where(burden.notna(), np.nan)
    return pd.DataFrame(
        {
            "ModelID": ids.tolist(),
            "sub_score": sub.to_numpy(),
            "evidence_tier": np.where(burden.notna(), "rna_only", "rna_only"),
            "source_detail": "checkpoint ligand RNA logistic on the frozen cohort reference; baseline, not IFN-induced protein",
            "biological_score": sub.to_numpy(),
        }
    )
