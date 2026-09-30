"""F4 checkpoint ligand load.

Biological claim: ligands for inhibitory receptors reduce a PBMC readout. Lower
ligand RNA scores higher.

Data behind it: DepMap log2(TPM+1) z-scores for CD274, PDCD1LG2, CD276, VTCN1,
LGALS9, PVR, NECTIN2, and HLA-E. Weights are in gene_sets.yaml. JAK mutation
status is applied later as a flag only when a mutation table is present.

Limitation: effector gating (parental Jurkat, HEK, Ramos) is not done here.
scoring.py sets this factor to not_applicable when the effector has no inhibitory
receptors. Baseline RNA is not IFN-gamma-induced protein. Heuristic, not a predictor.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from immunoscore.src.factors.common import sigmoid_scaled, weighted_z


def score(cell_lines: pd.DataFrame, taa: str, arm: str, effector: str, cfg: dict) -> pd.DataFrame:
    del taa, arm, effector
    params = cfg["params"]
    expr = cfg["expression"]
    ids = cell_lines["ModelID"].astype(str)
    burden = weighted_z(expr, cfg["genes"]["checkpoint"], int(params["zscore_ddof"])).reindex(ids)
    spec = params["f4"]
    scaled = sigmoid_scaled(burden, float(spec["sigmoid_center"]), float(spec["sigmoid_width"]))
    sub = 100.0 - 100.0 * scaled
    sub = sub.where(burden.notna(), np.nan)
    return pd.DataFrame(
        {
            "ModelID": ids.tolist(),
            "sub_score": sub.to_numpy(),
            "evidence_tier": np.where(burden.notna(), "rna_only", "rna_only"),
            "source_detail": "checkpoint ligand RNA z-score; baseline, not IFN-induced protein",
            "biological_score": sub.to_numpy(),
        }
    )
