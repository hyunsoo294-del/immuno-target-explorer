"""F6 glycocalyx and mucin barrier.

Biological claim: long mucins and sialylation physically obstruct crosslinking.
The barrier is effector-independent, so this weight stays high in every arm.
A lower barrier scores higher.

Data behind it: DepMap log2(TPM+1). Mucin genes are z-scored, then averaged with
weights proportional to extracellular_length_nm in mucins.yaml, so MUC16 counts
more than MUC13. Sialylation genes are an unweighted z mean. The 0.65/0.35 mix
is in scoring_params.yaml.

Limitation: length_nm values are order-of-magnitude constants, not assay measurements.
Bulk RNA is not the glycocalyx thickness. Heuristic for panel selection, not a predictor.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from immunoscore.src.factors.common import cohort_z, sigmoid_scaled


def _length_weighted_z(expr: pd.DataFrame, genes: list, lengths: dict, ddof: int) -> pd.Series:
    present = []
    for gene in genes:
        if gene in expr.columns and gene in lengths:
            present.append(gene)
    if not present:
        return pd.Series(np.nan, index=expr.index)
    total = float(sum(lengths[gene] for gene in present))
    acc = pd.Series(0.0, index=expr.index)
    coverage = pd.Series(0.0, index=expr.index)
    for gene in present:
        share = float(lengths[gene]) / total
        z = cohort_z(expr[gene], ddof)
        valid = z.notna()
        acc = acc.add(z.fillna(0.0) * share, fill_value=0.0)
        coverage = coverage.add(valid.astype(float) * share, fill_value=0.0)
    return (acc / coverage.where(coverage > 0, np.nan)).where(coverage > 0, np.nan)


def score(cell_lines: pd.DataFrame, taa: str, arm: str, effector: str, cfg: dict) -> pd.DataFrame:
    del taa, arm, effector
    params = cfg["params"]
    expr = cfg["expression"]
    ids = cell_lines["ModelID"].astype(str)
    ddof = int(params["zscore_ddof"])
    lengths = {
        gene: float(entry["extracellular_length_nm"])
        for gene, entry in (cfg["mucins"].get("mucins") or {}).items()
    }
    mucin = _length_weighted_z(expr, cfg["genes"]["mucin_bulk"], lengths, ddof)
    sial_genes = [gene for gene in cfg["genes"]["sialylation"] if gene in expr.columns]
    if sial_genes:
        sial = pd.concat({gene: cohort_z(expr[gene], ddof) for gene in sial_genes}, axis=1).mean(axis=1, skipna=True)
    else:
        sial = pd.Series(np.nan, index=expr.index)
    spec = params["f6"]
    barrier = float(spec["mucin_weight"]) * mucin + float(spec["sialylation_weight"]) * sial
    scaled = sigmoid_scaled(barrier, float(spec["sigmoid_center"]), float(spec["sigmoid_width"]))
    sub = (100.0 - 100.0 * scaled).reindex(ids)
    sub = sub.where(barrier.reindex(ids).notna(), np.nan)
    return pd.DataFrame(
        {
            "ModelID": ids.tolist(),
            "sub_score": sub.to_numpy(),
            "evidence_tier": "rna_only",
            "source_detail": "length-weighted mucin RNA plus sialylation RNA; lengths are unvalidated constants",
        }
    )
