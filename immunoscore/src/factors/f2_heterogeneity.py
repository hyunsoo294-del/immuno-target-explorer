"""F2 TAA expression heterogeneity.

Biological claim: an antigen-negative subpopulation caps killing (PBMC) and only
dilutes a bulk reporter. Uniform expression scores higher.

Data behind it: scRNA pos_frac and CV when data/processed/scrna_heterogeneity.parquet
exists. Otherwise a single default pos_frac and CV from scoring_params.yaml, marked
f2_source=default. Nothing is imputed silently.

Limitation: the default is the same number on every line, so it does not rank lines.
It keeps the weight from being dropped while telling the user the factor is not measured.
Heuristic for panel selection, not a predictor.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from immunoscore.src.config_loader import DATA_DIR


def _apply_mode(pos_frac: pd.Series, cv: pd.Series, mode: str, f2: dict) -> pd.Series:
    cv_ref = float(f2["cv_ref"])
    if mode == "cap":
        return 100.0 * pos_frac * (1.0 / (1.0 + cv / cv_ref))
    floor = float(f2["dilute_floor"])
    slope = float(f2["dilute_slope"])
    cv_scale = float(f2["dilute_cv_scale"])
    return 100.0 * (floor + slope * pos_frac) * (1.0 / (1.0 + cv_scale * cv / cv_ref))


def score(cell_lines: pd.DataFrame, taa: str, arm: str, effector: str, cfg: dict) -> pd.DataFrame:
    del arm
    f2 = cfg["params"]["f2"]
    mode = cfg["effectors"]["effectors"][effector]["heterogeneity_mode"]
    ids = cell_lines["ModelID"].astype(str)
    path = Path(cfg.get("scrna_path") or (DATA_DIR / "processed" / "scrna_heterogeneity.parquet"))
    source = "default"
    tier = "default"
    if path.exists():
        table = pd.read_parquet(path)
        table["ModelID"] = table["ModelID"].astype(str)
        subset = table[table["gene_symbol"] == taa]
        if not subset.empty:
            merged = ids.to_frame(name="ModelID").merge(subset, on="ModelID", how="left")
            has = merged["pos_frac"].notna()
            if has.any():
                source = "scRNAseq"
                tier = "rna_only"
                pos = merged["pos_frac"]
                cv = merged["cv"]
                sub = _apply_mode(pos.fillna(np.nan), cv.fillna(np.nan), mode, f2)
                sub = sub.where(has, np.nan)
                detail = np.where(
                    has,
                    "f2_source=scRNAseq; mode=" + mode,
                    "f2_source=default; no scRNA row",
                )
                # Lines without scRNA fall through to the default below only when none match.
                # Mixed coverage: fill the gaps with the declared default and say so.
                default_pos = float(f2["default_pos_frac"])
                default_cv = float(f2["default_cv"])
                gap = ~has
                if gap.any():
                    pos = pos.where(~gap, default_pos)
                    cv = cv.where(~gap, default_cv)
                    sub = _apply_mode(pos, cv, mode, f2)
                    detail = np.where(has, detail, "f2_source=default; line absent from scRNA")
                return pd.DataFrame(
                    {
                        "ModelID": ids.tolist(),
                        "sub_score": sub.to_numpy(),
                        "evidence_tier": np.where(has, tier, "default").tolist(),
                        "source_detail": detail.tolist(),
                        "pos_frac": pos.to_numpy(),
                        "f2_source": np.where(has, "scRNAseq", "default").tolist(),
                    }
                )
    default_pos = float(f2["default_pos_frac"])
    default_cv = float(f2["default_cv"])
    pos = pd.Series(default_pos, index=range(len(ids)))
    cv = pd.Series(default_cv, index=range(len(ids)))
    sub = _apply_mode(pos, cv, mode, f2)
    return pd.DataFrame(
        {
            "ModelID": ids.tolist(),
            "sub_score": sub.to_numpy(),
            "evidence_tier": "default",
            "source_detail": f"f2_source=default; pos_frac={default_pos}; cv={default_cv}; mode={mode}",
            "pos_frac": default_pos,
            "f2_source": "default",
        }
    )
