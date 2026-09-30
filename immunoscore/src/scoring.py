"""Assemble factor scores, apply effector gating, and attach flags.

Gating lives here so factor modules stay pure. Missing factors are NaN and their
weights are renormalized. A missing value is not a zero.

This is a panel-selection heuristic, not a validated predictor of assay Emax.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from immunoscore.src.confidence import completeness, row_confidence
from immunoscore.src.config_loader import combination_key
from immunoscore.src.factors import FACTORS


def _as_bool(value) -> bool:
    return bool(value) and value is not False


def apply_effector_gating(parts: dict[str, pd.DataFrame], effector: str, cfg: dict) -> list[str]:
    """Mutate F4 and F5 frames. Return factor ids that are inert for this effector."""
    eff = cfg["effectors"]["effectors"][effector]
    inert = []
    jurkat_pd1 = bool(cfg.get("jurkat_pd1"))
    checkpoint_inert = _as_bool(eff.get("checkpoint_inert")) and not (
        jurkat_pd1 and str(effector).startswith("Jurkat")
    )
    f4 = parts["F4_checkpoint"]
    if checkpoint_inert:
        f4["sub_score"] = 100.0
        f4["evidence_tier"] = "not_applicable"
        f4["source_detail"] = (
            "이 효과기는 억제 수용체가 없어 체크포인트 항을 쓰지 않습니다. "
            "생물학적 burden은 점수에 넣지 않았습니다."
        )
        inert.append("F4_checkpoint")
    f5 = parts["F5_adhesion"]
    competence = float(eff["synapse_competence"])
    raw = pd.to_numeric(f5["sub_score"], errors="coerce")
    f5["sub_score"] = raw * competence
    note = str(eff.get("adhesion_note") or "")
    f5["source_detail"] = (
        f5["source_detail"].astype(str) + f"; synapse_competence={competence:.2f}"
        + (f"; {note}" if note else "")
    )
    if _as_bool(eff.get("adhesion_inert")):
        inert.append("F5_adhesion")
    return inert


def _flags_for_row(
    row: pd.Series,
    *,
    effector: str,
    cfg: dict,
    inert_factors: list[str],
    taa_log: float,
    hla_low: bool,
    growth: str,
    ramos_positive: bool,
    surfaceome_missing: bool,
) -> list[str]:
    params = cfg["params"]
    gates = params["gates"]
    flags = []
    s1 = row.get("F1_expression")
    if pd.notna(s1) and float(s1) < float(gates["antigen_low_s1"]):
        flags.append("ANTIGEN_LOW")
    if row.get("f2_source") not in (None, "default") and pd.notna(row.get("pos_frac")):
        if float(row["pos_frac"]) < float(gates["high_heterogeneity_pos_frac"]):
            flags.append("HIGH_HETEROGENEITY")
    s6 = row.get("F6_glycocalyx")
    if pd.notna(s6) and float(s6) < float(gates["mucin_high_s6"]):
        flags.append("MUCIN_HIGH")
    if bool(row.get("jak_loss")):
        flags.append("JAK_LOSS")
    if float(row["data_completeness"]) < float(gates["low_confidence_completeness"]):
        flags.append("LOW_CONFIDENCE")
    if surfaceome_missing:
        flags.append("NO_SURFACE_EVIDENCE")
    if effector == "PBMC" and hla_low:
        flags.append("NK_BACKGROUND_RISK")
    eff = cfg["effectors"]["effectors"][effector]
    if str(growth).lower() == "suspension" and _as_bool(eff.get("suspension")):
        flags.append("SUSPENSION_PAIR")
    present = params["gates"]["taa_present_log2tpm"]
    if effector == "Ramos" and ramos_positive and pd.notna(taa_log) and float(taa_log) >= float(present):
        flags.append("RAMOS_CONFOUND")
    if inert_factors:
        flags.append("EFFECTOR_FACTOR_INERT")
    return flags


def score_models(cell_lines: pd.DataFrame, taa: str, arm: str, effector: str, cfg: dict) -> pd.DataFrame:
    if cell_lines.empty:
        return pd.DataFrame()
    lines = cell_lines.copy()
    lines["ModelID"] = lines["ModelID"].astype(str)
    parts = {}
    for factor_id, function in FACTORS:
        part = function(lines, taa, arm, effector, cfg).set_index("ModelID")
        parts[factor_id] = part
    inert = apply_effector_gating(parts, effector, cfg)
    key = combination_key(arm, effector, bool(cfg.get("jurkat_pd1")))
    weights = dict(cfg["weights"]["combinations"][key]["weights"])
    expr = cfg["expression"]
    taa_series = expr[taa] if taa in expr.columns else pd.Series(dtype=float)
    hla_genes = [gene for gene in cfg["params"]["hla_genes"] if gene in expr.columns]
    hla_threshold = float(cfg["params"]["gates"]["hla_low_log2tpm"])
    ramos_positive = _ramos_expresses(taa, expr, lines, cfg)
    surfaceome = cfg.get("surfaceome")
    surfaceome_missing = surfaceome is not None and taa not in surfaceome

    rows = []
    for model_id in lines["ModelID"]:
        scores = {}
        tiers = {}
        details = {}
        for factor_id, _function in FACTORS:
            part = parts[factor_id]
            if model_id not in part.index:
                scores[factor_id] = np.nan
                tiers[factor_id] = "default"
                details[factor_id] = "missing"
                continue
            record = part.loc[model_id]
            if isinstance(record, pd.DataFrame):
                record = record.iloc[0]
            scores[factor_id] = record["sub_score"]
            tiers[factor_id] = record["evidence_tier"]
            details[factor_id] = record["source_detail"]
        available = {factor: weight for factor, weight in weights.items() if pd.notna(scores.get(factor))}
        weight_sum = float(sum(available.values()))
        if weight_sum <= 0:
            total = np.nan
        else:
            total = sum(float(scores[factor]) * float(weight) / weight_sum for factor, weight in available.items())
        done = completeness(weights, tiers, scores, cfg["params"])
        f2 = parts["F2_heterogeneity"]
        f2_source = "default"
        pos_frac = np.nan
        if model_id in f2.index:
            record = f2.loc[model_id]
            if isinstance(record, pd.DataFrame):
                record = record.iloc[0]
            f2_source = record.get("f2_source", "default")
            pos_frac = record.get("pos_frac", np.nan)
        taa_log = taa_series.get(model_id, np.nan) if len(taa_series) else np.nan
        hla_low = False
        if hla_genes:
            if model_id in expr.index:
                hla_values = expr.loc[model_id, hla_genes]
                if not isinstance(hla_values, pd.Series):
                    hla_values = pd.Series({hla_genes[0]: hla_values})
            else:
                hla_values = pd.Series(dtype=float)
            classic = [gene for gene in ("HLA-A", "HLA-B", "HLA-C") if gene in hla_genes]
            if classic and pd.to_numeric(hla_values[classic], errors="coerce").min(skipna=True) < hla_threshold:
                hla_low = True
            if "B2M" in hla_genes and pd.notna(hla_values.get("B2M")) and float(hla_values["B2M"]) < hla_threshold:
                hla_low = True
        meta = lines.loc[lines["ModelID"] == model_id].iloc[0]
        row = {
            "ModelID": model_id,
            "cell_line": meta.get("CellLineName") or meta.get("StrippedCellLineName") or model_id,
            "cancer": meta.get("cancer") or meta.get("OncotreePrimaryDisease") or "",
            "lineage": meta.get("OncotreeLineage") or "",
            "growth": meta.get("GrowthPattern") or "",
            "score": total,
            "data_completeness": done,
            "confidence": row_confidence(done, str(tiers.get("F1_expression")), cfg["params"]),
            "combination": key,
            "f2_source": f2_source,
            "pos_frac": pos_frac,
            "taa_log2tpm": taa_log,
            "inert_factors": ",".join(inert),
        }
        for factor_id, _function in FACTORS:
            row[factor_id] = scores[factor_id]
            row[factor_id + "_tier"] = tiers[factor_id]
            row[factor_id + "_detail"] = details[factor_id]
        row["jak_loss"] = False
        flags = _flags_for_row(
            pd.Series(row),
            effector=effector,
            cfg=cfg,
            inert_factors=inert,
            taa_log=taa_log,
            hla_low=hla_low,
            growth=str(meta.get("GrowthPattern") or ""),
            ramos_positive=ramos_positive,
            surfaceome_missing=surfaceome_missing,
        )
        labels = cfg["params"]["flag_labels"]
        row["flags"] = flags
        row["flag_text"] = " ".join(labels.get(flag, flag) for flag in flags)
        rows.append(row)
    result = pd.DataFrame(rows)
    result = result.sort_values(["score", "taa_log2tpm"], ascending=[False, False], na_position="last")
    return result.reset_index(drop=True)


def _ramos_expresses(taa: str, expr: pd.DataFrame, lines: pd.DataFrame, cfg: dict) -> bool:
    if taa not in expr.columns:
        return False
    needle = str(cfg["params"].get("ramos_name_contains") or "RAMOS").upper()
    names = lines["StrippedCellLineName"].astype(str).str.upper() if "StrippedCellLineName" in lines.columns else pd.Series(dtype=str)
    # Ramos expression is a property of the effector line, looked up in the full expression index.
    # The scored panel may not contain Ramos. Search the model table stored on cfg when present.
    models = cfg.get("models")
    if models is None or "StrippedCellLineName" not in models.columns:
        hit = lines[names.str.contains(needle, na=False)]
        if hit.empty:
            return False
        value = expr[taa].get(str(hit.iloc[0]["ModelID"]), np.nan)
    else:
        hit = models[models["StrippedCellLineName"].astype(str).str.upper().str.contains(needle, na=False)]
        if hit.empty:
            return False
        value = expr[taa].get(str(hit.iloc[0]["ModelID"]), np.nan)
    threshold = float(cfg["params"]["gates"]["taa_present_log2tpm"])
    return pd.notna(value) and float(value) >= threshold


def build_cfg(overrides: dict | None = None) -> dict:
    """Runtime config: YAML plus the DepMap cohort and any UI overrides."""
    from immunoscore.src.config_loader import load_configs
    from immunoscore.src.data.depmap import cohort_expression, load_models, manifest

    configs = load_configs()
    expression = cohort_expression()
    models = load_models()
    expression = expression.loc[expression.index.intersection(models["ModelID"])]
    models = models[models["ModelID"].isin(expression.index)].reset_index(drop=True)
    cfg = {
        "weights": configs["weights"],
        "effectors": configs["effectors"],
        "params": configs["params"],
        "genes": configs["genes"],
        "mucins": configs["mucins"],
        "expression": expression,
        "models": models,
        "manifest": manifest(),
        "assay_duration_h": configs["params"]["assay_duration_h_default"],
        "abc50_base": configs["params"]["ABC50_base"],
        "hill_h_override": None,
        "jurkat_pd1": False,
        "taa_entry": {},
        "surfaceome": None,
    }
    if overrides:
        cfg.update(overrides)
    return cfg
