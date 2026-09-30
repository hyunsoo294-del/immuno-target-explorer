"""Turn the 0-100 score into a predicted Emax.

log10(Emax) = a + b * score, fit by ordinary least squares for one
(taa, arm, effector) combination.

The mapping is monotone. It assigns a magnitude and must not be used to
re-rank lines. Changing a or b cannot repair a wrong order.

A cached fit is valid only for the scores and YAML it was built from.
weights.yaml, scoring_params.yaml, morphology.yaml, and the calibration-line
scores are all part of that key. n < 4 is refused.
"""

from __future__ import annotations

import hashlib
import json
import math
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

from immunoscore.src.config_loader import CONFIG_DIR, DATA_DIR, load_configs
from immunoscore.src.scoring import format_flag_text, score_models

MIN_N = 4
_ABS_SCORE_TOL = 1e-6


def format_emax_interval(pred: float, lo: float, hi: float, extrapolated: bool) -> str:
    if pred is None or not math.isfinite(float(pred)):
        return ""
    text = f"{float(pred):,.0f} ({float(lo):,.0f}–{float(hi):,.0f})"
    if extrapolated:
        text += " 외삽"
    return text


def fit_log_emax(scores, emax_values) -> dict:
    """OLS of log10(Emax) on score. Raises ValueError when n < 4."""
    x = np.asarray(scores, dtype=float)
    y_lin = np.asarray(emax_values, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y_lin) & (y_lin > 0)
    x = x[mask]
    y = np.log10(y_lin[mask])
    n = int(x.size)
    if n < MIN_N:
        raise ValueError(f"n={n} < {MIN_N}")
    x_mean = float(x.mean())
    y_mean = float(y.mean())
    sxx = float(np.sum((x - x_mean) ** 2))
    if sxx <= 0:
        raise ValueError("score variance is 0; the slope is undefined")
    b = float(np.sum((x - x_mean) * (y - y_mean)) / sxx)
    a = float(y_mean - b * x_mean)
    fitted = a + b * x
    residual = y - fitted
    sst = float(np.sum((y - y_mean) ** 2))
    sse = float(np.sum(residual ** 2))
    r2 = float(1.0 - sse / sst) if sst > 0 else float("nan")
    dof = n - 2
    return {
        "a": a,
        "b": b,
        "r2": r2,
        "n": n,
        "dof": dof,
        "s2": float(sse / dof),
        "x_mean": x_mean,
        "sxx": sxx,
        "score_min": float(x.min()),
        "score_max": float(x.max()),
    }


def predict_emax(fit: dict, score: float) -> dict:
    """Point estimate and 95% prediction interval on the Emax scale."""
    from scipy.stats import t as student_t

    yhat = float(fit["a"]) + float(fit["b"]) * float(score)
    n = int(fit["n"])
    dof = int(fit["dof"])
    se = math.sqrt(
        float(fit["s2"])
        * (1.0 + 1.0 / n + (float(score) - float(fit["x_mean"])) ** 2 / float(fit["sxx"]))
    )
    tcrit = float(student_t.ppf(0.975, dof))
    lo = yhat - tcrit * se
    hi = yhat + tcrit * se
    extrapolated = float(score) < float(fit["score_min"]) - _ABS_SCORE_TOL or float(score) > float(fit["score_max"]) + _ABS_SCORE_TOL
    return {
        "emax_pred": 10.0 ** yhat,
        "emax_lo": 10.0 ** lo,
        "emax_hi": 10.0 ** hi,
        "extrapolated": extrapolated,
    }


def yaml_fingerprint() -> str:
    digest = hashlib.sha256()
    for path in (
        CONFIG_DIR / "weights.yaml",
        CONFIG_DIR / "scoring_params.yaml",
        DATA_DIR / "curated" / "morphology.yaml",
    ):
        digest.update(path.name.encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def _canon(name: str) -> str:
    token = str(name).strip()
    upper = token.upper()
    table = load_configs()["taa"].get("TAAs") or {}
    for symbol, entry in table.items():
        if symbol.upper() == upper:
            return symbol
        for alias in entry.get("aliases") or []:
            if str(alias).upper() == upper:
                return symbol
    return token


def _cache_path(symbol: str, arm: str, effector: str) -> Path:
    return DATA_DIR / "processed" / f"calibration_{symbol}_{arm}_{effector}.json"


def _scores_match(cached: dict, current: dict) -> bool:
    if set(cached) != set(current):
        return False
    for name, value in current.items():
        if not math.isclose(float(cached[name]), float(value), rel_tol=0, abs_tol=_ABS_SCORE_TOL):
            return False
    return True


def _caption(fit: dict) -> str:
    return (
        f"log10(Emax) = {float(fit['a']):.4f} + {float(fit['b']):.5f} × 점수"
        f" · n={int(fit['n'])} · R²={float(fit['r2']):.3f} · {fit['fit_date']}"
    )


def _empty_columns(result: pd.DataFrame) -> pd.DataFrame:
    out = result.copy()
    out["emax_pred"] = np.nan
    out["emax_lo"] = np.nan
    out["emax_hi"] = np.nan
    out["emax_extrapolated"] = False
    out["emax_text"] = ""
    return out


def _refusal(result: pd.DataFrame, reason: str) -> tuple[pd.DataFrame, dict]:
    print("calibration refused: " + reason)
    meta = {
        "ok": False,
        "caption": "캘리브레이션 없음 — 순위만 유효",
        "reason": reason,
        "a": None,
        "b": None,
        "r2": None,
        "n": 0,
        "fit_date": "",
    }
    return _empty_columns(result), meta


def _load_measured() -> pd.DataFrame:
    path = DATA_DIR / "curated" / "measured_emax.csv"
    if not path.exists():
        return pd.DataFrame(columns=["taa", "cell_line", "arm", "effector", "emax_value"])
    return pd.read_csv(path)


def _matching_rows(table: pd.DataFrame, taa: str, arm: str, effector: str) -> pd.DataFrame:
    if table.empty:
        return table
    symbol = _canon(taa)
    mask = (
        table["taa"].map(lambda value: _canon(value) == symbol)
        & table["arm"].astype(str).eq(str(arm))
        & table["effector"].astype(str).eq(str(effector))
    )
    return table.loc[mask].copy()


def attach_calibration(
    result: pd.DataFrame,
    taa: str,
    arm: str,
    effector: str,
    cfg: dict,
) -> tuple[pd.DataFrame, dict]:
    """Fit or reuse the combination's Emax line and attach intervals.

    Does not sort the frame. Predicted Emax stays empty when the fit is refused.
    """
    symbol = _canon(taa)
    measured = _matching_rows(_load_measured(), symbol, arm, effector)
    if len(measured) < MIN_N:
        have = measured["cell_line"].astype(str).tolist() if len(measured) else []
        return _refusal(
            result,
            f"{symbol} × {arm} × {effector}: measured n={len(measured)}, need {MIN_N}. rows={have or 'none'}.",
        )
    names = measured["cell_line"].astype(str).tolist()
    models = cfg.get("models")
    if models is None or "CellLineName" not in models.columns:
        return _refusal(result, f"{symbol} × {arm} × {effector}: model table is missing, so the lines cannot be scored.")
    chosen = models[models["CellLineName"].astype(str).isin(names)]
    found = set(chosen["CellLineName"].astype(str))
    missing = [name for name in names if name not in found]
    scored = score_models(chosen, symbol, arm, effector, cfg)
    paired = measured.merge(scored[["cell_line", "score"]], on="cell_line", how="left")
    paired = paired[pd.to_numeric(paired["score"], errors="coerce").notna()].copy()
    if len(paired) < MIN_N:
        unscored = [name for name in names if name not in set(paired["cell_line"].astype(str))]
        return _refusal(
            result,
            f"{symbol} × {arm} × {effector}: scored n={len(paired)}, need {MIN_N}. "
            f"not in DepMap: {missing or 'none'}. unscored: {unscored or 'none'}.",
        )
    score_map = {
        str(row.cell_line): float(row.score)
        for row in paired.itertuples(index=False)
    }
    fingerprint = yaml_fingerprint()
    path = _cache_path(symbol, arm, effector)
    fit = _fit_from_cache(path, fingerprint, score_map)
    if fit is None:
        try:
            fit = fit_log_emax(paired["score"].to_numpy(), paired["emax_value"].to_numpy())
        except ValueError as exc:
            return _refusal(result, f"{symbol} × {arm} × {effector}: {exc}")
        fit["fit_date"] = date.today().isoformat()
        fit["taa"] = symbol
        fit["arm"] = arm
        fit["effector"] = effector
        fit["scores"] = score_map
        fit["yaml_hash"] = fingerprint
        fit["emax"] = {
            str(row.cell_line): float(row.emax_value)
            for row in paired.itertuples(index=False)
        }
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(fit, ensure_ascii=False, indent=2), encoding="utf-8")
    meta = {
        "ok": True,
        "caption": _caption(fit),
        "reason": "",
        "a": fit["a"],
        "b": fit["b"],
        "r2": fit["r2"],
        "n": fit["n"],
        "fit_date": fit["fit_date"],
        "scores": fit["scores"],
        "taa": symbol,
        "arm": arm,
        "effector": effector,
    }
    return _apply_predictions(result, fit, cfg), meta


def _fit_from_cache(path: Path, fingerprint: str, score_map: dict) -> dict | None:
    if not path.exists():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    required = ("a", "b", "r2", "n", "dof", "s2", "x_mean", "sxx", "score_min", "score_max", "fit_date", "scores", "yaml_hash")
    if any(key not in payload for key in required):
        return None
    if payload["yaml_hash"] != fingerprint or not _scores_match(payload["scores"], score_map):
        return None
    return payload


def _apply_predictions(result: pd.DataFrame, fit: dict, cfg: dict) -> pd.DataFrame:
    out = _empty_columns(result)
    if out.empty or "score" not in out.columns:
        return out
    labels = (cfg.get("params") or {}).get("flag_labels") or {}
    factor_labels = (cfg.get("weights") or {}).get("factor_labels") or {}
    preds = []
    los = []
    his = []
    extras = []
    texts = []
    flag_lists = []
    flag_texts = []
    for record in out.itertuples(index=False):
        score = getattr(record, "score")
        flags = list(getattr(record, "flags") or [])
        zero = [item for item in str(getattr(record, "zero_variance_factors", "") or "").split(",") if item]
        names = [factor_labels.get(item, item) for item in zero]
        if pd.isna(score):
            preds.append(np.nan)
            los.append(np.nan)
            his.append(np.nan)
            extras.append(False)
            texts.append("")
            flag_lists.append(flags)
            flag_texts.append(format_flag_text(flags, labels, names))
            continue
        pred = predict_emax(fit, float(score))
        preds.append(pred["emax_pred"])
        los.append(pred["emax_lo"])
        his.append(pred["emax_hi"])
        extras.append(pred["extrapolated"])
        texts.append(format_emax_interval(pred["emax_pred"], pred["emax_lo"], pred["emax_hi"], pred["extrapolated"]))
        if pred["extrapolated"] and "EXTRAPOLATED" not in flags:
            flags.append("EXTRAPOLATED")
        flag_lists.append(flags)
        flag_texts.append(format_flag_text(flags, labels, names))
    out["emax_pred"] = preds
    out["emax_lo"] = los
    out["emax_hi"] = his
    out["emax_extrapolated"] = extras
    out["emax_text"] = texts
    out["flags"] = flag_lists
    out["flag_text"] = flag_texts
    return out
