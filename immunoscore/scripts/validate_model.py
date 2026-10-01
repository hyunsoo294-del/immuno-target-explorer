"""Join predicted scores to measured Emax and write a short validation note.

Spearman and Pearson are reported against log10(Emax). A p-value is printed only
when n >= 8. A correlation on a handful of lines is not evidence.
"""

from __future__ import annotations

import sys
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from immunoscore.src.calibration import attach_calibration  # noqa: E402
from immunoscore.src.config_loader import DATA_DIR  # noqa: E402
from immunoscore.src.genes import gene_identity  # noqa: E402
from immunoscore.src.scoring import build_cfg, score_models  # noqa: E402

MIN_N_FOR_P = 8


def _correlation(predicted: np.ndarray, measured_log: np.ndarray) -> dict:
    from scipy.stats import pearsonr, spearmanr

    mask = np.isfinite(predicted) & np.isfinite(measured_log)
    predicted = np.asarray(predicted, dtype=float)[mask]
    measured_log = np.asarray(measured_log, dtype=float)[mask]
    n = int(len(predicted))
    out = {"n": n, "spearman": None, "pearson": None, "p_value": None}
    if n < 3 or np.nanstd(predicted) == 0 or np.nanstd(measured_log) == 0:
        return out
    spearman = spearmanr(predicted, measured_log)
    pearson = pearsonr(predicted, measured_log)
    out["spearman"] = float(spearman.statistic)
    out["pearson"] = float(pearson.statistic)
    if n >= MIN_N_FOR_P:
        out["p_value"] = float(spearman.pvalue)
    return out


def evaluate(measured: pd.DataFrame | None = None) -> str:
    path = DATA_DIR / "curated" / "measured_emax.csv"
    table = measured if measured is not None else pd.read_csv(path)
    lines = []
    lines.append(f"# Validation {date.today().isoformat()}")
    lines.append("")
    lines.append("Predicted score versus dose-response Emax. This is a check, not a claim of validity.")
    lines.append("")
    for keys, group in table.groupby(["taa", "arm", "effector"], dropna=False):
        taa, arm, effector = keys
        identity = gene_identity(str(taa))
        cfg = build_cfg({"taa_entry": identity["curated_entry"]})
        names = group["cell_line"].astype(str).tolist()
        chosen = cfg["models"][cfg["models"]["CellLineName"].isin(names)]
        predicted = score_models(chosen, identity["gene_symbol"], str(arm), str(effector), cfg)
        predicted, fit = attach_calibration(predicted, identity["gene_symbol"], str(arm), str(effector), cfg)
        merged = group.merge(predicted, left_on="cell_line", right_on="cell_line", how="left")
        n = int(merged["score"].notna().sum())
        lines.append(f"## {taa} × {arm} × {effector}")
        lines.append("")
        lines.append(f"**n = {n}**")
        lines.append("")
        if n < MIN_N_FOR_P:
            lines.append(f"n is below {MIN_N_FOR_P}, so no p-value is reported.")
            lines.append("")
        measured_log = np.log10(merged["emax_value"].astype(float).to_numpy())
        scores = merged["score"].astype(float).to_numpy()
        stats = _correlation(scores, measured_log)
        above = merged[pd.to_numeric(merged["F1_expression"], errors="coerce") >= 30]
        lines.append(f"- Spearman vs log10(Emax): {_fmt(stats['spearman'])}")
        lines.append(f"- Pearson vs log10(Emax): {_fmt(stats['pearson'])}")
        pred_order = merged.sort_values("score", ascending=False)["cell_line"].tolist()
        meas_order = merged.sort_values("emax_value", ascending=False)["cell_line"].tolist()
        paired = list(zip(pred_order, meas_order))
        agree = sum(1 for left, right in paired if left == right)
        lines.append(f"- Rank agreement: {agree}/{n}")
        lines.append(f"- Predicted order: {', '.join(pred_order)}")
        lines.append(f"- Measured order: {', '.join(meas_order)}")
        if n:
            pred_range = float(np.nanmax(scores) / np.nanmin(scores)) if np.nanmin(scores) > 0 else float("nan")
            meas_range = float(merged["emax_value"].max() / merged["emax_value"].min())
            lines.append(f"- Predicted dynamic range: {pred_range:.2f}x")
            lines.append(f"- Measured dynamic range: {meas_range:.2f}x")
        lines.append("")
        lines.append("| cell line | score | predicted Emax | measured Emax | ratio |")
        lines.append("|---|---:|---:|---:|---:|")
        for record in merged.sort_values("emax_value", ascending=False).itertuples(index=False):
            score = record.score
            score_text = "" if pd.isna(score) else f"{float(score):.1f}"
            pred = getattr(record, "emax_pred", float("nan"))
            if pd.isna(pred):
                pred_text = ""
                ratio_text = ""
            else:
                pred_text = f"{float(pred):,.0f}"
                ratio_text = f"{float(pred) / float(record.emax_value):.2f}"
            lines.append(
                f"| {record.cell_line} | {score_text} | {pred_text} | {int(record.emax_value)} | {ratio_text} |"
            )
        lines.append("")
        if fit.get("ok"):
            lines.append(fit["caption"])
            lines.append("")
            lines.append("Calibration is a monotone transform of the score. It does not re-rank.")
        else:
            lines.append("캘리브레이션 없음 — 순위만 유효")
            lines.append("")
            lines.append(fit.get("reason") or "")
        lines.append("")
        lines.append(
            "Per-factor Spearman with log10(Emax). The antigen-positive column keeps lines with "
            "F1 gate ≥ 30, which is the set where adhesion was expected to carry the rank. "
            "Undefined when the factor does not vary. No p-values at this n."
        )
        lines.append("")
        lines.append("| factor | all lines | antigen-positive |")
        lines.append("|---|---:|---:|")
        above_log = np.log10(above["emax_value"].astype(float).to_numpy()) if len(above) else np.array([])
        for factor in (
            "F1_expression",
            "F2_heterogeneity",
            "F3_internalization",
            "F4_checkpoint",
            "F5_adhesion",
            "F6_glycocalyx",
            "F7_epitope_proximity",
            "accessibility",
        ):
            if factor not in merged.columns:
                continue
            values = pd.to_numeric(merged[factor], errors="coerce").to_numpy()
            factor_stats = _correlation(values, measured_log)
            if len(above):
                above_stats = _correlation(pd.to_numeric(above[factor], errors="coerce").to_numpy(), above_log)
                above_text = _fmt(above_stats["spearman"])
            else:
                above_text = "undefined"
            lines.append(f"| {factor} | {_fmt(factor_stats['spearman'])} | {above_text} |")
        lines.append("")
    return "\n".join(lines) + "\n"


def _fmt(value) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return "undefined"
    return f"{value:.3f}"


def main() -> None:
    text = evaluate()
    dest = ROOT / "docs" / f"validation_{date.today().isoformat()}.md"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(text, encoding="utf-8")
    print(text)
    print(f"wrote {dest}")


if __name__ == "__main__":
    main()
