# -*- coding: utf-8 -*-
"""Recompute RNA fold, drop IHC pseudocount folds, and flag bad CPTAC scales."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from csv_io import read_csv_safe
from data_qc import drop_ihc_pseudocount_fold, remap_rna_arrays, to_cdap_log_ratio

CSV_DIR = Path(__file__).resolve().parent / "expression_csv"


def recalc_comparison(path: Path) -> tuple[Path, int]:
    df = read_csv_safe(path)
    if df.empty:
        return path, 0
    repaired = drop_ihc_pseudocount_fold(remap_rna_arrays(df))
    repaired.to_csv(path, index=False, encoding="utf-8-sig")
    ok = int(pd.Series(repaired.get("tcga_inverse_sanity_ok")).fillna(False).sum()) if "tcga_inverse_sanity_ok" in repaired.columns else 0
    return path, ok


def recalc_cptac_summary(path: Path) -> tuple[Path, int]:
    df = read_csv_safe(path)
    if df.empty or "median_log2_ratio" not in df.columns:
        return path, 0
    converted, scale = to_cdap_log_ratio(df["median_log2_ratio"])
    # Summary-only files cannot recover sample-level ratios. Reject out-of-range medians.
    med = pd.to_numeric(df["median_log2_ratio"], errors="coerce")
    bad = med.abs() > 8
    df["cptac_scale"] = "cdap_log_ratio"
    df.loc[bad, "cptac_scale"] = "rejected_raw_intensity_refetch_required"
    df.loc[bad, "median_log2_ratio"] = pd.NA
    df.loc[bad, "mean_log2_ratio"] = pd.NA
    df.loc[bad, "fold_vs_reference"] = pd.NA
    df.to_csv(path, index=False, encoding="utf-8-sig")
    return path, int((~bad).sum())


def main() -> None:
    for path in sorted(CSV_DIR.glob("*_rna_protein_comparison_*.csv")):
        _, ok = recalc_comparison(path)
        print(f"{path.name}: RNA remapped; TCGA inverse OK rows={ok}")
    for path in sorted(CSV_DIR.glob("*_cptac_protein_*.csv")):
        _, kept = recalc_cptac_summary(path)
        print(f"{path.name}: CPTAC in-range rows kept={kept}")


if __name__ == "__main__":
    main()
