"""Fit an RNA log2(TPM+1) to ABC curve from data/curated/measured_abc.csv.

This is the hook for QIFIKIT or flow ABC. With no measured rows it exits without
writing a curve. The scorer keeps using the unvalidated fallback in scoring_params.yaml.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent))

from immunoscore.src.config_loader import DATA_DIR  # noqa: E402
from immunoscore.src.data.depmap import cohort_expression  # noqa: E402


def main() -> None:
    measured_path = DATA_DIR / "curated" / "measured_abc.csv"
    measured = pd.read_csv(measured_path)
    if measured.empty or measured["abc"].dropna().empty:
        print("measured_abc.csv has no ABC rows. Nothing to fit.")
        print("Add columns gene_symbol, ModelID, abc and run this script again.")
        return
    expr = cohort_expression()
    rows = []
    for record in measured.itertuples(index=False):
        gene = str(record.gene_symbol)
        model_id = str(record.ModelID)
        if gene not in expr.columns or model_id not in expr.index:
            continue
        log_value = expr.at[model_id, gene]
        if pd.isna(log_value) or float(record.abc) <= 0:
            continue
        rows.append({"gene_symbol": gene, "log2tpm": float(log_value), "abc": float(record.abc)})
    if len(rows) < 3:
        print(f"Only {len(rows)} paired rows. Need at least 3 to fit a line.")
        return
    import numpy as np

    frame = pd.DataFrame(rows)
    # log(ABC) = intercept + slope * log2(TPM+1). A global fit until a gene has enough points of its own.
    x = frame["log2tpm"].to_numpy()
    y = np.log(frame["abc"].to_numpy())
    slope, intercept = np.polyfit(x, y, 1)
    out = {
        "model": "log(abc) = intercept + slope * log2(TPM+1)",
        "slope": float(slope),
        "intercept": float(intercept),
        "n": int(len(frame)),
        "note": "Global fit. Not used automatically until estimate_surface_abc is pointed at this file.",
    }
    dest = DATA_DIR / "processed" / "abc_calibration.json"
    dest.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"wrote {dest}")


if __name__ == "__main__":
    main()
