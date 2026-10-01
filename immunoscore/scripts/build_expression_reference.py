"""Freeze cohort medians and logistic slopes for expression factors.

Run this when the DepMap matrix or cohort_logistic.span changes.
Scoring does not recompute these numbers from the user's selection.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from immunoscore.src.config_loader import load_configs  # noqa: E402
from immunoscore.src.data.depmap import cohort_expression, manifest  # noqa: E402
from immunoscore.src.data.expression_reference import build_reference_frame, write_reference  # noqa: E402


def main() -> None:
    load_configs.cache_clear()
    params = load_configs()["params"]["cohort_logistic"]
    span = float(params["span"])
    expression = cohort_expression()
    frame = build_reference_frame(expression, span)
    info = manifest()
    write_reference(
        frame,
        {
            "release": info.get("release") or "",
            "release_id": info.get("release_id") or "",
            "doi": info.get("doi") or "",
            "value": info.get("value") or "log2(TPM+1)",
            "cohort": "cancer lines; Normal and Fibroblast dropped",
            "n_models": int(expression.shape[0]),
            "n_genes": int(expression.shape[1]),
            "span": span,
            "quantile": "linear",
        },
    )
    finite = int(frame["k"].notna().sum())
    print(f"wrote {frame.shape[0]} genes, {finite} with a finite k, from {expression.shape[0]} lines")


if __name__ == "__main__":
    main()
