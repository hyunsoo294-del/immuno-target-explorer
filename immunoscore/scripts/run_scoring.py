"""CLI for the same engine the UI uses.

python immunoscore/scripts/run_scoring.py --taa HER2 --arm 4-1BB --effector Jurkat_NFkB \
    --diseases "Invasive Breast Carcinoma" --top 30 --out results/her2.xlsx
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from immunoscore.src.genes import gene_identity  # noqa: E402
from immunoscore.src.report import to_excel  # noqa: E402
from immunoscore.src.scoring import build_cfg, score_models  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--taa", required=True)
    parser.add_argument("--arm", required=True)
    parser.add_argument("--effector", required=True)
    parser.add_argument("--diseases", default="", help="Comma-separated OncotreePrimaryDisease values")
    parser.add_argument("--lines", default="", help="Comma-separated CellLineName or StrippedCellLineName")
    parser.add_argument("--top", type=int, default=0)
    parser.add_argument("--out", default="")
    parser.add_argument("--jurkat-pd1", action="store_true")
    parser.add_argument("--assay-hours", type=float, default=None)
    args = parser.parse_args()
    identity = gene_identity(args.taa)
    if identity.get("refused"):
        print(identity.get("warning") or "refused")
        genes = identity.get("biosynthesis_genes") or []
        if genes:
            print("biosynthesis genes (not scored): " + ", ".join(genes))
        return
    overrides = {
        "jurkat_pd1": args.jurkat_pd1,
        "taa_entry": identity["curated_entry"],
        "resolution_flags": list(identity.get("flags") or []),
    }
    if identity.get("confidence_override"):
        overrides["confidence_override"] = identity["confidence_override"]
    if args.assay_hours is not None:
        overrides["assay_duration_h"] = args.assay_hours
    cfg = build_cfg(overrides)
    models = cfg["models"]
    if args.lines:
        wanted = {item.strip().upper() for item in args.lines.split(",") if item.strip()}
        mask = models["CellLineName"].astype(str).str.upper().isin(wanted) | models[
            "StrippedCellLineName"
        ].astype(str).str.upper().isin(wanted)
        chosen = models.loc[mask]
    else:
        diseases = [item.strip() for item in args.diseases.split(",") if item.strip()]
        chosen = models[models["cancer"].isin(diseases)] if diseases else models
    if args.top and args.top > 0:
        expr = cfg["expression"]
        gene = identity["gene_symbol"]
        order = expr[gene].reindex(chosen["ModelID"]).sort_values(ascending=False)
        keep = list(order.head(args.top).index)
        chosen = chosen[chosen["ModelID"].isin(keep)]
    result = score_models(chosen, identity["gene_symbol"], args.arm, args.effector, cfg)
    print(result[["cell_line", "cancer", "score", "confidence", "flag_text"]].head(30).to_string(index=False))
    if args.out:
        dest = Path(args.out)
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(
            to_excel(
                result,
                {
                    "taa": identity["gene_symbol"],
                    "arm": args.arm,
                    "effector": args.effector,
                    "manifest": cfg.get("manifest"),
                },
            )
        )
        print(f"wrote {dest}")


if __name__ == "__main__":
    main()
