"""Kinker et al. 2020 heterogeneity (Broad SCP542).

The portal requires a login, so this script never pretends the data were fetched.
Drop the export in data/raw/scp542/ and re-run. Expected columns after you save a
table as data/processed/scrna_heterogeneity.parquet:
    ModelID, gene_symbol, pos_frac, cv, gini
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent))

from immunoscore.src.config_loader import DATA_DIR  # noqa: E402

RAW = DATA_DIR / "raw" / "scp542"
PROCESSED = DATA_DIR / "processed" / "scrna_heterogeneity.parquet"


def main() -> None:
    if PROCESSED.exists():
        print(f"scRNA heterogeneity table is present: {PROCESSED}")
        return
    print("SCP542 is not in this environment. F2 will use the declared default and mark f2_source=default.")
    print("Manual download (login required):")
    print("  https://singlecell.broadinstitute.org/single_cell/study/SCP542")
    print(f"  Save the study files under {RAW}")
    print("  Then write data/processed/scrna_heterogeneity.parquet with columns")
    print("  ModelID, gene_symbol, pos_frac, cv, gini")
    print("The scoring run will not fail because this file is missing.")
    RAW.mkdir(parents=True, exist_ok=True)


if __name__ == "__main__":
    main()
