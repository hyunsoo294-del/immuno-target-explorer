"""Load YAML configs. Paths stay inside the package."""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import yaml

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
CONFIG_DIR = PACKAGE_ROOT / "config"
DATA_DIR = PACKAGE_ROOT / "data"


def load_yaml(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    return data or {}


@lru_cache(maxsize=1)
def load_configs() -> dict:
    return {
        "weights": load_yaml(CONFIG_DIR / "weights.yaml"),
        "effectors": load_yaml(CONFIG_DIR / "effector_params.yaml"),
        "params": load_yaml(CONFIG_DIR / "scoring_params.yaml"),
        "genes": load_yaml(CONFIG_DIR / "gene_sets.yaml"),
        "mucins": load_yaml(DATA_DIR / "reference" / "mucins.yaml"),
        "taa": load_yaml(DATA_DIR / "curated" / "taa_properties.yaml"),
        "morphology": load_yaml(DATA_DIR / "curated" / "morphology.yaml"),
        "depmap_sources": load_yaml(CONFIG_DIR / "depmap_sources.yaml"),
    }


def combination_key(arm: str, effector: str, jurkat_pd1: bool) -> str:
    key = f"{arm}__{effector}"
    if jurkat_pd1 and effector.startswith("Jurkat"):
        variant = key + "_PD1"
        weights = load_configs()["weights"]["combinations"]
        if variant in weights:
            return variant
    return key
