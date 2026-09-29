# -*- coding: utf-8 -*-
"""Runtime settings.

Public dataset URLs are code defaults. Optional environment variables can
override them. Values are never invented. Secret values are not logged.
"""

from __future__ import annotations

import os
from pathlib import Path

# Keys this app knows how to read. A TAA-platform env file may contain more;
# unknown keys are ignored and never written into docs.
KNOWN_ENV_KEYS = (
    "TAA_XENA_TOIL_HUB",
    "TAA_XENA_TCGA_HUB",
    "TAA_XENA_PANCAN_HUB",
    "TAA_HPA_IMMUNE_CELL_URL",
    "TAA_CACHE_DIR",
    "TAA_ENV_FILE",
)

DEFAULTS = {
    "TAA_XENA_TOIL_HUB": "https://toil.xenahubs.net",
    "TAA_XENA_TCGA_HUB": "https://tcga.xenahubs.net",
    "TAA_XENA_PANCAN_HUB": "https://pancanatlas.xenahubs.net",
    "TAA_HPA_IMMUNE_CELL_URL": "https://www.proteinatlas.org/download/tsv/rna_immune_cell.tsv.zip",
}


def _parse_env_file(path: Path) -> dict[str, str]:
    found: dict[str, str] = {}
    if not path.is_file():
        return found
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        key = key.strip()
        value = value.strip().strip('"').strip("'")
        if key:
            found[key] = value
    return found


def load_settings(root: Path | None = None) -> dict[str, str]:
    """Return resolved settings and record which known keys came from the environment."""
    root = root or Path(__file__).resolve().parent.parent
    file_keys: dict[str, str] = {}
    env_file = os.environ.get("TAA_ENV_FILE", "").strip()
    candidates = []
    if env_file:
        candidates.append(Path(env_file))
    candidates.append(root / ".env")
    for candidate in candidates:
        file_keys.update(_parse_env_file(candidate))

    resolved = dict(DEFAULTS)
    used_from_env: list[str] = []
    for key in KNOWN_ENV_KEYS:
        if key == "TAA_ENV_FILE":
            continue
        if key in os.environ and os.environ[key].strip():
            resolved[key] = os.environ[key].strip()
            used_from_env.append(key)
        elif key in file_keys and file_keys[key].strip():
            resolved[key] = file_keys[key].strip()
            used_from_env.append(key)

    cache = resolved.get("TAA_CACHE_DIR") or str(root / "cache" / "taa")
    resolved["TAA_CACHE_DIR"] = cache
    resolved["_used_env_keys"] = ",".join(used_from_env)
    resolved["_env_file_present"] = "yes" if any(path.is_file() for path in candidates) else "no"
    return resolved
