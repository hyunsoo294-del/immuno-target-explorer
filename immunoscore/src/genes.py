"""Resolve a TAA query to one DepMap gene symbol."""

from __future__ import annotations

import re

from immunoscore.src.config_loader import load_configs
from immunoscore.src.data.depmap import load_expression, load_gene_ids

_EXTRA_ALIASES = {
    "PD-L1": "CD274",
    "PDL1": "CD274",
    "B7-H1": "CD274",
    "B7H1": "CD274",
    "B7-H3": "CD276",
    "PD-L2": "PDCD1LG2",
    "EPCAM": "EPCAM",
    "CD326": "EPCAM",
    "STEAP1": "STEAP1",
    "MSLN": "MSLN",
    "MESOTHELIN": "MSLN",
    "CD20": "MS4A1",
    "BCMA": "TNFRSF17",
}


def _alias_map() -> dict[str, str]:
    mapping = {key.upper(): value for key, value in _EXTRA_ALIASES.items()}
    for symbol, entry in (load_configs()["taa"].get("TAAs") or {}).items():
        mapping[symbol.upper()] = symbol
        for alias in entry.get("aliases") or []:
            mapping[str(alias).upper()] = symbol
    return mapping


def canonical_symbol(query: str) -> str:
    text = (query or "").strip()
    if not text:
        raise ValueError("empty gene query")
    token = re.sub(r"[^A-Za-z0-9-]", "", text).upper()
    if not token:
        raise ValueError("gene symbol not found")
    mapped = _alias_map().get(token, token)
    columns = set(load_expression().columns)
    if mapped in columns:
        return mapped
    # Case-insensitive match against the matrix.
    upper_to_symbol = {column.upper(): column for column in columns}
    if mapped.upper() in upper_to_symbol:
        return upper_to_symbol[mapped.upper()]
    if token in upper_to_symbol:
        return upper_to_symbol[token]
    raise ValueError(f"{text} is not in the DepMap expression matrix")


def gene_identity(query: str) -> dict:
    symbol = canonical_symbol(query)
    ids = load_gene_ids()
    ensembl = ""
    if symbol in ids.index:
        ensembl = str(ids.loc[symbol])
    curated = (load_configs()["taa"].get("TAAs") or {}).get(symbol)
    return {
        "query": query.strip(),
        "gene_symbol": symbol,
        "gene_id": ensembl,
        "curated": curated is not None,
        "curated_entry": curated or {},
    }


def taa_entry(symbol: str) -> dict:
    return (load_configs()["taa"].get("TAAs") or {}).get(symbol) or {}
