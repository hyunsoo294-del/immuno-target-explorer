"""Resolve a TAA query to one DepMap gene symbol."""

from __future__ import annotations

from immunoscore.src.config_loader import load_configs
from immunoscore.src.data.depmap import load_gene_ids
from immunoscore.src.data.gene_resolver import Resolution, resolve

_ROUTE_LABELS = {
    "exact": "정확 일치",
    "curated": "별칭 해석(큐레이션)",
    "hgnc": "별칭 해석(HGNC)",
}


def canonical_symbol(query: str) -> str:
    result = resolve(query)
    if result.refused:
        raise ValueError(result.warning or f"{query} is not a gene product this model can score")
    if not result.symbol:
        raise ValueError(result.failure_message())
    return result.symbol


def gene_identity(query: str) -> dict:
    result = resolve(query)
    if result.refused:
        return _identity(query, result, "", "")
    if not result.symbol:
        raise ValueError(result.failure_message())
    ids = load_gene_ids()
    ensembl = str(ids.loc[result.symbol]) if result.symbol in ids.index else ""
    return _identity(query, result, result.symbol, ensembl)


def _identity(query: str, result: Resolution, symbol: str, ensembl: str) -> dict:
    curated = (load_configs()["taa"].get("TAAs") or {}).get(symbol) if symbol else None
    return {
        "query": query.strip(),
        "gene_symbol": symbol,
        "gene_id": ensembl,
        "curated": curated is not None,
        "curated_entry": curated or {},
        "route": result.route,
        "route_label": _ROUTE_LABELS.get(result.route, ""),
        "antigen_class": result.antigen_class,
        "warning": result.warning,
        "biosynthesis_genes": list(result.biosynthesis_genes),
        "refused": result.refused,
        "suggestions": list(result.suggestions),
        "flags": list(result.flags),
        "confidence_override": result.confidence_override,
    }


def taa_entry(symbol: str) -> dict:
    return (load_configs()["taa"].get("TAAs") or {}).get(symbol) or {}
