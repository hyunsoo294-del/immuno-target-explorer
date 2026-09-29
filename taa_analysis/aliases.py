# -*- coding: utf-8 -*-
"""Gene aliases. HER2 and ERBB2 share one target id and one result path."""

from __future__ import annotations

import re

from data_sources import GENE_ALIASES

ERBB2_ENSEMBL = "ENSG00000141736"


def canonical_symbol(query: str) -> str:
    """Return the approved symbol used for retrieval. HER2 -> ERBB2."""
    text = (query or "").strip()
    if not text:
        raise ValueError("empty gene query")
    ensembl = re.search(r"ENSG\d{11}", text, flags=re.I)
    if ensembl and ensembl.group(0).upper() == ERBB2_ENSEMBL:
        return "ERBB2"
    compact = re.sub(r"\s+", " ", text.upper())
    for alias in sorted(GENE_ALIASES, key=len, reverse=True):
        pattern = r"(?<![A-Z0-9])" + re.escape(alias.upper()) + r"(?![A-Z0-9])"
        if re.search(pattern, compact):
            return GENE_ALIASES[alias]
    token = re.search(r"\b[A-Z][A-Z0-9-]{1,14}\b", text.upper())
    if not token:
        raise ValueError("gene symbol not found")
    raw = token.group(0)
    return GENE_ALIASES.get(raw, raw)


def target_identity(query: str) -> dict[str, str]:
    symbol = canonical_symbol(query)
    identity = {
        "query": query.strip(),
        "gene_symbol": symbol,
        "species": "Homo sapiens",
    }
    if symbol == "ERBB2":
        identity["gene_id"] = ERBB2_ENSEMBL
    return identity
