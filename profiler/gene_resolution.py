# -*- coding: utf-8 -*-
"""TAA / immune-marker alias resolution with official gene confirmation."""

from __future__ import annotations

import re
from dataclasses import asdict, dataclass, field
from typing import Any

from profiler.http_client import NotFoundError, request_json

CBIO_GENE = "https://www.cbioportal.org/api/genes/{symbol}"
CBIO_GENES_FETCH = "https://www.cbioportal.org/api/genes/fetch"

# Official Hugo symbols for the 17 priority TAAs + immune markers.
# Failed aliases from v13 handoff are included intentionally.
TAA_ALIASES: dict[str, str] = {
    # HER2
    "HER2": "ERBB2",
    "HER-2": "ERBB2",
    "HER 2": "ERBB2",
    "ERBB2": "ERBB2",
    "NEU": "ERBB2",
    # EGFR
    "EGFR": "EGFR",
    "HER1": "EGFR",
    "ERBB1": "EGFR",
    # PSMA
    "PSMA": "FOLH1",
    "FOLH1": "FOLH1",
    "GCPII": "FOLH1",
    # STEAP1
    "STEAP1": "STEAP1",
    "STEAP-1": "STEAP1",
    # Mesothelin
    "MSLN": "MSLN",
    "MESOTHELIN": "MSLN",
    "MESOTHELIN ": "MSLN",
    # B7-H3
    "B7-H3": "CD276",
    "B7H3": "CD276",
    "B7 H3": "CD276",
    "CD276": "CD276",
    # TROP2
    "TROP2": "TACSTD2",
    "TROP-2": "TACSTD2",
    "TROP 2": "TACSTD2",
    "TACSTD2": "TACSTD2",
    "EGP1": "TACSTD2",
    # Nectin-4
    "NECTIN-4": "NECTIN4",
    "NECTIN4": "NECTIN4",
    "NECTIN 4": "NECTIN4",
    "PVRL4": "NECTIN4",
    "PVRL-4": "NECTIN4",
    # CD19
    "CD19": "CD19",
    # CD20
    "CD20": "MS4A1",
    "MS4A1": "MS4A1",
    # BCMA
    "BCMA": "TNFRSF17",
    "CD269": "TNFRSF17",
    "TNFRSF17": "TNFRSF17",
    # GPRC5D
    "GPRC5D": "GPRC5D",
    # DLL3
    "DLL3": "DLL3",
    # CEA
    "CEA": "CEACAM5",
    "CD66E": "CEACAM5",
    "CD66e": "CEACAM5",
    "CEACAM5": "CEACAM5",
    # FAP
    "FAP": "FAP",
    "FAPα": "FAP",
    "FAPALPHA": "FAP",
    # CLDN18 (isoform limit)
    "CLDN18": "CLDN18",
    "CLAUDIN18": "CLDN18",
    "CLAUDIN-18": "CLDN18",
    "CLDN18.2": "CLDN18",
    "CLAUDIN18.2": "CLDN18",
    "CLAUDIN-18.2": "CLDN18",
    "CLDN18.1": "CLDN18",
}

IMMUNE_ALIASES: dict[str, str] = {
    "PD-1": "PDCD1",
    "PD1": "PDCD1",
    "PD 1": "PDCD1",
    "PDCD1": "PDCD1",
    "4-1BB": "TNFRSF9",
    "41BB": "TNFRSF9",
    "CD137": "TNFRSF9",
    "TNFRSF9": "TNFRSF9",
    "PD-L1": "CD274",
    "PDL1": "CD274",
    "CD274": "CD274",
    "CD40": "CD40",
    "CLEC9A": "CLEC9A",
    "CD3E": "CD3E",
    "CD8A": "CD8A",
    "LAG3": "LAG3",
    "LAG-3": "LAG3",
    "TIGIT": "TIGIT",
    "ICOS": "ICOS",
    "CD80": "CD80",
    "CD86": "CD86",
    "B2M": "B2M",
    "PTPRC": "PTPRC",
    "CD45": "PTPRC",
}

# MHCII is not a single gene — callers must pick one locus.
MHCII_LOCI = ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1"]

SPECIAL_ANTIGENS: dict[str, dict[str, Any]] = {
    "GD2": {
        "kind": "glycolipid",
        "mappable_gene": None,
        "proxy_enzymes": ["B4GALNT1", "ST8SIA1"],
        "message": (
            "GD2 is a glycolipid antigen, not a single gene product. "
            "RNA/proteomics cannot directly quantify GD2. "
            "B4GALNT1/ST8SIA1 are biosynthetic enzyme proxies only, not GD2 itself."
        ),
    },
    "CLDN18.2": {
        "kind": "isoform_limited",
        "mappable_gene": "CLDN18",
        "message": (
            "CLDN18.2 maps to gene CLDN18 only. Gene-collapsed RNA/protein cannot "
            "distinguish CLDN18.2 from CLDN18.1. Do not treat results as isoform-positive."
        ),
    },
    "CLAUDIN18.2": {
        "kind": "isoform_limited",
        "mappable_gene": "CLDN18",
        "message": (
            "CLDN18.2 maps to gene CLDN18 only. Gene-collapsed RNA/protein cannot "
            "distinguish CLDN18.2 from CLDN18.1. Do not treat results as isoform-positive."
        ),
    },
}

PRIORITY_TAAS = [
    {"input": "HER2", "symbol": "ERBB2", "category": "approved_adc_context", "therapy_types": ["ADC", "mAb", "CAR-T"]},
    {"input": "EGFR", "symbol": "EGFR", "category": "approved_mab_context", "therapy_types": ["mAb", "ADC"]},
    {"input": "PSMA", "symbol": "FOLH1", "category": "approved_radioligand_context", "therapy_types": ["radioligand", "CAR-T", "TCE"]},
    {"input": "STEAP1", "symbol": "STEAP1", "category": "clinical_development", "therapy_types": ["TCE", "ADC"]},
    {"input": "MSLN", "symbol": "MSLN", "category": "clinical_development", "therapy_types": ["CAR-T", "ADC", "TCE"]},
    {"input": "B7-H3", "symbol": "CD276", "category": "clinical_development", "therapy_types": ["ADC", "CAR-T"]},
    {"input": "TROP2", "symbol": "TACSTD2", "category": "approved_adc_context", "therapy_types": ["ADC"]},
    {"input": "Nectin-4", "symbol": "NECTIN4", "category": "approved_adc_context", "therapy_types": ["ADC"]},
    {"input": "CD19", "symbol": "CD19", "category": "approved_cell_therapy", "therapy_types": ["CAR-T", "TCE", "mAb"]},
    {"input": "CD20", "symbol": "MS4A1", "category": "approved_mab_context", "therapy_types": ["mAb", "CAR-T", "TCE"]},
    {"input": "BCMA", "symbol": "TNFRSF17", "category": "approved_cell_therapy", "therapy_types": ["CAR-T", "TCE", "ADC"]},
    {"input": "GPRC5D", "symbol": "GPRC5D", "category": "approved_cell_therapy", "therapy_types": ["CAR-T", "TCE"]},
    {"input": "DLL3", "symbol": "DLL3", "category": "clinical_development", "therapy_types": ["TCE", "ADC", "CAR-T"]},
    {"input": "CEA", "symbol": "CEACAM5", "category": "clinical_development", "therapy_types": ["TCE", "ADC", "CAR-T"]},
    {"input": "FAP", "symbol": "FAP", "category": "research_stromal", "therapy_types": ["radioligand", "CAR-T", "TCE"]},
    {"input": "CLDN18.2", "symbol": "CLDN18", "category": "approved_mab_context", "therapy_types": ["mAb", "ADC", "CAR-T"]},
    {"input": "GD2", "symbol": None, "category": "approved_cell_therapy_glycolipid", "therapy_types": ["mAb", "CAR-T"]},
]

TAA_PANEL_FOR_RANKING = [
    "ERBB2", "EGFR", "FOLH1", "STEAP1", "MSLN", "CD276", "TACSTD2", "NECTIN4",
    "CD19", "MS4A1", "TNFRSF17", "GPRC5D", "DLL3", "CEACAM5", "FAP", "CLDN18",
    "MET", "AXL", "CDH3", "GPC3", "ROR1", "ROR2", "CD70", "CD33", "CD123",
    "IL3RA", "SLAMF7", "FCRL5", "LY6E", "CDH6", "ENPP3", "TPBG", "MUC16",
    "MUC1", "EPCAM", "KDR", "FLT1", "PDGFRA", "KIT", "RET", "ALK", "ROS1",
]


@dataclass
class ResolvedTarget:
    query: str
    normalized_alias: str | None
    official_symbol: str | None
    entrez_gene_id: int | None
    mapping_source: str
    special_kind: str | None = None
    warnings: list[str] = field(default_factory=list)
    unresolved: bool = False
    raw_provider: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _normalize_key(text: str) -> str:
    raw = (text or "").strip()
    # Preserve HLA-DRA style hyphens in official symbols; only collapse spaces.
    compact = re.sub(r"\s+", " ", raw)
    return compact.upper()


def alias_lookup(query: str) -> tuple[str | None, str | None, list[str]]:
    """Return (official_symbol_or_None, special_kind, warnings)."""
    q = _normalize_key(query)
    warnings: list[str] = []

    # Exact special antigen first (GD2)
    for key, meta in SPECIAL_ANTIGENS.items():
        if q == key.upper() or q.replace(" ", "") == key.upper().replace(" ", ""):
            if meta["kind"] == "glycolipid":
                return None, "glycolipid", [meta["message"]]
            if meta["kind"] == "isoform_limited":
                return meta["mappable_gene"], "isoform_limited", [meta["message"]]

    if q in ("MHCII", "MHC II", "MHC-II", "HLA-DR"):
        return None, "multi_gene", [
            "MHCII is not a single gene. Choose one of: " + ", ".join(MHCII_LOCI)
        ]

    # Alias tables (case-insensitive keys stored upper)
    alias_map = {k.upper(): v for k, v in {**TAA_ALIASES, **IMMUNE_ALIASES}.items()}
    # Also try without spaces / with hyphen variants
    candidates = [q, q.replace(" ", "-"), q.replace("-", ""), q.replace(" ", ""), q.replace("-", " ")]
    for cand in candidates:
        if cand in alias_map:
            symbol = alias_map[cand]
            if q in ("CLDN18.2", "CLAUDIN18.2", "CLAUDIN-18.2") or "18.2" in q:
                warnings.append(SPECIAL_ANTIGENS["CLDN18.2"]["message"])
                return symbol, "isoform_limited", warnings
            if symbol == "FAP":
                warnings.append(
                    "FAP is often stromal/CAF-associated. Cancer cell-line lists alone "
                    "do not establish CAF expression or stromal fraction."
                )
            return symbol, None, warnings

    # If query already looks like a gene symbol, pass through without silent fuzzy match
    if re.fullmatch(r"[A-Z][A-Z0-9-]{0,14}", q.replace(" ", "")):
        return q.replace(" ", ""), None, warnings

    return None, None, [f"Unrecognized target name: {query}. Will not silently map to a similar gene."]


def fetch_gene_by_symbol(symbol: str) -> dict[str, Any]:
    return request_json("GET", CBIO_GENE.format(symbol=symbol), cache_ttl_s=86400)


def fetch_genes_batch(symbols: list[str]) -> list[dict[str, Any]]:
    if not symbols:
        return []
    return request_json(
        "POST",
        CBIO_GENES_FETCH,
        params={"geneIdType": "HUGO_GENE_SYMBOL"},
        json_body=symbols,
        cache_ttl_s=86400,
    )


def resolve_target(query: str, *, confirm_provider: bool = True) -> ResolvedTarget:
    q = (query or "").strip()
    if not q:
        return ResolvedTarget(
            query=q,
            normalized_alias=None,
            official_symbol=None,
            entrez_gene_id=None,
            mapping_source="empty",
            unresolved=True,
            warnings=["Empty query"],
        )

    symbol, special, warnings = alias_lookup(q)
    if special == "glycolipid" or (symbol is None and special == "multi_gene"):
        return ResolvedTarget(
            query=q,
            normalized_alias=_normalize_key(q),
            official_symbol=None,
            entrez_gene_id=None,
            mapping_source="special_antigen",
            special_kind=special,
            warnings=warnings,
            unresolved=True,
        )
    if symbol is None:
        return ResolvedTarget(
            query=q,
            normalized_alias=_normalize_key(q),
            official_symbol=None,
            entrez_gene_id=None,
            mapping_source="unresolved",
            unresolved=True,
            warnings=warnings or [f"Could not resolve: {q}"],
        )

    if not confirm_provider:
        return ResolvedTarget(
            query=q,
            normalized_alias=_normalize_key(q),
            official_symbol=symbol,
            entrez_gene_id=None,
            mapping_source="alias_table",
            special_kind=special,
            warnings=warnings,
        )

    try:
        gene = fetch_gene_by_symbol(symbol)
        return ResolvedTarget(
            query=q,
            normalized_alias=_normalize_key(q),
            official_symbol=gene.get("hugoGeneSymbol") or symbol,
            entrez_gene_id=int(gene["entrezGeneId"]),
            mapping_source="cbioportal_genes",
            special_kind=special,
            warnings=warnings,
            raw_provider=gene,
        )
    except NotFoundError:
        return ResolvedTarget(
            query=q,
            normalized_alias=_normalize_key(q),
            official_symbol=symbol,
            entrez_gene_id=None,
            mapping_source="alias_table_provider_missing",
            special_kind=special,
            warnings=warnings + [f"cBioPortal has no gene record for {symbol}"],
            unresolved=True,
        )
    except Exception as exc:
        return ResolvedTarget(
            query=q,
            normalized_alias=_normalize_key(q),
            official_symbol=symbol,
            entrez_gene_id=None,
            mapping_source="alias_table_provider_error",
            special_kind=special,
            warnings=warnings + [f"Provider confirmation failed: {exc}"],
            unresolved=False,
        )


FAILED_V13_ALIASES = [
    "Mesothelin",
    "TROP-2",
    "Nectin-4",
    "PVRL4",
    "CD20",
    "BCMA",
    "CD269",
    "CEA",
    "CD66e",
    "CLDN18.2",
    "Claudin18.2",
]
