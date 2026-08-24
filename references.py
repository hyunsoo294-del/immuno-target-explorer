# -*- coding: utf-8 -*-
"""Numbered scientific citations for chat answers and PowerPoint slides."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Any

import pandas as pd


@dataclass(frozen=True)
class Citation:
    number: int
    key: str
    text: str


def _as_df(value) -> pd.DataFrame:
    if isinstance(value, pd.DataFrame):
        return value
    return pd.DataFrame()


def _identity(bundle: dict[str, Any]) -> tuple[str, str, str]:
    identity = bundle.get("identity") or {}
    symbol = str(identity.get("symbol") or "NA")
    ensembl = str(identity.get("ensembl_id") or "NA")
    retrieved = bundle.get("retrieved_at") or datetime.now(timezone.utc).strftime("%Y-%m-%d")
    return symbol, ensembl, str(retrieved)[:10]


def _used_sources(bundle: dict[str, Any]) -> set[str]:
    comparison = _as_df(bundle.get("comparison_df"))
    tcga = _as_df(bundle.get("tcga_df"))
    tcga_summary = _as_df(bundle.get("tcga_summary_df"))
    ihc = _as_df(bundle.get("analytics_df"))
    cptac = _as_df(bundle.get("cptac_df"))
    cptac_summary = _as_df(bundle.get("cptac_summary_df"))
    expr = _as_df(bundle.get("expression_df"))

    used: set[str] = {"opentargets"}
    if not comparison.empty:
        if comparison.get("normal_rna_hpa_ntpm") is not None and comparison["normal_rna_hpa_ntpm"].notna().any():
            used.add("hpa")
        if comparison.get("normal_rna_gtex_ntpm") is not None and comparison["normal_rna_gtex_ntpm"].notna().any():
            used.add("gtex")
        if comparison.get("cancer_rna_tcga_median") is not None and comparison["cancer_rna_tcga_median"].notna().any():
            used.add("tcga")
            used.add("cbioportal")
        if comparison.get("protein_ihc_cancer_0_3") is not None and comparison["protein_ihc_cancer_0_3"].notna().any():
            used.add("hpa")
        if comparison.get("protein_ihc_normal_0_3") is not None and comparison["protein_ihc_normal_0_3"].notna().any():
            used.add("hpa")
    if not tcga.empty or not tcga_summary.empty:
        used.add("tcga")
        used.add("cbioportal")
    if not ihc.empty:
        used.add("hpa")
    if not cptac.empty or not cptac_summary.empty:
        used.add("cptac")
        used.add("cbioportal")
        used.add("pdc")
    if not expr.empty and "source" in expr.columns:
        sources = expr["source"].astype(str).str.lower()
        if sources.str.contains("gtex", na=False).any():
            used.add("gtex")
        if sources.str.contains("hpa|protein atlas", na=False).any():
            used.add("hpa")
    return used


def _catalog(symbol: str, ensembl: str, retrieved: str) -> dict[str, str]:
    hpa_url = f"https://www.proteinatlas.org/{ensembl}"
    gtex_url = f"https://gtexportal.org/home/gene/{ensembl}"
    ot_url = f"https://platform.opentargets.org/target/{ensembl}"
    ensembl_url = f"https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={ensembl}"
    return {
        "hpa": (
            f"Uhlen M, Fagerberg L, Hallstrom BM, et al. Tissue-based map of the human proteome. "
            f"Science. 2015;347(6220):1260419. doi:10.1126/science.1260419. "
            f"Human Protein Atlas (HPA), gene {symbol} ({ensembl}). "
            f"Accessed {retrieved}. Available from: {hpa_url}"
        ),
        "gtex": (
            f"GTEx Consortium. The GTEx Consortium atlas of genetic regulatory effects across human tissues. "
            f"Science. 2020;369(6509):1318-1330. doi:10.1126/science.aaz1776. "
            f"GTEx Portal v8, gene {symbol} ({ensembl}). "
            f"Accessed {retrieved}. Available from: {gtex_url}"
        ),
        "tcga": (
            f"Weinstein JN, Collisson EA, Mills GB, et al. The Cancer Genome Atlas Pan-Cancer analysis project. "
            f"Nat Genet. 2013;45(10):1113-1120. doi:10.1038/ng.2764. "
            f"TCGA Pan-Cancer Atlas RNA-seq for {symbol} ({ensembl}), retrieved via cBioPortal. "
            f"Accessed {retrieved}. Available from: https://www.cbioportal.org/"
        ),
        "cbioportal": (
            "Cerami E, Gao J, Dogrusoz U, et al. The cBio cancer genomics portal: an open platform for exploring "
            "multidimensional cancer genomics data. Cancer Discov. 2012;2(5):401-404. "
            "doi:10.1158/2159-8290.CD-12-0095. Gao J, Aksoy BA, Dogrusoz U, et al. Integrative analysis of complex "
            "cancer genomics and clinical profiles using the cBioPortal. Sci Signal. 2013;6(269):pl1. "
            f"Accessed {retrieved}. Available from: https://www.cbioportal.org/"
        ),
        "opentargets": (
            f"Ochoa D, Hercules A, Carmona M, et al. The next-generation Open Targets Platform: reimagined, "
            f"redesigned, rebuilt. Nucleic Acids Res. 2023;51(D1):D1353-D1359. doi:10.1093/nar/gkac1046. "
            f"Target {symbol} (Ensembl {ensembl}). Accessed {retrieved}. Available from: {ot_url}"
        ),
        "cptac": (
            f"Edwards NJ, Oberti M, Thangudu RR, et al. The CPTAC Data Portal: a resource for cancer proteomics "
            f"research. J Proteome Res. 2015;14(6):2707-2713. doi:10.1021/pr501254j. "
            f"CPTAC mass-spectrometry protein abundance for {symbol} ({ensembl}). "
            f"Accessed {retrieved}. Available from: https://proteomics.cancer.gov/programs/cptac"
        ),
        "pdc": (
            f"Thangudu RR, Rudnick PA, Holck M, et al. Proteomic Data Commons: a resource for proteogenomic analysis. "
            f"NCI Proteomic Data Commons (PDC). Gene {symbol} ({ensembl}). "
            f"Accessed {retrieved}. Available from: https://pdc.cancer.gov/"
        ),
        "ensembl": (
            f"Martin FJ, Amode MR, Aneja A, et al. Ensembl 2023. Nucleic Acids Res. 2023;51(D1):D933-D941. "
            f"doi:10.1093/nar/gkac958. Gene identifier {ensembl} ({symbol}). "
            f"Accessed {retrieved}. Available from: {ensembl_url}"
        ),
    }


SOURCE_ORDER = ("ensembl", "hpa", "gtex", "tcga", "cbioportal", "opentargets", "cptac", "pdc")


class ReferenceSet:
    def __init__(self, citations: list[Citation]):
        self.citations = citations
        self._by_key = {c.key: c for c in citations}

    def mark(self, *keys: str) -> str:
        nums = []
        seen = set()
        for key in keys:
            cite = self._by_key.get(key)
            if cite and cite.number not in seen:
                nums.append(cite.number)
                seen.add(cite.number)
        if not nums:
            return ""
        return "[" + ",".join(str(n) for n in nums) + "]"

    def markdown(self) -> str:
        if not self.citations:
            return ""
        lines = ["### References"]
        for cite in self.citations:
            lines.append(f"{cite.number}. {cite.text}")
        return "\n\n".join(lines)

    def numbered_texts(self) -> list[str]:
        return [f"{c.number}. {c.text}" for c in self.citations]


def build_references(bundle: dict[str, Any] | None) -> ReferenceSet:
    if not bundle:
        return ReferenceSet([])
    symbol, ensembl, retrieved = _identity(bundle)
    used = _used_sources(bundle)
    used.add("ensembl")
    catalog = _catalog(symbol, ensembl, retrieved)
    citations: list[Citation] = []
    number = 1
    for key in SOURCE_ORDER:
        if key in used and key in catalog:
            citations.append(Citation(number=number, key=key, text=catalog[key]))
            number += 1
    return ReferenceSet(citations)
