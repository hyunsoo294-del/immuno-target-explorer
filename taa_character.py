# -*- coding: utf-8 -*-
"""TAA character for bispecific / ADC fit: localization, size, internalization literature."""

from __future__ import annotations

import re
from typing import Any

import pandas as pd
import requests

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
EUROPEPMC_SEARCH = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
PUBMED_SEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
PUBMED_SUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
PUBMED_FETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

MEMBRANE_TERMS = (
    "cell membrane",
    "plasma membrane",
    "cell surface",
    "membrane",
    "single-pass",
    "multi-pass",
    "transmembrane",
    "secreted",
    "extracellular",
)
CYTOSOL_TERMS = ("cytoplasm", "cytosol", "cytoplasmic")
NUCLEAR_TERMS = ("nucleus", "nucleoplasm", "nuclear")
SECRETED_TERMS = ("secreted", "extracellular space", "extracellular region")

INTERNALIZE_POS = (
    "internaliz",
    "endocytos",
    "lysosom",
    "antibody-drug",
    "adc ",
    "rapid uptake",
    "efficient uptake",
)
DOMAIN_DEP = (
    "epitope-dependent",
    "domain-dependent",
    "epitope specific",
    "does not internalize",
    "poor internaliz",
    "weak internaliz",
    "recycling",
    "membrane persist",
    "surface persist",
    "not sufficient to induce",
)


def _get_json(session: requests.Session, url: str, params: dict[str, Any] | None = None) -> Any:
    from data_sources import _request

    response = _request(session, "GET", url, params=params or {})
    ctype = (response.headers.get("Content-Type") or "").lower()
    if "json" in ctype or response.text.lstrip().startswith("{") or response.text.lstrip().startswith("["):
        return response.json()
    return response.text


def _text_blob(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, str):
        return value
    if isinstance(value, dict):
        parts = []
        for key in ("value", "id", "name", "evidence"):
            if value.get(key):
                parts.append(_text_blob(value[key]))
        for item in value.get("texts") or value.get("locations") or []:
            parts.append(_text_blob(item))
        return " ".join(p for p in parts if p)
    if isinstance(value, list):
        return " ".join(_text_blob(item) for item in value)
    return str(value)


def _feature_span(feature: dict[str, Any]) -> tuple[int | None, int | None]:
    loc = feature.get("location") or {}
    start = (loc.get("start") or {}).get("value")
    end = (loc.get("end") or {}).get("value")
    try:
        start_i = int(start) if start is not None else None
    except (TypeError, ValueError):
        start_i = None
    try:
        end_i = int(end) if end is not None else None
    except (TypeError, ValueError):
        end_i = None
    return start_i, end_i


def fetch_uniprot(session: requests.Session, symbol: str, ensembl_id: str) -> dict[str, Any]:
    query = f"(gene_exact:{symbol}) AND (organism_id:9606) AND (reviewed:true)"
    try:
        payload = _get_json(
            session,
            UNIPROT_SEARCH,
            {
                "query": query,
                "format": "json",
                "size": 1,
                "fields": "accession,id,protein_name,gene_names,length,mass,cc_subcellular_location,cc_function,ft_signal,ft_transmem,ft_topo_dom,ft_domain,keyword,xref_ensembl",
            },
        )
    except Exception:
        payload = {}
    results = (payload or {}).get("results") or []
    if not results and ensembl_id:
        try:
            payload = _get_json(
                session,
                UNIPROT_SEARCH,
                {
                    "query": f"(xref:ensembl-{ensembl_id}) AND (organism_id:9606)",
                    "format": "json",
                    "size": 1,
                },
            )
            results = (payload or {}).get("results") or []
        except Exception:
            results = []
    if not results:
        return {}
    rec = results[0]
    accession = rec.get("primaryAccession") or ""
    rec_id = rec.get("uniProtkbId") or ""
    desc = rec.get("proteinDescription") or {}
    rec_name = (((desc.get("recommendedName") or {}).get("fullName") or {}).get("value")) or ""
    seq = rec.get("sequence") or {}
    length = seq.get("length")
    mass = seq.get("molWeight")
    mass_kda = round(float(mass) / 1000.0, 1) if mass else None

    locations: list[str] = []
    function_bits: list[str] = []
    for comment in rec.get("comments") or []:
        ctype = (comment.get("commentType") or "").upper()
        if ctype == "SUBCELLULAR LOCATION":
            locations.append(_text_blob(comment))
        elif ctype == "FUNCTION":
            function_bits.append(_text_blob(comment))

    keywords = [str((kw.get("name") or "")).strip() for kw in (rec.get("keywords") or []) if kw.get("name")]
    features: list[dict[str, Any]] = []
    tm_count = 0
    ecd_spans: list[str] = []
    cyto_spans: list[str] = []
    domains: list[str] = []
    signal = False
    for feat in rec.get("features") or []:
        ftype = str(feat.get("type") or "")
        start, end = _feature_span(feat)
        desc_f = str(feat.get("description") or "")
        span = f"{start}-{end}" if start and end else ""
        label = " ".join(p for p in (ftype, desc_f, span) if p).strip()
        if label:
            features.append({"type": ftype, "description": desc_f, "start": start, "end": end})
        low = (ftype + " " + desc_f).lower()
        if "transmem" in low:
            tm_count += 1
        if "signal" in low:
            signal = True
        if "extracellular" in low or "extracell" in low:
            ecd_spans.append(span or desc_f)
        if "cytoplasmic" in low:
            cyto_spans.append(span or desc_f)
        if ftype.lower() == "domain" and desc_f:
            domains.append(f"{desc_f} ({span})" if span else desc_f)

    loc_text = "; ".join(dict.fromkeys(x.strip() for x in locations if x.strip()))
    loc_low = loc_text.lower() + " " + " ".join(keywords).lower()
    if any(t in loc_low for t in ("plasma membrane", "cell membrane", "cell surface")) or tm_count:
        loc_class = "Cell surface"
    elif any(t in loc_low for t in SECRETED_TERMS) and tm_count == 0:
        loc_class = "Secreted / extracellular"
    elif any(t in loc_low for t in NUCLEAR_TERMS) and not tm_count:
        loc_class = "Nuclear"
    elif any(t in loc_low for t in CYTOSOL_TERMS) and not tm_count:
        loc_class = "Cytosol"
    elif "membrane" in loc_low:
        loc_class = "Membrane-associated"
    else:
        loc_class = "Unknown"

    return {
        "accession": accession,
        "uniprot_id": rec_id,
        "protein_name": rec_name,
        "length_aa": length,
        "mass_kda": mass_kda,
        "location_text": loc_text,
        "location_class": loc_class,
        "function": " ".join(function_bits)[:900],
        "keywords": keywords[:18],
        "tm_count": tm_count,
        "has_signal": signal,
        "ecd_spans": ecd_spans[:8],
        "cyto_spans": cyto_spans[:8],
        "domains": domains[:12],
        "features": features[:40],
        "uniprot_url": f"https://www.uniprot.org/uniprotkb/{accession}/entry" if accession else "",
    }


def _score_abstract(text: str) -> dict[str, Any]:
    blob = " " + re.sub(r"\s+", " ", (text or "").lower()) + " "
    pos_hits = [k.strip() for k in INTERNALIZE_POS if k in blob]
    dep_hits = [k.strip() for k in DOMAIN_DEP if k in blob]
    if dep_hits and not pos_hits:
        label = "Low / domain-restricted"
    elif dep_hits and pos_hits:
        label = "Domain-dependent"
    elif pos_hits:
        label = "Internalizes (literature)"
    else:
        label = "Mention only"
    return {"label": label, "pos_hits": pos_hits, "dep_hits": dep_hits}


def _paper_url(hit: dict[str, Any], pmid: str) -> str:
    if pmid:
        return f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
    doi = str(hit.get("doi") or "").strip()
    if doi:
        return f"https://doi.org/{doi}"
    urls = hit.get("fullTextUrlList")
    items = []
    if isinstance(urls, dict):
        raw = urls.get("fullTextUrl") or []
        items = raw if isinstance(raw, list) else [raw]
    elif isinstance(urls, list):
        items = urls
    for item in items:
        if isinstance(item, dict) and item.get("url"):
            return str(item["url"])
        if isinstance(item, str) and item.startswith("http"):
            return item
    return ""


def fetch_papers(session: requests.Session, symbol: str, protein_name: str) -> list[dict[str, Any]]:
    name_q = protein_name.split("(")[0].strip() if protein_name else ""
    query = (
        f'(TITLE_ABS:"{symbol}"'
        + (f' OR TITLE_ABS:"{name_q}"' if name_q and name_q.lower() != symbol.lower() else "")
        + ") AND (internalization OR endocytosis OR internalisation OR ADC OR \"antibody-drug\") "
        "AND (antibody OR bispecific OR engager)"
    )
    papers: list[dict[str, Any]] = []
    try:
        payload = _get_json(
            session,
            EUROPEPMC_SEARCH,
            {
                "query": query,
                "format": "json",
                "pageSize": 12,
                "resultType": "core",
                "sort": "CITED desc",
            },
        )
        hits = ((payload or {}).get("resultList") or {}).get("result") or []
        for hit in hits:
            abstract = str(hit.get("abstractText") or "")
            title = str(hit.get("title") or "").strip()
            pmid = str(hit.get("pmid") or "").strip()
            if not title:
                continue
            score = _score_abstract(title + " " + abstract)
            snippet = re.sub(r"<[^>]+>", "", abstract)
            snippet = re.sub(r"\s+", " ", snippet).strip()
            if len(snippet) > 420:
                snippet = snippet[:417] + "..."
            year = str(hit.get("pubYear") or "")
            journal = str(hit.get("journalTitle") or hit.get("infoSource") or "")
            papers.append(
                {
                    "title": title,
                    "year": year,
                    "journal": journal,
                    "pmid": pmid,
                    "url": _paper_url(hit, pmid),
                    "snippet": snippet,
                    "internalization_call": score["label"],
                    "source": "Europe PMC / PubMed",
                }
            )
    except Exception:
        papers = []

    if len(papers) >= 4:
        return papers[:8]

    # PubMed fallback when Europe PMC is thin.
    term = (
        f'({symbol}[Title/Abstract]) AND (internalization OR endocytosis OR "antibody-drug conjugate") '
        "AND (antibody OR bispecific)"
    )
    try:
        search = _get_json(
            session,
            PUBMED_SEARCH,
            {
                "db": "pubmed",
                "retmode": "json",
                "retmax": 8,
                "term": term,
                "tool": "ImmunoTargetExplorer",
                "email": "research@localhost",
            },
        )
        ids = (((search or {}).get("esearchresult") or {}).get("idlist") or [])
        if not ids:
            return papers[:8]
        summary = _get_json(
            session,
            PUBMED_SUMMARY,
            {
                "db": "pubmed",
                "retmode": "json",
                "id": ",".join(ids),
                "tool": "ImmunoTargetExplorer",
                "email": "research@localhost",
            },
        )
        result = (summary or {}).get("result") or {}
        for pmid in ids:
            rec = result.get(pmid) or {}
            title = str(rec.get("title") or "").strip()
            if not title:
                continue
            score = _score_abstract(title)
            papers.append(
                {
                    "title": title,
                    "year": str((rec.get("pubdate") or "")[:4]),
                    "journal": str(rec.get("fulljournalname") or rec.get("source") or ""),
                    "pmid": pmid,
                    "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                    "snippet": "",
                    "internalization_call": score["label"],
                    "source": "PubMed",
                }
            )
    except Exception:
        pass
    return papers[:8]


def _consensus_internalization(uniprot: dict[str, Any], papers: list[dict[str, Any]]) -> dict[str, str]:
    labels = [p.get("internalization_call") or "" for p in papers]
    loc = uniprot.get("location_class") or "Unknown"
    if loc in {"Cytosol", "Nuclear"}:
        fit = "Poor classical IgG TAA. Intracellular antigen is not freely accessible to intact antibodies."
        call = "Not a typical surface TAA"
    elif loc == "Secreted / extracellular":
        fit = "Soluble/secreted protein. Antibody can bind, but tumor-cell display and internalization need extra evidence."
        call = "Secreted - confirm membrane form"
    elif any("Domain-dependent" in x for x in labels) or any("Low" in x for x in labels):
        call = "Domain-dependent (literature)"
        fit = (
            "Binding an ECD epitope can change internalization. For TAA x immune-cell engagers, pick an epitope that keeps "
            "enough surface antigen. For ADC, prefer epitopes that drive lysosomal uptake."
        )
    elif any(x.startswith("Internalizes") for x in labels):
        call = "Internalizes after antibody binding"
        fit = (
            "Literature reports antibody-driven uptake. That can help ADC payload delivery, but may lower surface density "
            "for T-cell engagers if internalization is fast."
        )
    elif loc == "Cell surface" or uniprot.get("tm_count"):
        call = "Surface antigen; internalization not settled"
        fit = (
            "Membrane topology supports antibody access. Internalization after ECD engagement is not clearly settled in "
            "the retrieved papers - check epitope and domain before locking a binder."
        )
    else:
        call = "Insufficient literature"
        fit = "Need more papers or experimental internalization assays (flow, pHrodo, ADC kill)."
    return {"call": call, "fit": fit}


def collect_taa_character(
    session: requests.Session,
    symbol: str,
    ensembl_id: str,
    gene_name: str = "",
    hpa_summary: dict[str, Any] | None = None,
) -> dict[str, Any]:
    uniprot = fetch_uniprot(session, symbol, ensembl_id)
    protein_name = uniprot.get("protein_name") or gene_name or symbol
    papers = fetch_papers(session, symbol, protein_name)
    consensus = _consensus_internalization(uniprot, papers)
    hpa_loc = ""
    if isinstance(hpa_summary, dict):
        hpa_loc = str(hpa_summary.get("Subcellular location") or "")
        if isinstance(hpa_summary.get("Subcellular location"), (list, dict)):
            hpa_loc = _text_blob(hpa_summary.get("Subcellular location"))
    return {
        "symbol": symbol,
        "ensembl_id": ensembl_id,
        "uniprot": uniprot,
        "hpa_location": hpa_loc,
        "papers": papers,
        "internalization_call": consensus["call"],
        "modality_note": consensus["fit"],
        "paper_count": len(papers),
    }
