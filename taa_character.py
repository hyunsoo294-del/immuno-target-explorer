# -*- coding: utf-8 -*-
"""TAA character for bispecific / ADC fit: localization, size, internalization literature."""

from __future__ import annotations

import re
import statistics
from typing import Any

import pandas as pd
import requests

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
EUROPEPMC_SEARCH = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
PUBMED_SEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
PUBMED_SUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
PUBMED_FETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
CTGOV_STUDIES = "https://clinicaltrials.gov/api/v2/studies"

# Published acid-wash / confocal formulas (Rajan et al., Biomolecules 2020;10:955).
INT_FORMULA_ACID = "100 - ((MFI_37C / MFI_4C) * 100)"
INT_FORMULA_CONFOCAL = "(F_in / (F_in + F_out)) * 100"

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

    synonyms: list[str] = []
    for gene in rec.get("genes") or []:
        primary = ((gene.get("geneName") or {}).get("value") or "").strip()
        if primary:
            synonyms.append(primary)
        for syn in gene.get("synonyms") or []:
            value = str(syn.get("value") or "").strip()
            if value:
                synonyms.append(value)

    return {
        "accession": accession,
        "uniprot_id": rec_id,
        "protein_name": rec_name,
        "synonyms": list(dict.fromkeys(synonyms))[:12],
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


def _clean_text(text: str) -> str:
    blob = (text or "").replace("&lt;", "<").replace("&gt;", ">").replace("&amp;", "&")
    blob = blob.replace("–", "-").replace("—", "-")
    blob = re.sub(r"<[^>]+>", " ", blob)
    return re.sub(r"\s+", " ", blob).strip()


_INT_CONTEXT = ("internaliz", "endocytos", "acid-wash", "acid wash", "phrodo")
_FALSE_PCT = (
    "%id",
    "id/g",
    "ia/g",
    "id/ml",
    "injected dose",
    "homology",
    "sequence ident",
    "overall survival",
    "progression-free",
    "response rate",
    "objective response",
    "confidence interval",
    "tumor volume",
    "cell viabil",
    "growth inhib",
    "yield",
    "purity",
    "homology",
)
_AB_CTX = (
    "antibody",
    "mab",
    "adc",
    "trastuzumab",
    "bispecific",
    "engager",
    "immunoglobulin",
    "fragment",
    "receptor",
    "antigen",
)
_EXCLUDE_CTX = (
    "protab",
    "protac",
    "knockdown",
    "sirna",
    "crispr",
    "yield",
    "purity",
)
_REMAIN_CTX = (
    "remain",
    "residual surface",
    "surface persist",
    "still on the surface",
    "surface-bound remaining",
    "not internalized",
)


def _accept_pct_window(before: str, after: str, value: float) -> float | None:
    if value < 0 or value > 100:
        return None
    window = (before + after).lower()
    if not any(key in window for key in _INT_CONTEXT):
        return None
    if any(bad in window for bad in _FALSE_PCT) or any(bad in window for bad in _EXCLUDE_CTX):
        return None
    if not any(tok in window for tok in _AB_CTX):
        return None
    if any(tok in window for tok in ("of patient", "of tumor", "of case", "monocyte", "macrophage")):
        return None
    local = (before[-36:] + after[:36]).lower()
    if any(tok in local for tok in _REMAIN_CTX):
        value = round(100.0 - value, 1)
    return round(value, 1)


def extract_internalization_pcts(text: str) -> list[float]:
    """Pull % internalization from assay wording in abstracts or full text."""
    blob = _clean_text(text)
    if not blob:
        return []
    found: list[float] = []
    used: list[tuple[int, int]] = []
    for match in re.finditer(
        r"(\d{1,3}(?:\.\d+)?)\s*[-~]\s*(\d{1,3}(?:\.\d+)?)\s*(?:%|percent)",
        blob,
        flags=re.I,
    ):
        before = blob[max(0, match.start() - 90) : match.start()]
        after = blob[match.end() : match.end() + 90]
        value = _accept_pct_window(before, after, (float(match.group(1)) + float(match.group(2))) / 2.0)
        if value is not None:
            found.append(value)
            used.append(match.span())
    for match in re.finditer(r"(\d{1,3}(?:\.\d+)?)\s*(?:%|percent)", blob, flags=re.I):
        if any(match.start() < end and match.end() > start for start, end in used):
            continue
        before = blob[max(0, match.start() - 90) : match.start()]
        after = blob[match.end() : match.end() + 90]
        value = _accept_pct_window(before, after, float(match.group(1)))
        if value is not None:
            found.append(value)
    return list(dict.fromkeys(found))[:12]


def _is_degrader_paper(title: str, abstract: str = "") -> bool:
    blob = f"{title} {abstract}".lower()
    return any(tok in blob for tok in ("protab", "protac", "lytac", "abtac"))


def _paper_weight(title: str, abstract: str, is_trial: bool = False) -> float:
    if is_trial:
        return 3.0
    blob = f"{title} {abstract}".lower()
    if any(token in blob for token in ("phase i", "phase ii", "phase iii", "clinical trial", "patients with")):
        return 3.0
    return 2.0


def _hit_to_paper(hit: dict[str, Any], source: str) -> dict[str, Any] | None:
    abstract = str(hit.get("abstractText") or "")
    title = str(hit.get("title") or "").strip()
    pmid = str(hit.get("pmid") or "").strip()
    if not title:
        return None
    score = _score_abstract(title + " " + abstract)
    snippet = _clean_text(abstract)
    pcts = [] if _is_degrader_paper(title, abstract) else extract_internalization_pcts((title.rstrip(".") + ". ") + abstract)
    rec = {
        "title": title,
        "year": str(hit.get("pubYear") or ""),
        "journal": str(hit.get("journalTitle") or hit.get("infoSource") or ""),
        "pmid": pmid,
        "pmcid": str(hit.get("pmcid") or "").strip(),
        "url": _paper_url(hit, pmid),
        "snippet": snippet[:417] + "..." if len(snippet) > 420 else snippet,
        "abstract": snippet,
        "internalization_call": score["label"],
        "source": source,
        "pct_values": pcts,
        "pct_internalization": round(statistics.median(pcts), 1) if pcts else None,
        "weight": _paper_weight(title, abstract),
    }
    return rec


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


def _europepmc_hits(session: requests.Session, query: str, page_size: int = 20) -> list[dict[str, Any]]:
    try:
        payload = _get_json(
            session,
            EUROPEPMC_SEARCH,
            {
                "query": query,
                "format": "json",
                "pageSize": page_size,
                "resultType": "core",
                "sort": "CITED desc",
            },
        )
        return ((payload or {}).get("resultList") or {}).get("result") or []
    except Exception:
        return []


def _fill_pubmed_abstracts(session: requests.Session, papers: list[dict[str, Any]]) -> None:
    missing = [p["pmid"] for p in papers if p.get("pmid") and not p.get("abstract")]
    if not missing:
        return
    try:
        response = _get_json(
            session,
            PUBMED_FETCH,
            {
                "db": "pubmed",
                "id": ",".join(missing[:16]),
                "retmode": "xml",
                "rettype": "abstract",
                "tool": "ImmunoTargetExplorer",
                "email": "research@localhost",
            },
        )
    except Exception:
        return
    xml = response if isinstance(response, str) else str(response)
    blocks = re.findall(
        r"<PubmedArticle>(.*?)</PubmedArticle>",
        xml,
        flags=re.I | re.S,
    )
    by_pmid: dict[str, str] = {}
    for block in blocks:
        pmid_m = re.search(r"<PMID[^>]*>(\d+)</PMID>", block, flags=re.I)
        texts = re.findall(r"<AbstractText[^>]*>(.*?)</AbstractText>", block, flags=re.I | re.S)
        if pmid_m and texts:
            by_pmid[pmid_m.group(1)] = _clean_text(" ".join(texts))
    for rec in papers:
        abstract = by_pmid.get(str(rec.get("pmid") or ""))
        if not abstract:
            continue
        rec["abstract"] = abstract
        rec["snippet"] = abstract[:417] + "..." if len(abstract) > 420 else abstract
        if _is_degrader_paper(rec.get("title") or "", abstract):
            rec["pct_values"] = []
            rec["pct_internalization"] = None
            continue
        pcts = extract_internalization_pcts((rec.get("title") or "") + " " + abstract)
        rec["pct_values"] = pcts
        rec["pct_internalization"] = round(statistics.median(pcts), 1) if pcts else rec.get("pct_internalization")
        rec["weight"] = _paper_weight(rec.get("title") or "", abstract)
        rec["internalization_call"] = _score_abstract((rec.get("title") or "") + " " + abstract)["label"]


def _title_abs_clause(symbol: str, protein_name: str, synonyms: list[str] | None = None) -> str:
    terms = [symbol]
    name_q = protein_name.split("(")[0].strip() if protein_name else ""
    if name_q and name_q.lower() != symbol.lower() and len(name_q.split()) <= 4:
        terms.append(name_q)
    for syn in synonyms or []:
        if syn and syn.lower() not in {t.lower() for t in terms}:
            terms.append(syn)
    quoted = " OR ".join(f'TITLE_ABS:"{t}"' for t in terms[:6])
    return f"({quoted})"


def fetch_papers(
    session: requests.Session,
    symbol: str,
    protein_name: str,
    synonyms: list[str] | None = None,
) -> list[dict[str, Any]]:
    gene_q = _title_abs_clause(symbol, protein_name, synonyms)
    queries = [
        (
            f"{gene_q} AND (internalization OR internalisation OR endocytosis) AND "
            '(antibody OR mAb OR ADC OR "antibody-drug" OR bispecific OR engager)'
        ),
        (
            f"{gene_q} AND (internalization OR internalisation OR endocytosis) AND "
            '("internalization rate" OR "% internalized" OR "percent internalized" OR pHrodo OR "acid wash")'
        ),
    ]
    papers: list[dict[str, Any]] = []
    seen: set[str] = set()
    for query in queries:
        for hit in _europepmc_hits(session, query, page_size=20):
            rec = _hit_to_paper(hit, "Europe PMC / PubMed")
            if not rec:
                continue
            key = rec["pmid"] or rec["title"].lower()
            if key in seen:
                continue
            seen.add(key)
            papers.append(rec)

    alias_terms = " OR ".join(
        f"{t}[Title/Abstract]" for t in [symbol] + [s for s in (synonyms or [])[:4] if s]
    )
    quant_term = (
        f"({alias_terms}) AND internalization[Title/Abstract] AND "
        '(percent OR pHrodo OR "acid-wash" OR "acid wash" OR internalized) AND '
        "(antibody OR ADC OR mAb OR bispecific)"
    )
    try:
        search = _get_json(
            session,
            PUBMED_SEARCH,
            {
                "db": "pubmed",
                "retmode": "json",
                "retmax": 10,
                "term": quant_term,
                "tool": "ImmunoTargetExplorer",
                "email": "research@localhost",
            },
        )
        ids = (((search or {}).get("esearchresult") or {}).get("idlist") or [])
        if ids:
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
                if pmid in seen:
                    continue
                rec = result.get(pmid) or {}
                title = str(rec.get("title") or "").strip()
                if not title:
                    continue
                score = _score_abstract(title)
                pmcid = ""
                for item in rec.get("articleids") or []:
                    if str(item.get("idtype") or "").lower() == "pmc":
                        pmcid = str(item.get("value") or "")
                        break
                papers.append(
                    {
                        "title": title,
                        "year": str((rec.get("pubdate") or "")[:4]),
                        "journal": str(rec.get("fulljournalname") or rec.get("source") or ""),
                        "pmid": pmid,
                        "pmcid": pmcid,
                        "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                        "snippet": "",
                        "abstract": "",
                        "internalization_call": score["label"],
                        "source": "PubMed quantitative",
                        "pct_values": [],
                        "pct_internalization": None,
                        "weight": _paper_weight(title, ""),
                    }
                )
                seen.add(pmid)
    except Exception:
        pass

    if len(papers) < 4:
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
                    "retmax": 12,
                    "term": term,
                    "tool": "ImmunoTargetExplorer",
                    "email": "research@localhost",
                },
            )
            ids = (((search or {}).get("esearchresult") or {}).get("idlist") or [])
            if ids:
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
                    if pmid in seen:
                        continue
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
                            "abstract": "",
                            "internalization_call": score["label"],
                            "source": "PubMed",
                            "pct_values": [],
                            "pct_internalization": None,
                            "weight": _paper_weight(title, ""),
                        }
                    )
                    seen.add(pmid)
        except Exception:
            pass

    _fill_pubmed_abstracts(session, papers)

    def _rank(rec: dict[str, Any]) -> tuple:
        blob = f"{rec.get('title') or ''} {rec.get('snippet') or ''}".lower()
        has_term = any(tok in blob for tok in ("internaliz", "endocytos", "uptake"))
        has_pct = rec.get("pct_internalization") is not None
        return (not has_pct, not has_term, -(rec.get("weight") or 0))

    papers.sort(key=_rank)
    _enrich_pmc_pcts(session, papers)
    papers.sort(key=_rank)
    return papers[:16]


def _enrich_pmc_pcts(session: requests.Session, papers: list[dict[str, Any]], limit: int = 6) -> None:
    """If the abstract has no %, try NCBI PMC full text for a few OA papers."""
    from data_sources import _request

    used = 0
    for rec in papers:
        if used >= limit:
            break
        if rec.get("pct_values"):
            continue
        pmcid = str(rec.get("pmcid") or "").strip().upper().replace("PMC", "")
        if not pmcid.isdigit():
            continue
        try:
            xml = _request(
                session,
                "GET",
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                params={
                    "db": "pmc",
                    "id": pmcid,
                    "retmode": "xml",
                    "tool": "ImmunoTargetExplorer",
                    "email": "research@localhost",
                },
            ).text
        except Exception:
            continue
        used += 1
        pcts = extract_internalization_pcts(xml)
        if not pcts:
            continue
        rec["pct_values"] = pcts
        rec["pct_internalization"] = round(statistics.median(pcts), 1)
        rec["pct_from"] = "pmc_fulltext"
        rec["weight"] = _paper_weight(rec.get("title") or "", rec.get("snippet") or "", is_trial=False)


def fetch_clinical_trials(
    session: requests.Session,
    symbol: str,
    protein_name: str,
    synonyms: list[str] | None = None,
) -> list[dict[str, Any]]:
    terms = [symbol]
    for syn in synonyms or []:
        if syn and len(syn) <= 24 and syn.lower() not in {t.lower() for t in terms}:
            terms.append(syn)
    name_q = " OR ".join(terms[:4])
    query = (
        f'({name_q}) AND ("antibody-drug conjugate" OR ADC OR bispecific '
        "OR internalization OR endocytosis)"
    )
    try:
        payload = _get_json(
            session,
            CTGOV_STUDIES,
            {
                "query.term": query,
                "pageSize": 12,
                "countTotal": "true",
                "format": "json",
            },
        )
    except Exception:
        return []
    trials: list[dict[str, Any]] = []
    for study in (payload or {}).get("studies") or []:
        proto = study.get("protocolSection") or {}
        ident = proto.get("identificationModule") or {}
        desc = proto.get("descriptionModule") or {}
        status = proto.get("statusModule") or {}
        design = proto.get("designModule") or {}
        nct = ident.get("nctId") or ""
        title = ident.get("briefTitle") or ident.get("officialTitle") or ""
        summary = desc.get("briefSummary") or desc.get("detailedDescription") or ""
        phases = design.get("phases") or []
        blob = f"{title} {summary}"
        pcts = extract_internalization_pcts(blob)
        trials.append(
            {
                "nct_id": nct,
                "title": title,
                "status": status.get("overallStatus") or "",
                "phase": ", ".join(phases) if isinstance(phases, list) else str(phases or ""),
                "url": f"https://clinicaltrials.gov/study/{nct}" if nct else "https://clinicaltrials.gov/",
                "snippet": _clean_text(summary)[:280],
                "pct_values": pcts,
                "pct_internalization": round(statistics.median(pcts), 1) if pcts else None,
                "is_adc": any(tok in blob.lower() for tok in ("adc", "antibody-drug", "antibody drug")),
                "is_bispecific": "bispecific" in blob.lower(),
            }
        )
    return trials


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


ORDINAL_PCT = {
    "Not a typical surface TAA": 0.0,
    "Low / domain-restricted": 15.0,
    "Domain-dependent (literature)": 30.0,
    "Domain-dependent": 30.0,
    "Secreted - confirm membrane form": None,
    "Surface antigen; internalization not settled": None,
    "Insufficient literature": None,
    "Internalizes after antibody binding": 55.0,
    "Internalizes (literature)": 55.0,
}


def score_internalization_pct(
    uniprot: dict[str, Any],
    papers: list[dict[str, Any]],
    trials: list[dict[str, Any]],
    call: str,
) -> dict[str, Any]:
    """Weighted mean of literature/clinical % using the published acid-wash formula scale."""
    observations: list[dict[str, Any]] = []
    for rec in papers:
        if _is_degrader_paper(rec.get("title") or "", rec.get("snippet") or rec.get("abstract") or ""):
            rec["pct_values"] = []
            rec["pct_internalization"] = None
            continue
        pcts = rec.get("pct_values") or []
        if not pcts:
            continue
        median = round(float(statistics.median(pcts)), 1)
        rec["pct_internalization"] = median
        observations.append(
            {
                "source": "paper",
                "id": rec.get("pmid") or rec.get("title"),
                "pct": median,
                "weight": float(rec.get("weight") or 2.0),
            }
        )
    for trial in trials:
        pcts = trial.get("pct_values") or []
        if not pcts:
            continue
        median = round(float(statistics.median(pcts)), 1)
        trial["pct_internalization"] = median
        observations.append(
            {
                "source": "clinical_trial",
                "id": trial.get("nct_id") or trial.get("title"),
                "pct": median,
                "weight": 3.0,
            }
        )

    loc = uniprot.get("location_class") or "Unknown"
    if loc in {"Cytosol", "Nuclear"}:
        return {
            "internalization_pct": 0.0,
            "pct_n": 0,
            "pct_min": 0.0,
            "pct_max": 0.0,
            "pct_source": "topology",
            "pct_note": "Intracellular antigen: intact IgG internalization is not applicable.",
            "observations": observations,
            "formula_acid": INT_FORMULA_ACID,
            "formula_confocal": INT_FORMULA_CONFOCAL,
            "formula_mean": "sum(w_i * pct_i) / sum(w_i)",
        }

    if observations:
        num = sum(item["pct"] * item["weight"] for item in observations)
        den = sum(item["weight"] for item in observations) or 1.0
        mean_pct = round(num / den, 1)
        values = [item["pct"] for item in observations]
        return {
            "internalization_pct": mean_pct,
            "pct_n": len(observations),
            "pct_min": min(values),
            "pct_max": max(values),
            "pct_source": "literature_weighted_mean",
            "pct_note": (
                f"Weighted mean of {len(observations)} extracted % values "
                "(paper weight 2, clinical paper/trial weight 3)."
            ),
            "observations": observations,
            "formula_acid": INT_FORMULA_ACID,
            "formula_confocal": INT_FORMULA_CONFOCAL,
            "formula_mean": "sum(w_i * pct_i) / sum(w_i)",
        }

    mapped = ORDINAL_PCT.get(call)
    return {
        "internalization_pct": mapped,
        "pct_n": 0,
        "pct_min": mapped,
        "pct_max": mapped,
        "pct_source": "ordinal_map" if mapped is not None else "insufficient",
        "pct_note": (
            "No numeric % in retrieved abstracts/trials. "
            "Shown value is an ADC literature bin, not a measured assay."
            if mapped is not None
            else "No numeric internalization % found in papers or ClinicalTrials.gov records."
        ),
        "observations": [],
        "formula_acid": INT_FORMULA_ACID,
        "formula_confocal": INT_FORMULA_CONFOCAL,
        "formula_mean": "sum(w_i * pct_i) / sum(w_i)",
    }


def collect_taa_character(
    session: requests.Session,
    symbol: str,
    ensembl_id: str,
    gene_name: str = "",
    hpa_summary: dict[str, Any] | None = None,
) -> dict[str, Any]:
    uniprot = fetch_uniprot(session, symbol, ensembl_id)
    protein_name = uniprot.get("protein_name") or gene_name or symbol
    synonyms = uniprot.get("synonyms") or []
    papers = fetch_papers(session, symbol, protein_name, synonyms)
    trials = fetch_clinical_trials(session, symbol, protein_name, synonyms)
    consensus = _consensus_internalization(uniprot, papers)
    scored = score_internalization_pct(uniprot, papers, trials, consensus["call"])
    hpa_loc = ""
    if isinstance(hpa_summary, dict):
        hpa_loc = str(hpa_summary.get("Subcellular location") or "")
        if isinstance(hpa_summary.get("Subcellular location"), (list, dict)):
            hpa_loc = _text_blob(hpa_summary.get("Subcellular location"))
    adc_trials = sum(1 for t in trials if t.get("is_adc") or t.get("is_bispecific"))
    for rec in papers:
        rec.pop("abstract", None)
    return {
        "symbol": symbol,
        "ensembl_id": ensembl_id,
        "uniprot": uniprot,
        "hpa_location": hpa_loc,
        "papers": papers,
        "trials": trials,
        "internalization_call": consensus["call"],
        "modality_note": consensus["fit"],
        "paper_count": len(papers),
        "trial_count": len(trials),
        "adc_trial_count": adc_trials,
        "pct_schema": 2,
        **scored,
    }
