# -*- coding: utf-8 -*-
"""Fetch normal/cancer expression from public APIs and save lightweight CSVs."""

from __future__ import annotations

import re
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

OT_GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"
HPA_JSON = "https://www.proteinatlas.org/{ensembl}.json"
GTEX_GENE = "https://gtexportal.org/api/v2/reference/gene"
GTEX_MEDIAN = "https://gtexportal.org/api/v2/expression/medianGeneExpression"

LOCAL_SAVE_DIR = Path(__file__).resolve().parent / "expression_csv"
SAVE_DIR = LOCAL_SAVE_DIR
DRIVE_SAVE_DIR = Path(r"G:\내 드라이브\Novel Target Database For cursor") / "expression_csv"
BUNDLE_VERSION = 3
TIMEOUT = (12, 90)
HEADERS = {
    "User-Agent": "Mozilla/5.0 (compatible; ImmunoTargetExplorer/1.0; research use)",
    "Accept": "application/json, text/plain, */*",
    "Connection": "close",
}

GENE_ALIASES = {
    "PD-1": "PDCD1",
    "PD1": "PDCD1",
    "PD 1": "PDCD1",
    "PDCD1": "PDCD1",
    "PD-L1": "CD274",
    "PDL1": "CD274",
    "PD L1": "CD274",
    "B7-H1": "CD274",
    "B7H1": "CD274",
    "CD274": "CD274",
    "PD-L2": "PDCD1LG2",
    "PDL2": "PDCD1LG2",
    "CTLA-4": "CTLA4",
    "CTLA4": "CTLA4",
    "LAG-3": "LAG3",
    "LAG3": "LAG3",
    "TIM-3": "HAVCR2",
    "TIM3": "HAVCR2",
    "HAVCR2": "HAVCR2",
    "TIGIT": "TIGIT",
    "OX40": "TNFRSF4",
    "4-1BB": "TNFRSF9",
    "41BB": "TNFRSF9",
    "CD137": "TNFRSF9",
    "ICOS": "ICOS",
    "GITR": "TNFRSF18",
    "CD47": "CD47",
    "SIRPA": "SIRPA",
    "SIRP-ALPHA": "SIRPA",
    "CD39": "ENTPD1",
    "CD73": "NT5E",
    "TGF-B": "TGFB1",
    "TGFB": "TGFB1",
    "VEGF": "VEGFA",
    "HER2": "ERBB2",
    "HER-2": "ERBB2",
    "ERBB2": "ERBB2",
    "NEU": "ERBB2",
    "P53": "TP53",
    "TP53": "TP53",
    "STING": "STING1",
    "CGAS": "CGAS",
    "NKG2A": "KLRC1",
    "CD94": "KLRD1",
    "VSIG4": "VSIG4",
    "LILRB1": "LILRB1",
    "LILRB2": "LILRB2",
    "SIGLEC15": "SIGLEC15",
    "SIGLEC-15": "SIGLEC15",
    "CD24": "CD24",
    "PVR": "PVR",
    "CD155": "PVR",
    "CD226": "CD226",
    "DNAM-1": "CD226",
    "CD96": "CD96",
    "CD112": "NECTIN2",
    "NECTIN2": "NECTIN2",
}

CANCER_WORDS = (
    "cancer",
    "carcinoma",
    "melanoma",
    "lymphoma",
    "leukemia",
    "leukaemia",
    "tumor",
    "tumour",
    "neoplasm",
    "sarcoma",
    "glioma",
    "myeloma",
    "malignan",
    "blastoma",
    "adenocarcinoma",
)


def save_dirs() -> list[Path]:
    LOCAL_SAVE_DIR.mkdir(parents=True, exist_ok=True)
    dirs = [LOCAL_SAVE_DIR]
    try:
        if DRIVE_SAVE_DIR.parent.exists():
            dirs.append(DRIVE_SAVE_DIR)
    except OSError:
        pass
    return dirs


def _session() -> requests.Session:
    s = requests.Session()
    s.headers.update(HEADERS)
    retry_kw = dict(
        total=5,
        connect=5,
        read=5,
        backoff_factor=0.8,
        status_forcelist=(429, 500, 502, 503, 504),
        raise_on_status=False,
    )
    try:
        retry = Retry(allowed_methods=frozenset(["GET", "POST"]), **retry_kw)
    except TypeError:
        retry = Retry(method_whitelist=frozenset(["GET", "POST"]), **retry_kw)
    adapter = HTTPAdapter(max_retries=retry, pool_connections=8, pool_maxsize=8)
    s.mount("https://", adapter)
    s.mount("http://", adapter)
    return s


def _request(session: requests.Session, method: str, url: str, **kwargs) -> requests.Response:
    kwargs.setdefault("timeout", TIMEOUT)
    last: Exception | None = None
    for attempt in range(5):
        try:
            response = session.request(method, url, **kwargs)
            response.raise_for_status()
            return response
        except (requests.ConnectionError, requests.Timeout, requests.exceptions.ChunkedEncodingError) as exc:
            last = exc
            time.sleep(0.6 * (2 ** attempt))
    raise ConnectionError(
        "A public database closed the connection. Wait a few seconds and search again."
    ) from last


def graphql(session: requests.Session, query: str, variables: dict[str, Any]) -> dict[str, Any]:
    response = _request(
        session,
        "POST",
        OT_GRAPHQL,
        json={"query": query, "variables": variables},
        headers={"Content-Type": "application/json"},
    )
    payload = response.json()
    if payload.get("errors"):
        messages = "; ".join(err.get("message", str(err)) for err in payload["errors"])
        raise RuntimeError("Open Targets error: " + messages)
    return payload.get("data") or {}


def extract_gene_query(text: str) -> str | None:
    if not text:
        return None
    raw = text.strip()
    ensembl = re.search(r"ENSG\d{11}", raw, flags=re.I)
    if ensembl:
        return ensembl.group(0).upper()

    compact = re.sub(r"\s+", " ", raw.upper())
    for alias in sorted(GENE_ALIASES, key=len, reverse=True):
        pattern = r"(?<![A-Z0-9])" + re.escape(alias.upper()) + r"(?![A-Z0-9])"
        if re.search(pattern, compact):
            return GENE_ALIASES[alias]

    for token in re.findall(r"\b[A-Z][A-Z0-9-]{1,14}\b", raw.upper()):
        if token in {"CSV", "API", "RNA", "TPM", "HPA", "OT", "GTEX", "TCGA", "PDF"}:
            continue
        return GENE_ALIASES.get(token, token)
    return None


extract_gene_query = extract_gene_query


def resolve_target(session: requests.Session, query: str) -> dict[str, str]:
    query = query.strip()
    if query.upper().startswith("ENSG"):
        data = graphql(
            session,
            """
            query GetTarget($ensemblId: String!) {
              target(ensemblId: $ensemblId) {
                id
                approvedSymbol
                approvedName
              }
            }
            """,
            {"ensemblId": query.upper()},
        )
        target = data.get("target")
        if not target:
            raise ValueError("Ensembl ID not found: " + query)
        return {
            "ensembl_id": target["id"],
            "symbol": target["approvedSymbol"],
            "name": target.get("approvedName") or "",
        }

    data = graphql(
        session,
        """
        query SearchTarget($queryString: String!) {
          search(queryString: $queryString, entityNames: ["target"], page: {index: 0, size: 5}) {
            hits { id name entity description }
          }
        }
        """,
        {"queryString": query},
    )
    hits = (data.get("search") or {}).get("hits") or []
    hits = [h for h in hits if h.get("entity") == "target" and str(h.get("id", "")).startswith("ENSG")]
    if not hits:
        raise ValueError("Gene not found: " + query)

    symbol_upper = query.upper()
    chosen = hits[0]
    for hit in hits:
        if str(hit.get("name", "")).upper() == symbol_upper:
            chosen = hit
            break
    return {
        "ensembl_id": chosen["id"],
        "symbol": chosen.get("name") or query.upper(),
        "name": chosen.get("description") or "",
    }


def fetch_open_targets(session: requests.Session, ensembl_id: str) -> dict[str, Any]:
    data = graphql(
        session,
        """
        query TargetBundle($ensemblId: String!) {
          target(ensemblId: $ensemblId) {
            id
            approvedSymbol
            approvedName
            biotype
            functionDescriptions
            isEssential
            tractability { label modality value }
            baselineExpression(page: {index: 0, size: 250}) {
              count
              rows {
                unit median min q1 q3 max datasourceId datatypeId
                tissueBiosample { biosampleId biosampleName }
                celltypeBiosample { biosampleId biosampleName }
              }
            }
            associatedDiseases(page: {index: 0, size: 50}) {
              count
              rows {
                score
                disease {
                  id
                  name
                  therapeuticAreas { id name }
                }
              }
            }
          }
        }
        """,
        {"ensemblId": ensembl_id},
    )
    target = data.get("target")
    if not target:
        raise ValueError("Could not load target from Open Targets.")
    return target


def _is_cancer_disease(row: dict[str, Any]) -> bool:
    disease = row.get("disease") or {}
    names = [str(disease.get("name") or "")]
    for area in disease.get("therapeuticAreas") or []:
        names.append(str(area.get("name") or ""))
    blob = " ".join(names).lower()
    return any(word in blob for word in CANCER_WORDS)


def expression_rows_from_open_targets(target: dict[str, Any]) -> list[dict[str, Any]]:
    symbol = target.get("approvedSymbol") or ""
    ensembl_id = target.get("id") or ""
    rows: list[dict[str, Any]] = []
    for item in (target.get("baselineExpression") or {}).get("rows") or []:
        tissue = (item.get("tissueBiosample") or {}).get("biosampleName") or ""
        cell = (item.get("celltypeBiosample") or {}).get("biosampleName") or ""
        source = item.get("datasourceId") or "open_targets"
        condition = "normal" if "gtex" in source.lower() or not cell else "unspecified"
        if cell and "tabula" in source.lower():
            condition = "normal_single_cell"
        rows.append(
            {
                "gene_symbol": symbol,
                "ensembl_id": ensembl_id,
                "condition": condition,
                "tissue_or_cancer": tissue or cell or "unknown",
                "cell_type": cell,
                "expression_value": item.get("median"),
                "unit": item.get("unit") or "",
                "source": "Open Targets (" + source + ")",
                "extra": "",
            }
        )
    return rows


def association_rows_from_open_targets(target: dict[str, Any]) -> list[dict[str, Any]]:
    symbol = target.get("approvedSymbol") or ""
    ensembl_id = target.get("id") or ""
    rows: list[dict[str, Any]] = []
    for item in (target.get("associatedDiseases") or {}).get("rows") or []:
        if not _is_cancer_disease(item):
            continue
        disease = item.get("disease") or {}
        areas = ", ".join(
            area.get("name") or "" for area in (disease.get("therapeuticAreas") or []) if area.get("name")
        )
        rows.append(
            {
                "gene_symbol": symbol,
                "ensembl_id": ensembl_id,
                "disease": disease.get("name") or "",
                "disease_id": disease.get("id") or "",
                "association_score": item.get("score"),
                "therapeutic_areas": areas,
                "source": "Open Targets",
            }
        )
    rows.sort(key=lambda r: float(r.get("association_score") or 0), reverse=True)
    return rows


def fetch_hpa(session: requests.Session, ensembl_id: str, symbol: str) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    try:
        response = _request(session, "GET", HPA_JSON.format(ensembl=ensembl_id))
        info = response.json()
    except Exception:
        return [], {}

    rows: list[dict[str, Any]] = []

    def add_map(mapping: Any, condition: str, unit: str, extra: str) -> None:
        if isinstance(mapping, dict):
            for tissue, value in mapping.items():
                rows.append(
                    {
                        "gene_symbol": symbol,
                        "ensembl_id": ensembl_id,
                        "condition": condition,
                        "tissue_or_cancer": str(tissue),
                        "cell_type": "",
                        "expression_value": value,
                        "unit": unit,
                        "source": "Human Protein Atlas",
                        "extra": extra,
                    }
                )

    add_map(info.get("RNA tissue specific nTPM"), "normal", "nTPM", "tissue-specific")
    add_map(info.get("RNA cancer specific pTPM"), "cancer", "pTPM", "cancer-specific")

    summary: dict[str, Any] = {
        "RNA tissue specificity": info.get("RNA tissue specificity"),
        "RNA cancer specificity": info.get("RNA cancer specificity"),
        "RNA tissue distribution": info.get("RNA tissue distribution"),
        "RNA cancer distribution": info.get("RNA cancer distribution"),
        "RNA cell line specific nTPM": info.get("RNA cell line specific nTPM"),
        "Gene description": info.get("Gene description"),
        "Disease involvement": info.get("Disease involvement"),
        "Protein class": info.get("Protein class"),
    }

    prognostics = []
    for key, value in info.items():
        if not str(key).startswith("Cancer prognostics") or not isinstance(value, dict):
            continue
        cancer_name = str(key).replace("Cancer prognostics - ", "")
        prognostics.append(
            {
                "cancer": cancer_name,
                "prognostic": value.get("prognostic"),
                "prognostic_type": value.get("prognostic type"),
                "is_prognostic": value.get("is_prognostic"),
                "p_val": value.get("p_val"),
            }
        )
        if value.get("is_prognostic"):
            rows.append(
                {
                    "gene_symbol": symbol,
                    "ensembl_id": ensembl_id,
                    "condition": "cancer",
                    "tissue_or_cancer": cancer_name,
                    "cell_type": "",
                    "expression_value": "",
                    "unit": "",
                    "source": "Human Protein Atlas",
                    "extra": "prognostic=%s; type=%s; p=%s"
                    % (value.get("prognostic"), value.get("prognostic type"), value.get("p_val")),
                }
            )
    summary["prognostics"] = prognostics
    return rows, summary


def fetch_gtex_normal(session: requests.Session, ensembl_id: str, symbol: str) -> list[dict[str, Any]]:
    try:
        gene_res = _request(session, "GET", GTEX_GENE, params={"geneId": ensembl_id, "pageSize": 5})
        genes = (gene_res.json() or {}).get("data") or []
        if not genes:
            return []
        gencode_id = genes[0].get("gencodeId")
        if not gencode_id:
            return []
        exp_res = _request(
            session,
            "GET",
            GTEX_MEDIAN,
            params={"gencodeId": gencode_id, "datasetId": "gtex_v8"},
        )
        items = (exp_res.json() or {}).get("data") or exp_res.json()
        if isinstance(items, dict):
            items = items.get("medianGeneExpression") or items.get("data") or []
        rows: list[dict[str, Any]] = []
        for item in items or []:
            rows.append(
                {
                    "gene_symbol": symbol,
                    "ensembl_id": ensembl_id,
                    "condition": "normal",
                    "tissue_or_cancer": item.get("tissueSiteDetailId") or item.get("tissue") or "",
                    "cell_type": "",
                    "expression_value": item.get("median") if item.get("median") is not None else item.get("value"),
                    "unit": "TPM",
                    "source": "GTEx v8",
                    "extra": "",
                }
            )
        return rows
    except Exception:
        return []


def _safe_filename(symbol: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_-]+", "_", symbol).strip("_")
    return cleaned or "gene"


def _write_csv(path: Path, df: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, encoding="utf-8-sig")


def collect_and_save(gene_query: str) -> dict[str, Any]:
    from analysis import build_cancer_analytics_table
    from cptac import collect_cptac_for_gene
    from expression_compare import build_expression_comparison, build_tcga_summary

    session = _session()
    identity = resolve_target(session, gene_query)
    ensembl_id = identity["ensembl_id"]
    symbol = identity["symbol"]

    ot_target = fetch_open_targets(session, ensembl_id)
    if ot_target.get("approvedSymbol"):
        symbol = ot_target["approvedSymbol"]
        identity["symbol"] = symbol
        identity["name"] = ot_target.get("approvedName") or identity.get("name") or ""

    expression_rows = expression_rows_from_open_targets(ot_target)
    hpa_rows, hpa_summary = fetch_hpa(session, ensembl_id, symbol)
    expression_rows.extend(hpa_rows)

    has_numeric = any(
        row.get("expression_value") not in (None, "") and str(row.get("condition", "")).startswith("normal")
        for row in expression_rows
    )
    if not has_numeric:
        expression_rows.extend(fetch_gtex_normal(session, ensembl_id, symbol))

    associations = association_rows_from_open_targets(ot_target)

    expr_df = pd.DataFrame(expression_rows)
    if not expr_df.empty:
        expr_df = expr_df.drop_duplicates()
    else:
        expr_df = pd.DataFrame(
            columns=[
                "gene_symbol",
                "ensembl_id",
                "condition",
                "tissue_or_cancer",
                "cell_type",
                "expression_value",
                "unit",
                "source",
                "extra",
            ]
        )

    assoc_df = pd.DataFrame(associations)
    if assoc_df.empty:
        assoc_df = pd.DataFrame(
            columns=[
                "gene_symbol",
                "ensembl_id",
                "disease",
                "disease_id",
                "association_score",
                "therapeutic_areas",
                "source",
            ]
        )

    analytics_df = build_cancer_analytics_table(
        symbol=symbol,
        ensembl_id=ensembl_id,
        expr_df=expr_df,
        assoc_df=assoc_df,
        hpa_summary=hpa_summary,
        session=session,
    )
    comparison_df, tcga_df, normal_rna_df = build_expression_comparison(
        symbol=symbol,
        ensembl_id=ensembl_id,
        expr_df=expr_df,
        ihc_df=analytics_df,
        session=session,
    )
    tcga_summary_df = build_tcga_summary(tcga_df)
    cptac_df, cptac_summary_df, cptac_error = collect_cptac_for_gene(session, symbol, ensembl_id)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d")
    stem = _safe_filename(symbol)

    saved_expr = []
    saved_assoc = []
    saved_analytics = []
    saved_compare = []
    saved_cptac = []
    for folder in save_dirs():
        expr_path = folder / f"{stem}_normal_cancer_expression_{stamp}.csv"
        assoc_path = folder / f"{stem}_cancer_associations_{stamp}.csv"
        analytics_path = folder / f"{stem}_cancer_analytics_{stamp}.csv"
        compare_path = folder / f"{stem}_rna_protein_comparison_{stamp}.csv"
        cptac_path = folder / f"{stem}_cptac_protein_{stamp}.csv"
        _write_csv(expr_path, expr_df)
        _write_csv(assoc_path, assoc_df)
        _write_csv(analytics_path, analytics_df if not analytics_df.empty else pd.DataFrame())
        _write_csv(compare_path, comparison_df if not comparison_df.empty else pd.DataFrame())
        _write_csv(cptac_path, cptac_summary_df if not cptac_summary_df.empty else pd.DataFrame())
        saved_expr.append(str(expr_path))
        saved_assoc.append(str(assoc_path))
        saved_analytics.append(str(analytics_path))
        saved_compare.append(str(compare_path))
        saved_cptac.append(str(cptac_path))

    return {
        "identity": identity,
        "target": ot_target,
        "hpa_summary": hpa_summary,
        "expression_df": expr_df,
        "association_df": assoc_df,
        "analytics_df": analytics_df,
        "comparison_df": comparison_df,
        "tcga_df": tcga_df,
        "tcga_summary_df": tcga_summary_df,
        "cptac_df": cptac_df,
        "cptac_summary_df": cptac_summary_df,
        "cptac_error": cptac_error,
        "normal_rna_df": normal_rna_df,
        "expression_csv": saved_expr[0],
        "association_csv": saved_assoc[0],
        "analytics_csv": saved_analytics[0] if saved_analytics else "",
        "comparison_csv": saved_compare[0] if saved_compare else "",
        "cptac_csv": saved_cptac[0] if saved_cptac else "",
        "retrieved_at": datetime.now(timezone.utc).strftime("%Y-%m-%d"),
        "bundle_version": BUNDLE_VERSION,
        "all_expression_csv": saved_expr,
        "all_association_csv": saved_assoc,
        "all_analytics_csv": saved_analytics,
        "save_dir": str(save_dirs()[-1]),
    }
