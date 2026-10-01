# -*- coding: utf-8 -*-
"""Cancer cell-line RNA from the HPA per-line matrix, labeled with Cell Model Passports.

Toil's "Cell Line" rows are GTEx transformed fibroblasts, EBV-transformed lymphocytes,
and one CML line. They are not this cohort. HPA's "RNA cell line specific nTPM"
field is an enrichment summary and is not read here.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests

from taa_analysis import ANALYSIS_VERSION
from taa_analysis.aliases import canonical_symbol, target_identity
from taa_analysis.settings import load_settings
from taa_analysis.xena import RetrievalFailed

HPA_CELLINE_URL = "https://www.proteinatlas.org/download/tsv/rna_celline.tsv.zip"
HPA_CELLINE_VERSION = "2025-11-05"
HPA_CELLINE_PAGE = "https://www.proteinatlas.org/humanproteome/cell+line"
CMP_API = "https://api.cellmodelpassports.sanger.ac.uk"
CMP_PAGE = "https://cellmodelpassports.sanger.ac.uk/"

SOURCE_COLUMNS = ("Gene", "Gene name", "Cell line", "TPM", "pTPM", "nTPM")
UNKNOWN_CANCER = "Unknown"
_ENSEMBL = re.compile(r"ENSG\d{11}", re.I)


def _now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _cache_dir() -> Path:
    path = Path(load_settings()["TAA_CACHE_DIR"])
    path.mkdir(parents=True, exist_ok=True)
    return path


def normalize_cell_line_name(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(name or "").lower())


def log2_ntpm_plus_1(value) -> float | None:
    """Display transform for HPA nTPM. This does not convert nTPM into TPM."""
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return None
    try:
        if pd.isna(value):
            return None
    except TypeError:
        pass
    number = float(value)
    if number < 0:
        raise ValueError(f"negative nTPM {number}")
    return math.log2(number + 1.0)


def _request_json(session: requests.Session, url: str) -> dict:
    response = session.get(url, headers={"Accept": "application/vnd.api+json"}, timeout=90)
    response.raise_for_status()
    payload = response.json()
    if not isinstance(payload, dict):
        raise RetrievalFailed(url, "Cell Model Passports response was not an object")
    return payload


def fetch_cmp_annotations(session: requests.Session | None = None) -> pd.DataFrame:
    """Model name, cancer type, and sample site. One passport row per model."""
    own_session = session is None
    session = session or requests.Session()
    try:
        cancer_payload = _request_json(session, f"{CMP_API}/cancer_types?page[size]=100")
        cancer_names = {
            str(item.get("id")): str((item.get("attributes") or {}).get("name") or "")
            for item in cancer_payload.get("data") or []
        }
        rows: list[dict] = []
        page = 1
        while page <= 40:
            payload = _request_json(
                session,
                f"{CMP_API}/models?page[size]=100&page[number]={page}&include=sample",
            )
            samples = {
                item.get("id"): item
                for item in payload.get("included") or []
                if item.get("type") == "sample"
            }
            for model in payload.get("data") or []:
                attrs = model.get("attributes") or {}
                names = [str(name) for name in (attrs.get("names") or []) if str(name).strip()]
                sample_rel = ((model.get("relationships") or {}).get("sample") or {}).get("data") or {}
                sample = samples.get(sample_rel.get("id")) or {}
                sample_attrs = sample.get("attributes") or {}
                cancer_rel = ((sample.get("relationships") or {}).get("cancer_type") or {}).get("data") or {}
                cancer_id = str(cancer_rel.get("id") or "")
                rows.append(
                    {
                        "model_id": str(model.get("id") or ""),
                        "model_name": names[0] if names else "",
                        "names": "|".join(names),
                        "model_type": str(attrs.get("model_type") or ""),
                        "cancer_type": cancer_names.get(cancer_id, ""),
                        "cancer_type_id": cancer_id,
                        "sample_site": str(sample_attrs.get("sample_site") or ""),
                        "tissue_status": str(sample_attrs.get("tissue_status") or ""),
                    }
                )
            if not (payload.get("links") or {}).get("next"):
                break
            page += 1
        if not rows:
            raise RetrievalFailed(CMP_API, "Cell Model Passports returned no models")
        return pd.DataFrame(rows)
    finally:
        if own_session:
            session.close()


def load_cmp_annotations(cache_dir: Path | None = None, force: bool = False) -> pd.DataFrame:
    cache_dir = cache_dir or _cache_dir()
    dest = cache_dir / "cmp_models.csv"
    if dest.exists() and not force and dest.stat().st_size > 1000:
        frame = pd.read_csv(dest, dtype=str).fillna("")
    else:
        frame = fetch_cmp_annotations()
        frame.to_csv(dest, index=False)
    required = {"names", "cancer_type", "model_id"}
    missing = required.difference(frame.columns)
    if missing:
        raise RetrievalFailed(str(dest), f"cell-line annotation is missing {sorted(missing)}")
    return frame


def annotation_lookup(annotations: pd.DataFrame) -> dict[str, dict]:
    """Map a normalized cell-line name to one cancer label.

    Several passport synonyms may share a name. Disagreeing cancer types do not
    expand the expression table and do not pick a winner.
    """
    grouped: dict[str, dict] = {}
    for record in annotations.itertuples(index=False):
        cancer = str(getattr(record, "cancer_type", "") or "")
        model_id = str(getattr(record, "model_id", "") or "")
        sample_site = str(getattr(record, "sample_site", "") or "")
        tissue_status = str(getattr(record, "tissue_status", "") or "")
        for name in str(getattr(record, "names", "") or "").split("|"):
            key = normalize_cell_line_name(name)
            if not key:
                continue
            bucket = grouped.setdefault(
                key,
                {"cancers": set(), "model_ids": set(), "sites": set(), "tissues": set()},
            )
            if cancer:
                bucket["cancers"].add(cancer)
            if model_id:
                bucket["model_ids"].add(model_id)
            if sample_site:
                bucket["sites"].add(sample_site)
            if tissue_status:
                bucket["tissues"].add(tissue_status)
    lookup: dict[str, dict] = {}
    for key, bucket in grouped.items():
        cancers = bucket["cancers"]
        if len(cancers) == 1:
            cancer = next(iter(cancers))
            status = "original"
        elif len(cancers) > 1:
            cancer = ""
            status = "conflict"
        else:
            cancer = ""
            status = "unavailable"
        lookup[key] = {
            "cancer_type": cancer,
            "annotation_status": status,
            "model_id": "|".join(sorted(bucket["model_ids"])),
            "sample_site": next(iter(bucket["sites"])) if len(bucket["sites"]) == 1 else "",
            "tissue_status": next(iter(bucket["tissues"])) if len(bucket["tissues"]) == 1 else "",
        }
    return lookup


def ensure_hpa_celline_zip(cache_dir: Path | None = None, force: bool = False) -> Path:
    cache_dir = cache_dir or _cache_dir()
    dest = cache_dir / "rna_celline.tsv.zip"
    if dest.exists() and not force and dest.stat().st_size > 1_000_000:
        return dest
    response = requests.get(HPA_CELLINE_URL, timeout=180, stream=True)
    response.raise_for_status()
    temporary = dest.with_suffix(".zip.partial")
    with temporary.open("wb") as handle:
        for chunk in response.iter_content(chunk_size=1024 * 1024):
            if chunk:
                handle.write(chunk)
    temporary.replace(dest)
    return dest


def ensembl_base(gene_id: str) -> str:
    match = _ENSEMBL.search(str(gene_id or ""))
    return match.group(0).upper() if match else ""


def resolve_catalog_symbol(catalog: pd.DataFrame, gene_query: str, required: bool = True) -> str:
    """Map a symbol, alias, or Ensembl id onto one HPA gene name."""
    requested = canonical_symbol(gene_query)
    symbols = {
        str(symbol).upper(): str(symbol)
        for symbol in catalog.get("gene_symbol", pd.Series(dtype=str)).tolist()
    }
    if requested.upper() in symbols:
        return symbols[requested.upper()]
    ensembl = ensembl_base(gene_query) or (ensembl_base(requested) if requested.upper().startswith("ENSG") else "")
    if ensembl and "gene_id" in catalog.columns:
        hits = catalog[catalog["gene_id"].astype(str).str.upper() == ensembl]
        found = sorted({str(symbol) for symbol in hits["gene_symbol"].tolist() if str(symbol)})
        if len(found) == 1:
            return found[0]
        if len(found) > 1:
            raise RetrievalFailed(HPA_CELLINE_URL, f"multiple HPA symbols for {ensembl}: {found}")
    if required:
        raise RetrievalFailed(HPA_CELLINE_URL, f"no HPA cell-line rows for {requested}")
    return requested


def extract_hpa_symbols(zip_path: Path, symbols: list[str] | set[str], catalog_path: Path | None = None) -> dict[str, pd.DataFrame]:
    """One pass over the HPA matrix. Keys are the gene names found in the file."""
    wanted_names = set()
    wanted_ensembl = set()
    for symbol in symbols:
        key = str(symbol or "").upper()
        if key.startswith("ENSG"):
            base = ensembl_base(key)
            if base:
                wanted_ensembl.add(base)
        elif key:
            wanted_names.add(key)
    write_catalog = catalog_path is not None and not catalog_path.exists()
    catalog: dict[str, str] = {}
    buckets: dict[str, list[dict]] = {}
    with zipfile.ZipFile(zip_path) as archive:
        names = archive.namelist()
        if not names:
            raise RetrievalFailed(str(zip_path), "HPA cell-line zip is empty")
        with archive.open(names[0]) as handle:
            header = handle.readline().decode("utf-8").rstrip("\n").split("\t")
            index = {column: position for position, column in enumerate(header)}
            missing = [column for column in ("Gene", "Gene name", "Cell line", "nTPM") if column not in index]
            if missing:
                raise RetrievalFailed(names[0], f"HPA cell-line header missing {missing}")
            for raw in handle:
                parts = raw.decode("utf-8").rstrip("\n").split("\t")
                if len(parts) <= index["nTPM"]:
                    continue
                gene_name = parts[index["Gene name"]]
                gene_id = ensembl_base(parts[index["Gene"]]) or parts[index["Gene"]]
                if write_catalog:
                    catalog.setdefault(gene_name, gene_id)
                if gene_name.upper() not in wanted_names and gene_id not in wanted_ensembl:
                    continue
                buckets.setdefault(gene_name, []).append(
                    {
                        "gene_id": gene_id,
                        "gene_symbol": gene_name,
                        "cell_line": parts[index["Cell line"]],
                        "nTPM": parts[index["nTPM"]],
                    }
                )
    if write_catalog and catalog:
        pd.DataFrame(
            [{"gene_symbol": name, "gene_id": gene_id} for name, gene_id in sorted(catalog.items())]
        ).to_csv(catalog_path, index=False)
    return {name: pd.DataFrame(rows) for name, rows in buckets.items()}


def extract_hpa_gene_rows(zip_path: Path, symbol: str, catalog_path: Path | None = None) -> pd.DataFrame:
    """Rows for one gene. Symbol and Ensembl id both select that gene only."""
    found = extract_hpa_symbols(zip_path, [symbol], catalog_path=catalog_path)
    if len(found) != 1:
        raise RetrievalFailed(HPA_CELLINE_URL, f"no HPA cell-line rows for {symbol}")
    frame = next(iter(found.values()))
    if frame.empty:
        raise RetrievalFailed(HPA_CELLINE_URL, f"no HPA cell-line rows for {symbol}")
    return frame


def _with_display(frame: pd.DataFrame) -> pd.DataFrame:
    display_values = []
    errors = []
    for value in frame["nTPM"]:
        try:
            display_values.append(log2_ntpm_plus_1(value))
            errors.append("")
        except (TypeError, ValueError) as exc:
            display_values.append(None)
            errors.append(str(exc))
    out = frame.copy()
    out["display_value"] = display_values
    out["transform_error"] = errors
    return out


def attach_cancer_types(expression: pd.DataFrame, annotations: pd.DataFrame) -> pd.DataFrame:
    """Left-join cancer labels. Row count stays equal to the expression matrix."""
    lookup = annotation_lookup(annotations)
    before = len(expression)
    rows = []
    for record in expression.itertuples(index=False):
        found = lookup.get(normalize_cell_line_name(record.cell_line), {})
        cancer = found.get("cancer_type") or ""
        status = found.get("annotation_status") or "unavailable"
        rows.append(
            {
                "gene_id": record.gene_id,
                "gene_symbol": record.gene_symbol,
                "cell_line": record.cell_line,
                "sample_id": record.cell_line,
                "patient_id": record.cell_line,
                "nTPM": pd.to_numeric(record.nTPM, errors="coerce"),
                "cancer_type": cancer,
                "group_label": cancer or UNKNOWN_CANCER,
                "annotation_status": status,
                "model_id": found.get("model_id", ""),
                "sample_site": found.get("sample_site", ""),
                "tissue_status": found.get("tissue_status", ""),
            }
        )
    labeled = pd.DataFrame(rows)
    if len(labeled) != before:
        raise RetrievalFailed(HPA_CELLINE_URL, "cell-line annotation changed the row count")
    return _with_display(labeled)


def summarize_cell_line_groups(lines: pd.DataFrame) -> pd.DataFrame:
    """Same factors as the tumor view: n, median, Q1, Q3, unit, source, missing."""
    rows = []
    if lines.empty:
        return pd.DataFrame(rows)
    for label, group in lines.groupby("group_label", dropna=False):
        display = pd.to_numeric(group["display_value"], errors="coerce")
        raw = pd.to_numeric(group["nTPM"], errors="coerce")
        present = display.dropna()
        raw_present = raw.dropna()
        rows.append(
            {
                "cancer": label,
                "cancer_code": label,
                "group_label": label,
                "n_cell_lines": int(group["cell_line"].nunique()),
                "n_samples": int(len(group)),
                "n_patients": None,
                "median": float(present.median()) if not present.empty else None,
                "q1": float(present.quantile(0.25)) if not present.empty else None,
                "q3": float(present.quantile(0.75)) if not present.empty else None,
                "median_nTPM": float(raw_present.median()) if not raw_present.empty else None,
                "q1_nTPM": float(raw_present.quantile(0.25)) if not raw_present.empty else None,
                "q3_nTPM": float(raw_present.quantile(0.75)) if not raw_present.empty else None,
                "unit": "log2(nTPM+1)",
                "source_unit": "nTPM",
                "source": "Human Protein Atlas cell-line RNA",
                "dataset_id": "rna_celline.tsv",
                "dataset_version": HPA_CELLINE_VERSION,
                "missing_n": int(display.isna().sum()),
                "transform_error_n": int((group["transform_error"].fillna("") != "").sum()),
                "distribution": "boxplot" if len(present) >= 2 else ("point" if len(present) == 1 else "unavailable"),
                "evidence_status": "Observed",
                "specimen_mode": "cell_line",
                "denominator": "HPA cancer cell line RNA, one cell line",
                "annotation_source": CMP_PAGE,
            }
        )
    summary = pd.DataFrame(rows)
    if not summary.empty and summary["median"].notna().any():
        summary = summary.sort_values("median", ascending=False, na_position="last")
    return summary.reset_index(drop=True)


def load_cell_line_result(gene_query: str, force: bool = False) -> dict:
    """Cell-line nTPM and summary factors for the queried TAA, not a fixed gene."""
    identity = target_identity(gene_query)
    cache_dir = _cache_dir()
    catalog_path = cache_dir / "hpa_celline_genes.csv"
    symbol = ""
    lines = pd.DataFrame()
    meta: dict = {}
    if catalog_path.exists() and not force:
        catalog = pd.read_csv(catalog_path, dtype=str).fillna("")
        symbol = resolve_catalog_symbol(catalog, gene_query, required=True)
        slice_path = cache_dir / f"hpa_celline_{symbol}.csv"
        meta_path = cache_dir / f"hpa_celline_{symbol}.json"
        if slice_path.exists() and meta_path.exists():
            cached = pd.read_csv(slice_path)
            cached_symbols = set(cached["gene_symbol"].astype(str)) if "gene_symbol" in cached.columns else set()
            if cached_symbols == {symbol}:
                lines = cached
                meta = json.loads(meta_path.read_text(encoding="utf-8"))
    if lines.empty:
        requested = symbol or canonical_symbol(gene_query)
        zip_path = ensure_hpa_celline_zip(cache_dir, force=False)
        expression = extract_hpa_gene_rows(zip_path, requested, catalog_path=catalog_path)
        symbol = str(expression["gene_symbol"].iloc[0])
        annotations = load_cmp_annotations(cache_dir, force=False)
        lines = attach_cancer_types(expression, annotations)
        digest = hashlib.sha256(lines["nTPM"].astype(str).str.cat(sep=",").encode()).hexdigest()[:16]
        matched = int((lines["annotation_status"] == "original").sum())
        meta = {
            "gene_symbol": symbol,
            "gene_id": str(lines["gene_id"].iloc[0]) if not lines.empty else "",
            "probe_id": "",
            "source_name": "Human Protein Atlas cell-line RNA",
            "source_url": HPA_CELLINE_URL,
            "dataset_id": "rna_celline.tsv",
            "dataset_version": HPA_CELLINE_VERSION,
            "unit_source": "nTPM",
            "unit_display": "log2(nTPM+1)",
            "formula": "display = log2(nTPM + 1). nTPM is not relabeled as TPM.",
            "input_checksum": digest,
            "retrieved_at": _now(),
            "analysis_version": ANALYSIS_VERSION,
            "n_rows": int(len(lines)),
            "n_matched_annotations": matched,
            "n_unknown_cancer": int((lines["group_label"] == UNKNOWN_CANCER).sum()),
            "annotation_source": "Cell Model Passports",
            "annotation_url": CMP_PAGE,
            "cohort_note": (
                "Per-cell-line HPA nTPM. Cancer labels come from Cell Model Passports name match. "
                "Unmatched or conflicting names stay Unknown. Toil GTEx transformed cell lines are not included. "
                "The HPA enrichment field RNA cell line specific nTPM is not used."
            ),
        }
        slice_path = cache_dir / f"hpa_celline_{symbol}.csv"
        meta_path = cache_dir / f"hpa_celline_{symbol}.json"
        lines.to_csv(slice_path, index=False)
        meta_path.write_text(json.dumps(meta), encoding="utf-8")
    identity["gene_symbol"] = str(meta.get("gene_symbol") or symbol)
    identity["gene_id"] = str(meta.get("gene_id") or "")
    if "group_label" not in lines.columns and "cancer_type" in lines.columns:
        lines["group_label"] = lines["cancer_type"].replace("", UNKNOWN_CANCER)
    summary = summarize_cell_line_groups(lines)
    coverage_rows = []
    for label, group in lines.groupby("group_label", dropna=False):
        coverage_rows.append(
            {
                "cancer": label,
                "n_cell_lines": int(group["cell_line"].nunique()),
                "rna_status": "available",
                "rna_reason": "",
                "annotation_status": "original" if label != UNKNOWN_CANCER else "unavailable",
                "subtype_status": "unavailable",
                "subtype_reason": "subtype_missing",
                "immune_status": "unavailable",
                "immune_reason": "unsupported",
            }
        )
    coverage = pd.DataFrame(coverage_rows).sort_values("cancer").reset_index(drop=True) if coverage_rows else pd.DataFrame()
    return {"identity": identity, "meta": meta, "lines": lines, "summary": summary, "coverage": coverage}
