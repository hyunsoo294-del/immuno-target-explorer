# -*- coding: utf-8 -*-
"""RNA and protein expression comparison: HPA, GTEx, TCGA, Open Targets."""

from __future__ import annotations

import io
import math
import re
import zipfile
from pathlib import Path
from typing import Any

import pandas as pd
import requests

from csv_io import read_csv_safe
from data_qc import (
    RNA_PSEUDOCOUNT,
    drop_ihc_pseudocount_fold,
    expand_ihc_counts,
    mannwhitney_ihc,
    remap_rna_arrays,
    sanity_check_tcga_inverse,
    tcga_log2_plus1,
)

CACHE_DIR = Path(__file__).resolve().parent / "cache"
TIMEOUT = 120
PSEUDOCOUNT = RNA_PSEUDOCOUNT
IHC_MAP = {"not detected": 0, "low": 1, "medium": 2, "high": 3}

HPA_FILES = {
    "consensus_rna": (
        "https://www.proteinatlas.org/download/tsv/rna_tissue_consensus.tsv.zip",
        "rna_tissue_consensus.tsv",
    ),
    "gtex_rna": (
        "https://www.proteinatlas.org/download/tsv/rna_tissue_gtex.tsv.zip",
        "rna_tissue_gtex.tsv",
    ),
    "normal_ihc": (
        "https://www.proteinatlas.org/download/tsv/normal_ihc_data.tsv.zip",
        "normal_ihc_data.tsv",
    ),
}

# HPA pathology cancer type -> matched normal organ + TCGA pan-cancer atlas studies
CANCER_CONTEXT = {
    "breast cancer": {"organ": "breast", "studies": ["brca_tcga_pan_can_atlas_2018"]},
    "lung cancer": {"organ": "lung", "studies": ["luad_tcga_pan_can_atlas_2018", "lusc_tcga_pan_can_atlas_2018"]},
    "colorectal cancer": {"organ": "colon", "studies": ["coadread_tcga_pan_can_atlas_2018"]},
    "liver cancer": {"organ": "liver", "studies": ["lihc_tcga_pan_can_atlas_2018"]},
    "prostate cancer": {"organ": "prostate", "studies": ["prad_tcga_pan_can_atlas_2018"]},
    "ovarian cancer": {"organ": "ovary", "studies": ["ov_tcga_pan_can_atlas_2018"]},
    "pancreatic cancer": {"organ": "pancreas", "studies": ["paad_tcga_pan_can_atlas_2018"]},
    "renal cancer": {"organ": "kidney", "studies": ["kirc_tcga_pan_can_atlas_2018"]},
    "stomach cancer": {"organ": "stomach", "studies": ["stad_tcga_pan_can_atlas_2018"]},
    "thyroid cancer": {"organ": "thyroid gland", "studies": ["thca_tcga_pan_can_atlas_2018"]},
    "cervical cancer": {"organ": "cervix", "studies": ["cesc_tcga_pan_can_atlas_2018"]},
    "endometrial cancer": {"organ": "endometrium", "studies": ["ucec_tcga_pan_can_atlas_2018"]},
    "glioma": {"organ": "cerebral cortex", "studies": ["gbm_tcga_pan_can_atlas_2018", "lgg_tcga_pan_can_atlas_2018"]},
    "head and neck cancer": {"organ": "salivary gland", "studies": ["hnsc_tcga_pan_can_atlas_2018"]},
    "melanoma": {"organ": "skin", "studies": ["skcm_tcga_pan_can_atlas_2018"]},
    "skin cancer": {"organ": "skin", "studies": ["skcm_tcga_pan_can_atlas_2018"]},
    "testis cancer": {"organ": "testis", "studies": ["tgct_tcga_pan_can_atlas_2018"]},
    "urothelial cancer": {"organ": "urinary bladder", "studies": ["blca_tcga_pan_can_atlas_2018"]},
    "lymphoma": {"organ": "lymph node", "studies": ["dlbc_tcga_pan_can_atlas_2018"]},
}


def _norm(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", str(text or "").lower()).strip()


def _num(value: Any) -> float | None:
    try:
        if value is None or value == "":
            return None
        out = float(value)
        if math.isnan(out) or math.isinf(out):
            return None
        return out
    except (TypeError, ValueError):
        return None


def _ihc_score(level: Any) -> float | None:
    if level is None:
        return None
    return IHC_MAP.get(str(level).strip().lower())


def _ensure_hpa_file(key: str, session: requests.Session) -> Path:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    url, inner_name = HPA_FILES[key]
    dest = CACHE_DIR / inner_name
    if dest.exists() and dest.stat().st_size > 1000:
        return dest
    response = session.get(url, timeout=TIMEOUT)
    response.raise_for_status()
    with zipfile.ZipFile(io.BytesIO(response.content)) as zf:
        dest.write_bytes(zf.read(zf.namelist()[0]))
    return dest


def _filter_gene_tsv(path: Path, ensembl_id: str) -> pd.DataFrame:
    chunks = []
    for chunk in read_csv_safe(path, sep="\t", dtype=str, chunksize=200000):
        hit = chunk[chunk["Gene"].astype(str).str.upper() == ensembl_id.upper()]
        if not hit.empty:
            chunks.append(hit)
    if not chunks:
        return pd.DataFrame()
    return pd.concat(chunks, ignore_index=True)


def _match_organ_mean(df: pd.DataFrame, organ: str, value_col: str) -> float | None:
    if df is None or df.empty or value_col not in df.columns:
        return None
    work = df.copy()
    work["_tissue"] = work["Tissue"].map(_norm)
    organ_n = _norm(organ)
    matched = work[work["_tissue"].str.contains(organ_n, na=False) | work["_tissue"].eq(organ_n)]
    if matched.empty:
        # allow organ contained in longer GTEx names, e.g. "lung" vs "lung"
        matched = work[work["_tissue"].apply(lambda t: organ_n in t or t in organ_n)]
    values = pd.to_numeric(matched[value_col], errors="coerce").dropna()
    if values.empty:
        return None
    return float(values.mean())


def _normal_ihc_organ_mean(df: pd.DataFrame, organ: str) -> float | None:
    scores = _normal_ihc_score_array(df, organ)
    if scores.size == 0:
        return None
    return float(scores.mean())


def _normal_ihc_score_array(df: pd.DataFrame, organ: str):
    import numpy as np

    if df is None or df.empty:
        return np.asarray([], dtype=float)
    work = df.copy()
    work["_score"] = work["Level"].map(_ihc_score)
    work["_tissue"] = work["Tissue"].map(_norm)
    organ_n = _norm(organ)
    matched = work[work["_tissue"].str.contains(organ_n, na=False) | work["_tissue"].eq(organ_n)]
    scores = pd.to_numeric(matched["_score"], errors="coerce").dropna()
    return scores.to_numpy(dtype=float)


def _gtex_from_expression(expr_df: pd.DataFrame, organ: str) -> float | None:
    if expr_df is None or expr_df.empty:
        return None
    work = expr_df.copy()
    work["value_num"] = pd.to_numeric(work["expression_value"], errors="coerce")
    work = work.dropna(subset=["value_num"])
    src = work["source"].astype(str).str.lower()
    cond = work["condition"].astype(str).str.lower()
    normal = work[cond.str.contains("normal", na=False) | src.str.contains("gtex|open targets", na=False)]
    if normal.empty:
        return None
    organ_n = _norm(organ)
    tissue = normal["tissue_or_cancer"].astype(str).map(_norm)
    matched = normal[tissue.str.contains(organ_n, na=False) | tissue.eq(organ_n)]
    if matched.empty:
        return float(normal["value_num"].mean())
    return float(matched["value_num"].mean())


def fetch_entrez_id(session: requests.Session, symbol: str, ensembl_id: str | None = None) -> int | None:
    for query in (symbol, ensembl_id):
        if not query:
            continue
        try:
            response = session.get(f"https://www.cbioportal.org/api/genes/{query}", timeout=45)
            if response.status_code == 404:
                continue
            response.raise_for_status()
            payload = response.json()
            if isinstance(payload, list):
                payload = payload[0] if payload else {}
            entrez = payload.get("entrezGeneId")
            if entrez is not None:
                return int(entrez)
        except Exception:
            continue
    try:
        response = session.get(
            "https://www.cbioportal.org/api/genes",
            params={"keyword": symbol, "pageSize": 5, "projection": "SUMMARY"},
            timeout=45,
        )
        response.raise_for_status()
        hits = response.json() or []
        sym_u = symbol.upper()
        for hit in hits:
            if str(hit.get("hugoGeneSymbol", "")).upper() == sym_u:
                return int(hit["entrezGeneId"])
        if hits:
            return int(hits[0]["entrezGeneId"])
    except Exception:
        pass
    return None


def fetch_tcga_rna(session: requests.Session, entrez_id: int) -> pd.DataFrame:
    profile_ids = []
    for item in CANCER_CONTEXT.values():
        for study in item["studies"]:
            profile_ids.append(f"{study}_rna_seq_v2_mrna")
    profile_ids = sorted(set(profile_ids))
    try:
        response = session.post(
            "https://www.cbioportal.org/api/molecular-data/fetch?projection=SUMMARY",
            json={"entrezGeneIds": [entrez_id], "molecularProfileIds": profile_ids},
            timeout=90,
            headers={"Content-Type": "application/json"},
        )
        response.raise_for_status()
        rows = response.json() or []
    except Exception:
        return pd.DataFrame()

    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    df = df.dropna(subset=["value"])
    df["sample_code"] = df["sampleId"].astype(str).str.extract(r"-(\d{2})$", expand=False)
    df["sample_class"] = df["sample_code"].map(lambda x: "tumor" if x == "01" else ("normal" if x == "11" else "other"))
    study_to_cancer = {}
    for cancer, meta in CANCER_CONTEXT.items():
        for study in meta["studies"]:
            study_to_cancer[study] = cancer
    df["cancer_type"] = df["studyId"].map(study_to_cancer)
    return df


def _minmax(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    valid = values.dropna()
    if valid.empty:
        return pd.Series([float("nan")] * len(values), index=series.index)
    lo = float(valid.min())
    hi = float(valid.max())
    if hi == lo:
        return values.map(lambda v: 0.0 if pd.notna(v) else float("nan"))
    return (values - lo) / (hi - lo)


def _zscore(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    std = float(values.std(ddof=0))
    if std == 0 or math.isnan(std):
        return pd.Series([0.0] * len(values), index=series.index)
    return (values - float(values.mean())) / std


def build_expression_comparison(
    *,
    symbol: str,
    ensembl_id: str,
    expr_df: pd.DataFrame,
    ihc_df: pd.DataFrame,
    session: requests.Session | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return (comparison table, TCGA sample table, normal RNA table)."""
    sess = session or requests.Session()
    sess.headers.setdefault("User-Agent", "NovelTargetExplorer/1.0")

    consensus = _filter_gene_tsv(_ensure_hpa_file("consensus_rna", sess), ensembl_id)
    gtex = _filter_gene_tsv(_ensure_hpa_file("gtex_rna", sess), ensembl_id)
    normal_ihc = _filter_gene_tsv(_ensure_hpa_file("normal_ihc", sess), ensembl_id)

    for frame in (consensus, gtex):
        if not frame.empty:
            for col in ("nTPM", "TPM", "pTPM"):
                if col in frame.columns:
                    frame[col] = pd.to_numeric(frame[col], errors="coerce")

    entrez = fetch_entrez_id(sess, symbol, ensembl_id)
    tcga = fetch_tcga_rna(sess, entrez) if entrez else pd.DataFrame()

    rows: list[dict[str, Any]] = []
    for cancer_type, meta in CANCER_CONTEXT.items():
        organ = meta["organ"]
        hpa_normal_ntpm = _match_organ_mean(consensus, organ, "nTPM")
        gtex_ntpm = _match_organ_mean(gtex, organ, "nTPM")
        gtex_tpm = _match_organ_mean(gtex, organ, "TPM")
        ot_gtex = _gtex_from_expression(expr_df, organ)
        normal_rna = hpa_normal_ntpm if hpa_normal_ntpm is not None else gtex_ntpm
        if normal_rna is None:
            normal_rna = gtex_tpm if gtex_tpm is not None else ot_gtex
        protein_normal = _normal_ihc_organ_mean(normal_ihc, organ)
        normal_ihc_scores = _normal_ihc_score_array(normal_ihc, organ)

        ihc_cancer = None
        protein_cancer = None
        tumor_ihc_scores = expand_ihc_counts(0, 0, 0, 0)
        if ihc_df is not None and not ihc_df.empty:
            hit = ihc_df[ihc_df["cancer_type"].astype(str).str.lower() == cancer_type]
            if not hit.empty:
                rec = hit.iloc[0]
                protein_cancer = (
                    int(rec["ihc_high"]) * 3
                    + int(rec["ihc_medium"]) * 2
                    + int(rec["ihc_low"]) * 1
                ) / max(float(rec["patients_total"]), 1.0)
                ihc_cancer = protein_cancer
                tumor_ihc_scores = expand_ihc_counts(
                    rec.get("ihc_high"),
                    rec.get("ihc_medium"),
                    rec.get("ihc_low"),
                    rec.get("ihc_not_detected"),
                )

        tcga_tumor = None
        tcga_adj = None
        rna_n_tumor = 0
        rna_n_normal = 0
        if not tcga.empty:
            subset = tcga[tcga["cancer_type"] == cancer_type]
            tumor = subset[subset["sample_class"] == "tumor"]["value"]
            adj = subset[subset["sample_class"] == "normal"]["value"]
            rna_n_tumor = int(tumor.count())
            rna_n_normal = int(adj.count())
            if not tumor.empty:
                tcga_tumor = float(tumor.median())
            if not adj.empty:
                tcga_adj = float(adj.mean())

        cancer_rna_abs = tcga_tumor
        normal_rna_abs = normal_rna
        mwu = mannwhitney_ihc(tumor_ihc_scores, normal_ihc_scores)
        check = sanity_check_tcga_inverse(tcga_tumor)

        rows.append(
            {
                "gene_symbol": symbol,
                "ensembl_id": ensembl_id,
                "organ": organ,
                "cancer_type": cancer_type,
                "normal_rna_hpa_ntpm": hpa_normal_ntpm,
                "normal_rna_gtex_ntpm": gtex_ntpm,
                "normal_rna_ot_gtex": ot_gtex,
                "normal_rna_abs": normal_rna_abs,
                "cancer_rna_tcga_median": tcga_tumor,
                "cancer_rna_tcga_adjacent_mean": tcga_adj,
                "cancer_rna_abs": cancer_rna_abs,
                "rna_n_tumor": rna_n_tumor,
                "rna_n_adjacent_normal": rna_n_normal,
                "rna_log2_normal": math.log2(normal_rna_abs + PSEUDOCOUNT) if normal_rna_abs is not None else None,
                "rna_log2_cancer": tcga_log2_plus1(tcga_tumor),
                "tcga_inverse_2pow_minus1": check["reconstructed_from_2pow_minus1"],
                "tcga_inverse_sanity_ok": check["sanity_ok"],
                "protein_ihc_normal_0_3": protein_normal,
                "protein_ihc_cancer_0_3": ihc_cancer,
                "ihc_mwu_u": mwu["ihc_mwu_u"],
                "ihc_mwu_pvalue": mwu["ihc_mwu_pvalue"],
                "ihc_mwu_effect": mwu["ihc_mwu_effect"],
                "ihc_mwu_n_tumor": mwu["ihc_mwu_n_tumor"],
                "ihc_mwu_n_normal": mwu["ihc_mwu_n_normal"],
                "ihc_mwu_significant": mwu["ihc_mwu_significant"],
            }
        )

    comparison = pd.DataFrame(rows)
    if comparison.empty:
        return comparison, tcga, consensus

    comparison = remap_rna_arrays(comparison)
    comparison = drop_ihc_pseudocount_fold(comparison)
    comparison["rna_minmax"] = _minmax(comparison["cancer_rna_abs"])
    comparison["normal_rna_minmax"] = _minmax(comparison["normal_rna_abs"])
    comparison["rna_zscore_cancer"] = _zscore(comparison["cancer_rna_tcga_median"])
    comparison["protein_zscore_cancer"] = _zscore(comparison["protein_ihc_cancer_0_3"])
    comparison = comparison.sort_values(
        ["rna_log2_fold", "ihc_mwu_effect"],
        ascending=False,
        na_position="last",
    ).reset_index(drop=True)

    normal_rna_table = pd.DataFrame()
    if not consensus.empty:
        normal_rna_table = consensus.rename(columns={"Tissue": "organ", "nTPM": "normal_rna_ntpm"})[
            ["Gene", "Gene name", "organ", "normal_rna_ntpm"]
        ]
        normal_rna_table["source"] = "HPA consensus / GTEx"
    return comparison, tcga, normal_rna_table


def build_tcga_summary(tcga: pd.DataFrame) -> pd.DataFrame:
    """Per-cancer TCGA summary for tables and bar charts."""
    if tcga is None or tcga.empty:
        return pd.DataFrame()

    rows: list[dict[str, Any]] = []
    for cancer_type, subset in tcga.groupby("cancer_type"):
        if pd.isna(cancer_type) or not str(cancer_type).strip():
            continue
        tumor = subset[subset["sample_class"] == "tumor"]["value"]
        normal = subset[subset["sample_class"] == "normal"]["value"]
        tumor_med = float(tumor.median()) if not tumor.empty else None
        normal_mean = float(normal.mean()) if not normal.empty else None
        rows.append(
            {
                "cancer_type": cancer_type,
                "n_tumor": int(tumor.count()),
                "n_adjacent_normal": int(normal.count()),
                "tumor_median_rsem": tumor_med,
                "normal_mean_rsem": normal_mean,
                "tumor_median_log2": math.log2(tumor_med + 1) if tumor_med is not None else None,
                "normal_mean_log2": math.log2(normal_mean + 1) if normal_mean is not None else None,
                "tcga_inverse_2pow_minus1": (
                    (2 ** math.log2(tumor_med + 1)) - 1 if tumor_med is not None else None
                ),
                "tcga_inverse_sanity_ok": (
                    sanity_check_tcga_inverse(tumor_med, math.log2(tumor_med + 1))["sanity_ok"]
                    if tumor_med is not None
                    else False
                ),
                "fold_tumor_vs_adjacent": (
                    (tumor_med + PSEUDOCOUNT) / (normal_mean + PSEUDOCOUNT)
                    if tumor_med is not None and normal_mean is not None
                    else None
                ),
                "log2_fold_tumor_vs_adjacent": (
                    math.log2((tumor_med + PSEUDOCOUNT) / (normal_mean + PSEUDOCOUNT))
                    if tumor_med is not None and normal_mean is not None
                    else None
                ),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values("tumor_median_log2", ascending=False, na_position="last").reset_index(drop=True)
