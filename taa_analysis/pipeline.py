# -*- coding: utf-8 -*-
"""On-demand TAA views. Numbers without a retrieved source are not returned."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from taa_analysis import ANALYSIS_VERSION
from taa_analysis.aliases import canonical_symbol, target_identity
from taa_analysis.cohort import (
    DISEASE_CODE,
    HEMATOLOGIC_DISEASES,
    MODES,
    PATIENT_RULE,
    patient_id,
    select_one_sample_per_patient,
    specimen_class,
)
from taa_analysis.composition import all_cell_percentage_status
from taa_analysis.joins import safe_left_join
from taa_analysis.molecules import MOLECULE_PANEL, load_hpa_immune_cell_table
from taa_analysis.settings import load_settings
from taa_analysis.transforms import SOURCE_LOG2_TPM_0001, convert_stored_expression
from taa_analysis.xena import RetrievalFailed, gene_probe, xena_columns

TOIL_TPM_DATASET = "TcgaTargetGtex_rsem_gene_tpm"
TOIL_TPM_VERSION = "2016-09-03"
TOIL_TPM_UNIT = SOURCE_LOG2_TPM_0001
TOIL_PROBEMAP = "probeMap/gencode.v23.annotation.gene.probemap"
TOIL_PHENOTYPE = "TcgaTargetGTEX_phenotype.txt"
TOIL_PHENOTYPE_VERSION = "2016-09-15"
BRCA_CLINICAL = "TCGA.BRCA.sampleMap/BRCA_clinicalMatrix"
PAM50_FIELD = "PAM50Call_RNAseq"
CLINICAL_HER2_FIELD = "HER2_Final_Status_nature2012"

PAM50_LABELS = {
    "LumA": "Luminal A",
    "LumB": "Luminal B",
    "Her2": "HER2-enriched",
    "Basal": "Basal-like",
    "Normal": "Normal-like",
}

SUBTYPE_SYSTEM_PAM50 = "PAM50"
SUBTYPE_SYSTEM_CLINICAL_HER2 = "Clinical HER2 final status (nature2012)"


def _now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _cache_dir() -> Path:
    path = Path(load_settings()["TAA_CACHE_DIR"])
    path.mkdir(parents=True, exist_ok=True)
    return path


def _settings() -> dict[str, str]:
    return load_settings()


def load_phenotype(force: bool = False) -> pd.DataFrame:
    settings = _settings()
    dest = _cache_dir() / "toil_phenotype.csv"
    if dest.exists() and not force:
        frame = pd.read_csv(dest)
    else:
        payload = xena_columns(
            settings["TAA_XENA_TOIL_HUB"],
            TOIL_PHENOTYPE,
            [
                "sampleID",
                "_study",
                "_sample_type",
                "_primary_site",
                "primary disease or tissue",
                "detailed_category",
            ],
        )
        frame = pd.DataFrame(payload)
        frame.to_csv(dest, index=False)
    frame = frame.rename(
        columns={
            "sampleID": "sample_id",
            "_study": "study",
            "_sample_type": "sample_type",
            "_primary_site": "primary_site",
            "primary disease or tissue": "disease",
            "detailed_category": "detailed_category",
        }
    )
    frame["sample_id"] = frame["sample_id"].astype(str)
    frame["study"] = frame["study"].fillna("").astype(str)
    frame["sample_type"] = frame["sample_type"].fillna("").astype(str)
    frame["disease"] = frame["disease"].fillna("").astype(str)
    frame["patient_id"] = [
        patient_id(sample, study) for sample, study in zip(frame["sample_id"], frame["study"])
    ]
    frame["specimen_class"] = [
        specimen_class(sample_type, disease)
        for sample_type, disease in zip(frame["sample_type"], frame["disease"])
    ]
    frame["cancer_code"] = frame["disease"].map(DISEASE_CODE)
    frame["lineage"] = frame["disease"].map(
        lambda disease: "hematologic" if disease in HEMATOLOGIC_DISEASES else "solid"
    )
    return frame


def load_gene_slice(gene_query: str, force: bool = False) -> dict:
    """Bulk expression slice plus phenotype. Immune and single-cell matrices are not downloaded."""
    symbol = canonical_symbol(gene_query)
    identity = target_identity(gene_query)
    identity["gene_symbol"] = symbol
    settings = _settings()
    dest = _cache_dir() / f"expr_{symbol}.csv"
    meta_path = _cache_dir() / f"expr_{symbol}.json"
    if dest.exists() and meta_path.exists() and not force:
        values = pd.read_csv(dest)
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
    else:
        hub = settings["TAA_XENA_TOIL_HUB"]
        probe = gene_probe(hub, TOIL_PROBEMAP, symbol)
        payload = xena_columns(hub, TOIL_TPM_DATASET, ["sampleID", probe])
        values = pd.DataFrame(
            {
                "sample_id": payload.get("sampleID") or [],
                "source_value": payload.get(probe) or [],
            }
        )
        if values.empty:
            raise RetrievalFailed(TOIL_TPM_DATASET, f"no expression rows for {symbol}")
        values["sample_id"] = values["sample_id"].astype(str)
        values["probe_id"] = probe
        converted = values["source_value"].map(
            lambda value: _convert_or_missing(value)
        )
        values["tpm"] = [item["tpm"] for item in converted]
        values["display_value"] = [item["display_value"] for item in converted]
        values["transform_error"] = [item["error"] for item in converted]
        values.to_csv(dest, index=False)
        digest = hashlib.sha256(values["source_value"].astype(str).str.cat(sep=",").encode()).hexdigest()[:16]
        meta = {
            "gene_symbol": symbol,
            "gene_id": identity.get("gene_id", ""),
            "probe_id": probe,
            "source_name": "UCSC Xena Toil",
            "source_url": hub,
            "dataset_id": TOIL_TPM_DATASET,
            "dataset_version": TOIL_TPM_VERSION,
            "unit_source": TOIL_TPM_UNIT,
            "unit_display": "log2(TPM+1)",
            "formula": "TPM = 2^x - 0.001; display = log2(TPM + 1)",
            "input_checksum": digest,
            "retrieved_at": _now(),
            "analysis_version": ANALYSIS_VERSION,
            "n_rows": int(len(values)),
        }
        meta_path.write_text(json.dumps(meta), encoding="utf-8")
    phenotype = load_phenotype(force=force)
    merged = values.merge(phenotype, on="sample_id", how="left", validate="one_to_one")
    if len(merged) != len(values):
        raise RetrievalFailed(TOIL_TPM_DATASET, "phenotype join changed expression row count")
    return {"identity": identity, "meta": meta, "samples": merged}


def _convert_or_missing(value) -> dict:
    try:
        if pd.isna(value):
            return {"tpm": None, "display_value": None, "error": ""}
    except TypeError:
        pass
    try:
        converted = convert_stored_expression(value, TOIL_TPM_UNIT)
    except Exception as exc:
        return {"tpm": None, "display_value": None, "error": str(exc)}
    return {"tpm": converted["tpm"], "display_value": converted["display_value"], "error": ""}


def filter_samples(samples: pd.DataFrame, mode: str) -> pd.DataFrame:
    if mode not in MODES:
        raise ValueError(f"unknown specimen mode {mode}")
    work = samples[samples["specimen_class"] == mode].copy()
    if mode in {"solid_primary", "metastatic", "hematologic"}:
        work = work[work["study"] == "TCGA"]
        work = work[work["cancer_code"].notna()]
    elif mode == "normal_reference":
        work = work[work["study"] == "GTEX"]
    elif mode == "cell_line":
        work = work[work["study"] == "TCGA"]
    return work


def _group_key(mode: str) -> str:
    if mode == "normal_reference":
        return "primary_site"
    if mode == "cell_line":
        return "disease"
    return "cancer_code"


def patient_expression(samples: pd.DataFrame, mode: str, annotation: pd.DataFrame | None = None) -> pd.DataFrame:
    grouped = filter_samples(samples, mode)
    key = _group_key(mode)
    if mode != "normal_reference":
        grouped = grouped[grouped[key].astype(str).str.len() > 0]
    else:
        grouped = grouped[grouped[key].fillna("").astype(str).str.len() > 0]
        grouped["cancer_code"] = grouped[key]
    annotation_column = None
    if annotation is not None and not annotation.empty:
        grouped = safe_left_join(grouped, annotation, "sample_id")
        grouped["subtype_label"] = grouped["subtype_label"].fillna("Unknown")
        grouped["subtype_system"] = grouped["subtype_system"].fillna("")
        grouped["annotation_status"] = grouped["annotation_status"].fillna("original")
        annotation_column = "subtype_label"
    before_counts = grouped.groupby(key)["sample_id"].nunique().to_dict()
    patients = select_one_sample_per_patient(grouped, annotation_column)
    patients["group_label"] = patients[key]
    patients["n_samples_in_group"] = patients[key].map(before_counts)
    patients["patient_rule"] = PATIENT_RULE
    return patients


def summarize_expression(patients: pd.DataFrame, mode: str, meta: dict) -> pd.DataFrame:
    key = "cancer_code" if mode != "cell_line" else "disease"
    if mode == "cell_line":
        key = "disease"
    rows = []
    if patients.empty:
        return pd.DataFrame(rows)
    group_column = "primary_site" if mode == "normal_reference" else ("disease" if mode == "cell_line" else "cancer_code")
    for label, group in patients.groupby(group_column, dropna=False):
        values = pd.to_numeric(group["display_value"], errors="coerce").dropna()
        errors = int((group["transform_error"].fillna("") != "").sum()) if "transform_error" in group else 0
        missing = int(pd.to_numeric(group["display_value"], errors="coerce").isna().sum())
        row = {
            "cancer": label if mode != "solid_primary" else _disease_for_code(group, label),
            "cancer_code": label if mode not in {"normal_reference", "cell_line"} else "",
            "group_label": label,
            "n_patients": int(group["patient_id"].nunique()),
            "n_samples": int(group["samples_collapsed"].sum()) if "samples_collapsed" in group else int(len(group)),
            "median": float(values.median()) if not values.empty else None,
            "q1": float(values.quantile(0.25)) if len(values) else None,
            "q3": float(values.quantile(0.75)) if len(values) else None,
            "unit": meta.get("unit_display"),
            "source": meta.get("source_name"),
            "dataset_id": meta.get("dataset_id"),
            "dataset_version": meta.get("dataset_version"),
            "missing_n": missing,
            "transform_error_n": errors,
            "distribution": "boxplot" if len(values) >= 2 else ("point" if len(values) == 1 else "unavailable"),
            "evidence_status": "Observed",
            "specimen_mode": mode,
            "denominator": "bulk tumor RNA, one patient" if mode.startswith("solid") or mode == "metastatic" else mode,
        }
        if mode == "normal_reference":
            row["cancer"] = label
            row["denominator"] = "GTEx normal tissue, one donor"
            row["evidence_status"] = "Reference only"
        if mode == "cell_line":
            row["denominator"] = "cell line RNA, not patient tumor"
        rows.append(row)
    summary = pd.DataFrame(rows)
    if not summary.empty and summary["median"].notna().any():
        summary = summary.sort_values("median", ascending=False, na_position="last")
    return summary.reset_index(drop=True)


def _disease_for_code(group: pd.DataFrame, code: str) -> str:
    if "disease" in group.columns and group["disease"].astype(str).str.len().gt(0).any():
        return str(group["disease"].dropna().astype(str).iloc[0])
    return str(code)


def load_brca_subtype_annotations(force: bool = False) -> pd.DataFrame:
    dest = _cache_dir() / "brca_subtype.csv"
    if dest.exists() and not force:
        frame = pd.read_csv(dest)
    else:
        settings = _settings()
        payload = xena_columns(
            settings["TAA_XENA_TCGA_HUB"],
            BRCA_CLINICAL,
            ["sampleID", PAM50_FIELD, CLINICAL_HER2_FIELD],
        )
        frame = pd.DataFrame(payload).rename(
            columns={
                "sampleID": "sample_id",
                PAM50_FIELD: "pam50_raw",
                CLINICAL_HER2_FIELD: "clinical_her2_raw",
            }
        )
        frame.to_csv(dest, index=False)
    frame["sample_id"] = frame["sample_id"].astype(str)
    return frame


def subtype_annotation_table(system: str) -> pd.DataFrame:
    """Original BRCA annotations. PAM50 is not mapped onto clinical HER2."""
    if system not in {SUBTYPE_SYSTEM_PAM50, SUBTYPE_SYSTEM_CLINICAL_HER2}:
        return pd.DataFrame(columns=["sample_id", "subtype_system", "subtype_label", "subtype_raw", "subtype_source", "annotation_status"])
    source = load_brca_subtype_annotations()
    if system == SUBTYPE_SYSTEM_PAM50:
        labels = source["pam50_raw"].map(lambda value: PAM50_LABELS.get(value, "Unknown") if pd.notna(value) and str(value) not in {"", "nan", "None"} else "Unknown")
        raw = source["pam50_raw"]
        field = PAM50_FIELD
    else:
        raw = source["clinical_her2_raw"]
        labels = raw.map(lambda value: str(value) if pd.notna(value) and str(value) not in {"", "nan", "None"} else "Unknown")
        field = CLINICAL_HER2_FIELD
    return pd.DataFrame(
        {
            "sample_id": source["sample_id"],
            "subtype_system": system,
            "subtype_label": labels,
            "subtype_raw": raw,
            "subtype_source": f"https://tcga.xenahubs.net {BRCA_CLINICAL} {field}",
            "annotation_status": "original",
        }
    )


def available_subtype_systems(cancer_code: str) -> list[str]:
    if cancer_code == "BRCA":
        return [SUBTYPE_SYSTEM_PAM50, SUBTYPE_SYSTEM_CLINICAL_HER2]
    return []


def summarize_subtypes(patients: pd.DataFrame, meta: dict) -> pd.DataFrame:
    if patients.empty or "subtype_label" not in patients.columns:
        return pd.DataFrame()
    rows = []
    unclassified = int((patients["subtype_label"] == "Unknown").sum())
    for (system, label), group in patients.groupby(["subtype_system", "subtype_label"], dropna=False):
        values = pd.to_numeric(group["display_value"], errors="coerce").dropna()
        q1 = float(values.quantile(0.25)) if len(values) else None
        q3 = float(values.quantile(0.75)) if len(values) else None
        rows.append(
            {
                "subtype_system": system,
                "subtype": label,
                "n": int(group["patient_id"].nunique()),
                "median": float(values.median()) if len(values) else None,
                "q1": q1,
                "q3": q3,
                "iqr": None if q1 is None or q3 is None else q3 - q1,
                "source": group["subtype_source"].iloc[0] if "subtype_source" in group else meta.get("source_name"),
                "unclassified_n": unclassified,
                "unit": meta.get("unit_display"),
                "distribution": "violin_boxplot" if len(values) >= 2 else ("point" if len(values) == 1 else "unavailable"),
                "annotation_status": "original",
                "evidence_status": "Observed",
            }
        )
    summary = pd.DataFrame(rows)
    if not summary.empty and summary["median"].notna().any():
        summary = summary.sort_values("median", ascending=False, na_position="last")
    return summary.reset_index(drop=True)


def immune_composition_result(scope: str) -> dict:
    """No verified fraction matrix is connected. Signature scores are not shown as percents."""
    status = all_cell_percentage_status("CIBERSORT_LM22_relative")
    return {
        "scope": scope,
        "status": "Unavailable",
        "missing_reason": "unsupported",
        "algorithm": "",
        "evidence_status": "Unavailable",
        "detail": (
            "Whole-tumor immune percent needs quanTIseq or a suitable EPIC reference, and within-immune "
            "composition needs a verified CIBERSORT LM22 (or equivalent) fraction table. Neither table is "
            "connected. PanCan Atlas immune signature scores are not converted into percentages. "
            "No residual is labeled as cancer cells, and displayed medians are not forced to 100%."
        ),
        "all_cell_note": status["detail"],
        "table": pd.DataFrame(
            columns=[
                "cancer",
                "subtype",
                "cell_type",
                "all_cell_percent",
                "within_immune_percent",
                "n",
                "algorithm",
                "denominator",
                "missing_reason",
            ]
        ),
    }


def tumor_molecule_result(cancer_code: str, subtype: str, cell_type: str) -> dict:
    return {
        "status": "Unavailable",
        "missing_reason": "unavailable",
        "evidence_status": "Unavailable",
        "detail": (
            f"No eligible tumor single-cell matrix is connected for {cancer_code or 'unselected cancer'} / "
            f"{subtype or 'no subtype'} / {cell_type or 'no cell type'}. "
            "Healthy-blood HPA values are not substituted. Protein evidence is unavailable."
        ),
        "table": pd.DataFrame(),
    }


def healthy_blood_reference(genes: list[str] | None = None) -> pd.DataFrame:
    settings = _settings()
    panel_genes = genes or [item["gene_symbol"] for item in MOLECULE_PANEL]
    table = load_hpa_immune_cell_table(settings["TAA_HPA_IMMUNE_CELL_URL"], _cache_dir(), panel_genes)
    molecule_by_gene = {}
    for item in MOLECULE_PANEL:
        molecule_by_gene.setdefault(item["gene_symbol"], item["molecule"])
    table["molecule"] = table["Gene name"].map(molecule_by_gene)
    table["gene_symbol"] = table["Gene name"]
    table["cell_type"] = table["Immune cell"]
    table["rna_summary_nTPM"] = pd.to_numeric(table["nTPM"], errors="coerce")
    table["distribution"] = "point"
    table["analysis_version"] = ANALYSIS_VERSION
    return table


def coverage_table(samples: pd.DataFrame) -> pd.DataFrame:
    rows = []
    tcga = samples[(samples["study"] == "TCGA") & (samples["cancer_code"].notna())]
    for disease, code in sorted(DISEASE_CODE.items(), key=lambda item: item[1]):
        subset = tcga[tcga["cancer_code"] == code]
        primary = subset[subset["specimen_class"] == "solid_primary"]
        if disease in HEMATOLOGIC_DISEASES:
            rna_status = "available" if not subset[subset["specimen_class"] == "hematologic"].empty else "unavailable"
            rna_reason = "hematologic mode only; excluded from the default solid-tumor view"
        elif primary.empty:
            rna_status = "unavailable"
            rna_reason = "no TCGA primary solid-tumor sample in the Toil matrix"
        else:
            rna_status = "available"
            rna_reason = ""
        subtype_status = "available" if code == "BRCA" else "unavailable"
        subtype_reason = "" if code == "BRCA" else "subtype_missing"
        rows.append(
            {
                "cancer_code": code,
                "cancer": disease,
                "rna_status": rna_status,
                "rna_reason": rna_reason,
                "n_primary_samples": int(len(primary)),
                "subtype_status": subtype_status,
                "subtype_reason": subtype_reason,
                "immune_status": "unavailable",
                "immune_reason": "unsupported",
                "tumor_scrna_status": "unavailable",
                "tumor_scrna_reason": "unavailable",
            }
        )
    return pd.DataFrame(rows)


def expression_bundle(gene_query: str, mode: str = "solid_primary", cancer_code: str | None = None, subtype_system: str | None = None, force: bool = False) -> dict:
    """One filtered object shared by the summary table, the figure, and CSV export."""
    if mode == "cell_line":
        from taa_analysis.cell_lines import load_cell_line_result

        loaded_lines = load_cell_line_result(gene_query, force=force)
        subtype_summary = pd.DataFrame()
        if cancer_code:
            subtype_summary = pd.DataFrame(
                [
                    {
                        "subtype_system": "",
                        "subtype": "",
                        "n": 0,
                        "median": None,
                        "q1": None,
                        "q3": None,
                        "iqr": None,
                        "source": "",
                        "unclassified_n": None,
                        "unit": loaded_lines["meta"].get("unit_display"),
                        "distribution": "unavailable",
                        "annotation_status": "",
                        "evidence_status": "Unavailable",
                        "missing_reason": "subtype_missing",
                    }
                ]
            )
        return {
            "identity": loaded_lines["identity"],
            "meta": loaded_lines["meta"],
            "mode": mode,
            "cancer_code": cancer_code,
            "subtype_system": subtype_system,
            "patients": loaded_lines["lines"],
            "summary": loaded_lines["summary"],
            "subtype_patients": pd.DataFrame(),
            "subtype_summary": subtype_summary,
            "immune": immune_composition_result(cancer_code or "cell lines"),
            "molecules": tumor_molecule_result(cancer_code or "", subtype_system or "", ""),
            "coverage": loaded_lines["coverage"],
        }
    loaded = load_gene_slice(gene_query, force=force)
    patients = patient_expression(loaded["samples"], mode)
    summary = summarize_expression(patients, mode, loaded["meta"])
    subtype_patients = pd.DataFrame()
    subtype_summary = pd.DataFrame()
    if cancer_code and subtype_system and cancer_code == "BRCA":
        focused = loaded["samples"]
        focused = focused[focused["cancer_code"] == cancer_code]
        subtype_patients = patient_expression(focused, mode, subtype_annotation_table(subtype_system))
        subtype_summary = summarize_subtypes(subtype_patients, loaded["meta"])
    elif cancer_code and cancer_code != "BRCA":
        subtype_summary = pd.DataFrame(
            [
                {
                    "subtype_system": "",
                    "subtype": "",
                    "n": 0,
                    "median": None,
                    "q1": None,
                    "q3": None,
                    "iqr": None,
                    "source": "",
                    "unclassified_n": None,
                    "unit": loaded["meta"].get("unit_display"),
                    "distribution": "unavailable",
                    "annotation_status": "",
                    "evidence_status": "Unavailable",
                    "missing_reason": "subtype_missing",
                }
            ]
        )
    return {
        "identity": loaded["identity"],
        "meta": loaded["meta"],
        "mode": mode,
        "cancer_code": cancer_code,
        "subtype_system": subtype_system,
        "patients": patients,
        "summary": summary,
        "subtype_patients": subtype_patients,
        "subtype_summary": subtype_summary,
        "immune": immune_composition_result(cancer_code or "all cancers"),
        "molecules": tumor_molecule_result(cancer_code or "", subtype_system or "", ""),
        "coverage": coverage_table(loaded["samples"]),
    }


def export_frame(frame: pd.DataFrame) -> str:
    return frame.to_csv(index=False)
