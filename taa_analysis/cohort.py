# -*- coding: utf-8 -*-
"""TCGA labels, specimen modes, and one-row-per-patient selection."""

from __future__ import annotations

import pandas as pd

DISEASE_CODE = {
    "Acute Myeloid Leukemia": "LAML",
    "Adrenocortical Cancer": "ACC",
    "Bladder Urothelial Carcinoma": "BLCA",
    "Brain Lower Grade Glioma": "LGG",
    "Breast Invasive Carcinoma": "BRCA",
    "Cervical & Endocervical Cancer": "CESC",
    "Cholangiocarcinoma": "CHOL",
    "Colon Adenocarcinoma": "COAD",
    "Diffuse Large B-Cell Lymphoma": "DLBC",
    "Esophageal Carcinoma": "ESCA",
    "Glioblastoma Multiforme": "GBM",
    "Head & Neck Squamous Cell Carcinoma": "HNSC",
    "Kidney Chromophobe": "KICH",
    "Kidney Clear Cell Carcinoma": "KIRC",
    "Kidney Papillary Cell Carcinoma": "KIRP",
    "Liver Hepatocellular Carcinoma": "LIHC",
    "Lung Adenocarcinoma": "LUAD",
    "Lung Squamous Cell Carcinoma": "LUSC",
    "Mesothelioma": "MESO",
    "Ovarian Serous Cystadenocarcinoma": "OV",
    "Pancreatic Adenocarcinoma": "PAAD",
    "Pheochromocytoma & Paraganglioma": "PCPG",
    "Prostate Adenocarcinoma": "PRAD",
    "Rectum Adenocarcinoma": "READ",
    "Sarcoma": "SARC",
    "Skin Cutaneous Melanoma": "SKCM",
    "Stomach Adenocarcinoma": "STAD",
    "Testicular Germ Cell Tumor": "TGCT",
    "Thymoma": "THYM",
    "Thyroid Carcinoma": "THCA",
    "Uterine Carcinosarcoma": "UCS",
    "Uterine Corpus Endometrioid Carcinoma": "UCEC",
    "Uveal Melanoma": "UVM",
}

HEMATOLOGIC_DISEASES = {
    "Acute Myeloid Leukemia",
    "Diffuse Large B-Cell Lymphoma",
}

SOLID_PRIMARY_TYPES = {"Primary Tumor", "Primary Solid Tumor"}
METASTATIC_TYPES = {"Metastatic", "Additional Metastatic"}
NORMAL_TYPES = {"Solid Tissue Normal", "Normal Tissue"}
CELL_LINE_TYPES = {"Cell Line"}
BLOOD_CANCER_TYPES = {
    "Primary Blood Derived Cancer - Peripheral Blood",
    "Primary Blood Derived Cancer - Bone Marrow",
    "Recurrent Blood Derived Cancer - Bone Marrow",
    "Recurrent Blood Derived Cancer - Peripheral Blood",
    "Post treatment Blood Cancer - Bone Marrow",
    "Post treatment Blood Cancer - Blood",
}

PATIENT_RULE = (
    "One row per patient. Within the selected specimen class, keep a single sample. "
    "Prefer a sample that has an original subtype annotation when a subtype column is supplied; "
    "otherwise keep the lexicographically smallest sample ID. Weight is 1. "
    "Annotations are not copied onto a different sample, site, or time point."
)

MODES = (
    "solid_primary",
    "metastatic",
    "hematologic",
    "normal_reference",
    "cell_line",
)


def patient_id(sample_id: str, study: str) -> str:
    text = str(sample_id or "")
    parts = text.split("-")
    if study == "TCGA" and len(parts) >= 3:
        return "-".join(parts[:3])
    if study == "GTEX" and len(parts) >= 2:
        return "-".join(parts[:2])
    return text


def specimen_class(sample_type: str, disease: str) -> str:
    if disease in HEMATOLOGIC_DISEASES or sample_type in BLOOD_CANCER_TYPES:
        return "hematologic"
    if sample_type in SOLID_PRIMARY_TYPES:
        return "solid_primary"
    if sample_type in METASTATIC_TYPES:
        return "metastatic"
    if sample_type in NORMAL_TYPES:
        return "normal_reference"
    if sample_type in CELL_LINE_TYPES:
        return "cell_line"
    return "other"


def select_one_sample_per_patient(frame: pd.DataFrame, annotation_column: str | None = None) -> pd.DataFrame:
    """Collapse repeated samples. Raises if the result has more rows than patients."""
    if frame.empty:
        out = frame.copy()
        out["aggregation_weight"] = []
        out["samples_collapsed"] = []
        return out
    if frame["patient_id"].duplicated().any() is False and annotation_column is None:
        out = frame.copy()
        out["aggregation_weight"] = 1.0
        out["samples_collapsed"] = 1
        return out

    chosen_rows = []
    collapsed = []
    for patient, group in frame.groupby("patient_id", sort=False):
        pool = group
        if annotation_column and annotation_column in group.columns:
            labels = group[annotation_column].astype("string")
            annotated = group[labels.notna() & ~labels.str.lower().isin(["", "nan", "none", "unknown"])]
            if not annotated.empty:
                pool = annotated
        pool = pool.sort_values("sample_id")
        chosen_rows.append(pool.iloc[0])
        collapsed.append(int(len(group)))
        _ = patient
    out = pd.DataFrame(chosen_rows).reset_index(drop=True)
    out["aggregation_weight"] = 1.0
    out["samples_collapsed"] = collapsed
    if len(out) > frame["patient_id"].nunique():
        raise ValueError("patient collapse increased rows")
    if out["patient_id"].duplicated().any():
        raise ValueError("duplicate patients remain after collapse")
    return out
