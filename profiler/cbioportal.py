# -*- coding: utf-8 -*-
"""cBioPortal adapter — studies, profiles, molecular data, clinical samples."""

from __future__ import annotations

import hashlib
import io
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from profiler.http_client import request_json, session

API = "https://www.cbioportal.org/api"
DATAHUB_SAMPLE = (
    "https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/"
    "{study}/data_clinical_sample.txt"
)
CACHE_DIR = Path(__file__).resolve().parent.parent / "cache" / "cbioportal"

# TCGA PanCancer Atlas 2018 studies used for bulk tumor expression
TCGA_STUDIES: dict[str, dict[str, str]] = {
    "BRCA": {"study": "brca_tcga_pan_can_atlas_2018", "label": "Breast invasive carcinoma"},
    "LUAD": {"study": "luad_tcga_pan_can_atlas_2018", "label": "Lung adenocarcinoma"},
    "LUSC": {"study": "lusc_tcga_pan_can_atlas_2018", "label": "Lung squamous cell carcinoma"},
    "COADREAD": {"study": "coadread_tcga_pan_can_atlas_2018", "label": "Colorectal adenocarcinoma"},
    "STAD": {"study": "stad_tcga_pan_can_atlas_2018", "label": "Stomach adenocarcinoma"},
    "OV": {"study": "ov_tcga_pan_can_atlas_2018", "label": "Ovarian serous cystadenocarcinoma"},
    "PRAD": {"study": "prad_tcga_pan_can_atlas_2018", "label": "Prostate adenocarcinoma"},
    "PAAD": {"study": "paad_tcga_pan_can_atlas_2018", "label": "Pancreatic adenocarcinoma"},
    "BLCA": {"study": "blca_tcga_pan_can_atlas_2018", "label": "Bladder urothelial carcinoma"},
    "HNSC": {"study": "hnsc_tcga_pan_can_atlas_2018", "label": "Head and neck squamous cell carcinoma"},
    "LIHC": {"study": "lihc_tcga_pan_can_atlas_2018", "label": "Liver hepatocellular carcinoma"},
    "KIRC": {"study": "kirc_tcga_pan_can_atlas_2018", "label": "Kidney renal clear cell carcinoma"},
    "GBM": {"study": "gbm_tcga_pan_can_atlas_2018", "label": "Glioblastoma multiforme"},
    "LGG": {"study": "lgg_tcga_pan_can_atlas_2018", "label": "Brain lower grade glioma"},
    "SKCM": {"study": "skcm_tcga_pan_can_atlas_2018", "label": "Skin cutaneous melanoma"},
    "UCEC": {"study": "ucec_tcga_pan_can_atlas_2018", "label": "Uterine corpus endometrial carcinoma"},
    "ESCA": {"study": "esca_tcga_pan_can_atlas_2018", "label": "Esophageal carcinoma"},
    "SARC": {"study": "sarc_tcga_pan_can_atlas_2018", "label": "Sarcoma"},
    "THCA": {"study": "thca_tcga_pan_can_atlas_2018", "label": "Thyroid carcinoma"},
    "DLBC": {"study": "dlbc_tcga_pan_can_atlas_2018", "label": "Diffuse large B-cell lymphoma"},
    "LAML": {"study": "laml_tcga_pan_can_atlas_2018", "label": "Acute myeloid leukemia"},
    "MESO": {"study": "meso_tcga_pan_can_atlas_2018", "label": "Mesothelioma"},
    "CESC": {"study": "cesc_tcga_pan_can_atlas_2018", "label": "Cervical squamous cell carcinoma"},
    "TGCT": {"study": "tgct_tcga_pan_can_atlas_2018", "label": "Testicular germ cell tumors"},
    "KIRP": {"study": "kirp_tcga_pan_can_atlas_2018", "label": "Kidney renal papillary cell carcinoma"},
    "ACC": {"study": "acc_tcga_pan_can_atlas_2018", "label": "Adrenocortical carcinoma"},
    "CHOL": {"study": "chol_tcga_pan_can_atlas_2018", "label": "Cholangiocarcinoma"},
    "UCS": {"study": "ucs_tcga_pan_can_atlas_2018", "label": "Uterine carcinosarcoma"},
    "UVM": {"study": "uvm_tcga_pan_can_atlas_2018", "label": "Uveal melanoma"},
}

CELL_LINE_ASSAYS: dict[str, dict[str, str]] = {
    "depmap_25q3_rna": {
        "study": "ccle_broad_2025",
        "profile": "ccle_broad_2025_rna_seq_mrna",
        "sample_list": "ccle_broad_2025_all",
        "unit": "log2(TPM+1)",
        "transform": "none_already_log2_tpm_plus1",
        "label": "DepMap Public 25Q3 RNA",
        "assay": "RNA-seq",
    },
    "ccle_2019_rna": {
        "study": "ccle_broad_2019",
        "profile": "ccle_broad_2019_rna_seq_mrna",
        "sample_list": "ccle_broad_2019_all",
        "unit": "RPKM",
        "plot_unit": "log2(RPKM+1)",
        "transform": "table_rpkm_plot_log2_plus1",
        "label": "CCLE 2019 RNA",
        "assay": "RNA-seq",
    },
    "ccle_2019_protein": {
        "study": "ccle_broad_2019",
        "profile": "ccle_broad_2019_protein_quantification",
        "sample_list": "ccle_broad_2019_all",
        "unit": "normalized relative log2 (whole-cell MS)",
        "transform": "none_keep_normalized",
        "label": "CCLE/Nusinow 2020 protein",
        "assay": "TMT proteomics",
    },
}


@dataclass
class MolecularFetchResult:
    study_id: str
    molecular_profile_id: str
    entrez_gene_id: int
    symbol: str
    assay_key: str
    unit: str
    transform: str
    values: list[dict[str, Any]]
    n_measured: int
    n_missing_registered: int | None
    retrieved_at: str
    status: str
    error: str | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def to_dataframe(self) -> pd.DataFrame:
        if not self.values:
            return pd.DataFrame(columns=["sampleId", "value", "entrezGeneId", "molecularProfileId", "studyId"])
        return pd.DataFrame(self.values)


def now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def get_study(study_id: str) -> dict[str, Any]:
    return request_json("GET", f"{API}/studies/{study_id}", cache_ttl_s=86400)


def list_molecular_profiles(study_id: str) -> list[dict[str, Any]]:
    return request_json("GET", f"{API}/studies/{study_id}/molecular-profiles", cache_ttl_s=86400)


def resolve_rna_profile(study_id: str) -> str | None:
    profiles = list_molecular_profiles(study_id)
    preferred = f"{study_id}_rna_seq_v2_mrna"
    ids = [p.get("molecularProfileId") for p in profiles]
    if preferred in ids:
        return preferred
    for p in profiles:
        mid = p.get("molecularProfileId") or ""
        dtype = (p.get("datatype") or "").upper()
        if "RNA" in mid.upper() and "MRNA" in mid.upper() and "ZSCORE" not in mid.upper():
            return mid
        if dtype == "CONTINUOUS" and "mrna" in mid.lower() and "zscores" not in mid.lower():
            return mid
    return None


def resolve_rppa_profile(study_id: str) -> str | None:
    profiles = list_molecular_profiles(study_id)
    for p in profiles:
        mid = (p.get("molecularProfileId") or "")
        mol = (p.get("molecularAlterationType") or "").upper()
        if "RPPA" in mid.upper() or mol == "PROTEIN_LEVEL":
            if "ZSCORE" in mid.upper() or "Z-SCORE" in mid.upper():
                continue
            return mid
    return None


def fetch_molecular_data(
    profile_id: str,
    entrez_gene_id: int,
    sample_list_id: str,
    *,
    study_id: str = "",
    symbol: str = "",
    assay_key: str = "",
    unit: str = "",
    transform: str = "",
) -> MolecularFetchResult:
    retrieved = now_iso()
    try:
        data = request_json(
            "GET",
            f"{API}/molecular-profiles/{profile_id}/molecular-data",
            params={
                "entrezGeneId": entrez_gene_id,
                "sampleListId": sample_list_id,
                "projection": "SUMMARY",
            },
            cache_ttl_s=6 * 3600,
            timeout=(12, 120),
            max_attempts=8,
        )
        if not isinstance(data, list):
            raise RuntimeError(f"Unexpected molecular-data payload type: {type(data)}")

        measured = []
        for row in data:
            val = row.get("value")
            # Preserve actual 0; treat missing/None separately. cBio usually omits NA.
            if val is None:
                continue
            measured.append(
                {
                    "sampleId": row.get("sampleId"),
                    "patientId": row.get("patientId"),
                    "value": float(val),
                    "entrezGeneId": row.get("entrezGeneId"),
                    "molecularProfileId": row.get("molecularProfileId") or profile_id,
                    "studyId": row.get("studyId") or study_id,
                    "uniqueSampleKey": row.get("uniqueSampleKey"),
                }
            )

        # Validate that returned profile/gene match request
        mismatches = [
            r for r in measured
            if r.get("entrezGeneId") not in (None, entrez_gene_id)
            or (r.get("molecularProfileId") and r["molecularProfileId"] != profile_id)
        ]
        status = "ok" if measured else "no_measured_values"
        meta = {
            "requested_profile": profile_id,
            "requested_entrez": entrez_gene_id,
            "returned_n": len(measured),
            "mismatch_n": len(mismatches),
            "checksum": hashlib.sha256(
                "|".join(f"{r['sampleId']}:{r['value']}" for r in sorted(measured, key=lambda x: x["sampleId"] or "")).encode()
            ).hexdigest()[:16],
        }
        return MolecularFetchResult(
            study_id=study_id,
            molecular_profile_id=profile_id,
            entrez_gene_id=entrez_gene_id,
            symbol=symbol,
            assay_key=assay_key,
            unit=unit,
            transform=transform,
            values=measured,
            n_measured=len(measured),
            n_missing_registered=None,
            retrieved_at=retrieved,
            status=status,
            metadata=meta,
        )
    except Exception as exc:
        return MolecularFetchResult(
            study_id=study_id,
            molecular_profile_id=profile_id,
            entrez_gene_id=entrez_gene_id,
            symbol=symbol,
            assay_key=assay_key,
            unit=unit,
            transform=transform,
            values=[],
            n_measured=0,
            n_missing_registered=None,
            retrieved_at=retrieved,
            status="fetch_failed",
            error=str(exc),
        )


def load_clinical_samples(study_id: str) -> pd.DataFrame:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    dest = CACHE_DIR / f"{study_id}_clinical_sample.tsv"
    if dest.exists() and dest.stat().st_size > 500:
        age = time.time() - dest.stat().st_mtime
        if age < 7 * 86400:
            return _read_clinical(dest)

    url = DATAHUB_SAMPLE.format(study=study_id)
    sess = session()
    resp = sess.get(url, timeout=(12, 120))
    if resp.status_code != 200:
        # Fallback: cBioPortal clinical API
        try:
            data = request_json("GET", f"{API}/studies/{study_id}/samples", cache_ttl_s=86400)
            df = pd.DataFrame(data)
            return df
        except Exception as exc:
            raise RuntimeError(f"Could not load clinical samples for {study_id}: {exc}") from exc
    dest.write_bytes(resp.content)
    return _read_clinical(dest)


def _read_clinical(path: Path) -> pd.DataFrame:
    text = path.read_text(encoding="utf-8", errors="replace")
    # Skip meta header lines starting with #
    lines = [ln for ln in text.splitlines() if ln and not ln.startswith("#")]
    if not lines:
        return pd.DataFrame()
    return pd.read_csv(io.StringIO("\n".join(lines)), sep="\t", dtype=str)


def sample_list_id_for_study(study_id: str) -> str:
    return f"{study_id}_all"
