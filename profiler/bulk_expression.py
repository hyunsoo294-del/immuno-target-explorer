# -*- coding: utf-8 -*-
"""Cancer-type / subtype bulk RNA & RPPA expression from TCGA PanCancer Atlas."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from profiler.cbioportal import (
    TCGA_STUDIES,
    fetch_molecular_data,
    load_clinical_samples,
    resolve_rna_profile,
    resolve_rppa_profile,
    sample_list_id_for_study,
)


def _stats(values: list[float]) -> dict[str, Any]:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return {
            "n": 0,
            "median": None,
            "mean": None,
            "q1": None,
            "q3": None,
            "min": None,
            "max": None,
        }
    return {
        "n": int(arr.size),
        "median": float(np.median(arr)),
        "mean": float(np.mean(arr)),
        "q1": float(np.percentile(arr, 25)),
        "q3": float(np.percentile(arr, 75)),
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
    }


def _subtype_column(clinical: pd.DataFrame) -> str | None:
    for col in clinical.columns:
        cu = col.upper()
        if cu in ("SUBTYPE", "CANCER_TYPE_DETAILED", "ONCOTREE_CODE", "HISTOLOGY", "TUMOR_TYPE"):
            return col
    return None


def fetch_cancer_expression(
    symbol: str,
    entrez_gene_id: int,
    *,
    cancer_keys: list[str] | None = None,
    include_rppa: bool = True,
) -> dict[str, Any]:
    """Fetch per-cancer RNA (and optional RPPA) distributions for one gene."""
    keys = cancer_keys or list(TCGA_STUDIES.keys())
    rows: list[dict[str, Any]] = []
    point_rows: list[dict[str, Any]] = []
    notes = [
        "Bulk tumor expression is mixed-tissue signal, not pure tumor-cell expression.",
        "Cross-cancer differences may reflect batch/composition; not absolute expression.",
        "RNA values from TCGA PanCancer Atlas RNA-seq v2 mRNA profiles.",
        "Displayed RNA uses original profile scale; log2(value+1) is also reported when raw > 0.",
        "Clinical HER2+/IHC score/HER2-enriched subtype are not interchangeable.",
    ]

    for key in keys:
        meta = TCGA_STUDIES.get(key)
        if not meta:
            continue
        study = meta["study"]
        rna_profile = resolve_rna_profile(study)
        if not rna_profile:
            rows.append(
                {
                    "cancer_key": key,
                    "cancer_label": meta["label"],
                    "study_id": study,
                    "assay": "RNA",
                    "status": "profile_missing",
                    "n": 0,
                }
            )
            continue

        rna = fetch_molecular_data(
            rna_profile,
            entrez_gene_id,
            sample_list_id_for_study(study),
            study_id=study,
            symbol=symbol,
            assay_key="tcga_rna",
            unit="RNA-seq RSEM (profile native)",
            transform="also_report_log2_plus1",
        )
        if rna.status == "fetch_failed":
            rows.append(
                {
                    "cancer_key": key,
                    "cancer_label": meta["label"],
                    "study_id": study,
                    "molecular_profile_id": rna_profile,
                    "assay": "RNA",
                    "status": "fetch_failed",
                    "error": rna.error,
                    "n": 0,
                }
            )
            continue

        values = [v["value"] for v in rna.values]
        # Detect whether already log-like (typical RSEM can be large)
        log2_vals = [float(np.log2(v + 1.0)) for v in values if v is not None and v >= 0]
        st = _stats(values)
        st_log = _stats(log2_vals)
        rows.append(
            {
                "cancer_key": key,
                "cancer_label": meta["label"],
                "study_id": study,
                "molecular_profile_id": rna_profile,
                "assay": "RNA",
                "status": rna.status,
                "unit": rna.unit,
                "n": st["n"],
                "median": st["median"],
                "mean": st["mean"],
                "q1": st["q1"],
                "q3": st["q3"],
                "median_log2_plus1": st_log["median"],
                "retrieved_at": rna.retrieved_at,
                "checksum": rna.metadata.get("checksum"),
            }
        )
        for v in rna.values:
            point_rows.append(
                {
                    "cancer_key": key,
                    "cancer_label": meta["label"],
                    "study_id": study,
                    "sampleId": v["sampleId"],
                    "assay": "RNA",
                    "value": v["value"],
                    "value_log2_plus1": float(np.log2(v["value"] + 1.0)) if v["value"] >= 0 else None,
                    "unit": rna.unit,
                }
            )

        # Subtype breakdown when clinical annotation exists
        try:
            clinical = load_clinical_samples(study)
            sub_col = _subtype_column(clinical) if not clinical.empty else None
            if sub_col and "SAMPLE_ID" in [c.upper() for c in clinical.columns]:
                sid_col = next(c for c in clinical.columns if c.upper() == "SAMPLE_ID")
                clin_map = dict(zip(clinical[sid_col].astype(str), clinical[sub_col].astype(str)))
                by_sub: dict[str, list[float]] = {}
                for v in rna.values:
                    sub = clin_map.get(str(v["sampleId"]), "Unknown / missing")
                    if not sub or sub.lower() in ("nan", "none", ""):
                        sub = "Unknown / missing"
                    by_sub.setdefault(sub, []).append(v["value"])
                for sub, vals in sorted(by_sub.items(), key=lambda x: -len(x[1])):
                    sst = _stats(vals)
                    rows.append(
                        {
                            "cancer_key": key,
                            "cancer_label": meta["label"],
                            "subtype": sub,
                            "study_id": study,
                            "molecular_profile_id": rna_profile,
                            "assay": "RNA",
                            "status": "ok",
                            "unit": rna.unit,
                            "n": sst["n"],
                            "median": sst["median"],
                            "q1": sst["q1"],
                            "q3": sst["q3"],
                            "level": "subtype",
                        }
                    )
        except Exception:
            pass

        if include_rppa:
            rppa_profile = resolve_rppa_profile(study)
            if not rppa_profile:
                rows.append(
                    {
                        "cancer_key": key,
                        "cancer_label": meta["label"],
                        "study_id": study,
                        "assay": "RPPA",
                        "status": "not_measured",
                        "n": 0,
                        "note": "No PROTEIN_LEVEL/RPPA profile for this study",
                    }
                )
            else:
                rppa = fetch_molecular_data(
                    rppa_profile,
                    entrez_gene_id,
                    sample_list_id_for_study(study),
                    study_id=study,
                    symbol=symbol,
                    assay_key="tcga_rppa",
                    unit="RPPA normalized",
                    transform="none",
                )
                if rppa.status == "fetch_failed":
                    rows.append(
                        {
                            "cancer_key": key,
                            "cancer_label": meta["label"],
                            "study_id": study,
                            "molecular_profile_id": rppa_profile,
                            "assay": "RPPA",
                            "status": "fetch_failed",
                            "error": rppa.error,
                            "n": 0,
                        }
                    )
                elif rppa.n_measured == 0:
                    rows.append(
                        {
                            "cancer_key": key,
                            "cancer_label": meta["label"],
                            "study_id": study,
                            "molecular_profile_id": rppa_profile,
                            "assay": "RPPA",
                            "status": "not_measured",
                            "n": 0,
                            "note": "Target absent from this protein profile (not equivalent to zero expression)",
                        }
                    )
                else:
                    st = _stats([v["value"] for v in rppa.values])
                    rows.append(
                        {
                            "cancer_key": key,
                            "cancer_label": meta["label"],
                            "study_id": study,
                            "molecular_profile_id": rppa_profile,
                            "assay": "RPPA",
                            "status": "ok",
                            "unit": rppa.unit,
                            "n": st["n"],
                            "median": st["median"],
                            "q1": st["q1"],
                            "q3": st["q3"],
                            "retrieved_at": rppa.retrieved_at,
                        }
                    )

    return {
        "symbol": symbol,
        "entrez_gene_id": entrez_gene_id,
        "summary_rows": rows,
        "points": point_rows,
        "notes": notes,
        "source": "cBioPortal TCGA PanCancer Atlas 2018",
    }
