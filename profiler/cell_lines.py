# -*- coding: utf-8 -*-
"""Human-derived cancer cell line expression (DepMap/CCLE via cBioPortal)."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from profiler.cbioportal import CELL_LINE_ASSAYS, fetch_molecular_data, load_clinical_samples

# Literature directionality checks (not absolute truth values).
LITERATURE_DIRECTIONS = [
    {
        "target": "HER2",
        "symbol": "ERBB2",
        "high": "ZR-75-1",
        "low": "MDA-MB-231",
        "high_sample": "ZR751_BREAST",
        "low_sample": "MDAMB231_BREAST",
        "assay_note": "surface flow receptors/cell; DOI 10.4061/2010/814176",
        "high_depmap": "ACH-000097",
        "low_depmap": "ACH-000768",
    },
    {
        "target": "EGFR",
        "symbol": "EGFR",
        "high": "A431",
        "low": "MCF-7",
        "high_sample": "A431_SKIN",
        "low_sample": "MCF7_BREAST",
        "assay_note": "flow; DOI 10.1177/39.8.1856460",
    },
    {
        "target": "PSMA",
        "symbol": "FOLH1",
        "high": "LNCaP",
        "low": "PC3",
        "high_sample": "LNCAPCLONEFGC_PROSTATE",
        "low_sample": "PC3_PROSTATE",
        "assay_note": "intact-cell flow parental; DOI 10.3390/ma12050756; exclude engineered PC3-PIP",
    },
    {
        "target": "STEAP1",
        "symbol": "STEAP1",
        "high": "LNCaP",
        "low": "PC-3",
        "high_sample": "LNCAPCLONEFGC_PROSTATE",
        "low_sample": "PC3_PROSTATE",
        "assay_note": "Oslo-1 flow; DOI 10.1016/j.omto.2022.06.007",
    },
    {
        "target": "MSLN",
        "symbol": "MSLN",
        "high": "EKVX",
        "low": "A549",
        "high_sample": "EKVX_LUNG",
        "low_sample": "A549_LUNG",
        "assay_note": "surface flow; DOI 10.1158/1078-0432.CCR-06-2161; A549 not absolute negative",
    },
    {
        "target": "B7-H3",
        "symbol": "CD276",
        "high": "MDA-MB-231",
        "low": "Raji",
        "high_sample": "MDAMB231_BREAST",
        "low_sample": "RAJI_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE",
        "assay_note": "flow; DOI 10.3390/life12020157",
    },
    {
        "target": "TROP2",
        "symbol": "TACSTD2",
        "high": "MDA-MB-468",
        "low": "BT-549",
        "high_sample": "MDAMB468_BREAST",
        "low_sample": "BT549_BREAST",
        "assay_note": "native baseline flow; DOI 10.1136/jitc-2025-012442",
    },
    {
        "target": "Nectin-4",
        "symbol": "NECTIN4",
        "high": "MCF-7",
        "low": "MDA-MB-231",
        "high_sample": "MCF7_BREAST",
        "low_sample": "MDAMB231_BREAST",
        "assay_note": "flow; DOI 10.1186/1471-2407-7-73",
    },
    {
        "target": "CD19",
        "symbol": "CD19",
        "high": "Raji",
        "low": "K562",
        "high_sample": "RAJI_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE",
        "low_sample": "K562_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE",
        "assay_note": "flow; DOI 10.3390/ijms18122773",
    },
    {
        "target": "CD20",
        "symbol": "MS4A1",
        "high": "Raji",
        "low": "Jurkat",
        "high_sample": "RAJI_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE",
        "low_sample": "JURKAT_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE",
        "assay_note": "flow; DOI 10.7717/peerj.31",
    },
    {
        "target": "BCMA",
        "symbol": "TNFRSF17",
        "high": "RPMI8226",
        "low": "CCRF-CEM",
        "high_sample": "RPMI8226_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE",
        "low_sample": "CCRCEM_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE",
        "assay_note": "live-cell flow; DOI 10.1158/1078-0432.CCR-12-2422",
    },
    {
        "target": "GPRC5D",
        "symbol": "GPRC5D",
        "high": "RPMI8226",
        "low": "Raji",
        "high_sample": "RPMI8226_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE",
        "low_sample": "RAJI_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE",
        "assay_note": "flow; DOI 10.1182/bloodadvances.2020003805",
    },
    {
        "target": "DLL3",
        "symbol": "DLL3",
        "high": "NCI-H82",
        "mid": "NCI-H69",
        "low": "A549",
        "high_sample": "NCIH82_LUNG",
        "mid_sample": "NCIH69_LUNG",
        "low_sample": "A549_LUNG",
        "assay_note": "antibody radiobinding; DOI 10.1158/0008-5472.CAN-17-0299",
    },
]


def _normalize_name(text: str) -> str:
    return "".join(ch for ch in str(text).upper() if ch.isalnum())


def _clinical_lookup(study_id: str) -> pd.DataFrame:
    try:
        return load_clinical_samples(study_id)
    except Exception:
        return pd.DataFrame()


def _annotate(df: pd.DataFrame, clinical: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    out = df.copy()
    if clinical.empty:
        out["cancer_type"] = None
        out["depmap_id"] = None
        out["model_type"] = "Cell Line"
        out["proteomics_10plex_id"] = None
        out["proteomics_tmt_label"] = None
        return out

    cols = {c.upper(): c for c in clinical.columns}
    sid = cols.get("SAMPLE_ID")
    if not sid:
        return out
    clin = clinical.set_index(clinical[sid].astype(str))

    def pick(row_id: str, *keys: str) -> Any:
        if row_id not in clin.index:
            return None
        for k in keys:
            col = cols.get(k)
            if col and col in clin.columns:
                val = clin.at[row_id, col]
                if isinstance(val, pd.Series):
                    val = val.iloc[0]
                if pd.notna(val) and str(val).strip():
                    return str(val)
        return None

    out["cancer_type"] = out["sampleId"].map(
        lambda s: pick(str(s), "CANCER_TYPE", "CANCER TYPE", "ONCOTREE_PRIMARY_DISEASE", "PRIMARY_DISEASE")
    )
    out["cancer_type_detailed"] = out["sampleId"].map(
        lambda s: pick(str(s), "CANCER_TYPE_DETAILED", "CANCER TYPE DETAILED", "HISTOLOGY")
    )
    out["depmap_id"] = out["sampleId"].map(
        lambda s: pick(str(s), "DEPMAP_ID", "DEPMAP ID", "ACH_ID", "BROAD_ID")
    )
    out["model_type"] = out["sampleId"].map(
        lambda s: pick(str(s), "MODEL_TYPE", "TYPE", "SAMPLE_TYPE") or "Cell Line"
    )
    out["cell_line_name"] = out["sampleId"].map(
        lambda s: pick(str(s), "CELL_LINE_NAME", "CELL LINE NAME", "CCLE_NAME") or str(s).split("_")[0]
    )
    out["proteomics_10plex_id"] = out["sampleId"].map(
        lambda s: pick(str(s), "PROTEOMICS_10PLEX_ID", "PROTEOMICS 10PLEX ID")
    )
    out["proteomics_tmt_label"] = out["sampleId"].map(
        lambda s: pick(str(s), "PROTEOMICS_TMT_LABEL", "PROTEOMICS TMT LABEL")
    )
    return out


def fetch_cell_line_expression(
    symbol: str,
    entrez_gene_id: int,
    assay_key: str = "depmap_25q3_rna",
) -> dict[str, Any]:
    meta = CELL_LINE_ASSAYS[assay_key]
    result = fetch_molecular_data(
        meta["profile"],
        entrez_gene_id,
        meta["sample_list"],
        study_id=meta["study"],
        symbol=symbol,
        assay_key=assay_key,
        unit=meta["unit"],
        transform=meta["transform"],
    )
    df = result.to_dataframe()
    clinical = _clinical_lookup(meta["study"])
    df = _annotate(df, clinical)

    if not df.empty:
        if assay_key == "ccle_2019_rna":
            df["value_table"] = df["value"]
            df["value_plot"] = np.log2(df["value"].astype(float) + 1.0)
            df["unit_table"] = "RPKM"
            df["unit_plot"] = "log2(RPKM+1)"
        else:
            df["value_table"] = df["value"]
            df["value_plot"] = df["value"]
            df["unit_table"] = meta["unit"]
            df["unit_plot"] = meta.get("plot_unit") or meta["unit"]

    notes = [
        "Human-derived cultured cancer models (not patient tumors).",
        f"Assay: {meta['label']}; unit: {meta['unit']}; transform: {meta['transform']}.",
        "Protein values are whole-cell MS relative log2, not surface receptors/cell or % positive.",
        "Negative protein values are valid relative abundances; do not drop them as missing.",
        "NA (absent from profile) is distinct from measured 0.",
        "FAP stromal/CAF context is not established by cell-line lists alone.",
    ]
    if assay_key == "ccle_2019_protein":
        notes.append(
            "Simple ranking/comparison charts are limited to the same TMT batch (PROTEOMICS_10PLEX_ID) "
            "as a conservative UI policy; browsing all batches remains allowed."
        )
        notes.append(
            "Source: Nusinow et al. Cell 2020 DOI 10.1016/j.cell.2019.12.023; "
            "normalized file https://gygi.hms.harvard.edu/data/ccle/protein_quant_current_normalized.csv.gz"
        )

    return {
        "symbol": symbol,
        "entrez_gene_id": entrez_gene_id,
        "assay_key": assay_key,
        "assay_meta": meta,
        "status": result.status,
        "error": result.error,
        "n_measured": result.n_measured,
        "retrieved_at": result.retrieved_at,
        "metadata": result.metadata,
        "dataframe": df,
        "notes": notes,
        "source": "cBioPortal " + meta["study"],
    }


def filter_cell_lines(
    df: pd.DataFrame,
    *,
    query: str = "",
    cancer_type: str = "All",
    model_type: str = "Cell Line",
    include_organoid: bool = False,
) -> pd.DataFrame:
    if df is None or df.empty:
        return pd.DataFrame()
    out = df.copy()
    if not include_organoid and "model_type" in out.columns:
        mt = out["model_type"].fillna("Cell Line").astype(str)
        out = out[~mt.str.contains("organoid", case=False, na=False)]
        if model_type and model_type != "All":
            out = out[mt.str.contains(model_type, case=False, na=True) | (model_type == "Cell Line")]
    if cancer_type and cancer_type != "All" and "cancer_type" in out.columns:
        out = out[out["cancer_type"].fillna("").str.contains(cancer_type, case=False, na=False)]
    if query:
        qn = _normalize_name(query)
        mask = (
            out.get("sampleId", pd.Series(dtype=str)).fillna("").map(_normalize_name).str.contains(qn, na=False)
            | out.get("cell_line_name", pd.Series(dtype=str)).fillna("").map(_normalize_name).str.contains(qn, na=False)
            | out.get("depmap_id", pd.Series(dtype=str)).fillna("").map(_normalize_name).str.contains(qn, na=False)
        )
        out = out[mask]
    return out.reset_index(drop=True)


def compare_selected(df: pd.DataFrame, sample_ids: list[str], assay_key: str) -> dict[str, Any]:
    if df is None or df.empty or not sample_ids:
        return {"rows": [], "warning": "No selection"}
    sub = df[df["sampleId"].isin(sample_ids)].copy()
    warning = None
    if assay_key == "ccle_2019_protein" and "proteomics_10plex_id" in sub.columns:
        batches = sorted({b for b in sub["proteomics_10plex_id"].dropna().astype(str) if b})
        if len(batches) > 1:
            warning = (
                f"Selected lines span TMT batches {batches}. "
                "Comparison chart is withheld; table values remain visible."
            )
    return {"rows": sub.to_dict(orient="records"), "warning": warning}


def literature_direction_check(df: pd.DataFrame, symbol: str) -> list[dict[str, Any]]:
    """Compare RNA/protein order vs literature surface-expression direction. Separate from fidelity."""
    results = []
    if df is None or df.empty:
        return results

    def _value_for(sample_id: str) -> float | None:
        if not sample_id or "sampleId" not in df.columns:
            return None
        hits = df.loc[df["sampleId"].astype(str) == str(sample_id), "value"]
        if hits.empty:
            return None
        # Duplicate sample rows: take first measured value; do not silently average.
        val = hits.iloc[0]
        if pd.isna(val):
            return None
        return float(val)

    for item in LITERATURE_DIRECTIONS:
        if item["symbol"] != symbol:
            continue
        hv = _value_for(item["high_sample"])
        lv = _value_for(item["low_sample"])
        mid_id = item.get("mid_sample") or ""
        mv = _value_for(mid_id) if mid_id else None
        record = {
            "target": item["target"],
            "symbol": symbol,
            "literature_direction": f"{item['high']} > {item['low']}",
            "assay_note": item["assay_note"],
            "high_sample": item["high_sample"],
            "low_sample": item["low_sample"],
            "high_value": hv,
            "low_value": lv,
            "data_fidelity": "n/a — biological direction check",
            "status": "pending",
        }
        if hv is None or lv is None:
            record["status"] = "data_insufficient"
            record["biology_match"] = None
            record["note"] = "One or both comparison lines missing/unmeasured in this assay"
        else:
            match = hv > lv
            record["biology_match"] = match
            record["status"] = "pass_direction" if match else "mismatch_direction"
            record["note"] = (
                "Literature surface assay direction vs this molecular assay. "
                "Mismatch does not imply data corruption; assays measure different quantities."
            )
            if mv is not None and mid_id:
                record["mid_sample"] = mid_id
                record["mid_value"] = mv
                record["biology_match_mid"] = hv > mv > lv
        results.append(record)
    return results
