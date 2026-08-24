# -*- coding: utf-8 -*-
"""Immuno-oncology target explorer - RNA, TCGA, Protein IHC tabs."""

from __future__ import annotations

import importlib
import sys
from pathlib import Path

import pandas as pd
import streamlit as st

import charts
import csv_io
import data_sources
import pptx_export
import references

# Streamlit can keep old modules after file updates.
for _mod_name in ("csv_io", "charts", "data_sources", "analysis", "expression_compare", "cptac", "data_qc", "references", "pptx_export", "pptx_assets"):
    if _mod_name in sys.modules:
        importlib.reload(sys.modules[_mod_name])

from csv_io import read_csv_safe
from data_sources import SAVE_DIR, collect_and_save, extract_gene_query
from pptx_export import build_analysis_pptx
from references import build_references

BUNDLE_VERSION = 3

st.set_page_config(
    page_title="Immuno Target Explorer",
    page_icon="DNA",
    layout="wide",
    initial_sidebar_state="expanded",
)

CUSTOM_CSS = """
<style>
    .stApp { background: linear-gradient(180deg, #f7f5f2 0%, #eef2f6 100%); }
    [data-testid="stSidebar"] { background: #0b1f33; }
    [data-testid="stSidebar"] * { color: #e8eef4 !important; }
    [data-testid="stSidebar"] .stButton > button {
        background: #1aa6a6; color: #0b1f33 !important; border: 0;
        border-radius: 999px; width: 100%; font-weight: 600;
    }
    .hero { max-width: 920px; margin: 1rem auto 0.5rem auto; text-align: center; }
    .hero h1 { font-size: 2rem; font-weight: 700; color: #0b1f33; }
    .hero p { color: #5b6c7d; font-size: 1.02rem; }
    .save-note {
        background: #ecfdf5; border: 1px solid #a7f3d0; color: #065f46;
        padding: 0.75rem 1rem; border-radius: 12px; margin: 0.6rem 0;
    }
</style>
"""
st.markdown(CUSTOM_CSS, unsafe_allow_html=True)

RNA_COLUMNS = [
    "organ",
    "cancer_type",
    "normal_rna_hpa_ntpm",
    "normal_rna_gtex_ntpm",
    "cancer_rna_tcga_median",
    "cancer_rna_tcga_adjacent_mean",
    "rna_fold_induction",
    "log2_fold_normalized",
    "rna_n_tumor",
    "rna_fold_source",
]

PROTEIN_COLUMNS = [
    "organ",
    "cancer_type",
    "protein_ihc_normal_0_3",
    "protein_ihc_cancer_0_3",
    "ihc_mwu_pvalue",
    "ihc_mwu_effect",
    "ihc_mwu_significant",
]

TCGA_COLUMNS = [
    "cancer_type",
    "n_tumor",
    "n_adjacent_normal",
    "tumor_median_rsem",
    "normal_mean_rsem",
    "tumor_median_log2",
    "normal_mean_log2",
    "fold_tumor_vs_adjacent",
    "log2_fold_tumor_vs_adjacent",
]

IHC_COLUMNS = [
    "cancer_type",
    "ihc_score_0_3",
    "ihc_weighted_0_3",
    "ihc_high",
    "ihc_medium",
    "ihc_low",
    "ihc_not_detected",
    "patients_total",
    "ihc_mwu_pvalue",
    "ihc_mwu_effect",
    "ihc_mwu_significant",
    "protein_ihc_normal_0_3",
    "protein_ihc_cancer_0_3",
]

CPTAC_COLUMNS = [
    "cancer_type",
    "source",
    "n_tumor",
    "median_log2_ratio",
    "mean_log2_ratio",
    "cptac_scale",
    "q1_log2_ratio",
    "q3_log2_ratio",
    "fold_vs_reference",
    "pct_above_reference",
]


def init_state() -> None:
    if "messages" not in st.session_state:
        st.session_state.messages = []
    if "last_result" not in st.session_state:
        st.session_state.last_result = None
    elif st.session_state.last_result and _bundle_needs_comparison(st.session_state.last_result):
        st.session_state.last_result = _ensure_bundle(st.session_state.last_result)


def list_saved_csvs() -> list[Path]:
    if not SAVE_DIR.exists():
        return []
    return sorted(SAVE_DIR.glob("*.csv"), key=lambda p: p.stat().st_mtime, reverse=True)


def _as_df(value) -> pd.DataFrame:
    if isinstance(value, pd.DataFrame):
        return value
    return pd.DataFrame()


def _bundle_needs_comparison(bundle: dict) -> bool:
    if not bundle:
        return False
    if bundle.get("bundle_version", 0) < BUNDLE_VERSION:
        return True
    comparison = _as_df(bundle.get("comparison_df"))
    if not comparison.empty:
        return False
    return "comparison_df" not in bundle or bundle.get("comparison_error")


def _load_comparison_csv(bundle: dict) -> pd.DataFrame:
    csv_path = bundle.get("comparison_csv")
    if csv_path:
        path = Path(str(csv_path))
        if path.exists():
            return read_csv_safe(path)

    identity = bundle.get("identity") or {}
    symbol = identity.get("symbol") or ""
    if not symbol or not SAVE_DIR.exists():
        return pd.DataFrame()

    from data_sources import _safe_filename

    for stem in (symbol, _safe_filename(symbol)):
        matches = sorted(SAVE_DIR.glob(f"{stem}_rna_protein_comparison_*.csv"), reverse=True)
        if matches:
            bundle["comparison_csv"] = str(matches[0])
            return read_csv_safe(matches[0])
    return pd.DataFrame()


def _ensure_comparison_data(bundle: dict) -> dict:
    """Rebuild RNA/protein comparison when an older Streamlit session is missing it."""
    comparison = _as_df(bundle.get("comparison_df"))
    if not comparison.empty:
        bundle["bundle_version"] = BUNDLE_VERSION
        return bundle

    loaded = _load_comparison_csv(bundle)
    if not loaded.empty:
        bundle["comparison_df"] = loaded
        bundle["bundle_version"] = BUNDLE_VERSION
        bundle["comparison_error"] = None
        return bundle

    identity = bundle.get("identity") or {}
    symbol = identity.get("symbol")
    ensembl_id = identity.get("ensembl_id")
    if not symbol:
        bundle["comparison_error"] = "Gene identity missing from session."
        return bundle

    try:
        expr_df = _as_df(bundle.get("expression_df"))
        analytics_df = _as_df(bundle.get("analytics_df"))

        if expr_df.empty:
            refreshed = collect_and_save(symbol)
            refreshed["bundle_version"] = BUNDLE_VERSION
            return refreshed

        from data_sources import _session
        from expression_compare import build_expression_comparison, build_tcga_summary

        sess = _session()
        comparison_df, tcga_df, normal_rna_df = build_expression_comparison(
            symbol=symbol,
            ensembl_id=ensembl_id or "",
            expr_df=expr_df,
            ihc_df=analytics_df,
            session=sess,
        )
        bundle["comparison_df"] = comparison_df
        bundle["normal_rna_df"] = normal_rna_df
        if _as_df(bundle.get("tcga_df")).empty and not tcga_df.empty:
            bundle["tcga_df"] = tcga_df
            bundle["tcga_summary_df"] = build_tcga_summary(tcga_df)
        bundle["bundle_version"] = BUNDLE_VERSION
        bundle["comparison_error"] = None if not comparison_df.empty else "Comparison table is empty for this gene."
    except Exception as exc:
        bundle["comparison_error"] = str(exc)
    return bundle


def _ensure_bundle(bundle: dict | None) -> dict | None:
    if not bundle:
        return bundle
    if _bundle_needs_comparison(bundle):
        bundle = _ensure_comparison_data(bundle)
    bundle = _ensure_tcga_data(bundle)
    bundle = _ensure_cptac_data(bundle)
    return bundle


def _ensure_cptac_data(bundle: dict) -> dict:
    cptac = _as_df(bundle.get("cptac_df"))
    if not cptac.empty:
        return bundle

    identity = bundle.get("identity") or {}
    symbol = identity.get("symbol")
    ensembl_id = identity.get("ensembl_id")
    if not symbol:
        return bundle

    try:
        from cptac import collect_cptac_for_gene
        from data_sources import _session

        sess = _session()
        cptac_df, cptac_summary_df, cptac_error = collect_cptac_for_gene(sess, symbol, ensembl_id)
        bundle["cptac_df"] = cptac_df
        bundle["cptac_summary_df"] = cptac_summary_df
        bundle["cptac_error"] = cptac_error
    except Exception as exc:
        bundle["cptac_error"] = str(exc)
    return bundle


def _ensure_tcga_data(bundle: dict) -> dict:
    """Re-fetch TCGA if an older session bundle is missing per-sample data."""
    tcga = _as_df(bundle.get("tcga_df"))
    if not tcga.empty:
        return bundle

    identity = bundle.get("identity") or {}
    symbol = identity.get("symbol")
    ensembl_id = identity.get("ensembl_id")
    if not symbol:
        return bundle

    try:
        from data_sources import _session
        from expression_compare import build_tcga_summary, fetch_entrez_id, fetch_tcga_rna

        sess = _session()
        entrez = fetch_entrez_id(sess, symbol, ensembl_id)
        if not entrez:
            bundle["tcga_error"] = f"Could not resolve TCGA gene ID for {symbol}."
            return bundle
        tcga = fetch_tcga_rna(sess, entrez)
        bundle["tcga_df"] = tcga
        bundle["tcga_summary_df"] = build_tcga_summary(tcga)
        bundle["tcga_error"] = None if not tcga.empty else "TCGA returned no samples for this gene."
    except Exception as exc:
        bundle["tcga_error"] = str(exc)
    return bundle


def _ihc_view(bundle: dict) -> pd.DataFrame:
    analytics = _as_df(bundle.get("analytics_df"))
    if not analytics.empty:
        out = analytics.copy()
        if "protein_ihc_normal_0_3" not in out.columns:
            comparison = _as_df(bundle.get("comparison_df"))
            if not comparison.empty:
                merge_cols = ["cancer_type", "protein_ihc_normal_0_3", "protein_ihc_cancer_0_3", "ihc_mwu_pvalue", "ihc_mwu_effect", "ihc_mwu_significant"]
                avail = [c for c in merge_cols if c in comparison.columns]
                if len(avail) > 1:
                    out = out.merge(comparison[avail], on="cancer_type", how="left")
        return out

    comparison = _as_df(bundle.get("comparison_df"))
    if comparison.empty:
        return pd.DataFrame()
    rows = []
    for rec in comparison.to_dict("records"):
        if pd.notna(rec.get("protein_ihc_cancer_0_3")):
            rows.append(
                {
                    "cancer_type": rec.get("cancer_type"),
                    "ihc_weighted_0_3": rec.get("protein_ihc_cancer_0_3"),
                    "ihc_fold_vs_negative": None,
                    "protein_ihc_normal_0_3": rec.get("protein_ihc_normal_0_3"),
                    "protein_ihc_cancer_0_3": rec.get("protein_ihc_cancer_0_3"),
                    "ihc_mwu_pvalue": rec.get("ihc_mwu_pvalue"),
                    "ihc_mwu_effect": rec.get("ihc_mwu_effect"),
                    "ihc_mwu_significant": rec.get("ihc_mwu_significant"),
                }
            )
    return pd.DataFrame(rows)


def render_plotly(fig, empty_msg: str, *, key: str) -> None:
    if fig is None:
        st.info(empty_msg)
        return
    st.plotly_chart(fig, use_container_width=True, config={"displayModeBar": False}, key=key)


def render_table(df: pd.DataFrame, columns: list[str], formats: dict[str, str], cite: str = "") -> None:
    if df.empty:
        st.info("No data available for this section.")
        return
    view = df[[c for c in columns if c in df.columns]].copy()
    st.dataframe(
        view.style.format(formats, na_rep="-"),
        use_container_width=True,
        hide_index=True,
    )
    if cite:
        st.caption(f"Table sources {cite}")


def summarize_result(bundle: dict) -> str:
    identity = bundle["identity"]
    symbol = identity["symbol"]
    ensembl = identity["ensembl_id"]
    name = identity.get("name") or ""
    comparison = _as_df(bundle.get("comparison_df"))
    tcga_summary = _as_df(bundle.get("tcga_summary_df"))
    cptac_summary = _as_df(bundle.get("cptac_summary_df"))
    ihc = _ihc_view(bundle)
    refs = build_references(bundle)

    lines = [
        f"**{symbol}** (`{ensembl}`) {refs.mark('ensembl')} - {name}",
        "",
        "Four analysis panels are ready: **RNA**, **Protein**, **TCGA**, and **IHC pathology (HPA + CPTAC)** "
        f"{refs.mark('hpa', 'gtex', 'tcga', 'opentargets')}.",
        "",
    ]
    if not comparison.empty:
        top = comparison.dropna(subset=["log2_fold_normalized"]).head(1)
        if not top.empty:
            r = top.iloc[0]
            lines.append(
                f"- RNA top induction {refs.mark('hpa', 'gtex', 'tcga')}: **{r['cancer_type']}** "
                f"(log2((Tumor+0.25)/(Normal+0.25)) = {r['log2_fold_normalized']:.2f})"
            )
    if not tcga_summary.empty:
        t = tcga_summary.iloc[0]
        lines.append(
            f"- TCGA highest tumor RNA {refs.mark('tcga', 'cbioportal')}: **{t['cancer_type']}** "
            f"(median log2 {t['tumor_median_log2']:.2f})"
        )
    if not ihc.empty:
        h = ihc.sort_values("ihc_weighted_0_3", ascending=False, na_position="last").iloc[0]
        lines.append(
            f"- HPA IHC top score {refs.mark('hpa')}: **{h['cancer_type']}** "
            f"(weighted {h['ihc_weighted_0_3']:.2f})"
        )
    if not cptac_summary.empty:
        c = cptac_summary.iloc[0]
        lines.append(
            f"- CPTAC top tumor protein {refs.mark('cptac', 'pdc', 'cbioportal')}: **{c['cancer_type']}** "
            f"(median log2 ratio {c['median_log2_ratio']:.2f})"
        )
    lines.extend(["", "Research use only. Not for clinical decisions."])
    return "\n".join(lines)


def render_rna_tab(bundle: dict, panel_key: str) -> None:
    comparison = _as_df(bundle.get("comparison_df"))
    normal_rna = _as_df(bundle.get("normal_rna_df"))

    if comparison.empty:
        detail = bundle.get("comparison_error") or "Data not available."
        st.warning(f"RNA comparison could not be loaded. {detail}")
        if st.button("Retry RNA / Protein load", key=f"{panel_key}_retry_comparison_rna"):
            bundle.pop("comparison_df", None)
            bundle.pop("comparison_error", None)
            st.session_state.last_result = _ensure_comparison_data(bundle)
            st.rerun()
        return

    refs = build_references(bundle)
    st.markdown(f"### RNA expression: normal organ vs cancer {refs.mark('hpa', 'gtex', 'tcga')}")
    st.caption(
        f"Sources: HPA consensus nTPM {refs.mark('hpa')} mapped as Normal, "
        f"TCGA median {refs.mark('tcga', 'cbioportal')} mapped as Tumor. "
        "RNA log2 fold = log2((Tumor + 0.25) / (Normal + 0.25))."
    )
    render_table(
        comparison,
        RNA_COLUMNS,
        {
            "normal_rna_hpa_ntpm": "{:.2f}",
            "normal_rna_gtex_ntpm": "{:.2f}",
            "cancer_rna_tcga_median": "{:.1f}",
            "cancer_rna_tcga_adjacent_mean": "{:.1f}",
            "rna_fold_induction": "{:.2f}",
            "log2_fold_normalized": "{:.2f}",
        },
        cite=refs.mark("hpa", "gtex", "tcga", "cbioportal"),
    )
    c1, c2 = st.columns(2)
    with c1:
        render_plotly(
            charts.grouped_rna_bar(comparison),
            "No RNA comparison data for this target.",
            key=f"{panel_key}_rna_grouped",
        )
    with c2:
        render_plotly(
            charts.fold_bar(comparison, "log2_fold_normalized", "RNA log2 fold induction", "log2((Tumor+0.25)/(Normal+0.25))"),
            "RNA log2 fold could not be calculated.",
            key=f"{panel_key}_rna_fold",
        )
    render_plotly(
        charts.normal_rna_bar(normal_rna),
        "Normal organ RNA profile not available.",
        key=f"{panel_key}_rna_normal",
    )


def render_protein_comparison_tab(bundle: dict, panel_key: str) -> None:
    comparison = _as_df(bundle.get("comparison_df"))

    if comparison.empty:
        detail = bundle.get("comparison_error") or "Data not available."
        st.warning(f"Protein comparison could not be loaded. {detail}")
        if st.button("Retry RNA / Protein load", key=f"{panel_key}_retry_comparison_protein"):
            bundle.pop("comparison_df", None)
            bundle.pop("comparison_error", None)
            st.session_state.last_result = _ensure_comparison_data(bundle)
            st.rerun()
        return

    protein = comparison.dropna(subset=["protein_ihc_cancer_0_3"], how="any").copy()
    if protein.empty:
        protein = comparison.dropna(subset=["protein_ihc_normal_0_3"], how="any").copy()

    refs = build_references(bundle)
    st.markdown(f"### Protein expression: normal organ vs cancer (IHC 0-3) {refs.mark('hpa')}")
    st.caption(
        f"HPA normal tissue IHC vs HPA tumor TMA IHC {refs.mark('hpa')}, matched by organ. "
        "Scores are weighted 0-3 (Not detected=0, Low=1, Medium=2, High=3)."
    )

    if protein.empty:
        st.info(
            "No matched normal-vs-cancer protein IHC comparison for this target. "
            "Try the **IHC pathology** tab for raw patient counts, or search a gene with HPA protein data (e.g. HER2, PD-L1)."
        )
        return

    render_table(
        protein,
        PROTEIN_COLUMNS,
        {
            "protein_ihc_normal_0_3": "{:.2f}",
            "protein_ihc_cancer_0_3": "{:.2f}",
            "ihc_mwu_pvalue": "{:.3g}",
            "ihc_mwu_effect": "{:.2f}",
        },
        cite=refs.mark("hpa"),
    )
    render_plotly(
        charts.protein_grouped_bar(protein),
        "Protein IHC normal vs cancer chart unavailable.",
        key=f"{panel_key}_protein_grouped",
    )
    render_plotly(
        charts.protein_fold_horizontal(protein),
        "Protein fold induction chart unavailable.",
        key=f"{panel_key}_protein_fold",
    )
    st.caption(
        f"**IHC statistics** {refs.mark('hpa')}: Scores stay on the ordinal 0-3 scale "
        "(Not detected=0, Low=1, Medium=2, High=3). "
        "Tumor vs matched-normal distributions are compared with a Mann-Whitney U test; "
        "no 0.25 pseudocount fold is applied to IHC."
    )


def render_tcga_tab(bundle: dict, panel_key: str) -> None:
    bundle = _ensure_tcga_data(bundle)
    st.session_state.last_result = bundle

    tcga = _as_df(bundle.get("tcga_df"))
    tcga_summary = _as_df(bundle.get("tcga_summary_df"))

    refs = build_references(bundle)
    st.markdown(f"### TCGA RNA distribution {refs.mark('tcga', 'cbioportal')}")
    st.caption(
        f"Per-sample TCGA pan-cancer atlas RNA-seq {refs.mark('tcga')}, accessed via cBioPortal {refs.mark('cbioportal')}. "
        "Box plot uses tumor samples; table includes adjacent normals when available."
    )

    if tcga.empty:
        detail = bundle.get("tcga_error") or "Unknown error."
        st.warning(
            "TCGA data could not be loaded. Check internet connection and **search the gene again**. "
            f"If the problem persists, cBioPortal may be temporarily unavailable.\n\nDetail: `{detail}`"
        )
        if st.button("Retry TCGA load", key=f"{panel_key}_retry_tcga"):
            bundle["tcga_df"] = pd.DataFrame()
            bundle = _ensure_tcga_data(bundle)
            st.session_state.last_result = bundle
            st.rerun()
        return

    render_table(
        tcga_summary,
        TCGA_COLUMNS,
        {
            "tumor_median_rsem": "{:.1f}",
            "normal_mean_rsem": "{:.1f}",
            "tumor_median_log2": "{:.2f}",
            "normal_mean_log2": "{:.2f}",
            "fold_tumor_vs_adjacent": "{:.2f}",
            "log2_fold_tumor_vs_adjacent": "{:.2f}",
        },
        cite=refs.mark("tcga", "cbioportal"),
    )
    c1, c2 = st.columns(2)
    with c1:
        render_plotly(
            charts.tcga_box(tcga),
            "No TCGA tumor samples found for this gene.",
            key=f"{panel_key}_tcga_box",
        )
    with c2:
        render_plotly(
            charts.tcga_tumor_normal_bar(tcga_summary),
            "TCGA tumor vs adjacent-normal summary not available for this gene.",
            key=f"{panel_key}_tcga_bar",
        )


def render_hpa_ihc_tab(bundle: dict, panel_key: str) -> None:
    ihc = _ihc_view(bundle)

    refs = build_references(bundle)
    st.markdown(f"### HPA IHC pathology {refs.mark('hpa')}")
    st.caption(
        f"Patient-level HPA tumor TMA staining counts (High / Medium / Low / Not detected) {refs.mark('hpa')}. "
        "All available HPA cancer categories are shown (typically up to 20 types). "
        "HPA provides cancer-type summaries only - histological subtypes (e.g. LUAD vs LUSC) are not in this dataset."
    )

    if ihc.empty:
        st.info(
            "No HPA pathology IHC data for this target. "
            "Some genes have limited protein staining records in Human Protein Atlas."
        )
        return

    render_table(
        ihc,
        IHC_COLUMNS,
        {
            "ihc_score_0_3": "{:.0f}",
            "ihc_weighted_0_3": "{:.2f}",
            "ihc_mwu_pvalue": "{:.3g}",
            "ihc_mwu_effect": "{:.2f}",
        },
        cite=refs.mark("hpa"),
    )
    render_plotly(
        charts.ihc_stacked_bar(ihc),
        "IHC patient distribution unavailable.",
        key=f"{panel_key}_ihc_stacked",
    )
    render_plotly(
        charts.ihc_weighted_bar(ihc),
        "IHC weighted score unavailable.",
        key=f"{panel_key}_ihc_weighted",
    )
    if hasattr(charts, "ihc_fold_bar"):
        render_plotly(
            charts.ihc_fold_bar(ihc),
            "IHC fold vs negative unavailable.",
            key=f"{panel_key}_ihc_fold",
        )
    else:
        render_plotly(
            charts.fold_bar(ihc, "ihc_fold_vs_negative", "IHC fold vs negative", "Fold induction"),
            "IHC fold vs negative unavailable.",
            key=f"{panel_key}_ihc_fold",
        )


def render_cptac_tab(bundle: dict, panel_key: str) -> None:
    bundle = _ensure_cptac_data(bundle)
    st.session_state.last_result = bundle

    cptac = _as_df(bundle.get("cptac_df"))
    cptac_summary = _as_df(bundle.get("cptac_summary_df"))

    refs = build_references(bundle)
    st.markdown(f"### CPTAC proteomics (NCI / PDC via cBioPortal) {refs.mark('cptac', 'pdc', 'cbioportal')}")
    st.caption(
        f"Mass-spectrometry protein abundance from CPTAC {refs.mark('cptac')} tumor cohorts integrated in cBioPortal "
        f"{refs.mark('cbioportal')} (linked to NCI Proteomic Data Commons {refs.mark('pdc')}). "
        "Values are CDAP-style log2 ratios vs the study reference/cohort median (typically -5 to +5). "
        "Raw TMT intensities that were mislabeled as log-ratios are median-centered before display."
    )

    if cptac.empty:
        detail = bundle.get("cptac_error") or "No CPTAC protein data returned."
        st.warning(f"CPTAC proteomics could not be loaded. {detail}")
        if st.button("Retry CPTAC load", key=f"{panel_key}_retry_cptac"):
            bundle.pop("cptac_df", None)
            bundle.pop("cptac_summary_df", None)
            bundle.pop("cptac_error", None)
            st.session_state.last_result = _ensure_cptac_data(bundle)
            st.rerun()
        return

    render_table(
        cptac_summary,
        CPTAC_COLUMNS,
        {
            "median_log2_ratio": "{:.2f}",
            "mean_log2_ratio": "{:.2f}",
            "q1_log2_ratio": "{:.2f}",
            "q3_log2_ratio": "{:.2f}",
            "fold_vs_reference": "{:.2f}",
            "pct_above_reference": "{:.1f}",
        },
        cite=refs.mark("cptac", "pdc", "cbioportal"),
    )
    render_plotly(
        charts.cptac_box(cptac),
        "CPTAC sample distribution unavailable.",
        key=f"{panel_key}_cptac_box",
    )
    render_plotly(
        charts.cptac_median_bar(cptac_summary),
        "CPTAC median protein summary unavailable.",
        key=f"{panel_key}_cptac_median",
    )
    st.caption(
        "**Interpretation:** Each row is a CPTAC cancer cohort. "
        "Median log2 ratio > 0 means higher tumor protein than the study reference. "
        "pct_above_reference = percent of tumor samples above reference (log2 > 0)."
    )


def render_ihc_tab(bundle: dict, panel_key: str) -> None:
    tab_hpa, tab_cptac = st.tabs(["4a. HPA IHC pathology", "4b. CPTAC proteomics (NCI/PDC)"])
    with tab_hpa:
        render_hpa_ihc_tab(bundle, panel_key)
    with tab_cptac:
        render_cptac_tab(bundle, panel_key)


def render_ppt_download(bundle: dict, panel_key: str) -> None:
    identity = bundle.get("identity") or {}
    symbol = identity.get("symbol") or "gene"
    pptx_key = f"pptx_bytes_{symbol}"
    pptx_name = f"{symbol}_expression_analysis.pptx"
    payload = st.session_state.get(pptx_key)
    left, right = st.columns([1, 3])
    with left:
        if st.button("Create analysis PPT", key=f"{panel_key}_make_pptx"):
            with st.spinner("Building PowerPoint (charts, fold values, references)..."):
                try:
                    payload = build_analysis_pptx(bundle)
                    st.session_state[pptx_key] = payload
                    st.session_state[f"{pptx_key}_name"] = pptx_name
                    st.success("PPT ready. Click Download on the right.")
                except Exception as exc:
                    st.session_state[pptx_key] = None
                    payload = None
                    st.exception(exc)
    with right:
        if payload:
            st.download_button(
                "Download analysis PPT (.pptx)",
                data=payload,
                file_name=st.session_state.get(f"{pptx_key}_name") or pptx_name,
                mime="application/vnd.openxmlformats-officedocument.presentationml.presentation",
                key=f"{panel_key}_download_pptx",
                type="primary",
            )
        else:
            st.caption("Click Create first. Four slides use the current browser charts; figure notes are on the last slide.")


def render_analytics_panel(bundle: dict, panel_key: str) -> None:
    if _bundle_needs_comparison(bundle):
        with st.spinner("Loading RNA and protein comparison data..."):
            bundle = _ensure_bundle(bundle)
            st.session_state.last_result = bundle

    render_ppt_download(bundle, panel_key)
    tab_rna, tab_protein, tab_tcga, tab_ihc = st.tabs(
        ["1. RNA expression", "2. Protein expression", "3. TCGA distribution", "4. IHC pathology"]
    )
    with tab_rna:
        render_rna_tab(bundle, panel_key)
    with tab_protein:
        render_protein_comparison_tab(bundle, panel_key)
    with tab_tcga:
        render_tcga_tab(bundle, panel_key)
    with tab_ihc:
        render_ihc_tab(bundle, panel_key)
    st.divider()
    st.markdown(build_references(bundle).markdown())


def run_lookup(user_text: str) -> None:
    gene = extract_gene_query(user_text)
    if not gene:
        st.session_state.messages.append(
            {"role": "assistant", "content": "Gene name not found. Try `HER2`, `PD-L1`, `TIGIT`, or an Ensembl ID."}
        )
        return

    with st.spinner(f"Loading HPA, GTEx, TCGA, CPTAC and Open Targets for {gene}..."):
        try:
            bundle = collect_and_save(gene)
        except Exception as exc:
            st.session_state.messages.append(
                {"role": "assistant", "content": f"Lookup failed for **{gene}**.\n\n`{exc}`"}
            )
            return

    st.session_state.last_result = bundle
    save_html = (
        '<div class="save-note">CSV saved<br>'
        f"Comparison: <code>{bundle.get('comparison_csv', '')}</code><br>"
        f"Analytics: <code>{bundle.get('analytics_csv', '')}</code><br>"
        f"CPTAC: <code>{bundle.get('cptac_csv', '')}</code></div>"
    )
    st.session_state.messages.append(
        {
            "role": "assistant",
            "content": summarize_result(bundle),
            "save_html": save_html,
            "bundle": True,
        }
    )


init_state()

with st.sidebar:
    st.markdown("### Immuno-oncology explorer")
    st.caption("HPA / GTEx / TCGA / Open Targets")
    st.markdown("Search a gene. Results split into RNA, Protein, TCGA, and IHC (HPA + CPTAC) tabs.")
    st.markdown(f"`{SAVE_DIR}`")
    st.divider()
    st.markdown("**Quick search**")
    for symbol in ["ERBB2", "CD274", "PDCD1", "LAG3", "TIGIT", "HAVCR2"]:
        if st.button(symbol, key=f"quick_{symbol}"):
            st.session_state.messages = []
            st.session_state.messages.append({"role": "user", "content": symbol})
            run_lookup(symbol)
            st.rerun()
    st.divider()
    if st.session_state.last_result:
        identity = (st.session_state.last_result.get("identity") or {})
        current = identity.get("symbol") or "gene"
        if st.button(f"Refresh {current} panels", key="sidebar_refresh"):
            with st.spinner(f"Refreshing {current}..."):
                gene = identity.get("symbol") or current
                st.session_state.last_result = collect_and_save(gene)
                st.session_state.last_result["bundle_version"] = BUNDLE_VERSION
            st.rerun()
    st.divider()
    st.markdown("**Saved CSV files**")
    saved = list_saved_csvs()
    if not saved:
        st.caption("No files saved yet.")
    for path in saved[:12]:
        st.caption(path.name)

if not st.session_state.messages:
    st.markdown(
        """
        <div class="hero">
            <h1>RNA, Protein, TCGA and IHC analytics</h1>
            <p>Search a target gene to see normalized tables and biotech-style Plotly charts in four panels.</p>
        </div>
        """,
        unsafe_allow_html=True,
    )

last_bundle_idx = None
for idx, message in enumerate(st.session_state.messages):
    if message.get("bundle"):
        last_bundle_idx = idx

for idx, message in enumerate(st.session_state.messages):
    with st.chat_message(message["role"]):
        st.markdown(message["content"])
        if message.get("save_html"):
            st.markdown(message["save_html"], unsafe_allow_html=True)
        if message.get("bundle") and st.session_state.last_result and idx == last_bundle_idx:
            symbol = (st.session_state.last_result.get("identity") or {}).get("symbol") or "gene"
            render_analytics_panel(st.session_state.last_result, panel_key=f"panel_{idx}_{symbol}")

user_input = st.chat_input("Enter a gene name, e.g. HER2 or PD-L1")
if user_input:
    st.session_state.messages.append({"role": "user", "content": user_input})
    run_lookup(user_input)
    st.rerun()
