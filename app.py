# -*- coding: utf-8 -*-
"""TAA Immune Profiler — private research web app."""

from __future__ import annotations

import html
import json
import os
import time
from pathlib import Path

import pandas as pd
import streamlit as st

from profiler import APP_NAME, __version__
from profiler.auth import auth_status, token_fingerprint
from profiler.bulk_expression import fetch_cancer_expression
from profiler.cell_lines import (
    compare_selected,
    fetch_cell_line_expression,
    filter_cell_lines,
    literature_direction_check,
)
from profiler.cbioportal import CELL_LINE_ASSAYS, TCGA_STUDIES
from profiler.charts_ui import bar_cell_lines, box_by_cancer
from profiler.exports import df_to_csv_bytes, to_json_bytes
from profiler.gene_resolution import MHCII_LOCI, PRIORITY_TAAS, resolve_target
from profiler.single_cell import (
    DEFAULT_MARKERS,
    LINEAGE_GROUPS,
    demo_analysis_bundle,
    list_datasets,
    try_ucsc_index,
)
from profiler.taa_ranking import rank_coelevated_taas
from profiler.validation import run_full_validation, summarize

st.set_page_config(
    page_title=APP_NAME,
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

QUICK_START = ["HER2", "EGFR", "PSMA", "STEAP1", "B7-H3", "FAP", "TROP2", "Nectin-4", "CD19", "BCMA"]

CSS = """
<style>
@import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;500;600;700&family=Fraunces:opsz,wght@9..144,600;9..144,700&display=swap');
html, body, [class*="css"]  { font-family: 'IBM Plex Sans', sans-serif; font-size: 17px; color: #0b1f33; }
.stApp { background: linear-gradient(165deg, #f4f8fb 0%, #e7f1f4 45%, #f7fafc 100%); }
h1, h2, h3 { font-family: 'Fraunces', Georgia, serif; color: #0b1f33; }
[data-testid="stSidebar"] { background: #0b1f33; }
[data-testid="stSidebar"] * { color: #e8eef4 !important; }
.hero { text-align: center; max-width: 920px; margin: 0.4rem auto 1rem auto; }
.hero .brand { font-family: 'Fraunces', Georgia, serif; font-size: 2.35rem; font-weight: 700; color: #0b1f33; letter-spacing: -0.02em; }
.hero .sub { color: #355066; font-size: 1.05rem; margin-top: 0.35rem; }
.flow { display:flex; flex-wrap:wrap; gap:8px; justify-content:center; margin: 0.8rem 0 0.2rem; }
.flow span { background:#fff; border:1px solid #b9d0dd; color:#0b1f33; padding:6px 10px; border-radius:8px; font-size:0.92rem; }
.warn { background:#fff7ed; border:1px solid #fdba74; color:#7c2d12; padding:10px 14px; border-radius:10px; margin:0.5rem 0; }
.okbox { background:#ecfdf5; border:1px solid #6ee7b7; color:#065f46; padding:10px 14px; border-radius:10px; }
.failbox { background:#fef2f2; border:1px solid #fca5a5; color:#7f1d1d; padding:10px 14px; border-radius:10px; }
.meta { color:#355066; font-size:0.95rem; }
.stDataFrame { font-size: 15px; }
</style>
"""
st.markdown(CSS, unsafe_allow_html=True)


def init_state() -> None:
    defaults = {
        "authed": False,
        "resolved": None,
        "bulk": None,
        "cell_lines": None,
        "cell_assay": "depmap_25q3_rna",
        "selected_lines": [],
        "top20": None,
        "sc_dataset": None,
        "sc_result": None,
        "validation_run": None,
        "busy_msg": "",
        "request_id": 0,
    }
    for k, v in defaults.items():
        if k not in st.session_state:
            st.session_state[k] = v


def gate() -> bool:
    status = auth_status(
        token=st.session_state.get("access_token_input"),
        email=st.session_state.get("access_email_input"),
    )
    if st.session_state.authed or status["ok"]:
        st.session_state.authed = True
        return True

    st.markdown(
        f"<div class='hero'><div class='brand'>{html.escape(APP_NAME)}</div>"
        f"<div class='sub'>Private research tool · v{__version__}</div></div>",
        unsafe_allow_html=True,
    )
    st.markdown(
        "<div class='warn'>기본값은 비공개입니다. 소유자 승인 없이 공개로 전환하지 않습니다. "
        "액세스 토큰 또는 allowlist 이메일이 필요합니다.</div>",
        unsafe_allow_html=True,
    )
    c1, c2 = st.columns(2)
    with c1:
        st.text_input("Access token", type="password", key="access_token_input")
    with c2:
        st.text_input("Allowlisted email (optional)", key="access_email_input")
    if st.button("Enter", type="primary"):
        again = auth_status(
            token=st.session_state.get("access_token_input"),
            email=st.session_state.get("access_email_input"),
        )
        if again["ok"]:
            st.session_state.authed = True
            st.rerun()
        st.error(again.get("message") or "Access denied")
    st.caption(f"Token fingerprint: {token_fingerprint()} · set TAA_PROFILER_ACCESS_TOKEN / ALLOWLIST / DEV_OPEN")
    return False


def resolve_and_store(query: str) -> None:
    st.session_state.request_id += 1
    rid = st.session_state.request_id
    t0 = time.time()
    with st.status(f"Resolving “{query}”…", expanded=True) as status:
        st.write("1. Alias / special-antigen rules")
        resolved = resolve_target(query)
        st.write(f"2. Provider confirmation ({resolved.mapping_source})")
        if rid != st.session_state.request_id:
            st.write("Cancelled by newer request")
            return
        st.session_state.resolved = resolved.to_dict()
        st.session_state.bulk = None
        st.session_state.cell_lines = None
        st.session_state.top20 = None
        st.session_state.selected_lines = []
        elapsed = time.time() - t0
        status.update(label=f"Resolved in {elapsed:.1f}s", state="complete")


def ensure_cell_lines(assay_key: str) -> None:
    resolved = st.session_state.resolved
    if not resolved or not resolved.get("entrez_gene_id"):
        return
    if (
        st.session_state.cell_lines
        and st.session_state.cell_lines.get("assay_key") == assay_key
        and st.session_state.cell_lines.get("entrez_gene_id") == resolved["entrez_gene_id"]
    ):
        return
    t0 = time.time()
    with st.status(f"Fetching cell-line assay {assay_key}…", expanded=True) as status:
        st.write("cBioPortal molecular-data query (gene-level, not full matrix)")
        bundle = fetch_cell_line_expression(
            resolved["official_symbol"],
            int(resolved["entrez_gene_id"]),
            assay_key=assay_key,
        )
        # DataFrame is not serializable cleanly in session — keep as dict records
        df = bundle.pop("dataframe")
        bundle["records"] = df.to_dict(orient="records") if df is not None and not df.empty else []
        st.session_state.cell_lines = bundle
        st.session_state.cell_assay = assay_key
        st.session_state.selected_lines = []
        status.update(label=f"Done in {time.time()-t0:.1f}s · n={bundle.get('n_measured')}", state="complete" if bundle.get("status") != "fetch_failed" else "error")


def main() -> None:
    init_state()
    # Local/cloud private default: allow DEV_OPEN for this workspace if no token set
    if os.environ.get("TAA_PROFILER_DEV_OPEN", "1") and not os.environ.get("TAA_PROFILER_ACCESS_TOKEN"):
        os.environ.setdefault("TAA_PROFILER_DEV_OPEN", "1")

    if not gate():
        return

    st.markdown(
        f"<div class='hero'><div class='brand'>{html.escape(APP_NAME)}</div>"
        f"<div class='sub'>TAA × immune-target public-data explorer · research use only · v{__version__}</div>"
        "<div class='flow'>"
        "<span>1 Search TAA</span><span>2 Cancer RNA/protein</span><span>3 Cell lines</span>"
        "<span>4 Single-cell</span><span>5 Export + limits</span>"
        "</div></div>",
        unsafe_allow_html=True,
    )

    q = st.text_input("TAA / gene search", placeholder="HER2, PSMA, TROP-2, CLDN18.2, GD2…", key="search_box")
    bcols = st.columns([1, 1, 6])
    with bcols[0]:
        if st.button("Search", type="primary"):
            if q.strip():
                resolve_and_store(q.strip())
    with bcols[1]:
        if st.button("Clear"):
            st.session_state.resolved = None
            st.session_state.bulk = None
            st.session_state.cell_lines = None
            st.session_state.top20 = None

    st.write("Quick start:")
    qcols = st.columns(len(QUICK_START))
    for i, name in enumerate(QUICK_START):
        if qcols[i].button(name, key=f"qs_{name}"):
            resolve_and_store(name)

    resolved = st.session_state.resolved
    if not resolved:
        st.info("Search a target to open expression, cell-line, single-cell, and validation views.")
        with st.expander("MHCII note"):
            st.write("MHCII is not a single gene. Choose one of: " + ", ".join(MHCII_LOCI))
        return

    # Resolution banner
    if resolved.get("unresolved") and resolved.get("special_kind") == "glycolipid":
        st.markdown("<div class='warn'><b>GD2</b> — " + html.escape("; ".join(resolved.get("warnings") or [])) + "</div>", unsafe_allow_html=True)
    elif resolved.get("unresolved"):
        st.markdown("<div class='failbox'>Unresolved: " + html.escape(str(resolved.get("warnings"))) + "</div>", unsafe_allow_html=True)
        return
    else:
        warn_html = ""
        if resolved.get("warnings"):
            warn_html = "<div class='warn'>" + html.escape(" | ".join(resolved["warnings"])) + "</div>"
        st.markdown(
            f"<div class='okbox'><b>{html.escape(resolved.get('query',''))}</b> → "
            f"<b>{html.escape(resolved.get('official_symbol') or '')}</b> · Entrez "
            f"{resolved.get('entrez_gene_id')} · source {html.escape(resolved.get('mapping_source') or '')}"
            f"</div>{warn_html}",
            unsafe_allow_html=True,
        )

    tabs = st.tabs(
        [
            "Target expression",
            "Cancer cell lines",
            "Single-cell lab",
            "Evidence & limits",
            "Validation",
        ]
    )

    # ----- Tab 1: Target expression -----
    with tabs[0]:
        st.subheader("Cancer-type RNA / RPPA (TCGA PanCancer Atlas 2018)")
        st.caption("Bulk tumor mixed-tissue signal. RNA and RPPA are never merged onto one axis.")
        cancer_options = ["All"] + list(TCGA_STUDIES.keys())
        pick = st.multiselect("Cancers", options=list(TCGA_STUDIES.keys()), default=["BRCA", "LUAD", "PRAD", "PAAD", "OV"])
        if st.button("Load cancer expression", key="load_bulk"):
            with st.status("Querying cBioPortal per study…", expanded=True) as status:
                t0 = time.time()
                bulk = fetch_cancer_expression(
                    resolved["official_symbol"],
                    int(resolved["entrez_gene_id"]),
                    cancer_keys=pick or None,
                    include_rppa=True,
                )
                st.session_state.bulk = bulk
                status.update(label=f"Finished in {time.time()-t0:.1f}s", state="complete")

        bulk = st.session_state.bulk
        if bulk:
            summary = pd.DataFrame(bulk["summary_rows"])
            rna = summary[(summary.get("assay") == "RNA") & (summary.get("level").isna() if "level" in summary.columns else True)]
            if "level" in summary.columns:
                rna = summary[(summary["assay"] == "RNA") & (summary["level"].isna() | (summary["level"] != "subtype"))]
            else:
                rna = summary[summary["assay"] == "RNA"]
            st.dataframe(rna, use_container_width=True, hide_index=True)
            points = pd.DataFrame(bulk.get("points") or [])
            if not points.empty:
                fig = box_by_cancer(points, f"{resolved['official_symbol']} RNA log2(value+1) by cancer", unit="log2(expression+1)")
                st.plotly_chart(fig, use_container_width=True)
            for note in bulk.get("notes") or []:
                st.caption("• " + note)
            st.download_button(
                "Download cancer summary CSV",
                data=df_to_csv_bytes(
                    summary,
                    {
                        "requested_target": resolved.get("query"),
                        "resolved_symbol": resolved.get("official_symbol"),
                        "gene_id": resolved.get("entrez_gene_id"),
                        "source": bulk.get("source"),
                    },
                ),
                file_name=f"{resolved['official_symbol']}_cancer_expression.csv",
                mime="text/csv",
            )

            st.markdown("#### Other TAA relative elevation Top 20")
            st.caption(
                "High vs low quartile of the reference TAA in the same cohort. "
                "Not treatment-induced change. Metric = difference of mean log2(expr+1)."
            )
            study_default = "brca_tcga_pan_can_atlas_2018"
            study = st.selectbox(
                "Cohort for Top 20",
                options=[TCGA_STUDIES[k]["study"] for k in pick] or [study_default],
            )
            if st.button("Compute Top 20", key="top20"):
                with st.status("Ranking panel TAAs (rate-limited)…", expanded=True) as status:
                    t0 = time.time()
                    top = rank_coelevated_taas(study, resolved["official_symbol"], int(resolved["entrez_gene_id"]))
                    st.session_state.top20 = top
                    status.update(label=f"{top.get('status')} in {time.time()-t0:.1f}s", state="complete")
            top = st.session_state.top20
            if top:
                if top.get("status") != "ok":
                    st.warning(f"{top.get('status')}: {top.get('error')}")
                else:
                    st.dataframe(pd.DataFrame(top["rows"]), use_container_width=True, hide_index=True)
                    for lim in top.get("interpretation_limits") or []:
                        st.caption("• " + lim)
                    st.download_button(
                        "Download Top 20 CSV",
                        data=df_to_csv_bytes(pd.DataFrame(top["rows"]), {"study": study, "reference": resolved["official_symbol"]}),
                        file_name=f"{resolved['official_symbol']}_top20.csv",
                        mime="text/csv",
                    )

    # ----- Tab 2: Cancer cell lines -----
    with tabs[1]:
        st.markdown(
            "<div class='okbox'><b>Human-derived cancer cell lines</b> — cultured models via cBioPortal "
            "DepMap/CCLE. Not patient tumors.</div>",
            unsafe_allow_html=True,
        )
        assay_labels = {k: v["label"] + f" [{v['unit']}]" for k, v in CELL_LINE_ASSAYS.items()}
        assay_key = st.selectbox(
            "Assay",
            options=list(CELL_LINE_ASSAYS.keys()),
            format_func=lambda k: assay_labels[k],
            key="assay_select",
        )
        if assay_key != st.session_state.cell_assay:
            st.session_state.selected_lines = []
        if st.button("Load cell-line expression", key="load_cl") or (
            st.session_state.cell_lines
            and st.session_state.cell_lines.get("entrez_gene_id") == resolved.get("entrez_gene_id")
            and st.session_state.cell_lines.get("assay_key") == assay_key
        ):
            ensure_cell_lines(assay_key)

        bundle = st.session_state.cell_lines
        if not bundle:
            st.info("Load an assay to list models.")
        elif bundle.get("status") == "fetch_failed":
            st.markdown(f"<div class='failbox'>조회 실패 (not empty expression): {html.escape(str(bundle.get('error')))}</div>", unsafe_allow_html=True)
        else:
            df = pd.DataFrame(bundle.get("records") or [])
            st.caption(
                f"n_measured={bundle.get('n_measured')} · retrieved {bundle.get('retrieved_at')} · "
                f"checksum {bundle.get('metadata',{}).get('checksum')} · "
                f"unit table={CELL_LINE_ASSAYS[assay_key]['unit']}"
            )
            for note in bundle.get("notes") or []:
                st.caption("• " + note)

            f1, f2, f3, f4 = st.columns(4)
            with f1:
                qline = st.text_input("Filter name / DepMap ID")
            with f2:
                cancers = ["All"] + sorted({c for c in df.get("cancer_type", pd.Series(dtype=str)).dropna().unique()})
                cancer_f = st.selectbox("Cancer/lineage", cancers)
            with f3:
                sort_by = st.selectbox("Sort", ["value high→low", "value low→high", "name"])
            with f4:
                top_n = st.selectbox("Top N chart", [10, 20, 50], index=1)

            filtered = filter_cell_lines(df, query=qline, cancer_type=cancer_f)
            if sort_by == "value high→low":
                filtered = filtered.sort_values("value", ascending=False)
            elif sort_by == "value low→high":
                filtered = filtered.sort_values("value", ascending=True)
            else:
                filtered = filtered.sort_values("sampleId")

            # pagination 25
            page_size = 25
            n_pages = max(1, (len(filtered) + page_size - 1) // page_size)
            page = st.number_input("Page", min_value=1, max_value=n_pages, value=1)
            start = (page - 1) * page_size
            page_df = filtered.iloc[start : start + page_size]
            show_cols = [c for c in [
                "sampleId", "cell_line_name", "depmap_id", "cancer_type", "cancer_type_detailed",
                "model_type", "value_table", "unit_table", "proteomics_10plex_id", "proteomics_tmt_label",
            ] if c in page_df.columns]
            st.dataframe(page_df[show_cols], use_container_width=True, hide_index=True)

            chart_df = filtered.head(int(top_n))
            fig = bar_cell_lines(
                chart_df,
                f"{resolved['official_symbol']} — top {top_n} ({CELL_LINE_ASSAYS[assay_key]['label']})",
                unit=str(chart_df["unit_plot"].iloc[0]) if not chart_df.empty and "unit_plot" in chart_df.columns else "",
            )
            st.plotly_chart(fig, use_container_width=True)

            options = filtered["sampleId"].tolist()
            selected = st.multiselect(
                "Compare up to 12 lines (selection kept across list filters; cleared on assay change)",
                options=options,
                default=[s for s in st.session_state.selected_lines if s in options][:12],
                max_selections=12,
            )
            st.session_state.selected_lines = selected
            if selected:
                cmp = compare_selected(df, selected, assay_key)
                st.dataframe(pd.DataFrame(cmp["rows"]), use_container_width=True, hide_index=True)
                if cmp.get("warning"):
                    st.warning(cmp["warning"])
                elif cmp["rows"]:
                    fig2 = bar_cell_lines(
                        pd.DataFrame(cmp["rows"]),
                        "Selected comparison",
                        unit=CELL_LINE_ASSAYS[assay_key]["unit"],
                    )
                    st.plotly_chart(fig2, use_container_width=True)

            bio = literature_direction_check(df, resolved["official_symbol"])
            if bio:
                st.markdown("##### Literature surface-direction vs this assay (biology column — not fidelity)")
                st.dataframe(pd.DataFrame(bio), use_container_width=True, hide_index=True)

            cdl1, cdl2 = st.columns(2)
            with cdl1:
                st.download_button(
                    "CSV (filtered)",
                    data=df_to_csv_bytes(
                        filtered,
                        {
                            "symbol": resolved["official_symbol"],
                            "assay": assay_key,
                            "unit": CELL_LINE_ASSAYS[assay_key]["unit"],
                            "retrieved_at": bundle.get("retrieved_at"),
                        },
                    ),
                    file_name=f"{resolved['official_symbol']}_{assay_key}_filtered.csv",
                    mime="text/csv",
                )
            with cdl2:
                st.download_button(
                    "JSON (full + provenance)",
                    data=to_json_bytes(
                        {"records": bundle.get("records"), "meta": {k: bundle.get(k) for k in ("status", "n_measured", "retrieved_at", "metadata", "assay_meta", "notes")}},
                        {"symbol": resolved["official_symbol"], "entrez": resolved["entrez_gene_id"]},
                    ),
                    file_name=f"{resolved['official_symbol']}_{assay_key}.json",
                    mime="application/json",
                )

    # ----- Tab 3: Single-cell -----
    with tabs[2]:
        st.subheader("Public single-cell studies")
        st.caption("UCSC Cell Browser · TISCH2 catalog. Full matrices are not fabricated.")
        ucsc = try_ucsc_index()
        st.write(f"UCSC index probe: **{ucsc.get('status')}**", ucsc)
        cancer_f = st.selectbox("Cancer filter", ["All"] + list(TCGA_STUDIES.keys()), key="sc_cancer")
        datasets = list_datasets(None if cancer_f == "All" else cancer_f)
        st.dataframe(pd.DataFrame(datasets), use_container_width=True, hide_index=True)
        ds_id = st.selectbox("Select dataset", [d["id"] for d in datasets], format_func=lambda i: next(d["title"] for d in datasets if d["id"] == i))
        lineages = st.multiselect("Lineages (raw annotation filter — no auto cluster mapping)", LINEAGE_GROUPS, default=["T", "B", "DC"])
        markers = st.multiselect(
            "Markers",
            options=list(DEFAULT_MARKERS.values()) + MHCII_LOCI,
            default=["PDCD1", "TNFRSF9", "CD274", "CD40", "CLEC9A"],
        )
        st.caption("MHCII: pick specific HLA loci. Threshold default 0 = RNA>0, not protein-positive.")
        if st.button("Prepare analysis", key="sc_run"):
            st.session_state.sc_result = demo_analysis_bundle(ds_id, markers)
            st.session_state.sc_result["requested_lineages"] = lineages
        if st.session_state.sc_result:
            res = st.session_state.sc_result
            if res.get("status") == "awaiting_matrix":
                st.warning(res.get("message"))
                st.json({k: res[k] for k in res if k != "dataset"})
            else:
                st.json(res)
            ds = next((d for d in datasets if d["id"] == ds_id), {})
            st.markdown("**Provenance**")
            st.write(ds)

    # ----- Tab 4: Evidence -----
    with tabs[3]:
        st.subheader("Evidence & interpretation limits")
        st.markdown(
            """
- Public research exploration only — not clinical response, affinity, or efficacy scores.
- RNA, RPPA, and MS protein are separate assays; never one common expression axis.
- CLDN18.2: gene-level data cannot resolve isoforms.
- GD2: glycolipid — not gene-mappable expression.
- FAP: check stromal/CAF context; cell lines alone are insufficient.
- Protein not measured ≠ antigen absent.
- Top 20 co-elevation ≠ post-treatment increase.
- Single-cell: processed matrices only; dropout/dissociation biases apply.
            """
        )
        st.markdown("##### Priority TAA panel (approval context is indicative, verify with regulators)")
        st.dataframe(pd.DataFrame(PRIORITY_TAAS), use_container_width=True, hide_index=True)
        if resolved:
            st.download_button(
                "Export resolution JSON",
                data=to_json_bytes(resolved, {"tab": "evidence"}),
                file_name="resolution.json",
                mime="application/json",
            )

    # ----- Tab 5: Validation -----
    with tabs[4]:
        st.subheader("Validation")
        st.caption("Real check log only. Data fidelity and biology concordance are separate columns. No fake PASS.")
        if st.button("Run validation suite (aliases + HER2 regression + 17-target profiles)", type="primary"):
            progress = st.empty()
            def _prog(msg: str) -> None:
                progress.info(msg)
            with st.spinner("Running — respects 429 backoff; may take several minutes"):
                run = run_full_validation(progress=_prog)
            st.session_state.validation_run = run
            progress.empty()

        run = st.session_state.validation_run
        if not run:
            st.info("Not run yet this session. Prior claims of PASS without execution are invalid.")
        else:
            summ = run.get("summary") or summarize(run)
            st.markdown(
                f"<div class='okbox'>run <b>{html.escape(run['run_id'])}</b> · "
                f"started {html.escape(run.get('started_at') or '')} · finished {html.escape(run.get('finished_at') or '')}<br/>"
                f"pass={summ.get('pass')} mismatch={summ.get('mismatch')} fetch_failed={summ.get('fetch_failed')} "
                f"data_insufficient={summ.get('data_insufficient')} held={summ.get('held')} · "
                f"fidelity_pass={summ.get('fidelity_pass')} fidelity_fail={summ.get('fidelity_fail')} · "
                f"biology_pass_direction={summ.get('biology_pass_direction')} "
                f"biology_mismatch_direction={summ.get('biology_mismatch_direction')}"
                f"</div>",
                unsafe_allow_html=True,
            )
            checks = pd.DataFrame(run.get("checks") or [])
            if not checks.empty:
                view = checks[["name", "status", "data_fidelity", "biology_concordance", "recorded_at"]]
                st.dataframe(view, use_container_width=True, hide_index=True)
            st.download_button(
                "Download validation JSON",
                data=json.dumps(run, indent=2, default=str).encode("utf-8"),
                file_name=f"validation_{run['run_id']}.json",
                mime="application/json",
            )
            if run.get("result_path"):
                st.caption("Saved: " + run["result_path"])


if __name__ == "__main__":
    main()
else:
    main()
