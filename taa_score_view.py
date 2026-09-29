# -*- coding: utf-8 -*-
"""Simple TAA score on human cancer cell lines."""

from __future__ import annotations

import pandas as pd
import streamlit as st

from taa_analysis.pipeline import export_frame
from taa_analysis.score import formula_text, score_gene
from taa_analysis.xena import RetrievalFailed

TABLE_COLUMNS = [
    "cancer",
    "n_cell_lines",
    "median_nTPM",
    "q1_nTPM",
    "q3_nTPM",
    "detected_fraction",
    "selectivity_log2",
    "score",
]

CSS = """
<style>
@import url('https://fonts.googleapis.com/css2?family=Noto+Sans+KR:wght@400;600&family=Noto+Sans:wght@400;600&display=swap');
.stApp { background: #ffffff; }
html, body, [class*="css"] { font-family: "Noto Sans KR", "Noto Sans", sans-serif; }
.taa-score {
    background: #0b1f33;
    color: #e8eef4;
    border-radius: 14px;
    padding: 16px 18px;
    margin: 0.4rem 0 0.8rem 0;
}
.taa-score .num { font-size: 2.2rem; font-weight: 700; color: #5eead4; line-height: 1.1; }
.taa-score .meta { color: #b7c6d4; margin-top: 0.35rem; }
</style>
"""


def _chart(scored: pd.DataFrame):
    import plotly.graph_objects as go

    shown = scored[scored["scored"] & scored["score"].notna()].head(15)
    if shown.empty:
        return None
    labels = [str(item) for item in shown["cancer"].tolist()]
    figure = go.Figure(
        go.Bar(
            x=shown["score"],
            y=labels,
            orientation="h",
            marker_color="#1aa6a6",
            hovertemplate="%{y}<br>score %{x:.1f}<extra></extra>",
        )
    )
    figure.update_layout(
        xaxis_title="score 0–100",
        yaxis={"categoryorder": "array", "categoryarray": list(reversed(labels))},
        margin={"l": 10, "r": 10, "t": 10, "b": 10},
        height=max(320, 28 * len(labels) + 80),
        plot_bgcolor="#ffffff",
        paper_bgcolor="#ffffff",
    )
    return figure


def render_cell_line_score() -> None:
    st.markdown(CSS, unsafe_allow_html=True)
    st.markdown("## TAA score")
    st.caption(
        "사람 암 세포주 기준입니다. TAA를 검색하면 그 유전자의 HPA nTPM으로 "
        "암종별 발현, 검출 비율, 선택성을 계산하고 0–100 점수를 만듭니다."
    )
    gene = st.text_input(
        "TAA / Gene",
        value=st.session_state.get("score_gene", ""),
        placeholder="HER2, STEAP1, PD-L1, EGFR",
        key="score_gene_input",
    )
    load_col, retry_col = st.columns(2)
    load = load_col.button("조회 / Load", type="primary", key="score_load")
    retry = retry_col.button("다시 받기 / Retry fetch", key="score_retry")
    state_key = gene.strip()
    if st.session_state.get("score_key") not in (None, state_key) and not load and not retry:
        st.session_state.score_bundle = None
        st.session_state.score_error = ""
        st.info("유전자가 바뀌었습니다. 조회를 다시 누르세요.")
        return
    if load or retry:
        try:
            if not gene.strip():
                raise ValueError("empty gene query")
            with st.spinner("암 세포주 nTPM을 모아 점수를 계산하는 중"):
                st.session_state.score_bundle = score_gene(gene.strip(), force=bool(retry))
                st.session_state.score_key = state_key
                st.session_state.score_error = ""
        except RetrievalFailed as exc:
            st.session_state.score_error = str(exc)
            st.session_state.score_bundle = None
        except Exception as exc:
            st.session_state.score_error = f"retrieval_failed: {exc}"
            st.session_state.score_bundle = None

    if st.session_state.get("score_error"):
        st.error(f"조회에 실패했습니다. {st.session_state.score_error}")
        return
    bundle = st.session_state.get("score_bundle")
    if not bundle:
        st.info("HER2, STEAP1, PD-L1처럼 TAA를 입력하고 조회를 누르세요. 별칭은 승인 심볼로 바뀝니다.")
        return

    identity = bundle["identity"]
    meta = bundle["meta"]
    scored = bundle["scored"]
    ranked = scored[scored["scored"] & scored["score"].notna()] if not scored.empty else scored
    top = ranked.iloc[0] if not ranked.empty else None
    if top is not None:
        st.markdown(
            f'<div class="taa-score"><div class="num">{float(top["score"]):.1f}</div>'
            f'<div>{identity.get("gene_symbol")} · {top["cancer"]} · 세포주 {int(top["n_cell_lines"])}개</div>'
            f'<div class="meta">{identity.get("gene_id") or ""} · '
            f'암 세포주 {int(bundle["n_scored_lines"])}개 · 암종 {int(bundle["n_cancer_groups"])}개 · '
            f'nTPM 중앙값 {float(top["median_nTPM"]):.2f}</div></div>',
            unsafe_allow_html=True,
        )
    st.caption(
        f"{meta.get('source_name')} | {meta.get('dataset_id')} {meta.get('dataset_version')} | "
        f"단위 nTPM | {meta.get('formula')} | {identity.get('gene_symbol')}"
    )
    figure = _chart(scored)
    if figure is not None:
        st.plotly_chart(figure, width="stretch")
    table = ranked[TABLE_COLUMNS] if not ranked.empty else scored
    st.dataframe(
        table,
        width="stretch",
        hide_index=True,
        column_config={
            "cancer": "암종",
            "n_cell_lines": "세포주 n",
            "median_nTPM": st.column_config.NumberColumn("nTPM 중앙값", format="%.2f"),
            "q1_nTPM": st.column_config.NumberColumn("Q1", format="%.2f"),
            "q3_nTPM": st.column_config.NumberColumn("Q3", format="%.2f"),
            "detected_fraction": st.column_config.NumberColumn("검출 비율", format="%.2f"),
            "selectivity_log2": st.column_config.NumberColumn("선택성 log2", format="%.2f"),
            "score": st.column_config.NumberColumn("점수", format="%.1f"),
        },
    )
    if not ranked.empty:
        st.download_button(
            "점수표 CSV",
            export_frame(ranked),
            file_name=f"{identity.get('gene_symbol')}_cell_line_score.csv",
            mime="text/csv",
            key="score_csv",
        )
    with st.expander("점수 계산"):
        st.write(formula_text())
        components = ["cancer", "level_0_1", "prevalence_0_1", "selectivity_0_1", "score", "other_median_nTPM"]
        if not ranked.empty:
            st.dataframe(ranked[components], width="stretch", hide_index=True)
        held = scored[~scored["scored"]] if not scored.empty else scored
        if not held.empty:
            st.caption("Unknown과 Non-Cancerous는 암 점수에서 빼 두었습니다.")
            st.dataframe(held[["cancer", "n_cell_lines", "median_nTPM"]], width="stretch", hide_index=True)
    lines = bundle["lines"]
    if lines is not None and not lines.empty:
        columns = [
            column
            for column in ("cell_line", "group_label", "nTPM", "annotation_status")
            if column in lines.columns
        ]
        with st.expander("세포주 목록"):
            listing = lines[columns]
            if "nTPM" in listing.columns:
                listing = listing.sort_values("nTPM", ascending=False)
            st.dataframe(listing, width="stretch", hide_index=True)
