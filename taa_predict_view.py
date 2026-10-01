# -*- coding: utf-8 -*-
"""Search a TAA, choose cancers or cell lines, then show the immune-activation factors."""

from __future__ import annotations

import pandas as pd
import streamlit as st

from taa_analysis.immune_score import formula_text, prepare_gene, score_lines, select_lines
from taa_analysis.pipeline import export_frame
from taa_analysis.xena import RetrievalFailed

CSS = """
<style>
@import url('https://fonts.googleapis.com/css2?family=Noto+Sans+KR:wght@400;600&family=Noto+Sans:wght@400;600&display=swap');
.stApp { background: #ffffff; }
html, body, [class*="css"] { font-family: "Noto Sans KR", "Noto Sans", sans-serif; }
.pred-score { background: #0b1f33; color: #e8eef4; border-radius: 14px; padding: 16px 18px; margin: 0.6rem 0; }
.pred-score .num { font-size: 2.2rem; font-weight: 700; color: #5eead4; }
.pred-score .meta { color: #b7c6d4; }
</style>
"""


def _chart(per_line: pd.DataFrame):
    import plotly.graph_objects as go

    shown = per_line.dropna(subset=["score"]).head(20)
    if shown.empty:
        return None
    labels = [f"{row.cell_line} · {row.cancer}" for row in shown.itertuples(index=False)]
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
        xaxis_title="세포주별 점수 0–100",
        yaxis={"categoryorder": "array", "categoryarray": list(reversed(labels))},
        margin={"l": 10, "r": 10, "t": 10, "b": 10},
        height=max(360, 26 * len(labels) + 80),
        plot_bgcolor="#ffffff",
        paper_bgcolor="#ffffff",
    )
    return figure


def _show_result(bundle: dict, result: dict, selection_label: str) -> None:
    identity = bundle["identity"]
    hindrance = bundle["hindrance"]
    per_line = result["per_line"]
    st.markdown(
        f'<div class="pred-score"><div>{identity.get("gene_symbol")} · 세포주별 면역 활성 점수 · {selection_label}</div>'
        f'<div class="meta">{identity.get("gene_id") or ""} · 세포주 {result.get("n_cell_lines") or 0}개 · '
        f"각 행이 그 세포주의 점수입니다.</div></div>",
        unsafe_allow_html=True,
    )
    detail = (
        f"결합 방해 {hindrance.get('hindrance_0_1') if hindrance.get('hindrance_0_1') is not None else '없음'} · "
        f"세포외 길이 {hindrance.get('ecd_aa') if hindrance.get('ecd_aa') is not None else '없음'} aa · "
        f"세포외 N-glycan {hindrance.get('n_glycan_extracellular') if hindrance.get('n_glycan_extracellular') is not None else '없음'} · "
        f"{hindrance.get('location_class') or 'location 없음'}. "
        "결합 방해는 단백질 값이라 모든 세포주 행에 같습니다."
    )
    if hindrance.get("reason"):
        detail = f"{detail} {hindrance['reason']}"
    st.caption(detail)
    if per_line is None or per_line.empty:
        st.info("선택한 범위에 점수를 낼 세포주가 없습니다.")
        return
    figure = _chart(per_line)
    if figure is not None:
        st.plotly_chart(figure, width="stretch")
    st.dataframe(
        per_line,
        width="stretch",
        hide_index=True,
        column_config={
            "cell_line": "세포주",
            "cancer": "암종",
            "taa_nTPM": st.column_config.NumberColumn("TAA nTPM", format="%.2f"),
            "contact_nTPM": st.column_config.NumberColumn("접촉 nTPM", format="%.2f"),
            "hindrance_0_1": st.column_config.NumberColumn("결합 방해", format="%.3f"),
            "score": st.column_config.NumberColumn("점수", format="%.1f"),
        },
    )
    st.download_button(
        "세포주 점수 CSV",
        export_frame(per_line),
        file_name=f"{identity.get('gene_symbol')}_cell_line_scores.csv",
        mime="text/csv",
        key="pred_csv",
    )
    with st.expander("점수 계산"):
        st.write(formula_text())


def render_prediction() -> None:
    st.markdown(CSS, unsafe_allow_html=True)
    st.markdown("## TAA 면역 활성 예측")
    st.caption(
        "TAA를 검색한 뒤, 암종 전체 또는 human cancer cell line을 여러 개 고릅니다. "
        "점수는 세포주마다 따로 나옵니다. TAA 발현과 세포 간 접촉은 그 세포주 값이고, 결합 방해는 같은 단백질 값입니다."
    )
    gene = st.text_input(
        "TAA / Gene",
        value="",
        placeholder="HER2, STEAP1, PD-L1",
        key="pred_gene_input",
    )
    if st.button("조회 / Load", type="primary", key="pred_load"):
        try:
            if not gene.strip():
                raise ValueError("empty gene query")
            with st.spinner("이 TAA의 암 세포주와 접촉 분자, UniProt 결합 정보를 가져오는 중"):
                st.session_state.pred_bundle = prepare_gene(gene.strip())
                st.session_state.pred_result = None
                st.session_state.pred_error = ""
                st.session_state.pred_key = gene.strip()
        except RetrievalFailed as exc:
            st.session_state.pred_error = str(exc)
            st.session_state.pred_bundle = None
        except Exception as exc:
            st.session_state.pred_error = f"retrieval_failed: {exc}"
            st.session_state.pred_bundle = None

    if st.session_state.get("pred_error"):
        st.error(f"조회에 실패했습니다. {st.session_state.pred_error}")
        return
    bundle = st.session_state.get("pred_bundle")
    if not bundle:
        st.info("TAA를 입력하고 조회를 누르세요. 그다음 암종 또는 세포주를 고릅니다.")
        return

    identity = bundle["identity"]
    lines = bundle["lines"]
    st.caption(f"{identity.get('gene_symbol')} · {identity.get('gene_id') or ''} · 암 세포주 {lines['cell_line'].nunique()}개")
    cancers = sorted(lines["group_label"].dropna().astype(str).unique())
    left, right = st.columns(2)
    with left:
        st.markdown("### 암종")
        picked_cancers = st.multiselect("암종 카테고리", cancers, key="pred_cancers")
        cancer_click = st.button("선택한 암종 전체 계산", key="pred_cancer_btn")
    with right:
        st.markdown("### Human cancer cell lines")
        pool = lines if not picked_cancers else lines[lines["group_label"].astype(str).isin(picked_cancers)]
        cell_names = sorted(pool["cell_line"].astype(str).unique())
        cell_key = "pred_cells_" + "|".join(picked_cancers)
        picked_cells = st.multiselect("세포주 복수 선택", cell_names, key=cell_key)
        cell_click = st.button("선택한 세포주 계산", key="pred_cell_btn")

    if cancer_click:
        chosen = select_lines(lines, cancers=picked_cancers or None)
        label = "선택 암종 전체" if picked_cancers else "암종 전체"
        st.session_state.pred_result = score_lines(chosen, bundle["contact"], bundle["hindrance"].get("hindrance_0_1"))
        st.session_state.pred_label = label
    if cell_click:
        if not picked_cells:
            st.warning("세포주를 한 개 이상 선택하세요.")
        else:
            chosen = select_lines(lines, cell_lines=picked_cells)
            st.session_state.pred_result = score_lines(chosen, bundle["contact"], bundle["hindrance"].get("hindrance_0_1"))
            st.session_state.pred_label = f"세포주 {len(picked_cells)}개"

    result = st.session_state.get("pred_result")
    if result:
        _show_result(bundle, result, st.session_state.get("pred_label") or "")
