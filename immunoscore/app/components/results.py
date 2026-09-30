"""Inline results table, row detail, and Excel download."""

from __future__ import annotations

import pandas as pd
import streamlit as st

from immunoscore.src.report import DISPLAY_COLUMNS, display_frame, to_excel

INERT_LABELS = {
    "F4_checkpoint": "체크포인트",
    "F5_adhesion": "접촉",
}


def render_results(result: pd.DataFrame, request: dict, inert_factors: list[str]) -> None:
    if result is None or result.empty:
        return
    header = (
        f"결과 · {request['display_gene']} × {request['arm']} × {request['effector']}"
        f" · {request['scope_label']} · {len(result)}개 세포주"
    )
    left, right = st.columns([4, 1])
    left.markdown(f"**{header}**")
    right.download_button(
        "Excel 다운로드",
        data=to_excel(result, request),
        file_name=f"{request['gene']}_{request['arm']}_{request['effector']}.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    )
    st.caption("예비 패널 선택용 휴리스틱입니다. 검증된 Emax 예측기가 아닙니다. 각 행은 그 세포주 하나의 점수입니다.")
    shown = display_frame(result)
    grey_labels = [INERT_LABELS[factor] for factor in inert_factors if factor in INERT_LABELS]
    styler = shown.style
    if grey_labels:
        styler = styler.map(lambda _value: "color: #9aa0a6", subset=grey_labels)
        reason = "이 효과기에서 정보가 없어 회색으로 둔 열: " + ", ".join(grey_labels)
        if "F5_adhesion" in inert_factors:
            reason += ". HEK293에는 종양 세포 ICAM1/CD58의 대응 수용체가 없습니다."
        if "F4_checkpoint" in inert_factors:
            reason += " 억제 수용체가 없는 효과기라 체크포인트는 100으로 두고 점수에 거의 반영하지 않습니다."
        st.caption(reason)
    styler = styler.background_gradient(subset=["점수"], cmap="YlGn", vmin=0, vmax=100)
    styler = styler.format(precision=1, na_rep="")
    st.dataframe(styler, use_container_width=True, hide_index=True, height=460)
    picked = st.selectbox(
        "행을 고르면 근거가 열립니다",
        options=result["cell_line"].tolist(),
        index=0,
    )
    detail = result.loc[result["cell_line"] == picked].iloc[0]
    with st.container(border=True):
        st.markdown(f"**{detail['cell_line']} · {detail['cancer']} · 점수 {detail['score']:.1f}**")
        st.caption(
            f"신뢰도 {request['confidence_labels'].get(detail['confidence'], detail['confidence'])}"
            f" · 데이터 완전도 {detail['data_completeness']:.2f}"
            f" · 플래그 {detail['flag_text'] or '없음'}"
        )
        for source, label in DISPLAY_COLUMNS:
            if not str(source).startswith("F"):
                continue
            value = detail[source]
            text = "" if pd.isna(value) else f"{float(value):.1f}"
            st.markdown(f"{label} ({source}) · {text} · {detail[source + '_tier']}")
            st.caption(str(detail[source + "_detail"]))
