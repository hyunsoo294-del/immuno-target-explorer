"""Inline results table, row detail, and Excel download."""

from __future__ import annotations

import pandas as pd
import streamlit as st

from immunoscore.src.report import DISPLAY_COLUMNS, display_frame, to_excel

INERT_LABELS = {
    "F4_checkpoint": "체크포인트",
    "F5_adhesion": "접촉",
}
STRING_COLUMNS = ("세포주", "암종", "체크포인트", "접근성", "신뢰도", "플래그", "예측 Emax")


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
    notice = request.get("jurkat_notice") or ""
    if notice:
        st.markdown(notice)
    calibration = request.get("calibration") or {}
    if calibration.get("ok") and calibration.get("caption"):
        st.caption(calibration["caption"])
    else:
        st.caption("캘리브레이션 없음 — 순위만 유효")
    st.caption("발현(게이트)는 가산 점수가 아닙니다. 100은 ×1.00, 50은 ×0.50입니다. 점수는 그 배수를 변별력이 있는 요인에 곱한 값입니다.")
    shown = display_frame(result)
    zero_ids = []
    if "zero_variance_factors" in result.columns and len(result):
        zero_ids = [item for item in str(result.iloc[0]["zero_variance_factors"]).split(",") if item]
    label_by_id = {source: label for source, label in DISPLAY_COLUMNS}
    zero_labels = [label_by_id[factor] for factor in zero_ids if factor in label_by_id and label_by_id[factor] in shown.columns]
    grey_labels = [INERT_LABELS[factor] for factor in inert_factors if factor in INERT_LABELS and INERT_LABELS[factor] in shown.columns]
    styler = shown.style
    muted = list(dict.fromkeys(grey_labels + zero_labels))
    if muted:
        styler = styler.map(lambda _value: "color: #9aa0a6", subset=muted)
    if "예측 Emax" in shown.columns:
        styler = styler.map(
            lambda value: "color: #9a6700" if isinstance(value, str) and "외삽" in value else "",
            subset=["예측 Emax"],
        )
    reason_bits = []
    if grey_labels:
        reason = "이 효과기에서 정보가 없어 회색으로 둔 열: " + ", ".join(grey_labels)
        if "F5_adhesion" in inert_factors:
            reason += ". HEK293에는 종양 세포 ICAM1/CD58의 대응 수용체가 없습니다."
        if "F4_checkpoint" in inert_factors:
            reason += " 억제 수용체가 없는 효과기라 체크포인트는 — 이고 가중치는 0입니다."
        reason_bits.append(reason)
    if zero_labels:
        reason_bits.append(
            "결과 집합에서 분산이 없어 점수에서 뺀 열: " + ", ".join(zero_labels) + ". 값은 그대로 보입니다."
        )
    if reason_bits:
        st.caption(" ".join(reason_bits))
    numeric = [column for column in shown.columns if column not in STRING_COLUMNS]
    styler = styler.format(precision=1, subset=numeric, na_rep="")
    styler = styler.background_gradient(subset=["점수"], cmap="YlGn", vmin=0, vmax=100)
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
        gate = detail["F1_expression"]
        gate_text = "" if pd.isna(gate) else f"×{float(gate) / 100:.2f}"
        st.markdown(f"발현(게이트) · {gate_text} · {'' if pd.isna(gate) else f'{float(gate):.1f}'}")
        st.caption(str(detail["F1_expression_detail"]))
        emax_text = detail.get("emax_text") or ""
        if emax_text:
            st.markdown(f"예측 Emax · {emax_text}")
        access = detail.get("accessibility")
        access_text = "" if pd.isna(access) else f"{float(access):.2f}"
        st.markdown(f"접근성 · {access_text} · {detail.get('f5_accessibility_source')}")
        for source, label in DISPLAY_COLUMNS:
            if not str(source).startswith("F") or source == "F1_expression":
                continue
            value = detail[source]
            if detail.get(source + "_tier") == "not_applicable" or pd.isna(value):
                text = "—"
            else:
                text = f"{float(value):.1f}"
            st.markdown(f"{label} ({source}) · {text} · {detail[source + '_tier']}")
            st.caption(str(detail[source + "_detail"]))
        if detail.get("zero_variance_factors"):
            named = ", ".join(label_by_id.get(item, item) for item in str(detail["zero_variance_factors"]).split(",") if item)
            st.caption("이 요인들은 결과 집합에서 분산이 없어 점수에서 뺐습니다: " + named)
