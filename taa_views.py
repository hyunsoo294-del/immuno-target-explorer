# -*- coding: utf-8 -*-
"""Five TAA views. Tables and figures use one filtered data object."""

from __future__ import annotations

import pandas as pd
import streamlit as st

from taa_analysis.charts import matplotlib_png, plotly_expression, plotly_subtype
from taa_analysis.cohort import MODES
from taa_analysis.molecules import MOLECULE_PANEL
from taa_analysis.pipeline import (
    available_subtype_systems,
    expression_bundle,
    export_frame,
    healthy_blood_reference,
    tumor_molecule_result,
)
from taa_analysis.xena import RetrievalFailed

MODE_LABELS = {
    "solid_primary": "고형암 원발 종양 / Solid primary tumor",
    "metastatic": "전이 / Metastatic",
    "hematologic": "혈액암 / Hematologic",
    "normal_reference": "정상조직 참고 / Normal tissue reference",
    "cell_line": "세포주 / Cell line",
}

CSS = """
<style>
@import url('https://fonts.googleapis.com/css2?family=Noto+Sans+KR:wght@400;600&family=Noto+Sans:wght@400;600&display=swap');
.stApp { background: #ffffff; }
html, body, [class*="css"] { font-family: "Noto Sans KR", "Noto Sans", sans-serif; }
.taa-context {
    background: #f4f7fb;
    border: 1px solid #d7e1ea;
    border-radius: 10px;
    padding: 0.7rem 0.9rem;
    margin: 0.4rem 0 0.8rem 0;
    color: #1f2933;
}
</style>
"""


def _caption(meta: dict, mode: str, extra: str = "") -> str:
    return (
        f"{meta.get('source_name')} | {meta.get('dataset_id')} {meta.get('dataset_version')} | "
        f"{meta.get('unit_display')} | {meta.get('formula')} | mode={MODE_LABELS.get(mode, mode)} | "
        f"retrieved {meta.get('retrieved_at')} | bulk RNA, not surface HER2 and not an IHC score. {extra}"
    )


def _axis_title(meta: dict) -> str:
    return str(meta.get("unit_display") or "log2(TPM+1)")


def _downloads(name: str, frame: pd.DataFrame, patients: pd.DataFrame, meta: dict, mode: str) -> None:
    if frame is None or frame.empty:
        return
    st.download_button(
        "표 CSV / Table CSV",
        export_frame(frame),
        file_name=f"{name}_summary.csv",
        mime="text/csv",
        key=f"csv_{name}",
    )
    if patients is not None and not patients.empty and "display_value" in patients.columns:
        st.download_button(
            "환자 값 CSV / Patient values CSV",
            export_frame(patients),
            file_name=f"{name}_patients.csv",
            mime="text/csv",
            key=f"csv_patients_{name}",
        )
        title = _caption(meta, mode)
        st.download_button(
            "그림 PNG / Figure PNG",
            matplotlib_png(patients, frame, title, x_title=_axis_title(meta)),
            file_name=f"{name}.png",
            mime="image/png",
            key=f"png_{name}",
        )


def _show_figure(figure, patients: pd.DataFrame, summary: pd.DataFrame) -> None:
    if summary.empty or "distribution" not in summary.columns:
        st.info("표시할 분포가 없습니다. / No distribution to draw.")
        return
    if (summary["distribution"] == "unavailable").all():
        st.info("요약값만 있거나 값이 없습니다. 분포 그래프를 만들지 않습니다. / No distribution is drawn without values.")
        return
    st.plotly_chart(figure, width="stretch")
    single = summary[summary["distribution"] == "point"]
    if not single.empty:
        st.caption(
            "n=1 그룹은 점만 표시합니다. violin이나 error bar를 만들지 않습니다. / "
            "Groups with one value are points, not fabricated distributions."
        )
    _ = patients


def render_taa_analysis() -> None:
    st.markdown(CSS, unsafe_allow_html=True)
    st.markdown("## TAA 분석 / TAA analysis")
    st.caption(
        "유전자 → 암종 → 아형 → 면역세포 → 분자표. "
        "환자 종양 bulk RNA와 정상 혈액 참고는 같은 축에 섞지 않습니다."
    )
    gene = st.text_input("유전자 / Gene", value=st.session_state.get("taa_gene", "ERBB2"), key="taa_gene_input")
    mode = st.selectbox(
        "공통 필터 / Shared specimen filter",
        list(MODES),
        format_func=lambda item: MODE_LABELS[item],
        key="taa_mode",
    )
    load_col, retry_col = st.columns(2)
    load = load_col.button("조회 / Load", type="primary")
    retry = retry_col.button("다시 받기 / Retry fetch")
    state_key = f"{gene.strip()}|{mode}"
    if st.session_state.get("taa_key") not in (None, state_key) and not load and not retry:
        st.session_state.taa_bundle = None
        st.session_state.taa_error = ""
        st.info("유전자 또는 필터가 바뀌었습니다. 하위 선택을 비우고 조회를 다시 누르세요. / Upstream selection changed. Load again.")
        return
    if load or retry:
        try:
            with st.spinner("Xena에서 이 유전자의 bulk RNA slice와 cohort metadata를 가져오는 중"):
                st.session_state.taa_bundle = expression_bundle(gene.strip(), mode, force=bool(retry))
                st.session_state.taa_key = state_key
                st.session_state.taa_error = ""
                st.session_state.taa_cancer_loaded = None
                st.session_state.taa_system_loaded = None
        except RetrievalFailed as exc:
            st.session_state.taa_error = str(exc)
            st.session_state.taa_bundle = None
        except Exception as exc:
            st.session_state.taa_error = f"retrieval_failed: {exc}"
            st.session_state.taa_bundle = None

    if st.session_state.get("taa_error"):
        st.error(
            "조회에 실패했습니다. 원인: "
            f"{st.session_state.taa_error}. "
            "정상 혈액값이나 예시 숫자로 대체하지 않았습니다. 다시 받기를 사용하세요. / "
            "Fetch failed. No healthy-blood or demo values were substituted."
        )
        return

    bundle = st.session_state.get("taa_bundle")
    if not bundle:
        st.info("유전자를 입력하고 조회를 누르세요. HER2와 ERBB2는 같은 target으로 처리됩니다.")
        return

    identity = bundle["identity"]
    meta = bundle["meta"]
    summary = bundle["summary"]
    codes = [""] + [code for code in summary.get("cancer_code", pd.Series(dtype=str)).dropna().astype(str).tolist() if code]
    if "BRCA" in codes:
        default_index = codes.index("BRCA")
    else:
        default_index = 0
    cancer_label = (
        "세포주 암종 / Cell-line cancer"
        if mode == "cell_line"
        else "암종 / Cancer (아형은 선택한 암종 안에서만)"
    )
    cancer = st.selectbox(cancer_label, codes, index=default_index, key="taa_cancer_select")
    systems = available_subtype_systems(cancer) if cancer else []
    subtype_system = ""
    if systems:
        subtype_system = st.selectbox("아형 분류 / Subtype system", systems, key="taa_subtype_system_select")
    elif cancer:
        st.caption(f"{cancer}: 아형 annotation이 없습니다. subtype_missing. BRCA 전체를 아형으로 대체하지 않습니다.")
    cell = st.selectbox(
        "면역세포 / Immune cell",
        [""] + ["T", "CD8", "CD4 conventional", "Treg", "B", "Plasma", "NK", "Monocyte/Macrophage", "DC", "Granulocyte", "Mast", "Other/Unassigned immune"],
        key="taa_cell_select",
    )
    if cancer != st.session_state.get("taa_cancer_loaded") or subtype_system != st.session_state.get("taa_system_loaded"):
        if cancer or subtype_system:
            try:
                bundle = expression_bundle(identity["gene_symbol"], mode, cancer or None, subtype_system or None)
                st.session_state.taa_bundle = bundle
                st.session_state.taa_cancer_loaded = cancer
                st.session_state.taa_system_loaded = subtype_system
            except RetrievalFailed as exc:
                st.error(str(exc))
                return
    bundle = st.session_state.taa_bundle
    subtype_label = ""
    if not bundle["subtype_summary"].empty and "subtype" in bundle["subtype_summary"].columns:
        labels = [item for item in bundle["subtype_summary"]["subtype"].tolist() if item and item != "Unknown"]
        if labels:
            subtype_label = st.selectbox("아형 라벨 / Subtype label", [""] + labels, key="taa_subtype_label")

    st.markdown(
        f'<div class="taa-context">유전자 {identity.get("gene_symbol")} '
        f'({identity.get("gene_id") or "gene id from probe"}) · '
        f'필터 {MODE_LABELS.get(mode, mode)} · 암종 {cancer or "전체 / all"} · '
        f'아형 {subtype_system or "적용 안 함"} {subtype_label} · 면역세포 {cell or "선택 안 함"}</div>',
        unsafe_allow_html=True,
    )
    st.caption(_caption(meta, mode))

    tab1, tab2, tab3, tab4, tab5 = st.tabs(
        [
            "1 암종별 TAA / TAA by cancer",
            "2 아형별 TAA / TAA by subtype",
            "3 암종별 면역구성 / Immune by cancer",
            "4 아형별 면역구성 / Immune by subtype",
            "5 면역세포 분자 / Molecules",
        ]
    )
    with tab1:
        if mode == "cell_line":
            st.markdown(
                "**HPA cell-line RNA (nTPM).** 환자 종양 bulk RNA가 아니고, "
                "HPA enrichment 요약 필드도 아닙니다. 화면 값은 log2(nTPM+1)입니다."
            )
            note = meta.get("cohort_note") or ""
            if note:
                st.caption(note)
        else:
            st.markdown("**Bulk RNA from tumor tissue.** 악성세포 표면 HER2 밀도, IHC 3+, 임상 HER2 양성률이 아닙니다.")
        _show_figure(
            plotly_expression(
                bundle["patients"],
                bundle["summary"],
                _caption(meta, mode, "one cell line" if mode == "cell_line" else "patient-level"),
                x_title=_axis_title(meta),
            ),
            bundle["patients"],
            bundle["summary"],
        )
        st.dataframe(bundle["summary"], width="stretch", hide_index=True)
        if mode == "cell_line" and not bundle["patients"].empty:
            line_columns = [
                column
                for column in (
                    "cell_line",
                    "group_label",
                    "nTPM",
                    "display_value",
                    "annotation_status",
                    "model_id",
                    "sample_site",
                    "tissue_status",
                    "transform_error",
                )
                if column in bundle["patients"].columns
            ]
            with st.expander("세포주 목록 / Cell lines"):
                st.dataframe(bundle["patients"][line_columns], width="stretch", hide_index=True)
        _downloads("cancer", bundle["summary"], bundle["patients"], meta, mode)
        with st.expander("커버리지 / Coverage"):
            st.dataframe(bundle["coverage"], width="stretch", hide_index=True)

    with tab2:
        if not cancer:
            st.info("암종을 선택하면 그 암종의 아형만 계산합니다. 전체 요약에 아형을 강제하지 않습니다.")
        elif bundle["subtype_summary"].empty:
            st.info("아형 결과가 없습니다.")
        else:
            missing = ""
            if "missing_reason" in bundle["subtype_summary"].columns:
                missing = str(bundle["subtype_summary"]["missing_reason"].iloc[0])
            if missing == "subtype_missing":
                st.warning(f"{cancer}: subtype_missing. 이 코호트에 연결된 annotation이 없어 아형 그림을 만들지 않습니다.")
            else:
                _show_figure(
                    plotly_subtype(
                        bundle["subtype_patients"],
                        bundle["subtype_summary"],
                        _caption(meta, mode, subtype_system or ""),
                        x_title=_axis_title(meta),
                    ),
                    bundle["subtype_patients"],
                    bundle["subtype_summary"],
                )
            st.dataframe(bundle["subtype_summary"], width="stretch", hide_index=True)
            _downloads("subtype", bundle["subtype_summary"], bundle["subtype_patients"], meta, mode)
            st.caption("PAM50 HER2-enriched와 임상 HER2-positive는 다른 분류입니다. Basal-like와 TNBC도 같습니다로 두지 않습니다.")

    immune = bundle["immune"]
    with tab3:
        st.warning(immune["detail"])
        st.caption(f"상태 {immune['status']} · 이유 {immune['missing_reason']} · 분모 all-cell 과 within-immune 은 분리됩니다.")
        st.dataframe(immune["table"], width="stretch", hide_index=True)

    with tab4:
        st.warning("아형별 면역구성은 같은 환자 immune table을 subtype metadata와 결합해야 합니다. 현재 immune table이 없어 암종 평균을 아형에 복사하지 않습니다.")
        st.dataframe(immune["table"], width="stretch", hide_index=True)

    with tab5:
        molecules = tumor_molecule_result(cancer, subtype_label or subtype_system, cell)
        st.warning(molecules["detail"])
        st.markdown("기본 분자 패널 / Default molecule panel")
        st.dataframe(pd.DataFrame(MOLECULE_PANEL), width="stretch", hide_index=True)
        st.caption("HLA-DRA 한 유전자를 MHC-II 단백질 총량으로 표기하지 않습니다. RNA 검출률을 FACS 양성률로 표기하지 않습니다.")
        show_ref = st.checkbox(
            "건강인 혈액 HPA 참고만 따로 보기 / Show healthy-blood HPA reference (not tumor TIL)",
            value=False,
            key="taa_blood_reference",
        )
        if show_ref:
            try:
                reference = healthy_blood_reference()
                if cell:
                    st.caption(
                        f"선택한 bulk/scRNA 세포명 `{cell}` 과 HPA 건강 혈액 세포명은 같은 taxonomy가 아닙니다. "
                        "자동으로 종양 TIL에 연결하지 않고 참고표를 옆에 둡니다."
                    )
                st.dataframe(
                    reference[
                        [
                            "molecule",
                            "gene_symbol",
                            "cell_type",
                            "rna_summary_nTPM",
                            "unit",
                            "evidence_status",
                            "specimen_type",
                            "cancer_type",
                            "protein_evidence",
                            "limitation",
                            "source_url",
                        ]
                    ],
                    width="stretch",
                    hide_index=True,
                )
                picked = reference[(reference["gene_symbol"] == "TNFRSF9") & (reference["cell_type"] == "memory CD8 T-cell")]
                pdcd1 = reference[(reference["gene_symbol"] == "PDCD1") & (reference["cell_type"] == "memory CD8 T-cell")]
                if not picked.empty:
                    st.caption(
                        f"memory CD8 T-cell TNFRSF9 {picked.iloc[0]['rna_summary_nTPM']} nTPM, "
                        f"PDCD1 {pdcd1.iloc[0]['rna_summary_nTPM'] if not pdcd1.empty else 'unavailable'} nTPM. "
                        "건강인 혈액 참고이며 BRCA TIL이 아닙니다. 분포가 아니라 점 요약입니다."
                    )
            except Exception as exc:
                st.error(f"HPA 참고 조회 실패 / reference fetch failed: {exc}. 종양 결과에 대체값을 넣지 않았습니다.")
