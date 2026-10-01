"""Single scrolling page. No tabs, no sidebar, no modals.

streamlit run app.py
"""

from __future__ import annotations

import html
import json

import numpy as np
import pandas as pd
import streamlit as st

from immunoscore.app.components.results import render_results
from immunoscore.src.calibration import attach_calibration
from immunoscore.src.config_loader import combination_key, load_configs
from immunoscore.src.genes import gene_identity
from immunoscore.src.scoring import build_cfg, score_models

PAGE_CSS = """
<style>
    [data-testid="stSidebar"] { display: none; }
    [data-testid="collapsedControl"] { display: none; }
    .block-container { max-width: 1180px; padding-top: 1.2rem; }
    h1 { font-size: 1.7rem; }
</style>
"""


def _tpm(log2_tpm_plus_1: float) -> str:
    if pd.isna(log2_tpm_plus_1):
        return "NA"
    tpm = float(np.exp2(log2_tpm_plus_1) - 1.0)
    if tpm >= 100:
        return f"{tpm:.0f}"
    return f"{tpm:.1f}"


@st.cache_data(show_spinner=False)
def _cached_cfg() -> dict:
    # The expression frame is not hashable for the cache key; this wrapper has no args.
    return build_cfg()


def _ranked_weights(weights: dict, labels: dict, order: list, digits: int | None) -> str:
    ranked = sorted(weights.items(), key=lambda item: (-item[1], order.index(item[0]) if item[0] in order else 99))
    parts = []
    for factor, weight in ranked:
        if float(weight) <= 0:
            continue
        label = labels.get(factor, factor)
        number = f"{float(weight):.{digits}f}" if digits is not None else f"{float(weight):g}"
        parts.append(f"{label} {number}")
    return " / ".join(parts)


def _weight_caption(arm: str, effector: str, jurkat_pd1: bool, effective: dict | None = None) -> tuple[str, str]:
    """Return (visible caption, nominal tooltip). Tooltip is empty before a result exists."""
    configs = load_configs()
    key = combination_key(arm, effector, jurkat_pd1)
    spec = configs["weights"]["combinations"][key]
    order = configs["weights"]["factor_display_order"]
    labels = configs["weights"]["factor_labels"]
    readout = configs["weights"]["readout_labels"].get(spec["readout"], spec["readout"])
    nominal = _ranked_weights(spec["weights"], labels, order, None)
    fitted = str(spec.get("fitted_on") or "")
    fit_tag = "미적합" if "not fitted" in fitted else "실측 적합"
    if effective:
        shown = _ranked_weights(effective, labels, order, 1)
        text = f"{arm} × {effector} — {readout} readout · 유효 가중치: {shown} · × 발현게이트 · {fit_tag}"
        return text, f"명목 가중치: {nominal}"
    text = f"{arm} × {effector} — {readout} readout · 가중치: {nominal} · × 발현게이트 · {fit_tag}"
    return text, ""


def _nominal_weights(arm: str, effector: str, jurkat_pd1: bool) -> dict:
    configs = load_configs()
    key = combination_key(arm, effector, jurkat_pd1)
    return dict(configs["weights"]["combinations"][key]["weights"])


def _jurkat_notice(arm: str, effector: str, jurkat_pd1: bool) -> str:
    if not str(effector).startswith("Jurkat"):
        return ""
    spec = load_configs()["weights"]["combinations"][combination_key(arm, effector, jurkat_pd1)]
    fitted = str(spec.get("fitted_on") or "")
    if "not fitted" in fitted:
        return (
            "이 Jurkat 조합은 실측으로 맞추지 않았습니다. 발현을 게이트로 바꾼 구조만 적용했습니다. "
            "4-1BB Jurkat NF-κB에 맞춘 가중치를 여기로 옮기지 않았습니다."
        )
    return (
        "이 가중치는 HER2 실측 4점으로 맞춘 작업 가설입니다. 검증된 모델이 아닙니다. "
        "이 패널에서는 항원 밀도와 세포 형태가 반대로 놓여 있어, 높은 밀도가 Emax를 낮추는지와 "
        "Calu-3의 접촉이 더 나은지를 구분할 수 없습니다."
    )


def _runtime(cfg: dict, identity: dict, settings: dict) -> dict:
    runtime = dict(cfg)
    runtime.update(settings)
    runtime["taa_entry"] = identity["curated_entry"]
    runtime["resolution_flags"] = list(identity.get("flags") or [])
    if identity.get("confidence_override"):
        runtime["confidence_override"] = identity["confidence_override"]
    return runtime


def _score(cfg: dict, models: pd.DataFrame, identity: dict, arm: str, effector: str, settings: dict):
    runtime = _runtime(cfg, identity, settings)
    result = score_models(models, identity["gene_symbol"], arm, effector, runtime)
    result, calibration = attach_calibration(result, identity["gene_symbol"], arm, effector, runtime)
    return result, calibration


def _same_request(request: dict, symbol: str, arm: str, effector: str, jurkat_pd1: bool) -> bool:
    return bool(
        request
        and request.get("gene") == symbol
        and request.get("arm") == arm
        and request.get("effector") == effector
        and bool(request.get("jurkat_pd1")) == bool(jurkat_pd1)
    )


def render() -> None:
    st.set_page_config(page_title="TAA 면역 활성 예측", layout="wide", initial_sidebar_state="collapsed")
    st.markdown(PAGE_CSS, unsafe_allow_html=True)
    st.title("TAA 면역 활성 예측")
    st.caption("TAA를 검색한 뒤, 암종 전체 또는 human cancer cell line을 여러 개 고릅니다. 점수는 세포주마다 따로 나옵니다.")

    try:
        cfg = _cached_cfg()
    except FileNotFoundError as exc:
        st.error(str(exc))
        return

    models = cfg["models"]
    expression = cfg["expression"]
    params = cfg["params"]
    arms = cfg["effectors"]["arms"]
    manifest = cfg.get("manifest") or {}

    query = st.text_input("TAA / Gene", value=st.session_state.get("gene_query", "HER2"))
    if st.button("조회 / Load", type="primary"):
        st.session_state["gene_query"] = query
        try:
            st.session_state["identity"] = gene_identity(query)
            st.session_state["result"] = None
        except ValueError as exc:
            st.session_state["identity"] = None
            st.error(str(exc))

    identity = st.session_state.get("identity")
    if not identity:
        st.caption(manifest.get("release") or "DepMap")
        return

    if identity.get("refused"):
        st.warning(identity.get("warning") or "이 항원은 유전자 점수로 계산하지 않습니다.")
        genes = identity.get("biosynthesis_genes") or []
        if genes:
            st.caption("생합성 유전자는 참고용입니다. 점수 표는 만들지 않습니다: " + ", ".join(genes))
        return

    symbol = identity["gene_symbol"]
    ensembl = identity["gene_id"] or "Ensembl 없음"
    n_lines = int(models["ModelID"].isin(expression.index).sum())
    query_text = str(identity.get("query") or symbol).strip()
    head = f"{query_text} → {symbol}" if query_text.upper() != symbol.upper() else symbol
    route = identity.get("route_label") or ""
    status = f"{head} · {ensembl} · 암 세포주 {n_lines:,}개"
    if route:
        status += f" · {route}"
    warnings = []
    if identity.get("warning"):
        warnings.append(identity["warning"])
    if not identity["curated"]:
        warnings.append("큐레이션 항목 없음")
    if cfg.get("surfaceome") is None:
        warnings.append("표면단백 참조 파일 없음")
    if warnings:
        status += " · " + " · ".join(warnings)
    st.markdown(f"**{status}**")
    release = manifest.get("release")
    if release:
        st.caption(f"발현은 {release} log2(TPM+1) 입니다. HPA는 이 화면의 1차 입력이 아닙니다.")

    st.divider()
    arm_row = st.container()
    caption_box = st.empty()
    scope_slot = st.container()
    advanced_slot = st.container()
    if "arm" not in st.session_state:
        st.session_state["arm"] = "CD3"
    with arm_row:
        arm_col, effector_col = st.columns(2)
        arm = arm_col.selectbox("면역 타겟", list(arms), key="arm")
        previous = st.session_state.get("_prev_arm")
        if previous != arm:
            st.session_state["effector"] = arms[arm]["default_effector"]
            st.session_state["_prev_arm"] = arm
        effector = effector_col.selectbox("효과기 세포", arms[arm]["effectors"], key="effector")
    with advanced_slot:
        with st.expander("고급 설정", expanded=False):
            duration = st.slider(
                "assay_duration_h",
                min_value=int(params["assay_duration_h_min"]),
                max_value=int(params["assay_duration_h_max"]),
                value=int(params["assay_duration_h_default"]),
            )
            abc50 = st.number_input(
                "ABC50_base",
                min_value=1.0,
                value=float(params["f1_expression"]["abc50_base"]),
            )
            hill_text = st.text_input("Hill h override", value="", help="비우면 효과기 YAML의 h를 씁니다.")
            jurkat_pd1 = st.checkbox("Jurkat PD-1 리포터", value=False)
    hill_override = None
    if str(hill_text).strip():
        hill_override = float(hill_text)
    settings = {
        "assay_duration_h": duration,
        "abc50_base": abc50,
        "hill_h_override": hill_override,
        "jurkat_pd1": jurkat_pd1,
    }
    pd1_on = bool(jurkat_pd1 and str(effector).startswith("Jurkat"))

    def _paint_caption(effective: dict | None) -> None:
        text, tip = _weight_caption(arm, effector, pd1_on, effective)
        if tip:
            caption_box.markdown(
                f'<p title="{html.escape(tip, quote=True)}" style="color:#6b7280;font-size:0.875rem;margin:0;">{html.escape(text)}</p>',
                unsafe_allow_html=True,
            )
        else:
            caption_box.caption(text)

    held = st.session_state.get("request") or {}
    held_result = st.session_state.get("result")
    held_effective = None
    if _same_request(held, symbol, arm, effector, jurkat_pd1) and held_result is not None and len(held_result):
        raw = held_result.iloc[0].get("effective_weights")
        if raw:
            held_effective = json.loads(raw)
    _paint_caption(held_effective)

    gene_values = expression[symbol]
    diseases = sorted(models["cancer"].dropna().astype(str).unique())
    with scope_slot:
        left, right = st.columns(2)
        with left:
            st.markdown("**암종**")
            chosen_diseases = st.multiselect("암종 카테고리", diseases, placeholder="Choose options")
            disease_models = models[models["cancer"].isin(chosen_diseases)]
            cap = int(params["cancer_line_cap_warn"])
            use_cap = False
            top_n = cap
            if len(disease_models) > cap:
                st.warning(f"선택한 암종에 세포주가 {len(disease_models)}개 있습니다. 약 {cap}개를 넘습니다.")
                use_cap = st.checkbox(f"TAA 발현 상위 {cap}개만 계산", value=True)
                if use_cap:
                    top_n = int(
                        st.number_input("상위 N", min_value=1, max_value=int(len(disease_models)), value=cap)
                    )
            cancer_clicked = st.button("선택한 암종 전체 계산")
        with right:
            st.markdown("**Human cancer cell lines**")
            order = gene_values.reindex(models["ModelID"]).sort_values(ascending=False)
            ordered = models.set_index("ModelID").loc[order.index]
            labels = {}
            for model_id, meta in ordered.iterrows():
                cancer = str(meta["cancer"])
                short = cancer.split(" ")[0]
                labels[model_id] = f"{meta['CellLineName']} — {short} — TPM {_tpm(order.loc[model_id])}"
            chosen_lines = st.multiselect(
                "세포주 복수 선택",
                options=list(ordered.index),
                format_func=lambda model_id: labels.get(model_id, model_id),
                placeholder="Choose options",
            )
            lines_clicked = st.button("선택한 세포주 계산")

    if cancer_clicked or lines_clicked:
        if cancer_clicked:
            chosen = disease_models
            scope = "선택한 암종"
            if chosen.empty:
                st.warning("암종을 먼저 고치세요.")
                return
            if use_cap:
                ranked_ids = gene_values.reindex(chosen["ModelID"]).sort_values(ascending=False).head(top_n).index
                chosen = chosen[chosen["ModelID"].isin(ranked_ids)]
                scope = f"선택한 암종 상위 {len(chosen)}"
        else:
            chosen = models[models["ModelID"].isin(chosen_lines)]
            scope = "선택한 세포주"
            if chosen.empty:
                st.warning("세포주를 먼저 고치세요.")
                return
        with st.spinner("세포주별 점수를 계산하는 중"):
            result, calibration = _score(cfg, chosen, identity, arm, effector, settings)
        effective = {}
        if len(result) and result.iloc[0].get("effective_weights"):
            effective = json.loads(result.iloc[0]["effective_weights"])
        st.session_state["result"] = result
        st.session_state["request"] = {
            "gene": symbol,
            "display_gene": identity["query"],
            "arm": arm,
            "effector": effector,
            "scope_label": scope,
            "jurkat_pd1": jurkat_pd1,
            "assay_duration_h": duration,
            "abc50_base": abc50,
            "manifest": {key: manifest.get(key) for key in ("release", "release_id", "doi", "downloaded_at", "note")},
            "confidence_labels": params["confidence_labels"],
            "jurkat_notice": _jurkat_notice(arm, effector, pd1_on),
            "effective_weights": effective,
            "nominal_weights": _nominal_weights(arm, effector, pd1_on),
            "calibration": {
                key: calibration.get(key)
                for key in ("ok", "caption", "reason", "a", "b", "r2", "n", "fit_date", "taa", "arm", "effector")
            },
        }
        _paint_caption(effective or None)

    result = st.session_state.get("result")
    request = st.session_state.get("request") or {}
    if result is not None and _same_request(request, symbol, arm, effector, jurkat_pd1):
        inert = []
        if "inert_factors" in result.columns and len(result):
            inert = [item for item in str(result.iloc[0]["inert_factors"]).split(",") if item]
        render_results(result, request, inert)


def main() -> None:
    render()


if __name__ == "__main__":
    main()
