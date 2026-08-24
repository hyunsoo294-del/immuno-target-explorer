# -*- coding: utf-8 -*-
"""Publication-style Plotly charts for biotech target expression."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import plotly.graph_objects as go

NAVY = "#0B1F33"
TEAL = "#1AA6A6"
CORAL = "#FF5A5F"
GOLD = "#E3B23C"
SLATE = "#5B6C7D"
MINT = "#7CDEB3"
LAVENDER = "#9B8AFB"
INK = "#1F2937"
PAPER = "#F7F5F2"
GRID = "#E6E2DC"

LAYOUT_BASE = dict(
    template="plotly_white",
    font=dict(family="Arial, Helvetica, sans-serif", size=13, color=INK),
    paper_bgcolor=PAPER,
    plot_bgcolor=PAPER,
    hoverlabel=dict(bgcolor="white", font_size=12, font_family="Arial"),
)


def _chart_height(n: int, per_row: int = 34, minimum: int = 380) -> int:
    return max(minimum, n * per_row + 180)


def _apply_axes(
    fig: go.Figure,
    *,
    horizontal: bool = False,
    height: int | None = None,
    legend: str = "none",
    left_margin: int | None = None,
) -> go.Figure:
    """Apply shared styling. legend: none | bottom | top."""
    margin = dict(
        l=left_margin if left_margin is not None else (150 if horizontal else 60),
        r=80 if horizontal else 24,
        t=72,
        b=80,
    )
    layout = dict(LAYOUT_BASE, margin=margin)

    if height:
        layout["height"] = height

    if legend == "bottom":
        layout["legend"] = dict(
            orientation="h",
            yanchor="top",
            y=-0.22,
            x=0.5,
            xanchor="center",
            bgcolor="rgba(0,0,0,0)",
            tracegroupgap=14,
            font=dict(size=12),
        )
        layout["margin"]["b"] = 130
    elif legend == "top":
        layout["margin"]["t"] = 108
        layout["title"] = dict(pad=dict(t=8, b=36))
        layout["legend"] = dict(
            orientation="h",
            yanchor="bottom",
            y=1.01,
            x=0,
            xanchor="left",
            bgcolor="rgba(0,0,0,0)",
            tracegroupgap=14,
            font=dict(size=12),
        )

    fig.update_layout(**layout)

    if horizontal:
        fig.update_xaxes(gridcolor=GRID, zeroline=True, zerolinecolor="#D6D1C9", title_standoff=10)
        fig.update_yaxes(showgrid=False, automargin=True, title_standoff=10)
    else:
        fig.update_xaxes(showgrid=False, tickangle=-35, title_standoff=10)
        fig.update_yaxes(gridcolor=GRID, zeroline=True, zerolinecolor="#D6D1C9", title_standoff=10)
    return fig


def _sort_cancers(df: pd.DataFrame, value_col: str) -> pd.DataFrame:
    if df.empty or value_col not in df.columns:
        return df
    return df.sort_values(value_col, ascending=False, na_position="last")


def _bar_text_positions(values: pd.Series) -> list[str]:
    return ["outside" if float(v) >= 0 else "inside" for v in values]


def _horizontal_value_axis(values: pd.Series, pad_ratio: float = 0.18) -> dict:
    vals = pd.to_numeric(values, errors="coerce").dropna()
    if vals.empty:
        return {}
    lo = float(vals.min())
    hi = float(vals.max())
    span = max(hi - lo, 0.5)
    pad = span * pad_ratio
    if lo < 0:
        lo -= pad
    if hi > 0:
        hi += pad
    if lo >= 0:
        lo = min(lo, 0)
    return dict(range=[lo, hi])


def grouped_rna_bar(comparison: pd.DataFrame) -> go.Figure | None:
    if comparison is None or comparison.empty:
        return None
    df = comparison.copy()
    for col in ("rna_log2_normal", "rna_log2_cancer"):
        if col not in df.columns:
            return None
    df = df.dropna(subset=["rna_log2_normal", "rna_log2_cancer"], how="all")
    if df.empty:
        return None

    sort_col = "log2_fold_normalized" if "log2_fold_normalized" in df.columns else "rna_log2_cancer"
    if sort_col in df.columns and df[sort_col].notna().any():
        df = _sort_cancers(df, sort_col).iloc[::-1]
    else:
        df = _sort_cancers(df, "rna_log2_normal").iloc[::-1]

    has_cancer = df["rna_log2_cancer"].notna().any()
    x_vals = pd.concat([df["rna_log2_normal"], df["rna_log2_cancer"]], ignore_index=True)
    x_range = _horizontal_value_axis(x_vals)
    if not x_range:
        x_range = dict(range=[0, 8])

    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            y=df["cancer_type"],
            x=df["rna_log2_normal"],
            name="Normal tissue (HPA/GTEx)",
            orientation="h",
            marker_color=TEAL,
            marker_line_width=0,
            customdata=df.get("normal_rna_abs"),
            hovertemplate=(
                "Cancer type: %{y}<br>"
                "Normal log2(TPM+1): %{x:.2f}<br>"
                "Normal TPM: %{customdata:.2f}<extra></extra>"
            ),
        )
    )
    if has_cancer:
        fig.add_trace(
            go.Bar(
                y=df["cancer_type"],
                x=df["rna_log2_cancer"],
                name="Tumor (TCGA median)",
                orientation="h",
                marker_color=CORAL,
                marker_line_width=0,
                customdata=df.get("cancer_rna_abs"),
                hovertemplate=(
                    "Cancer type: %{y}<br>"
                    "Tumor log2(RSEM+1): %{x:.2f}<br>"
                    "TCGA median: %{customdata:.1f}<extra></extra>"
                ),
            )
        )

    title = "RNA: normal tissue vs tumor (log2 scale)"
    if not has_cancer:
        title += " - tumor data unavailable"

    fig.update_layout(
        barmode="group",
        title=title,
        xaxis_title="log2(expression + 1)",
        xaxis=x_range,
        bargap=0.22,
        bargroupgap=0.08,
    )
    return _apply_axes(fig, horizontal=True, height=_chart_height(len(df)), legend="bottom")


def fold_bar(comparison: pd.DataFrame, value_col: str, title: str, y_title: str) -> go.Figure | None:
    if comparison is None or comparison.empty or value_col not in comparison.columns:
        return None
    df = comparison.dropna(subset=[value_col]).copy()
    if df.empty:
        return None
    df = _sort_cancers(df, value_col).iloc[::-1]
    fig = go.Figure(
        go.Bar(
            y=df["cancer_type"],
            x=df[value_col],
            orientation="h",
            marker_color=[CORAL if v >= 0 else TEAL for v in df[value_col]],
            marker_line_width=0,
            text=df[value_col].map(lambda v: f"{v:.2f}"),
            textposition=_bar_text_positions(df[value_col]),
            insidetextanchor="middle",
            cliponaxis=False,
        )
    )
    fig.add_vline(x=0, line_color=SLATE, line_width=1)
    fig.update_layout(title=title, xaxis_title=y_title, xaxis=_horizontal_value_axis(df[value_col]), showlegend=False)
    return _apply_axes(fig, horizontal=True, height=_chart_height(len(df)))


def protein_grouped_bar(comparison: pd.DataFrame) -> go.Figure | None:
    if comparison is None or comparison.empty:
        return None
    df = comparison.dropna(subset=["protein_ihc_cancer_0_3"]).copy()
    if df.empty:
        return None
    df = _sort_cancers(df, "protein_ihc_cancer_0_3").iloc[::-1]
    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            y=df["cancer_type"],
            x=df["protein_ihc_normal_0_3"],
            name="Normal IHC (0-3)",
            orientation="h",
            marker_color=NAVY,
            marker_line_width=0,
        )
    )
    fig.add_trace(
        go.Bar(
            y=df["cancer_type"],
            x=df["protein_ihc_cancer_0_3"],
            name="Cancer IHC (0-3)",
            orientation="h",
            marker_color=GOLD,
            marker_line_width=0,
        )
    )
    fig.update_layout(
        barmode="group",
        title="Protein IHC: matched normal organ vs tumor",
        xaxis_title="IHC score (0-3)",
        xaxis=dict(range=[0, 3.4]),
        bargap=0.22,
    )
    return _apply_axes(fig, horizontal=True, height=_chart_height(len(df)), legend="bottom")


def tcga_box(tcga: pd.DataFrame) -> go.Figure | None:
    if tcga is None or tcga.empty:
        return None
    tumor = tcga[tcga["sample_class"] == "tumor"].copy()
    tumor = tumor.dropna(subset=["cancer_type"])
    if tumor.empty:
        return None
    tumor["log2"] = (tumor["value"] + 1).map(lambda x: math.log2(x) if pd.notna(x) and x >= 0 else None)
    tumor = tumor.dropna(subset=["log2"])
    if tumor.empty:
        return None
    order = tumor.groupby("cancer_type")["log2"].median().sort_values(ascending=False).index.tolist()
    fig = go.Figure()
    fig.add_trace(
        go.Box(
            x=tumor["cancer_type"],
            y=tumor["log2"],
            marker_color=CORAL,
            line_color=NAVY,
            fillcolor="rgba(255,90,95,0.35)",
            boxpoints="outliers",
            whiskerwidth=0.6,
        )
    )
    fig.update_layout(title="TCGA tumor RNA distribution", yaxis_title="log2(RSEM + 1)", showlegend=False)
    fig.update_xaxes(categoryorder="array", categoryarray=order)
    fig = _apply_axes(fig)
    fig.update_layout(margin=dict(b=130))
    return fig


def tcga_tumor_normal_bar(tcga_summary: pd.DataFrame) -> go.Figure | None:
    if tcga_summary is None or tcga_summary.empty:
        return None
    df = tcga_summary.dropna(subset=["tumor_median_log2"]).copy()
    if df.empty:
        return None
    df = df.sort_values("tumor_median_log2", ascending=False).iloc[::-1]
    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            y=df["cancer_type"],
            x=df["tumor_median_log2"],
            name="Tumor median",
            orientation="h",
            marker_color=CORAL,
            marker_line_width=0,
        )
    )
    if df["normal_mean_log2"].notna().any():
        fig.add_trace(
            go.Bar(
                y=df["cancer_type"],
                x=df["normal_mean_log2"],
                name="Adjacent normal mean",
                orientation="h",
                marker_color=TEAL,
                marker_line_width=0,
            )
        )
    fig.update_layout(
        barmode="group",
        title="TCGA median RNA: tumor vs adjacent normal",
        xaxis_title="log2(RSEM + 1)",
        bargap=0.22,
    )
    return _apply_axes(fig, horizontal=True, height=_chart_height(len(df)), legend="bottom")


def normal_rna_bar(normal_rna: pd.DataFrame) -> go.Figure | None:
    if normal_rna is None or normal_rna.empty:
        return None
    df = normal_rna.copy()
    df["normal_rna_ntpm"] = pd.to_numeric(df["normal_rna_ntpm"], errors="coerce")
    df = df.dropna(subset=["normal_rna_ntpm"]).sort_values("normal_rna_ntpm", ascending=False)
    if df.empty:
        return None
    df["log2"] = (df["normal_rna_ntpm"] + 1).map(lambda x: math.log2(x) if pd.notna(x) and x >= 0 else None)
    df = df.iloc[::-1]
    fig = go.Figure(
        go.Bar(
            y=df["organ"],
            x=df["log2"],
            orientation="h",
            marker_color=MINT,
            marker_line_color=TEAL,
            marker_line_width=0.4,
        )
    )
    fig.update_layout(title="HPA/GTEx normal RNA by organ", xaxis_title="log2(nTPM + 1)", showlegend=False)
    return _apply_axes(fig, horizontal=True, height=_chart_height(min(len(df), 20)))


def ihc_stacked_bar(ihc: pd.DataFrame) -> go.Figure | None:
    if ihc is None or ihc.empty:
        return None
    df = _sort_cancers(ihc.copy(), "ihc_weighted_0_3").iloc[::-1]
    fig = go.Figure()
    layers = [
        ("ihc_not_detected", "Not detected (0)", SLATE),
        ("ihc_low", "Low (1)", TEAL),
        ("ihc_medium", "Medium (2)", GOLD),
        ("ihc_high", "High (3)", CORAL),
    ]
    for col, label, color in layers:
        if col in df.columns:
            fig.add_trace(
                go.Bar(
                    y=df["cancer_type"],
                    x=pd.to_numeric(df[col], errors="coerce").fillna(0),
                    name=label,
                    orientation="h",
                    marker_color=color,
                    marker_line_width=0,
                    hovertemplate=f"{label}: %{{x}} patients<extra></extra>",
                )
            )
    fig.update_layout(
        barmode="stack",
        title="IHC pathology: patient staining distribution",
        xaxis_title="Number of patients (HPA tumor TMA cohort per cancer type)",
        bargap=0.22,
    )
    fig = _apply_axes(fig, horizontal=True, height=_chart_height(len(df)) + 50, legend="bottom")
    fig.update_layout(margin=dict(b=145), legend=dict(y=-0.24))
    return fig


def ihc_weighted_bar(ihc: pd.DataFrame) -> go.Figure | None:
    if ihc is None or ihc.empty or "ihc_weighted_0_3" not in ihc.columns:
        return None
    df = _sort_cancers(ihc.dropna(subset=["ihc_weighted_0_3"]).copy(), "ihc_weighted_0_3").iloc[::-1]
    if df.empty:
        return None
    fig = go.Figure(
        go.Bar(
            y=df["cancer_type"],
            x=df["ihc_weighted_0_3"],
            orientation="h",
            marker_color=LAVENDER,
            marker_line_color=NAVY,
            marker_line_width=0.5,
            text=df["ihc_weighted_0_3"].map(lambda v: f"{v:.2f}"),
            textposition="outside",
            cliponaxis=False,
        )
    )
    fig.update_layout(
        title="IHC weighted score (0-3)",
        xaxis_title="Weighted IHC score",
        xaxis=dict(range=[0, min(3.4, float(df["ihc_weighted_0_3"].max()) * 1.3 + 0.2)]),
        showlegend=False,
    )
    return _apply_axes(fig, horizontal=True, height=_chart_height(len(df)))


def ihc_fold_bar(ihc: pd.DataFrame) -> go.Figure | None:
    if ihc is None or ihc.empty or "ihc_mwu_pvalue" not in ihc.columns:
        return None
    df = ihc.dropna(subset=["ihc_mwu_pvalue"]).copy()
    if df.empty:
        return None
    df["neglog10_p"] = (-np.log10(pd.to_numeric(df["ihc_mwu_pvalue"], errors="coerce").clip(lower=1e-300)))
    df = _sort_cancers(df, "neglog10_p").iloc[::-1]
    fig = go.Figure(
        go.Bar(
            y=df["cancer_type"],
            x=df["neglog10_p"],
            orientation="h",
            marker_color=[CORAL if bool(v) else TEAL for v in df.get("ihc_mwu_significant", [False] * len(df))],
            marker_line_width=0,
            text=df["neglog10_p"].map(lambda v: f"{v:.2f}"),
            textposition="outside",
            cliponaxis=False,
        )
    )
    fig.update_layout(
        title="IHC Mann-Whitney U: -log10(p), tumor vs matched normal",
        xaxis_title="-log10(p-value)",
        showlegend=False,
    )
    return _apply_axes(fig, horizontal=True, height=_chart_height(len(df)))


def protein_fold_horizontal(comparison: pd.DataFrame) -> go.Figure | None:
    if comparison is None or comparison.empty:
        return None
    value_col = "ihc_mwu_effect" if "ihc_mwu_effect" in comparison.columns else None
    if value_col is None or comparison[value_col].dropna().empty:
        return None
    df = comparison.dropna(subset=[value_col]).copy()
    df = _sort_cancers(df, value_col).iloc[::-1]
    fig = go.Figure(
        go.Bar(
            y=df["cancer_type"],
            x=df[value_col],
            orientation="h",
            marker_color=[CORAL if v >= 0 else TEAL for v in df[value_col]],
            marker_line_width=0,
            text=df[value_col].map(lambda v: f"{v:.2f}"),
            textposition=_bar_text_positions(df[value_col]),
            insidetextanchor="middle",
            cliponaxis=False,
        )
    )
    fig.add_vline(x=0, line_color=SLATE, line_width=1)
    fig.update_layout(
        title="IHC rank-biserial effect (tumor vs normal, Mann-Whitney)",
        xaxis_title="Effect size (>0 = higher tumor IHC rank)",
        xaxis=_horizontal_value_axis(df[value_col]),
        showlegend=False,
    )
    return _apply_axes(fig, horizontal=True, height=_chart_height(len(df)))


def cptac_box(cptac: pd.DataFrame) -> go.Figure | None:
    if cptac is None or cptac.empty:
        return None
    df = cptac.dropna(subset=["cancer_type", "value"]).copy()
    if df.empty:
        return None
    order = df.groupby("cancer_type")["value"].median().sort_values(ascending=False).index.tolist()
    fig = go.Figure()
    fig.add_trace(
        go.Box(
            x=df["cancer_type"],
            y=df["value"],
            marker_color=TEAL,
            line_color=NAVY,
            fillcolor="rgba(26,166,166,0.35)",
            boxpoints="outliers",
            whiskerwidth=0.6,
        )
    )
    fig.update_layout(
        title="CPTAC tumor protein distribution",
        yaxis_title="log2 protein abundance ratio",
        showlegend=False,
    )
    fig.update_xaxes(categoryorder="array", categoryarray=order)
    fig = _apply_axes(fig)
    fig.update_layout(margin=dict(b=130))
    return fig


def cptac_median_bar(cptac_summary: pd.DataFrame) -> go.Figure | None:
    if cptac_summary is None or cptac_summary.empty:
        return None
    df = cptac_summary.dropna(subset=["median_log2_ratio"]).copy()
    if df.empty:
        return None
    df = df.sort_values("median_log2_ratio", ascending=False).iloc[::-1]
    fig = go.Figure(
        go.Bar(
            y=df["cancer_type"],
            x=df["median_log2_ratio"],
            orientation="h",
            marker_color=TEAL,
            marker_line_color=NAVY,
            marker_line_width=0.4,
            text=df["median_log2_ratio"].map(lambda v: f"{v:.2f}"),
            textposition=_bar_text_positions(df["median_log2_ratio"]),
            insidetextanchor="middle",
            cliponaxis=False,
        )
    )
    fig.add_vline(x=0, line_color=SLATE, line_width=1)
    fig.update_layout(
        title="CPTAC median tumor protein by cancer type",
        xaxis_title="Median log2 abundance ratio",
        xaxis=_horizontal_value_axis(df["median_log2_ratio"]),
        showlegend=False,
    )
    return _apply_axes(fig, horizontal=True, height=_chart_height(len(df)))
