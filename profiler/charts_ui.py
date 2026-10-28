# -*- coding: utf-8 -*-
"""Plotly charts for TAA Immune Profiler (PPT-friendly white backgrounds)."""

from __future__ import annotations

from typing import Any

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

FONT = dict(family="Source Serif 4, Georgia, serif", size=15, color="#0b1f33")
LAYOUT = dict(
    paper_bgcolor="#ffffff",
    plot_bgcolor="#ffffff",
    font=FONT,
    margin=dict(l=80, r=30, t=50, b=80),
)


def box_by_cancer(points: pd.DataFrame, title: str, y_col: str = "value_log2_plus1", unit: str = "") -> go.Figure:
    if points is None or points.empty:
        fig = go.Figure()
        fig.update_layout(title="No data", **LAYOUT)
        return fig
    fig = px.box(
        points,
        x="cancer_key",
        y=y_col,
        points=False,
        color_discrete_sequence=["#1d6f9a"],
    )
    fig.update_layout(
        title=title,
        xaxis_title="Cancer study",
        yaxis_title=unit or y_col,
        **LAYOUT,
    )
    fig.update_xaxes(tickangle=-45)
    return fig


def bar_cell_lines(df: pd.DataFrame, title: str, value_col: str = "value_plot", name_col: str = "sampleId", unit: str = "") -> go.Figure:
    if df is None or df.empty:
        fig = go.Figure()
        fig.update_layout(title="No data", **LAYOUT)
        return fig
    plot_df = df.sort_values(value_col, ascending=True)
    fig = go.Figure(
        go.Bar(
            x=plot_df[value_col],
            y=plot_df[name_col].astype(str),
            orientation="h",
            marker_color="#0e8a7d",
            hovertemplate="%{y}<br>%{x}<extra></extra>",
        )
    )
    fig.update_layout(
        title=title,
        xaxis_title=unit or value_col,
        yaxis_title="",
        height=max(360, 22 * len(plot_df) + 120),
        **LAYOUT,
    )
    fig.update_yaxes(automargin=True)
    return fig


def heatmap_markers(stats_df: pd.DataFrame, title: str) -> go.Figure:
    if stats_df is None or stats_df.empty:
        fig = go.Figure()
        fig.update_layout(title="No marker stats", **LAYOUT)
        return fig
    fig = px.imshow(
        stats_df,
        aspect="auto",
        color_continuous_scale="Teal",
        labels=dict(color="value"),
    )
    fig.update_layout(title=title, **LAYOUT)
    return fig
