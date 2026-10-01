# -*- coding: utf-8 -*-
"""Figures drawn from the same patient-level table as the comparison grid."""

from __future__ import annotations

import io

import pandas as pd

IMMUNE_COLORS = {
    "T": "#4C78A8",
    "B": "#F58518",
    "Plasma": "#E45756",
    "NK": "#72B7B2",
    "Monocyte/Macrophage": "#54A24B",
    "DC": "#EECA3B",
    "Granulocyte": "#B279A2",
    "Mast": "#FF9DA6",
    "Other/Unassigned immune": "#9D755D",
}


def _groups(patients: pd.DataFrame, summary: pd.DataFrame) -> list[tuple[str, pd.Series]]:
    if patients.empty or summary.empty:
        return []
    label_column = "group_label" if "group_label" in patients.columns else "cancer_code"
    if label_column not in patients.columns:
        label_column = "disease" if "disease" in patients.columns else patients.columns[0]
    order = summary["group_label"].tolist() if "group_label" in summary.columns else summary.iloc[:, 0].tolist()
    groups = []
    for label in order:
        values = pd.to_numeric(patients.loc[patients[label_column].astype(str) == str(label), "display_value"], errors="coerce").dropna()
        groups.append((str(label), values))
    return groups


def plotly_expression(patients: pd.DataFrame, summary: pd.DataFrame, title: str, x_title: str = "log2(TPM+1)"):
    import plotly.graph_objects as go

    figure = go.Figure()
    for label, values in _groups(patients, summary):
        if len(values) >= 2:
            figure.add_trace(
                go.Box(
                    x=values,
                    y=[label] * len(values),
                    orientation="h",
                    name=label,
                    boxpoints=False,
                    marker_color="#1f4e79",
                    line_color="#1f4e79",
                    showlegend=False,
                )
            )
        elif len(values) == 1:
            figure.add_trace(
                go.Scatter(
                    x=values,
                    y=[label],
                    mode="markers",
                    name=label,
                    marker={"size": 10, "color": "#b45309"},
                    showlegend=False,
                    hovertemplate=f"{label}: %{{x}} (n=1, point only)<extra></extra>",
                )
            )
    figure.update_layout(
        title=title,
        paper_bgcolor="white",
        plot_bgcolor="white",
        xaxis_title=x_title,
        yaxis_title="",
        height=max(420, 28 * max(1, len(summary))),
        margin={"l": 160, "r": 24, "t": 60, "b": 48},
        font={"family": "Noto Sans, Noto Sans KR, sans-serif", "color": "#1f2933"},
    )
    labels = [label for label, _ in _groups(patients, summary)]
    figure.update_yaxes(categoryorder="array", categoryarray=list(reversed(labels)))
    return figure


def plotly_subtype(patients: pd.DataFrame, summary: pd.DataFrame, title: str, x_title: str = "log2(TPM+1)"):
    import plotly.graph_objects as go

    figure = go.Figure()
    if patients.empty or "subtype_label" not in patients.columns:
        figure.update_layout(title=title, paper_bgcolor="white", plot_bgcolor="white")
        return figure
    order = summary["subtype"].tolist() if "subtype" in summary.columns else []
    for label in order:
        values = pd.to_numeric(
            patients.loc[patients["subtype_label"].astype(str) == str(label), "display_value"],
            errors="coerce",
        ).dropna()
        if len(values) >= 2:
            figure.add_trace(
                go.Violin(
                    x=values,
                    y=[label] * len(values),
                    orientation="h",
                    name=label,
                    side="positive",
                    points=False,
                    box_visible=True,
                    meanline_visible=False,
                    fillcolor="rgba(31,78,121,0.25)",
                    line_color="#1f4e79",
                    showlegend=False,
                )
            )
        elif len(values) == 1:
            figure.add_trace(
                go.Scatter(
                    x=values,
                    y=[label],
                    mode="markers",
                    marker={"size": 10, "color": "#b45309"},
                    showlegend=False,
                    name=label,
                )
            )
    figure.update_layout(
        title=title,
        paper_bgcolor="white",
        plot_bgcolor="white",
        xaxis_title=x_title,
        height=max(420, 70 * max(1, len(order))),
        margin={"l": 160, "r": 24, "t": 60, "b": 48},
        font={"family": "Noto Sans, Noto Sans KR, sans-serif", "color": "#1f2933"},
    )
    figure.update_yaxes(categoryorder="array", categoryarray=list(reversed(order)))
    return figure


def matplotlib_png(patients: pd.DataFrame, summary: pd.DataFrame, title: str, x_title: str = "log2(TPM+1)") -> bytes:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = []
    series = []
    singles = []
    for label, values in _groups(patients, summary):
        if len(values) >= 2:
            labels.append(label)
            series.append(values.to_numpy())
        elif len(values) == 1:
            singles.append((label, float(values.iloc[0])))
    height = max(4.5, 0.38 * (len(labels) + len(singles) + 1))
    figure, axis = plt.subplots(figsize=(10.5, height))
    if series:
        box_kwargs = dict(
            orientation="horizontal",
            patch_artist=True,
            boxprops={"facecolor": "#d6e3f0", "edgecolor": "#1f4e79"},
            medianprops={"color": "#1f4e79"},
        )
        # First category sits at the bottom; reverse so the highest median is at the top.
        draw_labels = list(reversed(labels))
        draw_series = list(reversed(series))
        try:
            axis.boxplot(draw_series, tick_labels=draw_labels, **box_kwargs)
        except TypeError:
            axis.boxplot(draw_series, labels=draw_labels, **box_kwargs)
    for label, value in singles:
        axis.scatter([value], [label], color="#b45309", zorder=3)
    axis.set_xlabel(x_title)
    axis.set_title(title, loc="left", fontsize=10)
    figure.tight_layout()
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=120)
    plt.close(figure)
    return buffer.getvalue()
