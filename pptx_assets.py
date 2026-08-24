# -*- coding: utf-8 -*-
"""Biotech-themed slide imagery (protein structure, IHC microscopy style)."""

from __future__ import annotations

import io
import math
from typing import Literal

import numpy as np

Theme = Literal["cover", "rna", "protein", "tcga", "ihc", "cptac", "methods"]

NAVY = "#0B1F33"
TEAL = "#1AA6A6"
CORAL = "#FF5A5F"
MINT = "#7CDEB3"
SLATE = "#5B6C7D"
PAPER = "#F4F6F8"
WHITE = "#FFFFFF"
DAB = "#8B4A6B"
HEMATOXYLIN = "#4A5568"


def _canvas(w: int = 1920, h: int = 1080):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(w / 160, h / 160), dpi=160)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    fig.patch.set_facecolor(PAPER)
    return fig, ax, plt


def _save(fig, plt) -> bytes:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", pad_inches=0, facecolor=fig.get_facecolor())
    plt.close(fig)
    return buf.getvalue()


def _gradient(ax, c1: str, c2: str, alpha: float = 1.0) -> None:
    import matplotlib.colors as mcolors
    import matplotlib.pyplot as plt

    grad = np.linspace(0, 1, 256).reshape(256, 1)
    cmap = mcolors.LinearSegmentedColormap.from_list("bg", [c1, c2])
    ax.imshow(grad, aspect="auto", cmap=cmap, alpha=alpha, extent=(0, 1, 0, 1), origin="lower", zorder=0)


def _protein_ribbon(ax, cx: float, cy: float, scale: float, color: str, alpha: float = 0.85) -> None:
    t = np.linspace(0, 4 * math.pi, 240)
    x = cx + scale * 0.42 * np.cos(t) * (1 + 0.18 * np.sin(3 * t))
    y = cy + scale * 0.28 * np.sin(t) * (1 + 0.12 * np.cos(2 * t))
    ax.plot(x, y, color=color, linewidth=scale * 28, solid_capstyle="round", alpha=alpha, zorder=2)
    for phase in (0.6, 2.1, 3.8):
        px = cx + scale * 0.42 * math.cos(phase) * (1 + 0.18 * math.sin(3 * phase))
        py = cy + scale * 0.28 * math.sin(phase) * (1 + 0.12 * math.cos(2 * phase))
        ax.scatter([px], [py], s=scale * 520, color=WHITE, edgecolors=color, linewidths=2.2, zorder=3, alpha=0.95)


def _microscopy_texture(ax, seed: int = 7) -> None:
    rng = np.random.default_rng(seed)
    h, w = 480, 640
    base = rng.normal(0.72, 0.06, (h, w))
    for _ in range(140):
        cx, cy = rng.uniform(0, w), rng.uniform(0, h)
        r = rng.uniform(8, 34)
        yy, xx = np.ogrid[:h, :w]
        mask = (xx - cx) ** 2 + (yy - cy) ** 2 <= r**2
        base[mask] += rng.uniform(0.05, 0.22)
    dab = np.clip(base * 0.55 + rng.normal(0, 0.03, (h, w)), 0, 1)
    rgb = np.zeros((h, w, 3))
    rgb[..., 0] = 0.35 + dab * 0.45
    rgb[..., 1] = 0.28 + dab * 0.18
    rgb[..., 2] = 0.42 + dab * 0.25
    ax.imshow(rgb, extent=(0, 1, 0, 1), origin="lower", alpha=0.92, zorder=1, aspect="auto")


def _network_nodes(ax, accent: str) -> None:
    rng = np.random.default_rng(42)
    pts = rng.uniform(0.08, 0.92, (18, 2))
    for i, (x, y) in enumerate(pts):
        ax.scatter([x], [y], s=90 + (i % 3) * 40, color=accent, alpha=0.35, zorder=2)
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            if rng.random() > 0.72:
                ax.plot([pts[i, 0], pts[j, 0]], [pts[i, 1], pts[j, 1]], color=accent, alpha=0.12, linewidth=1.2, zorder=1)


def render_slide_image(theme: Theme, *, symbol: str = "") -> bytes:
    fig, ax, plt = _canvas()

    if theme == "cover":
        _gradient(ax, NAVY, "#123A5C")
        _protein_ribbon(ax, 0.68, 0.52, 1.05, TEAL, alpha=0.75)
        _protein_ribbon(ax, 0.72, 0.46, 0.82, CORAL, alpha=0.55)
        _network_nodes(ax, MINT)
        ax.add_patch(plt.Rectangle((0, 0), 0.58, 1, color=NAVY, alpha=0.82, zorder=4, transform=ax.transAxes))
        ax.text(
            0.06,
            0.82,
            "TARGET EXPRESSION",
            color=TEAL,
            fontsize=11,
            fontweight="bold",
            fontfamily="sans-serif",
            zorder=5,
        )
        ax.text(
            0.06,
            0.62,
            symbol or "GENE",
            color=WHITE,
            fontsize=34,
            fontweight="bold",
            fontfamily="sans-serif",
            zorder=5,
        )
        ax.text(
            0.06,
            0.12,
            "Modern biotech analytics  |  RNA  Protein  TCGA  IHC",
            color="#C8D6E5",
            fontsize=10,
            fontfamily="sans-serif",
            zorder=5,
        )
        return _save(fig, plt)

    if theme == "ihc":
        _gradient(ax, "#EFE8EE", PAPER)
        _microscopy_texture(ax, seed=11)
        ax.add_patch(plt.Rectangle((0, 0), 0.42, 1, color=NAVY, alpha=0.78, zorder=4, transform=ax.transAxes))
        return _save(fig, plt)

    if theme == "protein":
        _gradient(ax, PAPER, "#E8F4F4")
        _protein_ribbon(ax, 0.55, 0.5, 1.2, TEAL)
        _protein_ribbon(ax, 0.62, 0.55, 0.75, CORAL, alpha=0.6)
        ax.add_patch(plt.Rectangle((0, 0), 0.38, 1, color=NAVY, alpha=0.82, zorder=4, transform=ax.transAxes))
        return _save(fig, plt)

    accent = TEAL
    if theme == "tcga":
        accent = CORAL
    elif theme == "cptac":
        accent = MINT
    elif theme == "methods":
        accent = SLATE

    _gradient(ax, WHITE, PAPER)
    _network_nodes(ax, accent)
    ax.add_patch(plt.Rectangle((0, 0), 0.36, 1, color=NAVY, alpha=0.85, zorder=4, transform=ax.transAxes))
    if theme == "rna":
        _protein_ribbon(ax, 0.7, 0.55, 0.65, TEAL, alpha=0.45)
    return _save(fig, plt)


def render_side_accent(theme: Theme) -> bytes:
    """Narrow vertical strip for content slide margins."""
    fig, ax, plt = _canvas(w=420, h=1080)
    if theme == "ihc":
        _microscopy_texture(ax, seed=3)
    else:
        _gradient(ax, NAVY, "#16324A")
        if theme in ("protein", "cover"):
            _protein_ribbon(ax, 0.5, 0.55, 0.55, TEAL, alpha=0.7)
        else:
            _network_nodes(ax, TEAL if theme == "rna" else CORAL)
    return _save(fig, plt)
