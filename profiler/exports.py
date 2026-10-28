# -*- coding: utf-8 -*-
"""Export helpers with provenance metadata preserved."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from typing import Any

import pandas as pd


def provenance_block(meta: dict[str, Any]) -> dict[str, Any]:
    base = {
        "exported_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "tool": "TAA Immune Profiler",
        "disclaimer": (
            "Research exploration of public datasets only. "
            "Not a clinical decision score, antibody affinity, or treatment-response predictor."
        ),
    }
    base.update(meta or {})
    return base


def df_to_csv_bytes(df: pd.DataFrame, meta: dict[str, Any] | None = None) -> bytes:
    # Prefix provenance as comment lines for CSV readability
    lines = []
    prov = provenance_block(meta or {})
    for k, v in prov.items():
        lines.append(f"# {k}: {v}")
    body = df.to_csv(index=False)
    return ("\n".join(lines) + "\n" + body).encode("utf-8-sig")


def to_json_bytes(payload: dict[str, Any], meta: dict[str, Any] | None = None) -> bytes:
    out = {"provenance": provenance_block(meta or {}), "data": payload}
    return json.dumps(out, indent=2, default=str).encode("utf-8")


def fig_to_svg_bytes(fig) -> bytes:
    return fig.to_image(format="svg")


def fig_to_png_bytes(fig, scale: int = 2) -> bytes:
    return fig.to_image(format="png", scale=scale)
