# -*- coding: utf-8 -*-
"""Sample joins that refuse row expansion."""

from __future__ import annotations

import pandas as pd


class JoinError(ValueError):
    """Join would duplicate samples or apply an annotation to the wrong row."""


def safe_left_join(left: pd.DataFrame, right: pd.DataFrame, on: str) -> pd.DataFrame:
    if left.empty:
        return left.copy()
    if on not in left.columns or on not in right.columns:
        raise JoinError(f"join key {on} is missing")
    duplicated = right[on].astype(str).duplicated().any()
    if duplicated:
        raise JoinError("right table has duplicate keys; refusing many-to-many expansion")
    merged = left.merge(right, on=on, how="left", validate="many_to_one")
    if len(merged) != len(left):
        raise JoinError(
            f"join changed row count from {len(left)} to {len(merged)}"
        )
    return merged
