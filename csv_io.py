# -*- coding: utf-8 -*-
"""CSV/TSV readers that tolerate mixed encodings."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def read_csv_safe(path: Path, **kwargs):
    """Read CSV/TSV with UTF-8 first, then common Windows encodings."""
    chunksize = kwargs.get("chunksize")
    last_error: Exception | None = None
    for encoding in ("utf-8-sig", "utf-8", "cp949", "latin-1"):
        try:
            return pd.read_csv(path, encoding=encoding, **kwargs)
        except UnicodeDecodeError as exc:
            last_error = exc
        except pd.errors.EmptyDataError:
            return pd.DataFrame() if not chunksize else iter(())
    if last_error:
        raise last_error
    return pd.read_csv(path, **kwargs)
