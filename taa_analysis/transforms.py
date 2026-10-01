# -*- coding: utf-8 -*-
"""Expression unit conversions. Other units are never relabeled as TPM."""

from __future__ import annotations

import math

TOLERANCE = 1e-6

SOURCE_LOG2_TPM_0001 = "log2(tpm+0.001)"
DISPLAY_LOG2_TPM_1 = "log2(TPM+1)"
FORMULA_0001_TO_1 = "TPM = 2^x - 0.001; display = log2(TPM + 1)"


class TransformError(ValueError):
    """Stored value cannot be converted without inventing a scale."""


def _as_float(value) -> float:
    if value is None or value == "":
        raise TransformError("missing expression value")
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise TransformError("non-numeric expression value") from exc
    if math.isnan(number) or math.isinf(number):
        raise TransformError("non-finite expression value")
    return number


def recover_tpm_from_log2_tpm_plus_0001(stored: float, atol: float = TOLERANCE) -> float:
    """Invert log2(TPM+0.001). Negative TPM beyond tolerance is an error."""
    stored = _as_float(stored)
    tpm = math.pow(2.0, stored) - 0.001
    if tpm < -atol:
        raise TransformError(
            f"back-transformed TPM {tpm} is negative beyond tolerance {atol}"
        )
    if tpm < 0:
        return 0.0
    return float(tpm)


def log2_tpm_plus_1(tpm: float) -> float:
    tpm = _as_float(tpm)
    if tpm < 0:
        raise TransformError("TPM must be non-negative before log2(TPM+1)")
    return float(math.log2(tpm + 1.0))


def convert_stored_expression(stored: float, source_unit: str) -> dict:
    """Convert one stored value. Refuses double-logging and unit relabeling."""
    unit = " ".join(str(source_unit or "").strip().lower().split())
    stored_value = _as_float(stored)
    if unit == SOURCE_LOG2_TPM_0001:
        tpm = recover_tpm_from_log2_tpm_plus_0001(stored_value)
        return {
            "source_value": stored_value,
            "source_unit": SOURCE_LOG2_TPM_0001,
            "tpm": tpm,
            "display_value": log2_tpm_plus_1(tpm),
            "unit": DISPLAY_LOG2_TPM_1,
            "formula": FORMULA_0001_TO_1,
            "double_log_applied": False,
        }
    if unit in {DISPLAY_LOG2_TPM_1.lower(), "log2(tpm+1)", "log2(tpm + 1)"}:
        return {
            "source_value": stored_value,
            "source_unit": DISPLAY_LOG2_TPM_1,
            "tpm": None,
            "display_value": stored_value,
            "unit": DISPLAY_LOG2_TPM_1,
            "formula": "identity; source is already log2(TPM+1), so it is not logged again",
            "double_log_applied": False,
        }
    raise TransformError(
        f"refusing to relabel {source_unit!r} as TPM. "
        "nTPM, FPKM, RSEM expected count, and microarray values stay on their own scale."
    )
