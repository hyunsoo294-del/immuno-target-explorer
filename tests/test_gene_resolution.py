# -*- coding: utf-8 -*-
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from profiler.gene_resolution import FAILED_V13_ALIASES, alias_lookup, resolve_target


def test_v13_failed_aliases_map():
    expect = {
        "Mesothelin": "MSLN",
        "TROP-2": "TACSTD2",
        "Nectin-4": "NECTIN4",
        "PVRL4": "NECTIN4",
        "CD20": "MS4A1",
        "BCMA": "TNFRSF17",
        "CD269": "TNFRSF17",
        "CEA": "CEACAM5",
        "CD66e": "CEACAM5",
        "CLDN18.2": "CLDN18",
        "Claudin18.2": "CLDN18",
    }
    for q, symbol in expect.items():
        got, special, warnings = alias_lookup(q)
        assert got == symbol, f"{q} -> {got}, expected {symbol}"
        if "18.2" in q:
            assert special == "isoform_limited"
            assert warnings


def test_gd2_not_mapped_to_gene():
    got, special, warnings = alias_lookup("GD2")
    assert got is None
    assert special == "glycolipid"
    assert warnings


def test_no_silent_fuzzy():
    got, special, warnings = alias_lookup("NOTAREALGENEXYZ")
    assert got is None or got == "NOTAREALGENEXYZ"  # pass-through only if symbol-shaped
    # Unrecognized free text should not map
    got2, _, warnings2 = alias_lookup("some random phrase")
    assert got2 is None
    assert warnings2


def test_resolve_her2_live():
    r = resolve_target("HER2")
    assert r.official_symbol == "ERBB2"
    assert r.entrez_gene_id == 2064
