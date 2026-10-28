# -*- coding: utf-8 -*-
"""Validation suite — real checks with fidelity vs biology columns separated."""

from __future__ import annotations

import json
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

from profiler.cell_lines import LITERATURE_DIRECTIONS, fetch_cell_line_expression, literature_direction_check
from profiler.cbioportal import CELL_LINE_ASSAYS
from profiler.gene_resolution import (
    FAILED_V13_ALIASES,
    PRIORITY_TAAS,
    SPECIAL_ANTIGENS,
    resolve_target,
)

RESULTS_DIR = Path(__file__).resolve().parent.parent / "validation_results"

# HER2 regression fixtures from 2026-10-20 (release-tagged, not eternal truth)
HER2_REGRESSION = {
    "symbol": "ERBB2",
    "entrez": 2064,
    "depmap_25q3_rna": {
        "ZR751_BREAST": 6.7664,
        "MDAMB231_BREAST": 4.5967,
    },
    "ccle_2019_rna": {
        "ZR751_BREAST": 53.7175,
        "MDAMB231_BREAST": 12.8311,
    },
    "ccle_2019_protein": {
        "ZR751_BREAST": 0.8146,
        "MDAMB231_BREAST": 1.2413,
    },
    "tolerance_abs": 1e-4,  # source float reporting precision; not for covering real diffs
}


def _now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _new_run() -> dict[str, Any]:
    return {
        "run_id": str(uuid.uuid4())[:8],
        "started_at": _now(),
        "finished_at": None,
        "app_version": "0.2.0",
        "checks": [],
    }


def _add(
    run: dict[str, Any],
    *,
    name: str,
    status: str,
    data_fidelity: str,
    biology_concordance: str,
    detail: dict[str, Any],
) -> None:
    run["checks"].append(
        {
            "name": name,
            "status": status,  # pending|running|pass|mismatch|data_insufficient|fetch_failed|held
            "data_fidelity": data_fidelity,  # pass|fail|n/a|not_run
            "biology_concordance": biology_concordance,  # pass|fail|n/a|not_run|not_applicable
            "detail": detail,
            "recorded_at": _now(),
        }
    )


def check_alias_mappings(run: dict[str, Any] | None = None, progress: Callable | None = None) -> dict[str, Any]:
    run = run or _new_run()
    cases = []
    # All failed v13 aliases + core 17
    inputs = list(FAILED_V13_ALIASES) + [t["input"] for t in PRIORITY_TAAS] + ["HER2", "EGFR", "GD2", "PD-1", "4-1BB"]
    # unique preserve order
    seen = set()
    uniq = []
    for i in inputs:
        k = i.upper()
        if k not in seen:
            seen.add(k)
            uniq.append(i)

    for query in uniq:
        if progress:
            progress(f"Alias resolve: {query}")
        resolved = resolve_target(query, confirm_provider=True)
        alias_expect = {
            "MESOTHELIN": "MSLN",
            "TROP-2": "TACSTD2",
            "TROP2": "TACSTD2",
            "NECTIN-4": "NECTIN4",
            "PVRL4": "NECTIN4",
            "CD20": "MS4A1",
            "BCMA": "TNFRSF17",
            "CD269": "TNFRSF17",
            "CEA": "CEACAM5",
            "CD66E": "CEACAM5",
            "CLDN18.2": "CLDN18",
            "CLAUDIN18.2": "CLDN18",
            "HER2": "ERBB2",
            "EGFR": "EGFR",
            "PSMA": "FOLH1",
            "STEAP1": "STEAP1",
            "MSLN": "MSLN",
            "B7-H3": "CD276",
            "CD19": "CD19",
            "GPRC5D": "GPRC5D",
            "DLL3": "DLL3",
            "FAP": "FAP",
            "PD-1": "PDCD1",
            "4-1BB": "TNFRSF9",
        }
        expected_symbol = alias_expect.get(query.upper())
        if expected_symbol is None:
            for t in PRIORITY_TAAS:
                if t["input"].upper() == query.upper():
                    expected_symbol = t["symbol"]
                    break

        if query.upper() == "GD2":
            ok = resolved.unresolved and resolved.special_kind == "glycolipid"
            _add(
                run,
                name=f"alias:{query}",
                status="pass" if ok else "mismatch",
                data_fidelity="pass" if ok else "fail",
                biology_concordance="not_applicable",
                detail=resolved.to_dict(),
            )
            cases.append(ok)
            continue

        ok = (
            not resolved.unresolved
            and resolved.official_symbol == expected_symbol
            and resolved.entrez_gene_id is not None
        )
        # CLDN18.2 must carry isoform warning
        if "18.2" in query.upper().replace(" ", ""):
            ok = ok and resolved.special_kind == "isoform_limited" and any("CLDN18.2" in w or "isoform" in w.lower() for w in resolved.warnings)

        _add(
            run,
            name=f"alias:{query}",
            status="pass" if ok else "mismatch",
            data_fidelity="pass" if ok else "fail",
            biology_concordance="not_applicable",
            detail={
                "query": query,
                "expected_symbol": expected_symbol,
                "resolved": resolved.to_dict(),
            },
        )
        cases.append(ok)

    run["alias_summary"] = {
        "n": len(cases),
        "n_pass": sum(1 for c in cases if c),
        "n_fail": sum(1 for c in cases if not c),
    }
    return run


def check_her2_regression(run: dict[str, Any] | None = None, progress: Callable | None = None) -> dict[str, Any]:
    run = run or _new_run()
    resolved = resolve_target("HER2")
    if resolved.entrez_gene_id != HER2_REGRESSION["entrez"]:
        _add(
            run,
            name="her2_entrez",
            status="mismatch",
            data_fidelity="fail",
            biology_concordance="n/a",
            detail=resolved.to_dict(),
        )
        return run

    for assay_key, expected in [
        ("depmap_25q3_rna", HER2_REGRESSION["depmap_25q3_rna"]),
        ("ccle_2019_rna", HER2_REGRESSION["ccle_2019_rna"]),
        ("ccle_2019_protein", HER2_REGRESSION["ccle_2019_protein"]),
    ]:
        if progress:
            progress(f"HER2 regression: {assay_key}")
        bundle = fetch_cell_line_expression("ERBB2", resolved.entrez_gene_id, assay_key=assay_key)
        if bundle["status"] == "fetch_failed":
            _add(
                run,
                name=f"her2_regression:{assay_key}",
                status="fetch_failed",
                data_fidelity="fail",
                biology_concordance="not_run",
                detail={"error": bundle.get("error"), "retrieved_at": bundle.get("retrieved_at")},
            )
            continue

        df = bundle["dataframe"]
        by_id = {str(r["sampleId"]): float(r["value"]) for _, r in df.iterrows()}
        diffs = {}
        ok = True
        for sid, exp_val in expected.items():
            got = by_id.get(sid)
            if got is None:
                ok = False
                diffs[sid] = {"expected": exp_val, "got": None, "abs_diff": None}
            else:
                diff = abs(got - exp_val)
                diffs[sid] = {"expected": exp_val, "got": got, "abs_diff": diff}
                if diff > HER2_REGRESSION["tolerance_abs"]:
                    ok = False

        # Protein direction vs surface literature (biology column only)
        bio = "n/a"
        if assay_key == "ccle_2019_protein":
            z = by_id.get("ZR751_BREAST")
            m = by_id.get("MDAMB231_BREAST")
            if z is not None and m is not None:
                # Literature surface: ZR > MDA; protein relative: historically MDA > ZR
                bio = "mismatch_vs_surface_literature" if not (z > m) else "match_surface_literature"
                detail_note = (
                    "Whole-cell protein relative values need not match surface receptors/cell. "
                    "Do not swap values to force surface order."
                )
            else:
                detail_note = "Missing samples"
                bio = "not_run"
        else:
            detail_note = "RNA fidelity check against 2026-10-20 regression snapshot"
            if assay_key == "depmap_25q3_rna":
                z = by_id.get("ZR751_BREAST")
                m = by_id.get("MDAMB231_BREAST")
                if z is not None and m is not None:
                    bio = "pass_direction" if z > m else "mismatch_direction"

        _add(
            run,
            name=f"her2_regression:{assay_key}",
            status="pass" if ok else "mismatch",
            data_fidelity="pass" if ok else "fail",
            biology_concordance=bio,
            detail={
                "assay": assay_key,
                "diffs": diffs,
                "tolerance_abs": HER2_REGRESSION["tolerance_abs"],
                "n_measured": bundle["n_measured"],
                "retrieved_at": bundle["retrieved_at"],
                "checksum": bundle.get("metadata", {}).get("checksum"),
                "note": detail_note,
                "source_release_note": "Compare to 2026-10-20 snapshot; source may update — record reason if changed.",
            },
        )
    return run


def check_17_targets_cell_lines(
    run: dict[str, Any] | None = None,
    progress: Callable | None = None,
    assays: list[str] | None = None,
) -> dict[str, Any]:
    run = run or _new_run()
    assays = assays or list(CELL_LINE_ASSAYS.keys())

    for taa in PRIORITY_TAAS:
        query = taa["input"]
        if progress:
            progress(f"Target profiles: {query}")
        resolved = resolve_target(query)
        if query.upper() == "GD2":
            _add(
                run,
                name=f"profiles:{query}",
                status="held",
                data_fidelity="n/a",
                biology_concordance="not_applicable",
                detail={"reason": "glycolipid — no gene-level expression quantification"},
            )
            continue
        if resolved.unresolved or not resolved.entrez_gene_id:
            _add(
                run,
                name=f"profiles:{query}",
                status="mismatch",
                data_fidelity="fail",
                biology_concordance="not_run",
                detail=resolved.to_dict(),
            )
            continue

        for assay_key in assays:
            bundle = fetch_cell_line_expression(
                resolved.official_symbol or taa["symbol"],
                resolved.entrez_gene_id,
                assay_key=assay_key,
            )
            if bundle["status"] == "fetch_failed":
                _add(
                    run,
                    name=f"profile:{query}:{assay_key}",
                    status="fetch_failed",
                    data_fidelity="fail",
                    biology_concordance="not_run",
                    detail={"error": bundle.get("error"), "retrieved_at": bundle.get("retrieved_at")},
                )
                continue

            # Fidelity: response gene/profile match
            meta = bundle.get("metadata") or {}
            fidelity_ok = meta.get("mismatch_n", 0) == 0 and bundle["status"] in ("ok", "no_measured_values")
            status = "pass" if fidelity_ok else "mismatch"
            if bundle["n_measured"] == 0:
                status = "data_insufficient"
            bio_results = literature_direction_check(bundle["dataframe"], resolved.official_symbol or "")
            bio_status = "n/a"
            if bio_results:
                r0 = bio_results[0]
                bio_status = r0.get("status") or "n/a"

            _add(
                run,
                name=f"profile:{query}:{assay_key}",
                status=status,
                data_fidelity="pass" if fidelity_ok else "fail",
                biology_concordance=bio_status,
                detail={
                    "symbol": resolved.official_symbol,
                    "entrez": resolved.entrez_gene_id,
                    "n_measured": bundle["n_measured"],
                    "retrieved_at": bundle["retrieved_at"],
                    "unit": bundle["assay_meta"]["unit"],
                    "biology": bio_results,
                    "category": taa["category"],
                    "therapy_types": taa["therapy_types"],
                    "note": (
                        "not_measured ≠ antigen absent. "
                        "Biology column compares surface-literature direction to this assay only."
                    ),
                },
            )
    return run


def run_full_validation(progress: Callable | None = None) -> dict[str, Any]:
    run = _new_run()
    if progress:
        progress("Starting alias checks…")
    check_alias_mappings(run, progress=progress)
    if progress:
        progress("HER2 regression…")
    check_her2_regression(run, progress=progress)
    if progress:
        progress("17-target cell-line profiles…")
    check_17_targets_cell_lines(run, progress=progress)
    run["finished_at"] = _now()
    run["summary"] = summarize(run)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    path = RESULTS_DIR / f"validation_{run['run_id']}_{run['started_at'].replace(':', '')}.json"
    path.write_text(json.dumps(run, indent=2), encoding="utf-8")
    run["result_path"] = str(path)
    return run


def summarize(run: dict[str, Any]) -> dict[str, int]:
    counts = {
        "pass": 0,
        "mismatch": 0,
        "fetch_failed": 0,
        "data_insufficient": 0,
        "held": 0,
        "pending": 0,
        "total": 0,
        "fidelity_pass": 0,
        "fidelity_fail": 0,
        "biology_pass_direction": 0,
        "biology_mismatch_direction": 0,
    }
    for c in run.get("checks") or []:
        counts["total"] += 1
        st = c.get("status") or "pending"
        if st in counts:
            counts[st] += 1
        if c.get("data_fidelity") == "pass":
            counts["fidelity_pass"] += 1
        elif c.get("data_fidelity") == "fail":
            counts["fidelity_fail"] += 1
        bio = c.get("biology_concordance") or ""
        if bio == "pass_direction":
            counts["biology_pass_direction"] += 1
        elif bio in ("mismatch_direction", "mismatch_vs_surface_literature"):
            counts["biology_mismatch_direction"] += 1
    return counts
