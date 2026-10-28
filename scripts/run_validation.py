#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CLI validation runner for TAA Immune Profiler."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from profiler.validation import (  # noqa: E402
    check_17_targets_cell_lines,
    check_alias_mappings,
    check_her2_regression,
    run_full_validation,
    summarize,
)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run TAA Immune Profiler validation")
    parser.add_argument(
        "--suite",
        choices=["aliases", "her2", "profiles", "full"],
        default="full",
    )
    parser.add_argument("--out", type=str, default="")
    args = parser.parse_args()

    def progress(msg: str) -> None:
        print(f"[validation] {msg}", flush=True)

    if args.suite == "aliases":
        run = check_alias_mappings(progress=progress)
    elif args.suite == "her2":
        run = check_her2_regression(progress=progress)
    elif args.suite == "profiles":
        run = check_17_targets_cell_lines(progress=progress)
    else:
        run = run_full_validation(progress=progress)

    run["summary"] = summarize(run)
    text = json.dumps(run, indent=2, default=str)
    out = Path(args.out) if args.out else Path(run.get("result_path") or "validation_results/last.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text, encoding="utf-8")
    print(json.dumps(run["summary"], indent=2))
    print(f"Wrote {out}")

    # Exit non-zero only on fidelity failures / mismatches in critical suites — still report biology separately
    fidelity_fail = run["summary"].get("fidelity_fail", 0)
    mismatch = run["summary"].get("mismatch", 0)
    fetch_failed = run["summary"].get("fetch_failed", 0)
    if fidelity_fail or mismatch or fetch_failed:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
