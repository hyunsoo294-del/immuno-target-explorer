"""Download a DepMap release and write the parquet cache plus MANIFEST.json.

Do not hardcode a release in calling code. Pass --release.
`current` tries the portal file index, then the newest figshare mirror listed in
config/depmap_sources.yaml if the portal answers with a challenge page.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent))

from immunoscore.src.config_loader import DATA_DIR, load_configs  # noqa: E402

RAW = DATA_DIR / "raw"
PROCESSED = DATA_DIR / "processed"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def portal_table(api: str) -> pd.DataFrame | None:
    response = requests.get(api, timeout=60)
    text = response.text[:200].lstrip().lower()
    if "html" in response.headers.get("content-type", "") or text.startswith("<!doctype") or text.startswith("<html"):
        print("Portal file index returned a challenge page, not a file table.", file=sys.stderr)
        return None
    try:
        return pd.read_csv(pd.io.common.StringIO(response.text))
    except Exception as exc:  # noqa: BLE001
        print(f"Portal file index could not be parsed: {exc}", file=sys.stderr)
        return None


def resolve(release: str) -> tuple[str, dict, str]:
    sources = load_configs()["depmap_sources"]
    mirrors = sources.get("mirrors") or {}
    note = ""
    if release != "current":
        if release not in mirrors or not (mirrors[release].get("files") or {}):
            raise SystemExit(f"No public URL is configured for release {release}.")
        return release, mirrors[release], note
    table = portal_table(sources["portal_files_api"])
    if table is not None and {"release", "filename", "url"}.issubset(set(table.columns)):
        current = table.sort_values("release").iloc[-1]["release"]
        files = {}
        for name in sources["needed_files"]:
            hit = table[(table["release"] == current) & (table["filename"] == name)]
            if hit.empty:
                continue
            files[name] = str(hit.iloc[0]["url"])
        if files:
            return str(current), {"label": str(current), "files": files, "doi": ""}, note
    # Newest mirror key in the YAML, which is listed newest-first.
    for key, entry in mirrors.items():
        if entry.get("files"):
            note = (
                "Portal current release was not readable from this environment. "
                f"Using public mirror {entry['label']}."
            )
            print(note, file=sys.stderr)
            return key, entry, note
    raise SystemExit("No DepMap files could be resolved.")


def download(url: str, dest: Path) -> None:
    if dest.exists() and dest.stat().st_size > 0:
        print(f"exists {dest.name}")
        return
    print(f"download {dest.name}")
    with requests.get(url, stream=True, timeout=120) as response:
        response.raise_for_status()
        dest.parent.mkdir(parents=True, exist_ok=True)
        with dest.open("wb") as handle:
            for chunk in response.iter_content(1 << 20):
                if chunk:
                    handle.write(chunk)


def process(release_id: str, entry: dict, note: str) -> None:
    expr_path = RAW / "OmicsExpressionProteinCodingGenesTPMLogp1.csv"
    model_path = RAW / "Model.csv"
    print("reading expression")
    frame = pd.read_csv(expr_path, index_col=0)
    frame.index.name = "ModelID"
    frame.columns = [str(column).split(" (")[0].strip() for column in frame.columns]
    if frame.columns.duplicated().any():
        frame = frame.T.groupby(level=0).mean().T
    frame = frame.astype("float32")
    PROCESSED.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(PROCESSED / "expression_log2tpm.parquet")
    pd.read_csv(model_path).to_parquet(PROCESSED / "models.parquet", index=False)
    manifest = {
        "release": entry.get("label") or release_id,
        "release_id": release_id,
        "doi": entry.get("doi") or "",
        "downloaded_at": datetime.now(timezone.utc).strftime("%Y-%m-%d"),
        "note": note,
        "value": "log2(TPM+1)",
        "expression_shape": [int(frame.shape[0]), int(frame.shape[1])],
        "files": {
            path.name: {"sha256": sha256(path), "bytes": path.stat().st_size}
            for path in (model_path, expr_path)
            if path.exists()
        },
    }
    (RAW / "MANIFEST.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print("processed", frame.shape)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--release", default="current")
    parser.add_argument("--process-only", action="store_true")
    args = parser.parse_args()
    release_id, entry, note = resolve(args.release)
    RAW.mkdir(parents=True, exist_ok=True)
    if not args.process_only:
        for name, url in (entry.get("files") or {}).items():
            if name not in load_configs()["depmap_sources"]["needed_files"]:
                continue
            if not url:
                continue
            download(url, RAW / name)
    process(release_id, entry, note)


if __name__ == "__main__":
    main()
