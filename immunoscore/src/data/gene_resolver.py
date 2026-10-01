"""Resolve a TAA query to one DepMap gene symbol.

Order, first hit wins: exact column, curated alias, HGNC alias or previous
symbol. Comparisons are whole tokens. A prefix, substring, or regular
expression is not a match. Near-miss names are returned as unused suggestions.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache

import yaml

from immunoscore.src.config_loader import DATA_DIR, load_configs
from immunoscore.src.data.depmap import load_expression

HGNC_PATH = DATA_DIR / "reference" / "hgnc_complete_set.txt"
ALIAS_PATH = DATA_DIR / "reference" / "taa_aliases.yaml"
_SUGGESTION_LIMIT = 2
_SUGGESTION_COUNT = 5
_GLYCAN = {"glycan", "glycolipid"}


@dataclass
class Resolution:
    query: str
    symbol: str | None = None
    route: str = "none"
    antigen_class: str = ""
    warning: str = ""
    biosynthesis_genes: list[str] = field(default_factory=list)
    suggestions: list[str] = field(default_factory=list)
    refused: bool = False
    flags: list[str] = field(default_factory=list)
    confidence_override: str | None = None
    note: str = ""

    def failure_message(self) -> str:
        text = f"{self.query} 은(는) DepMap 발현 매트릭스에서 확정하지 못했습니다. 제안은 사용하지 않았습니다."
        if self.note:
            text = f"{text} {self.note}"
        if self.suggestions:
            text = f"{text} 가까운 이름: {', '.join(self.suggestions)}"
        return text


def resolve(
    query: str,
    *,
    symbols: tuple[str, ...] | list[str] | None = None,
    aliases: dict | None = None,
    hgnc_index: dict[str, tuple[str, ...]] | None = None,
) -> Resolution:
    text = (query or "").strip()
    if not text:
        raise ValueError("empty gene query")
    token = text.upper()
    chosen = tuple(symbols) if symbols is not None else _depmap_symbols()
    upper_map = _unique_upper_map(chosen)
    table = aliases if aliases is not None else _alias_table()
    hgnc = hgnc_index if hgnc_index is not None else _hgnc_index(chosen)

    exact = upper_map.get(token)
    if exact == "":
        return _unused(text, [name for name in chosen if name.upper() == token], "같은 이름의 유전자가 둘 이상입니다.")
    if exact:
        return Resolution(query=text, symbol=exact, route="exact", antigen_class="protein")

    curated = table.get(token)
    if curated is not None:
        return _from_curated(text, curated, upper_map)

    hits = hgnc.get(token, ())
    if len(hits) == 1:
        return Resolution(query=text, symbol=hits[0], route="hgnc", antigen_class="protein")
    if len(hits) > 1:
        return _unused(text, list(hits), "여러 유전자가 이 별칭을 씁니다.")
    return _unused(text, _suggestions(token, chosen, table))


def _from_curated(query: str, entry: dict, upper_map: dict[str, str]) -> Resolution:
    antigen_class = str(entry.get("antigen_class") or "protein")
    warning = " ".join(str(entry.get("warning") or "").split())
    genes = [str(gene) for gene in (entry.get("biosynthesis_genes") or [])]
    if antigen_class in _GLYCAN:
        return Resolution(
            query=query,
            route="curated",
            antigen_class=antigen_class,
            warning=warning,
            biosynthesis_genes=genes,
            refused=True,
            note=str(entry.get("note") or ""),
        )
    target = str(entry.get("symbol") or "")
    found = upper_map.get(target.upper()) if target else None
    if not found:
        return _unused(query, _suggestions(query.upper(), tuple(upper_map.values()), {}), f"{target or query} 가 발현 매트릭스에 없습니다.")
    flags: list[str] = []
    override = None
    if antigen_class == "isoform_variant":
        flags = ["ISOFORM_UNRESOLVED"]
        override = "LOW"
    return Resolution(
        query=query,
        symbol=found,
        route="curated",
        antigen_class=antigen_class,
        warning=warning,
        flags=flags,
        confidence_override=override,
        note=str(entry.get("note") or ""),
    )


def _unused(query: str, suggestions: list[str], note: str = "") -> Resolution:
    return Resolution(query=query, route="none", suggestions=suggestions, note=note)


def _unique_upper_map(symbols: tuple[str, ...] | list[str]) -> dict[str, str]:
    found: dict[str, str] = {}
    for symbol in symbols:
        key = str(symbol).upper()
        if key in found:
            found[key] = ""
        else:
            found[key] = str(symbol)
    return found


def _suggestions(token: str, symbols: tuple[str, ...] | list[str], aliases: dict) -> list[str]:
    ranked: list[tuple[int, str]] = []
    seen: set[str] = set()
    pool = [(str(symbol).upper(), str(symbol)) for symbol in symbols]
    for key, entry in aliases.items():
        label = str(key)
        target = entry.get("symbol")
        if target:
            label = f"{key} ({target})"
        pool.append((str(key).upper(), label))
    for upper, label in pool:
        if upper == token or label in seen:
            continue
        distance = _bounded_distance(token, upper, _SUGGESTION_LIMIT)
        if distance <= _SUGGESTION_LIMIT:
            ranked.append((distance, label))
            seen.add(label)
    ranked.sort(key=lambda item: (item[0], item[1]))
    return [label for _distance, label in ranked[:_SUGGESTION_COUNT]]


def _bounded_distance(left: str, right: str, limit: int) -> int:
    if abs(len(left) - len(right)) > limit:
        return limit + 1
    previous = list(range(len(right) + 1))
    for i, left_char in enumerate(left, start=1):
        current = [i]
        smallest = i
        for j, right_char in enumerate(right, start=1):
            insert = current[j - 1] + 1
            delete = previous[j] + 1
            replace = previous[j - 1] + (left_char != right_char)
            best = min(insert, delete, replace)
            current.append(best)
            if best < smallest:
                smallest = best
        if smallest > limit:
            return limit + 1
        previous = current
    return previous[-1]


@lru_cache(maxsize=1)
def _depmap_symbols() -> tuple[str, ...]:
    return tuple(str(column) for column in load_expression().columns)


@lru_cache(maxsize=1)
def _alias_table() -> dict[str, dict]:
    raw = yaml.safe_load(ALIAS_PATH.read_text(encoding="utf-8")) or {}
    table: dict[str, dict] = {}
    for key, entry in (raw.get("aliases") or {}).items():
        table[str(key).upper()] = dict(entry or {})
    for symbol, entry in (load_configs()["taa"].get("TAAs") or {}).items():
        table.setdefault(str(symbol).upper(), {"symbol": symbol, "antigen_class": "protein"})
        for alias in entry.get("aliases") or []:
            table.setdefault(str(alias).upper(), {"symbol": symbol, "antigen_class": "protein"})
    return table


@lru_cache(maxsize=1)
def _hgnc_index(symbols: tuple[str, ...]) -> dict[str, tuple[str, ...]]:
    if not HGNC_PATH.exists():
        return {}
    import csv

    upper_map = _unique_upper_map(symbols)
    found: dict[str, set[str]] = {}
    with HGNC_PATH.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("status") != "Approved":
                continue
            canonical = upper_map.get(str(row.get("symbol") or "").upper())
            if not canonical:
                continue
            blobs = (row.get("alias_symbol") or "", row.get("prev_symbol") or "")
            for blob in blobs:
                for part in blob.split("|"):
                    token = part.strip().upper()
                    if not token or token == canonical.upper():
                        continue
                    found.setdefault(token, set()).add(canonical)
    return {token: tuple(sorted(names)) for token, names in found.items()}
