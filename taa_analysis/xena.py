# -*- coding: utf-8 -*-
"""UCSC Xena scheme queries. Failures stay failures."""

from __future__ import annotations

import json
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


class RetrievalFailed(RuntimeError):
    def __init__(self, source: str, detail: str):
        self.source = source
        self.missing_reason = "retrieval_failed"
        super().__init__(f"{source}: {detail}")


def xena_post(hub: str, query: str, timeout: int = 180) -> object:
    url = hub.rstrip("/") + "/data/"
    request = Request(url, query.encode("utf-8"), {"Content-Type": "text/plain"})
    try:
        with urlopen(request, timeout=timeout) as response:
            raw = response.read().decode("utf-8")
    except (HTTPError, URLError, TimeoutError, OSError) as exc:
        raise RetrievalFailed(hub, str(exc)) from exc
    if raw.lstrip().startswith("<"):
        raise RetrievalFailed(hub, "Xena returned an error page instead of data")
    try:
        return json.loads(raw)
    except json.JSONDecodeError as exc:
        raise RetrievalFailed(hub, "Xena response was not JSON") from exc


def _quote(value: str) -> str:
    return json.dumps(value)


def xena_columns(hub: str, dataset: str, columns: list[str]) -> dict:
    selected = " ".join(_quote(column) for column in columns)
    query = f"(xena-query {{:select [{selected}] :from [{_quote(dataset)}]}})"
    payload = xena_post(hub, query)
    if not isinstance(payload, dict):
        raise RetrievalFailed(dataset, "unexpected Xena column payload")
    return payload


def gene_probe(hub: str, probemap: str, gene_symbol: str) -> str:
    query = (
        "(xena-query {:select [\"name\" \"genes\"] :from ["
        + _quote(probemap)
        + "] :where [:in :any \"genes\" ["
        + _quote(gene_symbol)
        + "]]})"
    )
    payload = xena_post(hub, query)
    names = payload.get("name") if isinstance(payload, dict) else None
    if not names:
        raise RetrievalFailed(probemap, f"no probe for {gene_symbol}")
    return str(names[0])
