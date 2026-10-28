# -*- coding: utf-8 -*-
"""Shared HTTP client with Retry-After / 429 backoff and concurrency limits."""

from __future__ import annotations

import hashlib
import json
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import requests
from requests.adapters import HTTPAdapter

CACHE_DIR = Path(__file__).resolve().parent.parent / "cache" / "http"
DEFAULT_TIMEOUT = (12, 90)
HEADERS = {
    "User-Agent": "TAAImmuneProfiler/0.2 (research; private)",
    "Accept": "application/json, text/plain, */*",
    "Connection": "close",
}

_lock = threading.Lock()
_inflight: dict[str, threading.Event] = {}
_results: dict[str, Any] = {}
_sem = threading.Semaphore(3)


@dataclass
class HttpError:
    status: int | None
    kind: str
    message: str
    url: str
    retry_after: float | None = None


class RateLimitedError(RuntimeError):
    def __init__(self, message: str, retry_after: float | None = None):
        super().__init__(message)
        self.retry_after = retry_after


class AuthError(RuntimeError):
    pass


class NotFoundError(RuntimeError):
    pass


def _cache_key(method: str, url: str, params: dict | None, body: Any) -> str:
    payload = json.dumps({"m": method, "u": url, "p": params or {}, "b": body}, sort_keys=True, default=str)
    return hashlib.sha256(payload.encode()).hexdigest()


def session() -> requests.Session:
    s = requests.Session()
    s.headers.update(HEADERS)
    adapter = HTTPAdapter(pool_connections=8, pool_maxsize=8, max_retries=0)
    s.mount("https://", adapter)
    s.mount("http://", adapter)
    return s


def request_json(
    method: str,
    url: str,
    *,
    params: dict | None = None,
    json_body: Any = None,
    timeout: tuple | float = DEFAULT_TIMEOUT,
    max_attempts: int = 6,
    cache_ttl_s: int = 3600,
    use_cache: bool = True,
    sess: requests.Session | None = None,
) -> Any:
    """GET/POST JSON with 429 Retry-After, bounded backoff, and TTL disk cache."""
    key = _cache_key(method, url, params, json_body)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_path = CACHE_DIR / f"{key}.json"

    if use_cache and cache_path.exists():
        age = time.time() - cache_path.stat().st_mtime
        if age <= cache_ttl_s:
            try:
                return json.loads(cache_path.read_text(encoding="utf-8"))
            except Exception:
                pass

    with _lock:
        if key in _inflight:
            event = _inflight[key]
            waiter = True
        else:
            event = threading.Event()
            _inflight[key] = event
            waiter = False

    if waiter:
        event.wait(timeout=180)
        if key in _results:
            return _results[key]
        if use_cache and cache_path.exists():
            return json.loads(cache_path.read_text(encoding="utf-8"))

    owned = sess or session()
    last_err: Exception | None = None
    try:
        for attempt in range(max_attempts):
            with _sem:
                try:
                    resp = owned.request(
                        method,
                        url,
                        params=params,
                        json=json_body,
                        timeout=timeout,
                    )
                except (requests.ConnectionError, requests.Timeout, requests.exceptions.ChunkedEncodingError) as exc:
                    last_err = exc
                    time.sleep(min(32.0, 0.8 * (2**attempt)))
                    continue

            if resp.status_code == 429:
                ra = resp.headers.get("Retry-After")
                wait = float(ra) if ra and str(ra).replace(".", "", 1).isdigit() else min(32.0, 1.5 * (2**attempt))
                last_err = RateLimitedError(f"HTTP 429 from {url}", retry_after=wait)
                time.sleep(wait)
                continue
            if resp.status_code in (401, 403):
                raise AuthError(f"HTTP {resp.status_code} auth failure for {url}")
            if resp.status_code == 404:
                raise NotFoundError(f"HTTP 404 not found: {url}")
            if resp.status_code >= 500:
                last_err = RuntimeError(f"HTTP {resp.status_code} server error for {url}")
                time.sleep(min(32.0, 0.8 * (2**attempt)))
                continue
            if resp.status_code >= 400:
                raise RuntimeError(f"HTTP {resp.status_code} for {url}: {resp.text[:300]}")

            try:
                data = resp.json()
            except Exception as exc:
                raise RuntimeError(f"Non-JSON response from {url}: {resp.text[:200]}") from exc

            if use_cache:
                try:
                    cache_path.write_text(json.dumps(data), encoding="utf-8")
                except Exception:
                    pass
            with _lock:
                _results[key] = data
            return data

        raise last_err or RuntimeError(f"Request failed for {url}")
    finally:
        with _lock:
            _inflight.pop(key, None)
            event.set()
            _results.pop(key, None)


def classify_error(exc: Exception) -> HttpError:
    if isinstance(exc, RateLimitedError):
        return HttpError(429, "rate_limited", str(exc), "", retry_after=exc.retry_after)
    if isinstance(exc, AuthError):
        return HttpError(401, "auth", str(exc), "")
    if isinstance(exc, NotFoundError):
        return HttpError(404, "not_found", str(exc), "")
    if isinstance(exc, requests.Timeout):
        return HttpError(None, "timeout", str(exc), "")
    return HttpError(None, "error", str(exc), "")
