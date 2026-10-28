# -*- coding: utf-8 -*-
"""Server-side private access gate (allowlist). Default: private."""

from __future__ import annotations

import hashlib
import hmac
import os
from typing import Any


def _env(name: str, default: str = "") -> str:
    return os.environ.get(name, default).strip()


def is_public_mode() -> bool:
    # Never default to public. Explicit opt-in only.
    return _env("TAA_PROFILER_PUBLIC", "0") in ("1", "true", "TRUE", "yes")


def allowlist_emails() -> set[str]:
    raw = _env("TAA_PROFILER_ALLOWLIST", "")
    return {e.strip().lower() for e in raw.split(",") if e.strip()}


def access_token() -> str:
    return _env("TAA_PROFILER_ACCESS_TOKEN", "")


def require_auth() -> bool:
    if is_public_mode():
        return False
    # Private by default even if token unset — then only local/dev bypass with explicit flag
    return True


def verify_token(provided: str | None) -> bool:
    expected = access_token()
    if not expected:
        # No token configured: allow only when explicitly enabled for local dev
        return _env("TAA_PROFILER_DEV_OPEN", "0") in ("1", "true", "TRUE")
    if not provided:
        return False
    return hmac.compare_digest(provided.strip(), expected)


def verify_email(email: str | None) -> bool:
    allowed = allowlist_emails()
    if not allowed:
        return False
    if not email:
        return False
    return email.strip().lower() in allowed


def auth_status(token: str | None = None, email: str | None = None) -> dict[str, Any]:
    if is_public_mode():
        return {"ok": True, "mode": "public", "warning": "PUBLIC MODE — requires explicit owner approval"}
    if verify_token(token) or verify_email(email):
        return {"ok": True, "mode": "private_allowlist"}
    if not access_token() and _env("TAA_PROFILER_DEV_OPEN", "0") in ("1", "true", "TRUE"):
        return {"ok": True, "mode": "dev_open"}
    return {
        "ok": False,
        "mode": "private",
        "message": "Private research tool. Provide an access token or allowlisted email.",
    }


def token_fingerprint() -> str:
    t = access_token()
    if not t:
        return "unset"
    return hashlib.sha256(t.encode()).hexdigest()[:10]
