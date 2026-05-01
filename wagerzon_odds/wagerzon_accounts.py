"""Registry of configured Wagerzon accounts.

Discovers accounts from environment variables matching the pattern
`WAGERZON{SUFFIX}_USERNAME` / `WAGERZON{SUFFIX}_PASSWORD`, where SUFFIX is
either empty (primary account) or a single uppercase ASCII letter (e.g. J, C).

Usage:
    from wagerzon_accounts import list_accounts, get_account
    for acct in list_accounts():
        print(acct.label, acct.username)
    j = get_account("WagerzonJ")
"""
from __future__ import annotations

import logging
import os
import re
from dataclasses import dataclass
from pathlib import Path

try:
    from dotenv import load_dotenv
    _ENV_PATH = Path(__file__).resolve().parent.parent / "bet_logger" / ".env"
    if _ENV_PATH.exists():
        load_dotenv(_ENV_PATH)
except ImportError:
    pass

logger = logging.getLogger(__name__)

_SUFFIX_RE = re.compile(r"^WAGERZON([A-Z]?)_USERNAME$")


class AccountNotFoundError(KeyError):
    """Raised by get_account when the label is not in the current registry."""


@dataclass(frozen=True)
class WagerzonAccount:
    label: str
    suffix: str
    username: str
    password: str


def list_accounts() -> list[WagerzonAccount]:
    """Discover WZ accounts from environment.

    Returns the primary account (suffix '') first if present, then alphabetical
    by suffix. Pairs missing one half (USERNAME or PASSWORD) are skipped with a
    warning.
    """
    found: dict[str, WagerzonAccount] = {}
    for key, val in os.environ.items():
        m = _SUFFIX_RE.match(key)
        if not m:
            continue
        suffix = m.group(1)
        pw_key = f"WAGERZON{suffix}_PASSWORD"
        password = os.environ.get(pw_key)
        if not val or not password:
            logger.warning(
                "Wagerzon account WAGERZON%s skipped: %s set but %s missing/empty",
                suffix or "(primary)", key,
                pw_key if not password else key,
            )
            continue
        label = "Wagerzon" + suffix
        found[suffix] = WagerzonAccount(
            label=label, suffix=suffix, username=val, password=password,
        )

    ordered: list[WagerzonAccount] = []
    if "" in found:
        ordered.append(found.pop(""))
    for suffix in sorted(found.keys()):
        ordered.append(found[suffix])
    return ordered


def get_account(label: str) -> WagerzonAccount:
    """Look up an account by label (e.g. 'WagerzonJ'). Raises AccountNotFoundError."""
    for acct in list_accounts():
        if acct.label == label:
            return acct
    raise AccountNotFoundError(label)
