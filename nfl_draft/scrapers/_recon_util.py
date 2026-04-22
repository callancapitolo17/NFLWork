"""Shared helpers for NFL Draft reconnaissance scripts.

Each recon script (recon_dk.py / recon_fd.py / recon_bm.py / recon_wz.py) captures
the raw API response that powers NFL Draft futures markets on its book and writes
it to nfl_draft/tests/fixtures/<book>/draft_markets.json. The parser subagents that
run afterwards work OFFLINE from these fixtures — they never hit the live book.

Keep this file tiny: just path resolution + a pretty-JSON writer. Anything
scraping-specific (Playwright, curl_cffi, auth) lives in each book's own script.
"""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path
from typing import Any

# Resolve repo root by walking up from this file. nfl_draft/ lives at repo root,
# so scrapers/_recon_util.py -> scrapers -> nfl_draft -> repo root.
_THIS_DIR = Path(__file__).resolve().parent
NFL_DRAFT_ROOT = _THIS_DIR.parent            # .../NFLWork/nfl_draft
FIXTURES_ROOT = NFL_DRAFT_ROOT / "tests" / "fixtures"


def _main_repo_root() -> Path:
    """Return the working tree root of the MAIN git repo, even from a worktree.

    In a worktree, NFL_DRAFT_ROOT.parent resolves to the worktree's checkout
    directory, not the main repo. Shared artifacts (.env files, session cookies)
    live in the main repo's working tree. This helper resolves to the main
    repo so callers can look up those shared paths regardless of where the
    script is invoked from.
    """
    try:
        # git-common-dir points at the shared .git directory (main repo's .git).
        # Its parent is the main repo's working tree root.
        common = subprocess.check_output(
            ["git", "-C", str(NFL_DRAFT_ROOT), "rev-parse", "--git-common-dir"],
            text=True, stderr=subprocess.DEVNULL,
        ).strip()
        if common:
            git_dir = Path(common)
            if not git_dir.is_absolute():
                git_dir = (NFL_DRAFT_ROOT / git_dir).resolve()
            return git_dir.parent
    except (subprocess.CalledProcessError, FileNotFoundError, OSError):
        pass
    # Fallback: if git isn't available, assume NFL_DRAFT_ROOT is inside the
    # main repo (the common non-worktree case).
    return NFL_DRAFT_ROOT.parent


def fixture_dir(book: str) -> Path:
    """Return (and create if needed) the fixture directory for a given book."""
    d = FIXTURES_ROOT / book
    d.mkdir(parents=True, exist_ok=True)
    return d


def fixture_path(book: str) -> Path:
    """Canonical fixture file path: nfl_draft/tests/fixtures/<book>/draft_markets.json."""
    return fixture_dir(book) / "draft_markets.json"


def save_fixture(book: str, data: Any, *, meta: dict | None = None) -> Path:
    """Write captured raw API data to the canonical fixture path.

    The payload is wrapped in a small envelope so a future parser can see WHERE
    the data came from. The envelope looks like:

        {
            "captured_from": "https://...",       # source URL if known
            "captured_at":   "2026-04-17T12:34:56Z",
            "book":          "draftkings",
            "data":          <raw book response, unchanged>,
        }

    Passing a dict for `meta` lets the caller stash extra diagnostics
    (e.g. league_id discovered, request body used). Returns the Path written.
    """
    from datetime import datetime, timezone

    envelope = {
        "book": book,
        "captured_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "meta": meta or {},
        "data": data,
    }
    path = fixture_path(book)
    with path.open("w", encoding="utf-8") as f:
        json.dump(envelope, f, indent=2, default=str)
    return path


def print_diagnostics(book: str, path: Path, url: str | None, data: Any) -> None:
    """Print a one-screen summary of what was captured so the user can sanity-check.

    Shows: file path, size on disk, source URL, and the top-level JSON keys (or
    list length). This is the minimum needed to confirm the fixture looks right
    before handing it off to a parser subagent.
    """
    size = path.stat().st_size if path.exists() else 0
    print(f"\n  {'=' * 58}")
    print(f"  [{book.upper()}] fixture saved")
    print(f"  {'=' * 58}")
    print(f"  path:    {path}")
    print(f"  size:    {size:,} bytes")
    if url:
        print(f"  source:  {url}")
    if isinstance(data, dict):
        keys = list(data.keys())[:15]
        print(f"  top-level keys ({len(data)}): {keys}")
    elif isinstance(data, list):
        print(f"  list of {len(data)} items")
        if data and isinstance(data[0], dict):
            print(f"  first item keys: {list(data[0].keys())[:10]}")
    print()


def ensure_fixture_dirs() -> None:
    """Create per-book fixture dirs + .gitkeep so the tree exists in git."""
    for book in ("betonline", "bookmaker", "draftkings", "fanduel", "hoop88", "kalshi", "wagerzon"):
        d = fixture_dir(book)
        gk = d / ".gitkeep"
        if not gk.exists():
            gk.touch()


def load_env(env_path: Path | None = None) -> None:
    """Load env vars from bet_logger/.env if python-dotenv is installed.

    The scripts fall back to os.getenv() regardless, so missing dotenv is not
    fatal — but if the lib is there (it is in the scraper venvs) we load the
    file so a user doesn't have to `export` credentials manually.

    Works in both main repo and worktree: resolves bet_logger/.env relative
    to the main repo's working tree, not NFL_DRAFT_ROOT.parent (which
    differs in a worktree).
    """
    if env_path is None:
        env_path = _main_repo_root() / "bet_logger" / ".env"
    if not env_path.exists():
        return
    try:
        from dotenv import load_dotenv
        load_dotenv(env_path, override=False)
    except ImportError:
        # Minimal fallback: naive KEY=VALUE parsing. Good enough for our 4 creds.
        for line in env_path.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            k, v = line.split("=", 1)
            os.environ.setdefault(k.strip(), v.strip().strip('"').strip("'"))
