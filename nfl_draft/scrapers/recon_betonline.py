#!/usr/bin/env python3
"""BetOnline NFL Draft markets recon.

Goal: capture the raw JSON that powers BetOnline's NFL Draft markets page
(https://www.betonline.ag/sportsbook/football/nfl-draft) and save it to
nfl_draft/tests/fixtures/betonline/draft_markets.json.

Discovery (2026-04-21 via live Playwright XHR capture + iterative probe):

  BetOnline's NFL Draft markets live under the `futures-and-props` sport
  tree — NOT `football`. Each bucket in the sidebar ("1st Round",
  "1st Round Props", "Draft Position", etc.) is its own contest page keyed
  by a URL slug. There are 11 buckets.

  Two-level API flow:
    1) GET https://www.betonline.ag/get-token
         Anonymous JWT used as `Authorization: Bearer`.
    2) POST https://api-offering.betonline.ag/api/offering/Sports/get-menu
         Returns the full sidebar, from which we extract the 11 NFL Draft
         sub-URLs. The slug after `nfl-draft/` is ContestType2.
    3) For each slug, POST
         https://api-offering.betonline.ag/api/offering/Sports/get-contests-by-contest-type2
         Body: {"ContestType":"nfl-draft","ContestType2":"<slug>","filterTime":0}
         Returns ContestOfferings.DateGroup[].DescriptionGroup[].TimeGroup[]
         .ContestExtended.ContestGroupLine[].Contestants[] — player names +
         American odds in .Line.MoneyLine.Line.

  Auth requires BetOnline-specific headers (gsetting, contests, gmt/utc-offset,
  actual-time, iso-time, utc-time) — same set as bet_logger/scraper_betonline.py.

Cloudflare is in front of api-offering.betonline.ag. We reuse the cookie jar
from `bet_logger/recon_betonline_cookies.json`, which carries `cf_clearance`
and is refreshed by `bet_logger/recon_betonline.py`.

If the script returns an HTTP 403 HTML error page, run:
    python bet_logger/recon_betonline.py
to refresh cookies, then re-run this.
"""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from nfl_draft.scrapers._recon_util import (
    _main_repo_root,
    print_diagnostics,
    save_fixture,
)

from curl_cffi import requests as cffi_requests

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

SITE_URL = "https://www.betonline.ag/"
TOKEN_URL = "https://www.betonline.ag/get-token"
MENU_URL = "https://api-offering.betonline.ag/api/offering/Sports/get-menu"
CONTESTS_URL = (
    "https://api-offering.betonline.ag/api/offering/Sports/get-contests-by-contest-type2"
)

CONTEST_TYPE = "nfl-draft"  # the slug under `/sportsbook/futures-and-props/`.

# Cookie jar maintained by bet_logger/recon_betonline.py. Carries cf_clearance.
# Resolved against main repo root so this runs from a worktree too.
COOKIE_PATH = _main_repo_root() / "bet_logger" / "recon_betonline_cookies.json"


# ---------------------------------------------------------------------------
# Session setup
# ---------------------------------------------------------------------------

def _load_cookies(session: cffi_requests.Session) -> int:
    """Load cookies from bet_logger's recon jar into the session.

    The jar is Playwright's list-of-dicts format: each entry has
    name/value/domain. We set every cookie with its stored domain so
    curl_cffi routes them correctly across .betonline.ag /
    api.betonline.ag / api-offering.betonline.ag / etc.
    """
    if not COOKIE_PATH.exists():
        return 0
    try:
        raw = json.loads(COOKIE_PATH.read_text())
    except Exception:
        return 0
    n = 0
    for c in raw:
        name = c.get("name")
        value = c.get("value")
        domain = c.get("domain")
        if not name or value is None:
            continue
        try:
            session.cookies.set(name, value, domain=domain)
            n += 1
        except Exception:
            continue
    return n


def _build_auth_headers(token: str) -> dict:
    """Build the full header set BetOnline's offering API expects.

    Mirrors bet_logger/scraper_betonline.py.build_api_headers. The time
    headers are re-computed per call; BetOnline appears not to enforce
    strict drift but sending them is cheap and matches the browser.
    """
    now = datetime.now(timezone.utc)
    ms = int(now.timestamp() * 1000)
    return {
        "Content-Type": "application/json",
        "Accept": "application/json, text/plain, */*",
        "Origin": "https://www.betonline.ag",
        "Referer": "https://www.betonline.ag/",
        "Authorization": f"Bearer {token}",
        "gsetting": "bolnasite",
        "contests": "na",
        "gmt-offset": "-8",
        "utc-offset": "480",
        "actual-time": str(ms),
        "iso-time": now.strftime("%Y-%m-%dT%H:%M:%S.") + f"{ms % 1000:03d}Z",
        "utc-time": now.strftime("%a, %d %b %Y %H:%M:%S GMT"),
    }


# ---------------------------------------------------------------------------
# Menu probe -> list of NFL Draft slugs
# ---------------------------------------------------------------------------

def _discover_slugs(session: cffi_requests.Session, headers: dict) -> list[str]:
    """Return the list of NFL Draft ContestType2 slugs from the sidebar menu.

    We walk the MenuItems tree and collect every child whose URL starts with
    `sportsbook/futures-and-props/nfl-draft/`. The final path segment is the
    slug (e.g. `1st-round`, `1st-round-props`, `draft-position`).

    Defined as live discovery rather than hardcoded so new buckets BetOnline
    adds between now and the draft flow through automatically.
    """
    resp = session.post(MENU_URL, data="", headers=headers, timeout=30)
    resp.raise_for_status()
    menu = resp.json()
    slugs: list[str] = []

    def walk(node):
        if isinstance(node, dict):
            url = node.get("URL") or ""
            if isinstance(url, str) and f"{CONTEST_TYPE}/" in url:
                slug = url.rstrip("/").rsplit("/", 1)[-1]
                # Root `nfl-draft` entry also has a URL pointing at 1st-round;
                # skip ones we already have so the list stays deduped.
                if slug and slug != CONTEST_TYPE and slug not in slugs:
                    slugs.append(slug)
            for v in node.values():
                walk(v)
        elif isinstance(node, list):
            for v in node:
                walk(v)

    walk(menu)
    return slugs


# ---------------------------------------------------------------------------
# Contest fetch per slug
# ---------------------------------------------------------------------------

def _fetch_contest(session: cffi_requests.Session, headers: dict,
                   slug: str) -> dict | None:
    """Fetch one contest-type2 bucket (e.g. `1st-round`). Returns parsed JSON."""
    body = {
        "ContestType": CONTEST_TYPE,
        "ContestType2": slug,
        "filterTime": 0,
    }
    try:
        resp = session.post(CONTESTS_URL, json=body, headers=headers, timeout=30)
    except Exception as exc:
        print(f"    slug={slug!r}: POST failed: {exc!r}")
        return None
    if resp.status_code != 200:
        snippet = resp.text[:200] if resp.text else ""
        print(f"    slug={slug!r}: HTTP {resp.status_code} — {snippet}")
        return None
    try:
        return resp.json()
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def _count_runners(contest_json: dict) -> int:
    """Sum total Contestants across all markets in one slug's response."""
    n = 0
    co = (contest_json or {}).get("ContestOfferings") or {}
    for dg in (co.get("DateGroup") or []):
        for desc in (dg.get("DescriptionGroup") or []):
            for tg in (desc.get("TimeGroup") or []):
                ce = tg.get("ContestExtended") or {}
                # ContestExtended can ship as dict OR list depending on count.
                ces = ce if isinstance(ce, list) else [ce]
                for c in ces:
                    for cgl in (c.get("ContestGroupLine") or []):
                        n += len(cgl.get("Contestants") or [])
    return n


if __name__ == "__main__":
    print("=" * 60)
    print("  BETONLINE NFL DRAFT RECON")
    print("=" * 60)

    session = cffi_requests.Session(impersonate="chrome")
    loaded = _load_cookies(session)
    print(f"  loaded {loaded} cookies from {COOKIE_PATH.name}")

    # Warm up Cloudflare.
    try:
        session.get(SITE_URL, timeout=20)
    except Exception as exc:
        print(f"  cloudflare warmup failed: {exc!r}")
        sys.exit(1)

    # Get anonymous JWT.
    try:
        token = session.get(TOKEN_URL, timeout=20).json()["token"]
    except Exception as exc:
        print(f"  get-token failed: {exc!r}")
        sys.exit(1)
    print(f"  got anonymous JWT ({len(token)} chars)")

    headers = _build_auth_headers(token)

    # Discover the sub-market slugs live.
    try:
        slugs = _discover_slugs(session, headers)
    except Exception as exc:
        print(f"  get-menu failed: {exc!r}")
        sys.exit(1)
    if not slugs:
        print("  no NFL Draft slugs found in the menu — cookies may be stale.")
        print("  Run: python bet_logger/recon_betonline.py")
        sys.exit(1)
    print(f"  discovered {len(slugs)} NFL Draft slugs: {slugs}")

    # Fetch each slug; bundle into one fixture keyed by slug.
    bundle: dict[str, dict] = {}
    total_runners = 0
    for slug in slugs:
        # refresh time headers per call (cheap, mirrors browser)
        headers = _build_auth_headers(token)
        data = _fetch_contest(session, headers, slug)
        if data is None:
            continue
        n = _count_runners(data)
        total_runners += n
        bundle[slug] = data
        print(f"    slug={slug!r}: {n} runners")

    if not bundle:
        print("  no buckets captured — aborting.")
        sys.exit(1)

    meta = {
        "method": "rest_with_cloudflare_cookies",
        "endpoint": CONTESTS_URL,
        "contest_type": CONTEST_TYPE,
        "slugs_captured": list(bundle.keys()),
        "total_runners": total_runners,
    }
    path = save_fixture("betonline", bundle, meta=meta)
    print_diagnostics("betonline", path, CONTESTS_URL, bundle)
    print(f"  total runners across all buckets: {total_runners}")
