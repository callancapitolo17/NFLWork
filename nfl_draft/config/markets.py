"""Per-book market label -> canonical market_id mappings.

MARKET_MAP is the authoritative list of (book, book_label, book_subject,
market_id) tuples consumed by lib/quarantine.py::write_or_quarantine: a row
whose (book, book_label, book_subject) key isn't in this list falls into
draft_odds_unmapped for later triage.

Population strategy
-------------------
Manually listing hundreds of player names gets stale the minute a book
re-shuffles its board. So instead, each `_<book>_entries()` helper:

  1. Loads that book's checked-in fixture
     (nfl_draft/tests/fixtures/<book>/draft_markets.json).
  2. Runs the same parse_response() the live scraper uses.
  3. For every OddsRow that comes back with a *structured* market_group
     (pick_outright / first_at_position / top_N_range), emits a MARKET_MAP
     tuple using build_market_id() with the row's own book_label + book_subject.
     Props are left unmapped on purpose - they land in draft_odds_unmapped
     for later triage.

This keeps the construction rules in lib/market_map.py as the single source
of truth: MARKET_MAP never drifts from build_market_id(). When a book rotates
its roster, re-run the recon capture, commit the fixture, and every new
player is automatically mapped at next seed.

If a fixture is missing (e.g. during a fresh clone before recon has run)
the helper returns [] so this module still imports cleanly.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

from nfl_draft.lib.market_map import build_market_id


_FIXTURES_ROOT = Path(__file__).resolve().parent.parent / "tests" / "fixtures"


def _load_fixture(book: str) -> dict | None:
    """Load a book's draft fixture JSON, or None if absent/unparseable."""
    path = _FIXTURES_ROOT / book / "draft_markets.json"
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text())
    except Exception:
        return None


# ---------------------------------------------------------------------------
# DraftKings
# ---------------------------------------------------------------------------

# "Number 1 Pick" -> 1.  Used at MARKET_MAP build time so each pick_outright
# row ends up with the right `pick_number` arg into build_market_id().
_DK_PICK_LABEL_RE = re.compile(r"Number\s+(\d+)\s+Pick", re.IGNORECASE)

# "To be Drafted Top 5" / "To be Drafted in Round 1"
_DK_TOPN_LABEL_RE = re.compile(r"Top\s+(\d+)|Round\s+(\d+)", re.IGNORECASE)

# "1st Wide Receiver Selected" -> WR.  Position codes intentionally match
# lib/market_map.py's build_market_id('first_at_position', position=...).
_DK_POSITION_MAP = {
    "quarterback": "QB", "qb": "QB",
    "running back": "RB", "rb": "RB",
    "wide receiver": "WR", "wr": "WR",
    "tight end": "TE", "te": "TE",
    "cornerback": "CB", "cb": "CB",
    "safety": "S",
    "linebacker": "LB", "lb": "LB",
    "defensive tackle": "DT",
    "defensive end": "DE",
    "edge": "EDGE",
    "offensive lineman": "OL", "offensive line": "OL", "ol": "OL",
}
_DK_FIRST_POS_RE = re.compile(
    r"^1st\s+(.+?)\s+Selected\s*$", re.IGNORECASE,
)


def _dk_entries() -> list[tuple[str, str, str, str]]:
    """Build DK MARKET_MAP rows from the committed fixture.

    Local import of the scraper avoids circular imports at module-load time
    (scrapers.draftkings has no deps on this module, but keeping the import
    lazy is cheap insurance).
    """
    raw = _load_fixture("draftkings")
    if raw is None:
        return []
    from nfl_draft.scrapers.draftkings import parse_response
    rows = parse_response(raw)

    entries: list[tuple[str, str, str, str]] = []
    seen: set[tuple[str, str, str]] = set()
    for r in rows:
        key = (r.book, r.book_label, r.book_subject)
        if key in seen:
            continue
        mid = _dk_market_id_for(r.market_group, r.book_label, r.book_subject)
        if mid is None:
            # prop / unknown -> leave in quarantine.
            continue
        entries.append((*key, mid))
        seen.add(key)
    return entries


def _dk_market_id_for(group: str, label: str, subject: str) -> str | None:
    """Return a canonical market_id for a DK row, or None for props."""
    if group == "pick_outright":
        m = _DK_PICK_LABEL_RE.search(label)
        if not m:
            return None
        return build_market_id(
            "pick_outright", pick_number=int(m.group(1)), player=subject,
        )
    if group.startswith("top_") and group.endswith("_range"):
        # group is 'top_5_range' / 'top_10_range' / 'top_32_range'
        try:
            high = int(group.split("_")[1])
        except (IndexError, ValueError):
            return None
        return build_market_id(
            "top_n_range", range_low=1, range_high=high, player=subject,
        )
    if group == "first_at_position":
        m = _DK_FIRST_POS_RE.match(label)
        if not m:
            return None
        pos = _DK_POSITION_MAP.get(m.group(1).strip().lower())
        if not pos:
            return None
        return build_market_id("first_at_position", position=pos, player=subject)
    # 'prop_*' markets intentionally unmapped.
    return None


# ---------------------------------------------------------------------------
# FanDuel
# ---------------------------------------------------------------------------

def _fd_entries() -> list[tuple[str, str, str, str]]:
    """Build FD MARKET_MAP rows from the committed fixture."""
    raw = _load_fixture("fanduel")
    if raw is None:
        return []
    from nfl_draft.scrapers.fanduel import parse_response, PICK_LABEL_RE, TOP_N_LABEL_RE, FIRST_POS_LABEL_RE, POSITION_MAP
    rows = parse_response(raw)

    entries: list[tuple[str, str, str, str]] = []
    seen: set[tuple[str, str, str]] = set()
    for r in rows:
        key = (r.book, r.book_label, r.book_subject)
        if key in seen:
            continue
        mid = _fd_market_id_for(
            r.market_group, r.book_label, r.book_subject,
            pick_re=PICK_LABEL_RE, topn_re=TOP_N_LABEL_RE,
            first_pos_re=FIRST_POS_LABEL_RE, position_map=POSITION_MAP,
        )
        if mid is None:
            continue  # leave props unmapped
        entries.append((*key, mid))
        seen.add(key)
    return entries


def _fd_market_id_for(group, label, subject, *, pick_re, topn_re, first_pos_re, position_map):
    """Return a canonical market_id for an FD row, or None for props."""
    if group == "pick_outright":
        m = pick_re.search(label)
        if not m:
            return None
        return build_market_id(
            "pick_outright", pick_number=int(m.group(1)), player=subject,
        )
    if group.startswith("top_") and group.endswith("_range"):
        try:
            high = int(group.split("_")[1])
        except (IndexError, ValueError):
            return None
        return build_market_id(
            "top_n_range", range_low=1, range_high=high, player=subject,
        )
    if group == "first_at_position":
        m = first_pos_re.search(label)
        if not m:
            return None
        pos_raw = m.group(1).strip().lower()
        pos = position_map.get(pos_raw)
        if not pos:
            return None
        return build_market_id("first_at_position", position=pos, player=subject)
    return None


# ---------------------------------------------------------------------------
# Aggregate MARKET_MAP (de-duped across all books)
# ---------------------------------------------------------------------------

MARKET_MAP: list[tuple[str, str, str, str]] = list(
    dict.fromkeys(
        _dk_entries()
        + _fd_entries()
        # BM / WZ entries appended as each book's parser lands.
    ).keys()
)
