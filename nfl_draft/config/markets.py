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


def _bm_entries() -> list[tuple[str, str, str, str]]:
    """Build Bookmaker MARKET_MAP rows from the committed fixture."""
    raw = _load_fixture("bookmaker")
    if raw is None:
        return []
    from nfl_draft.scrapers.bookmaker import (
        parse_response, PICK_HTM_RE, NTH_POS_HTM_RE, POSITION_MAP,
    )
    rows = parse_response(raw)

    entries: list[tuple[str, str, str, str]] = []
    seen: set[tuple[str, str, str]] = set()
    for r in rows:
        key = (r.book, r.book_label, r.book_subject)
        if key in seen:
            continue
        mid = _bm_market_id_for(
            r.market_group, r.book_label, r.book_subject,
            pick_re=PICK_HTM_RE, nth_pos_re=NTH_POS_HTM_RE,
            position_map=POSITION_MAP,
        )
        if mid is None:
            continue
        entries.append((*key, mid))
        seen.add(key)
    return entries


def _bm_market_id_for(group, label, subject, *, pick_re, nth_pos_re, position_map):
    """BM-specific market_id builder - mirrors DK/FD but keyed off BM's
    'htm' field (e.g. '1st Overall Pick 2026 NFL Draft')."""
    if group == "pick_outright":
        m = pick_re.search(label)
        if not m:
            return None
        return build_market_id(
            "pick_outright", pick_number=int(m.group(1)), player=subject,
        )
    if group == "first_at_position":
        m = nth_pos_re.match(label)
        if not m:
            return None
        pos = position_map.get(m.group(2).strip().lower())
        if not pos:
            return None
        return build_market_id("first_at_position", position=pos, player=subject)
    if group.startswith("nth_at_position_"):
        # e.g. group='nth_at_position_2'; the BM book_label carries the
        # position word ('2nd Wide Receiver Selected'). Re-parse to get
        # BOTH the nth and the position word.
        m = nth_pos_re.match(label)
        if not m:
            return None
        try:
            nth_val = int(m.group(1))
        except (ValueError, IndexError):
            return None
        pos = position_map.get(m.group(2).strip().lower())
        if not pos:
            return None
        return build_market_id(
            "nth_at_position", nth=nth_val, position=pos, player=subject,
        )
    return None


def _wz_entries() -> list[tuple[str, str, str, str]]:
    """Build Wagerzon MARKET_MAP rows from the committed fixture."""
    raw = _load_fixture("wagerzon")
    if raw is None:
        return []
    from nfl_draft.scrapers.wagerzon import (
        parse_response, PICK_DESC_RE, NTH_POS_HTM_RE, POSITION_MAP,
    )
    rows = parse_response(raw)

    entries: list[tuple[str, str, str, str]] = []
    seen: set[tuple[str, str, str]] = set()
    for r in rows:
        key = (r.book, r.book_label, r.book_subject)
        if key in seen:
            continue
        mid = _wz_market_id_for(
            r.market_group, r.book_label, r.book_subject,
            pick_desc_re=PICK_DESC_RE, nth_pos_re=NTH_POS_HTM_RE,
            position_map=POSITION_MAP,
        )
        if mid is None:
            continue
        entries.append((*key, mid))
        seen.add(key)
    return entries


def _h88_entries() -> list[tuple[str, str, str, str]]:
    """Build Hoop88 MARKET_MAP rows from the committed fixture.

    H88 ships one propDescription per market (book_label), and each
    contestant becomes a separate OddsRow. We map the structured groups
    (pick_outright / first_at_position) to canonical market_ids; the
    'Mr Irrelevant Position' prop stays unmapped (no canonical type).
    """
    raw = _load_fixture("hoop88")
    if raw is None:
        return []
    from nfl_draft.scrapers.hoop88 import (
        parse_response, PICK_LABEL_RE, FIRST_POS_LABEL_RE, POSITION_MAP,
    )
    rows = parse_response(raw)

    entries: list[tuple[str, str, str, str]] = []
    seen: set[tuple[str, str, str]] = set()
    for r in rows:
        key = (r.book, r.book_label, r.book_subject)
        if key in seen:
            continue
        mid = _h88_market_id_for(
            r.market_group, r.book_label, r.book_subject,
            pick_re=PICK_LABEL_RE, first_pos_re=FIRST_POS_LABEL_RE,
            position_map=POSITION_MAP,
        )
        if mid is None:
            continue  # props stay quarantined.
        entries.append((*key, mid))
        seen.add(key)
    return entries


def _h88_market_id_for(group, label, subject, *, pick_re, first_pos_re, position_map):
    """H88-specific market_id builder.

    book_label is the propDescription we queried with (e.g. '2nd Overall
    Pick', '1st Wide Receiver'). That string carries the pick_number /
    position, so we re-parse it here rather than threading extra state
    through the parser.
    """
    if group == "pick_outright":
        m = pick_re.match(label)
        if not m:
            return None
        return build_market_id(
            "pick_outright", pick_number=int(m.group(1)), player=subject,
        )
    if group == "first_at_position":
        m = first_pos_re.match(label)
        if not m:
            return None
        pos = position_map.get(m.group(1).strip().lower())
        if not pos:
            return None
        return build_market_id("first_at_position", position=pos, player=subject)
    return None


def _wz_market_id_for(group, label, subject, *, pick_desc_re, nth_pos_re, position_map):
    """WZ-specific market_id builder. book_label is typically the league
    Description (for pick_outright + top_n) or the game.htm (for
    first_at_position)."""
    if group == "pick_outright":
        m = pick_desc_re.search(label)
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
        m = nth_pos_re.match(label)
        if not m:
            return None
        pos_raw = m.group(2).strip().lower()
        pos = position_map.get(pos_raw)
        if not pos:
            return None
        return build_market_id("first_at_position", position=pos, player=subject)
    if group.startswith("nth_at_position_"):
        # Same logic as BM: re-parse the game.htm label to get nth + position.
        m = nth_pos_re.match(label)
        if not m:
            return None
        try:
            nth_val = int(m.group(1))
        except (ValueError, IndexError):
            return None
        pos = position_map.get(m.group(2).strip().lower())
        if not pos:
            return None
        return build_market_id(
            "nth_at_position", nth=nth_val, position=pos, player=subject,
        )
    return None


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
# Kalshi
# ---------------------------------------------------------------------------
#
# Unlike the sportsbooks, Kalshi doesn't expose a single "market group" field --
# we have to infer the canonical market_type from the series_ticker prefix and,
# for the PICK / TOP series, from the ticker's middle segment. Ticker anatomy:
#
#   KXNFLDRAFT1-26-TCHA                -> "who is 1st overall pick?" (subject=player)
#   KXNFLDRAFTPICK-26-10-TSIM          -> "who is 10th overall?" (pick_number in segment[2])
#   KXNFLDRAFTTOP-26-R1-GJAC           -> "1st round" (maps to top_32)
#   KXNFLDRAFTTOP-26-5-XXX             -> top_5 (range_high in segment[2])
#   KXNFLDRAFTQB-26P1-TSIM             -> "1st QB drafted" (only P1 maps; P2+ have no
#                                         canonical 'nth_at_position' type)
#
# Kalshi-only series (no structured canonical type yet) are intentionally left
# unmapped so they fall into draft_odds_unmapped:
#   - KXNFLDRAFT1ST     (team makes 1st overall)
#   - KXNFLDRAFTTEAM    (team drafts specific player)
#   - KXNFLDRAFTMATCHUP (player X drafted before player Y)
#   - KXNFLFIRSTPICK    (duplicate of KXNFLDRAFT1)
#   - KXNFLDRAFTOU      (over/under-pick lines)
#
# Position code for each series_ticker prefix. Used by _kalshi_entries for the
# first_at_position mapping.
_KALSHI_POSITION_SERIES = {
    "KXNFLDRAFTQB":   "QB",
    "KXNFLDRAFTRB":   "RB",
    "KXNFLDRAFTWR":   "WR",
    "KXNFLDRAFTTE":   "TE",
    "KXNFLDRAFTLB":   "LB",
    "KXNFLDRAFTDT":   "DT",
    "KXNFLDRAFTEDGE": "EDGE",
    "KXNFLDRAFTOL":   "OL",
    "KXNFLDRAFTDB":   "DB",
}


def _kalshi_entries() -> list[tuple[str, str, str, str]]:
    """Build Kalshi MARKET_MAP rows from the committed /markets fixture.

    The fixture is a per-series dict (captured at the recon step), because
    Kalshi's /markets endpoint requires a series_ticker query param -- unlike
    the sportsbook fixtures which are a single monolithic board response.

    book_label comes from scrapers.kalshi._kalshi_book_label: for PICK/TOP
    series it's a pick_number-scoped prefix so a single player appearing at
    multiple pick numbers doesn't collapse onto one MARKET_MAP key.

    Returns an empty list (not an error) if the fixture is missing.
    """
    path = _FIXTURES_ROOT / "kalshi" / "markets_response.json"
    if not path.exists():
        return []
    try:
        raw = json.loads(path.read_text())
    except Exception:
        return []

    # Fixture shape: {"series_responses": {series_ticker: {"markets": [...]}}}.
    # Older single-response shape ({"markets": [...]}) is ignored -- we log
    # empty and let the caller refresh the fixture.
    series_responses = raw.get("series_responses")
    if not isinstance(series_responses, dict):
        return []

    # Local import (mirrors DK/FD pattern; avoids circular import at module load).
    from nfl_draft.scrapers.kalshi import _kalshi_book_label

    entries: list[tuple[str, str, str, str]] = []
    seen: set[tuple[str, str, str]] = set()

    for series_ticker, resp in series_responses.items():
        for market in resp.get("markets", []) or []:
            ticker = market.get("ticker") or ""
            # Prefer yes_sub_title (what the parser emits as book_subject),
            # fall back to subtitle / custom_strike fields.
            subject = (
                market.get("yes_sub_title")
                or market.get("subtitle")
                or (market.get("custom_strike") or {}).get("Person")
                or (market.get("custom_strike") or {}).get("Team")
                or ""
            ).strip()
            if not subject:
                continue

            mid = _kalshi_market_id_for(series_ticker, ticker, subject)
            if mid is None:
                continue  # Kalshi-only / unmappable market -> quarantine at runtime.

            # Must match what the scraper emits at runtime.
            book_label = _kalshi_book_label(series_ticker, ticker)
            key = ("kalshi", book_label, subject)
            if key in seen:
                continue
            entries.append((*key, mid))
            seen.add(key)
    return entries


def _kalshi_market_id_for(series_ticker: str, ticker: str, subject: str) -> str | None:
    """Return a canonical market_id for a Kalshi market, or None if unmappable.

    Dispatches on series_ticker prefix, falling back to parsing the ticker's
    middle segment for series where pick_number / range_high varies per market.
    """
    # pick_outright, pick 1 only (series is dedicated to 1st overall).
    if series_ticker == "KXNFLDRAFT1":
        return build_market_id("pick_outright", pick_number=1, player=subject)

    # pick_outright, pick N extracted from ticker (e.g. KXNFLDRAFTPICK-26-10-TSIM).
    if series_ticker == "KXNFLDRAFTPICK":
        parts = ticker.split("-")
        if len(parts) >= 3 and parts[2].isdigit():
            return build_market_id(
                "pick_outright", pick_number=int(parts[2]), player=subject,
            )
        return None

    # top_n_range: segment[2] is either R1 (= top 32, 1st round) or a number.
    if series_ticker == "KXNFLDRAFTTOP":
        parts = ticker.split("-")
        if len(parts) < 3:
            return None
        tag = parts[2]
        if tag == "R1":
            high = 32
        elif tag.isdigit():
            high = int(tag)
        else:
            return None
        return build_market_id(
            "top_n_range", range_low=1, range_high=high, player=subject,
        )

    # first_at_position: only the P1 (1st-at-position) markets have a canonical
    # type. Nth-at-position (N>=2) is Kalshi-only.
    pos = _KALSHI_POSITION_SERIES.get(series_ticker)
    if pos is not None:
        parts = ticker.split("-")
        if len(parts) >= 2 and parts[1] == "26P1":
            return build_market_id(
                "first_at_position", position=pos, player=subject,
            )
        return None

    # Kalshi-only series (team_first_pick, matchup, team-drafts-player, etc.)
    # fall through unmapped on purpose.
    return None


# ---------------------------------------------------------------------------
# Aggregate MARKET_MAP (de-duped across all books)
# ---------------------------------------------------------------------------

MARKET_MAP: list[tuple[str, str, str, str]] = list(
    dict.fromkeys(
        _dk_entries()
        + _fd_entries()
        + _bm_entries()
        + _wz_entries()
        + _h88_entries()
        + _kalshi_entries()
    ).keys()
)
