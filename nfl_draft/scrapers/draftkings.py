"""DraftKings NFL Draft markets adapter.

Pulls NFL Draft futures from DK's public sportscontent REST API and turns them
into OddsRow records for the nfl_draft pipeline. The recon script
(nfl_draft/scrapers/recon_dk.py) handles live capture + fixture generation; the
parser here consumes that same shape, so it works identically online and from
the saved fixture.

DK fixture shape (from recon_dk.py::run_rest_phase envelope):
    {
      "book": "draftkings",
      "data": {
        "league_id": 88808,
        "draft_category_id": 1803,
        "draft_category_name": "2026 Draft",
        "draft_subcategories": {
          "<sub_id>": {
            "markets":    [ {id, name, marketType, subcategoryId, ...}, ... ],
            "selections": [ {marketId, label, displayOdds.american, ...}, ... ],
            ...
          }, ...
        },
        "subcategory_catalog": [ {id, name, categoryId}, ... ],
      }
    }

DK ships American odds as strings with a Unicode minus sign ('\u2212', not '-').
We normalize that before parsing.
"""

from __future__ import annotations

import re
import sys
from datetime import datetime
from pathlib import Path
from typing import List, Optional

from nfl_draft.scrapers._base import OddsRow
from nfl_draft.lib.market_map import build_market_id, slug


# ---------------------------------------------------------------------------
# Low-level helpers
# ---------------------------------------------------------------------------

def _american_to_int(raw: str) -> Optional[int]:
    """Convert DK's 'displayOdds.american' string (e.g. '\u2212150', '+180',
    'EVEN') into a signed int. Returns None if the string is unparseable.

    DK uses Unicode minus (U+2212) not ASCII '-'; handle both.
    """
    if raw is None:
        return None
    s = str(raw).strip().replace("\u2212", "-").replace("+", "")
    if not s:
        return None
    if s.upper() in ("EVEN", "EV", "PK"):
        return 100
    try:
        return int(s)
    except ValueError:
        # Decimal-looking odds? Only happens if DK ever rotates format.
        try:
            return int(round(float(s)))
        except ValueError:
            return None


# Maps "1st <PositionName> Selected" / "2nd <PositionName> Selected" strings
# to our abbreviated position codes (used in market_id construction).
_POSITION_ALIASES = {
    "quarterback": "QB", "qb": "QB",
    "running back": "RB", "rb": "RB",
    "wide receiver": "WR", "wr": "WR",
    "tight end": "TE", "te": "TE",
    "cornerback": "CB", "cb": "CB",
    "safety": "S", "s": "S",
    "linebacker": "LB", "lb": "LB",
    "defensive tackle": "DT", "dt": "DT",
    "defensive end": "DE", "de": "DE",
    "edge": "EDGE",
    "defensive line/edge": "DL",
    "defensive line": "DL", "dl": "DL",
    "offensive lineman": "OL", "offensive line": "OL", "ol": "OL",
    "kicker/punter/long snapper": "ST",
    "kicker": "K", "k": "K",
    "punter": "P",
}


def _normalize_position(raw: str) -> Optional[str]:
    """Strip '1st'/'2nd' + 'Selected'/'Drafted' fluff, return canonical abbr."""
    if not raw:
        return None
    s = raw.strip().lower()
    # Drop common words so "1st Wide Receiver Selected" -> "wide receiver"
    for junk in (
        "1st ", "2nd ", "3rd ", "first ", "second ", "third ",
        " selected", " drafted", " off the board", "the first ",
    ):
        s = s.replace(junk, "")
    s = s.strip()
    return _POSITION_ALIASES.get(s, raw.strip().upper() if len(raw.strip()) <= 5 else None)


# ---------------------------------------------------------------------------
# Subcategory classifiers
# ---------------------------------------------------------------------------

# Matches "Number 1 Pick", "Number 10 Pick", "Number 11 Pick [All Bets Action]".
_PICK_NUMBER_RE = re.compile(r"Number\s+(\d+)\s+Pick", re.IGNORECASE)
# Matches "To be Drafted Top 5", "To be Drafted Top 10", "To be Drafted in Round 1".
_TOP_N_RE = re.compile(r"Top\s+(\d+)|Round\s+(\d+)", re.IGNORECASE)
# Matches "Which Team will <player> be Drafted By?"
_PLAYER_DRAFTED_BY_RE = re.compile(
    r"Which Team will\s+(.+?)\s+be Drafted By", re.IGNORECASE,
)
# "1st Wide Receiver Selected", "2nd Quarterback Selected", "1st OL Selected"
_NTH_POSITION_RE = re.compile(
    r"^(1st|2nd|3rd|First|Second|Third)\s+(.+?)\s+Selected\s*$", re.IGNORECASE,
)


def _emit_pick_outright(market: dict, sels: list, now: datetime) -> list[OddsRow]:
    """Emit rows for a 'Number N Pick' market (market_type = pick_outright)."""
    m_name = (market.get("name") or "").strip()
    match = _PICK_NUMBER_RE.search(m_name)
    if not match:
        return []
    # We use the bracket-free market name as book_label so
    # 'Number 1 Pick [All Bets Action]' doesn't drift over time.
    book_label = f"Number {int(match.group(1))} Pick"
    rows: list[OddsRow] = []
    for s in sels:
        label = (s.get("label") or "").strip()
        american = _american_to_int((s.get("displayOdds") or {}).get("american"))
        if not label or american is None:
            continue
        rows.append(OddsRow(
            book="draftkings",
            book_label=book_label,
            book_subject=label,
            american_odds=american,
            fetched_at=now,
            market_group="pick_outright",
        ))
    return rows


def _emit_top_n(market: dict, sels: list, now: datetime) -> list[OddsRow]:
    """Emit rows for 'To be Drafted Top N' / 'Round 1 Pick' markets
    (market_type = top_n_range)."""
    m_name = (market.get("name") or "").strip()
    match = _TOP_N_RE.search(m_name)
    if not match:
        return []
    # Round 1 has 32 picks.
    high = int(match.group(1) or 0) or (32 if match.group(2) else 0)
    if not high:
        return []
    book_label = m_name  # Preserve the exact DK market name
    rows: list[OddsRow] = []
    for s in sels:
        label = (s.get("label") or "").strip()
        american = _american_to_int((s.get("displayOdds") or {}).get("american"))
        if not label or american is None:
            continue
        rows.append(OddsRow(
            book="draftkings",
            book_label=book_label,
            book_subject=label,
            american_odds=american,
            fetched_at=now,
            market_group=f"top_{high}_range",
        ))
    return rows


def _emit_first_at_position(market: dict, sels: list, now: datetime) -> list[OddsRow]:
    """Emit rows for '1st/2nd <Position> Selected' markets
    (market_type = first_at_position)."""
    m_name = (market.get("name") or "").strip()
    match = _NTH_POSITION_RE.match(m_name)
    if not match:
        return []
    book_label = m_name
    rows: list[OddsRow] = []
    for s in sels:
        label = (s.get("label") or "").strip()
        american = _american_to_int((s.get("displayOdds") or {}).get("american"))
        if not label or american is None:
            continue
        rows.append(OddsRow(
            book="draftkings",
            book_label=book_label,
            book_subject=label,
            american_odds=american,
            fetched_at=now,
            market_group="first_at_position",
        ))
    return rows


def _emit_prop(market: dict, sels: list, now: datetime, group: str) -> list[OddsRow]:
    """Fallback: emit rows for anything that doesn't fit the 5 core types."""
    m_name = (market.get("name") or "").strip()
    rows: list[OddsRow] = []
    for s in sels:
        label = (s.get("label") or "").strip()
        american = _american_to_int((s.get("displayOdds") or {}).get("american"))
        if not label or american is None:
            continue
        rows.append(OddsRow(
            book="draftkings",
            book_label=m_name,
            book_subject=label,
            american_odds=american,
            fetched_at=now,
            market_group=group,
        ))
    return rows


# Subcategory NAMES -> handler. Unknown subcategories fall through to _emit_prop.
# Using the subcategory_catalog `name` (not the numeric id) keeps us resilient
# when DK renumbers ids between seasons.
_SUB_HANDLERS: dict[str, tuple[str, callable]] = {
    "Pick Number":           ("pick_outright",     _emit_pick_outright),
    "Pick Number 6-10":      ("pick_outright",     _emit_pick_outright),
    "Pick Number 11-15":     ("pick_outright",     _emit_pick_outright),
    "Top 5 Pick":            ("top_n_range",       _emit_top_n),
    "Top 10 Pick":           ("top_n_range",       _emit_top_n),
    "Round 1 Pick":          ("top_n_range",       _emit_top_n),
    "1st Selected By Position": ("first_at_position", _emit_first_at_position),
    "2nd Selected Position":    ("first_at_position", _emit_first_at_position),
}


# ---------------------------------------------------------------------------
# Public parser
# ---------------------------------------------------------------------------

def parse_response(raw: dict) -> List[OddsRow]:
    """Walk a DK draft envelope (or its `.data` inner dict) and emit one
    OddsRow per selection.

    Accepts either the full envelope (with `book`/`data` keys from the recon
    script) or just the inner `data` dict. This keeps fixtures + live calls
    compatible.
    """
    # Unwrap recon envelope if present.
    data = raw.get("data") if isinstance(raw, dict) and "data" in raw else raw
    if not isinstance(data, dict):
        return []

    subs = data.get("draft_subcategories") or {}
    catalog = data.get("subcategory_catalog") or []
    # Build sub_id -> sub_name index so a subcategory w/o its own catalog entry
    # still has a human-readable name.
    id_to_name = {str(s["id"]): s.get("name", "") for s in catalog if "id" in s}

    now = datetime.now()
    rows: List[OddsRow] = []

    for sub_id, sub in subs.items():
        sub_name = id_to_name.get(str(sub_id), "")
        handler_tuple = _SUB_HANDLERS.get(sub_name)
        # Build a marketId -> selections index (selections carry their own
        # marketId; markets have no 'selections' array inline).
        sels_by_market: dict[str, list] = {}
        for s in sub.get("selections", []) or []:
            mid = str(s.get("marketId") or "")
            if mid:
                sels_by_market.setdefault(mid, []).append(s)

        for m in sub.get("markets", []) or []:
            mid = str(m.get("id") or "")
            sels = sels_by_market.get(mid, [])
            if not sels:
                continue
            if handler_tuple is not None:
                _, handler = handler_tuple
                rows.extend(handler(m, sels, now))
            else:
                # Subcategories not in _SUB_HANDLERS (Mr. Irrelevant, Draft
                # Position, Position Totals, College Props, Player H2Hs,
                # Player Drafted By) all become props - they don't fit the
                # 5 structured market_types cleanly.
                rows.extend(_emit_prop(m, sels, now, group=f"prop_{slug(sub_name)}"))

    return rows


# ---------------------------------------------------------------------------
# Live fetch (delegates to recon_dk.py for auth + REST orchestration)
# ---------------------------------------------------------------------------

def fetch_raw() -> dict:
    """Call the recon script's REST path to get a fresh DK envelope.

    We reuse the recon module's `run_rest_phase` instead of duplicating the
    Akamai / curl_cffi setup. Returns the same shape parse_response() expects.
    """
    # Local import so the scraper can be imported in test/fixture contexts
    # without requiring curl_cffi to be installed.
    _THIS = Path(__file__).resolve()
    sys.path.insert(0, str(_THIS.parent.parent.parent))
    from nfl_draft.scrapers import recon_dk  # type: ignore[import-not-found]

    data, _url, meta = recon_dk.run_rest_phase()
    if data is None:
        raise RuntimeError("DK REST fetch returned no data; see recon_dk logs.")
    # Mirror the envelope shape save_fixture() writes, so parse_response is
    # agnostic to origin (fixture vs live).
    return {"book": "draftkings", "meta": meta, "data": data}


def fetch_draft_odds() -> List[OddsRow]:
    """Thin wrapper: live fetch -> parse."""
    return parse_response(fetch_raw())
