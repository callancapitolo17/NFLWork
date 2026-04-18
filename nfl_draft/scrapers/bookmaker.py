"""Bookmaker.eu NFL Draft markets adapter.

BM's BetslipProxy.aspx/GetSchedule returns an ASP.NET envelope with
Schedule.Data.Leagues.League[].dateGroup[].game[].Derivatives.line[].

Each `game` is one market (e.g. '1st Overall Pick 2026 NFL Draft',
'1st Wide Receiver Selected'). Each `line` is one runner with
`tmname` (player) and `odds` (signed American int as a string).

Envelope shape (produced by recon_bm.py::run_browser_phase):
    {
      "book": "bookmaker",
      "meta": {"league_id": "13425", "league_name": "ODDS TO WIN", ...},
      "data": {
        "Schedule": {
          "Data": {
            "Leagues": {
              "League": [
                {
                  "Description": "ODDS TO WIN",
                  "dateGroup": [
                    {
                      "date": "Apr 23",
                      "game": [
                        {
                          "htm": "1st Overall Pick 2026 NFL Draft",
                          "gdesc": "1st Round",
                          "idgm": "46094039",
                          "Derivatives": {
                            "line": [
                              {"tmname": "Fernando Mendoza", "odds": "-30000"},
                              ...
                            ]
                          }
                        }, ...
                      ]
                    }
                  ]
                }
              ]
            }
          }
        }
      }
    }

BM lists may ship as a dict (single item) or a list (multiple items);
_as_list() normalizes that.
"""

from __future__ import annotations

import re
import sys
from datetime import datetime
from pathlib import Path
from typing import List, Optional

from nfl_draft.scrapers._base import OddsRow


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _as_list(x) -> list:
    """BM's ASP.NET XML-to-JSON bridge ships a singleton as a dict and a
    multi as a list. Normalize so callers always get a list."""
    if x is None:
        return []
    if isinstance(x, list):
        return x
    return [x]


def _american(raw: str) -> Optional[int]:
    """BM ships American odds as string ints like '-110', '15000'. 0 means
    'no line' - treat it the same as missing."""
    if raw is None:
        return None
    try:
        iv = int(str(raw).strip())
    except (TypeError, ValueError):
        return None
    if iv == 0:
        return None
    return iv


# ---------------------------------------------------------------------------
# Market classifier
# ---------------------------------------------------------------------------

# "1st Overall Pick" / "10th Overall Pick" / "2nd Overall Pick"
_PICK_HTM_RE = re.compile(
    r"(\d+)(?:st|nd|rd|th)\s+Overall\s+Pick", re.IGNORECASE,
)
# "1st Wide Receiver Selected" / "2nd Offensive Line Selected" / etc
_NTH_POS_HTM_RE = re.compile(
    r"^(\d+)(?:st|nd|rd|th)\s+(.+?)\s+Selected", re.IGNORECASE,
)


def _classify_market(htm: str, gdesc: str) -> tuple[str, dict]:
    """Return (market_group, hints_dict) for a BM 'game'.

    `hints_dict` carries the derived pick_number / position / nth so the
    MARKET_MAP builder can pass them to build_market_id() without
    re-parsing.
    """
    s = (htm or "").strip()

    m = _PICK_HTM_RE.search(s)
    if m:
        return "pick_outright", {"pick_number": int(m.group(1))}

    m = _NTH_POS_HTM_RE.match(s)
    if m:
        nth = int(m.group(1))
        pos_raw = m.group(2).strip()
        group = "first_at_position" if nth == 1 else f"nth_at_position_{nth}"
        return group, {"nth": nth, "position_raw": pos_raw}

    # Everything else (Mr. Irrelevant, etc) -> generic prop.
    return f"prop_{(gdesc or 'bm').lower().replace(' ', '_')}", {}


# ---------------------------------------------------------------------------
# Public parser
# ---------------------------------------------------------------------------

def parse_response(raw: dict) -> List[OddsRow]:
    """Walk BM's Schedule envelope and emit OddsRows for every line on every
    game under every draft league.

    Accepts either the full recon envelope or just the inner
    `Schedule`-bearing payload.
    """
    data = raw.get("data") if isinstance(raw, dict) and "data" in raw else raw
    if not isinstance(data, dict):
        return []

    sched = data.get("Schedule") or {}
    leagues = _as_list(((sched.get("Data") or {}).get("Leagues") or {}).get("League"))

    now = datetime.now()
    rows: List[OddsRow] = []

    for lg in leagues:
        for dg in _as_list(lg.get("dateGroup")):
            for game in _as_list(dg.get("game")):
                htm = (game.get("htm") or "").strip()
                gdesc = (game.get("gdesc") or "").strip()
                if not htm:
                    continue
                group, _hints = _classify_market(htm, gdesc)
                lines = _as_list((game.get("Derivatives") or {}).get("line"))
                for ln in lines:
                    name = (ln.get("tmname") or "").strip()
                    american = _american(ln.get("odds"))
                    if not name or american is None:
                        continue
                    rows.append(OddsRow(
                        book="bookmaker",
                        book_label=htm,
                        book_subject=name,
                        american_odds=american,
                        fetched_at=now,
                        market_group=group,
                    ))
    return rows


# ---------------------------------------------------------------------------
# Live fetch (delegates to recon_bm.py)
# ---------------------------------------------------------------------------

def fetch_raw() -> dict:
    """Refresh BM via the recon script's REST probe path.

    Note: BM's NFL Draft league ID is discovered live (recon script writes
    `meta.league_id` into the saved fixture). If that stale ID is rotated by
    BM, the recon probe will surface a new one and this adapter picks it up
    automatically.
    """
    _THIS = Path(__file__).resolve()
    sys.path.insert(0, str(_THIS.parent.parent.parent))
    from nfl_draft.scrapers import recon_bm  # type: ignore[import-not-found]

    data, _url, meta = recon_bm.run_rest_phase()
    if data is None:
        raise RuntimeError("BM REST fetch returned no data; see recon_bm logs.")
    return {"book": "bookmaker", "meta": meta, "data": data}


def fetch_draft_odds() -> List[OddsRow]:
    return parse_response(fetch_raw())


# ---------------------------------------------------------------------------
# MARKET_MAP helpers
# ---------------------------------------------------------------------------

PICK_HTM_RE = _PICK_HTM_RE
NTH_POS_HTM_RE = _NTH_POS_HTM_RE

# BM position words -> canonical abbr used in build_market_id('first_at_position').
POSITION_MAP = {
    "quarterback": "QB", "qb": "QB",
    "running back": "RB", "rb": "RB",
    "wide receiver": "WR", "wr": "WR",
    "tight end": "TE", "te": "TE",
    "cornerback": "CB",
    "safety": "S",
    "linebacker": "LB",
    "offensive line": "OL", "offensive lineman": "OL", "ol": "OL",
    "defensive back": "DB",
    "defensive tackle": "DT",
    "defensive end": "DE",
    "edge": "EDGE",
}
