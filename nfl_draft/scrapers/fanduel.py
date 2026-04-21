"""FanDuel NFL Draft markets adapter.

FD exposes NFL Draft futures inside the NFL content-managed-page bundle under
tab id 391 ('NFL Draft'). Each tab->card->coupon walk resolves to a market_id
inside `attachments.markets`; runners on each market carry both decimal and
*American* odds, already signed ints - much nicer than DK's string format.

Envelope shape (produced by recon_fd.py::run_rest_phase):
    {
      "book": "fanduel",
      "meta": { ... draft_tab_id / draft_market_count / etc ... },
      "data": {
        "layout": {
          "tabs":    { "391": {"cards":[{"id":...},...]}, ... },
          "cards":   { "<card_id>": {"title":..., "coupons":[{"id":...}]}, },
          "coupons": { "<coupon_id>": {"marketId":"734...", ...}, },
        },
        "attachments": {
          "markets": {
            "734.xxx": {
              "marketName":"2026 NFL Draft - Number 1 Overall Pick",
              "marketType":"NUMBER_X_OVERALL_PICK",
              "runners":[{"runnerName":"Fernando Mendoza",
                          "winRunnerOdds":{"americanDisplayOdds":{"americanOdds":-20000}}},
                         ...],
            },
          },
        },
      },
    }
"""

from __future__ import annotations

import re
import sys
from datetime import datetime
from pathlib import Path
from typing import List, Optional

from nfl_draft.scrapers._base import OddsRow


# ---------------------------------------------------------------------------
# Market-type classifier: FD's marketType strings -> internal market_group
# ---------------------------------------------------------------------------

# Each FD marketType maps to (market_group, any extra hints needed by the
# MARKET_MAP builder in config/markets.py).
_MARKET_TYPE_GROUPS = {
    "NUMBER_X_OVERALL_PICK":   "pick_outright",
    "TEAM_TO_DRAFT_PLAYER_X":  "prop_team_to_draft_player",
    "SPP_-_NFL_DRAFT_TOP_5":   "top_5_range",
    "SPP_-_NFL_DRAFT_TOP_10":  "top_10_range",
    "SPP_-_NFL_DRAFT_TOP_32":  "top_32_range",
    "FIRST_POSITION_X_DRAFTED": "first_at_position",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _american(runner: dict) -> Optional[int]:
    """Extract a signed American odd as int from a runner dict; None if absent."""
    odds = (runner.get("winRunnerOdds") or {}).get("americanDisplayOdds") or {}
    for key in ("americanOddsInt", "americanOdds"):
        val = odds.get(key)
        if val is None:
            continue
        try:
            iv = int(val)
            if iv == 0:
                return None  # FD sometimes ships 0 for suspended runners
            return iv
        except (TypeError, ValueError):
            continue
    return None


def _collect_market_ids(data: dict, tab_id: str) -> list[str]:
    """Walk layout.tabs[tab_id].cards[*].coupons[*] and collect every market_id
    that belongs to the NFL Draft tab.

    Handles both:
      - coupon.marketId (single-market coupons; most of the draft card)
      - coupon.display[*].rows[*].marketIds (multi-market coupon grids)
    """
    layout = data.get("layout") or {}
    tabs = layout.get("tabs") or {}
    cards_all = layout.get("cards") or {}
    coupons_all = layout.get("coupons") or {}

    tab = tabs.get(str(tab_id)) or {}
    out: list[str] = []
    for card_ref in tab.get("cards", []) or []:
        cid = card_ref.get("id") if isinstance(card_ref, dict) else card_ref
        card = cards_all.get(str(cid))
        if not card:
            continue
        # If the card has hasAttachments=false, its coupons don't have market
        # data in this bundle - skip entirely (FD lazy-loads some cards).
        if not card.get("hasAttachments", True):
            continue
        for coup_ref in card.get("coupons", []) or []:
            coup_id = coup_ref.get("id") if isinstance(coup_ref, dict) else coup_ref
            coup = coupons_all.get(str(coup_id))
            if not coup:
                continue
            if not coup.get("hasAttachments", True):
                continue
            mid = coup.get("marketId")
            if mid:
                out.append(str(mid))
            for disp in coup.get("display", []) or []:
                for row in disp.get("rows", []) or []:
                    for m in row.get("marketIds", []) or []:
                        out.append(str(m))
    # De-dup preserving order.
    seen: set[str] = set()
    uniq: list[str] = []
    for m in out:
        if m not in seen:
            seen.add(m)
            uniq.append(m)
    return uniq


def _find_draft_tab_id(data: dict) -> Optional[str]:
    """Find the tab whose title contains 'draft' (case-insensitive)."""
    tabs = (data.get("layout") or {}).get("tabs") or {}
    for tid, t in tabs.items():
        if "draft" in ((t or {}).get("title", "") or "").lower():
            return str(tid)
    return None


# ---------------------------------------------------------------------------
# Public parser
# ---------------------------------------------------------------------------

def parse_response(raw: dict) -> List[OddsRow]:
    """Walk FD's content-managed-page envelope and emit OddsRows for every
    runner on every market under the NFL Draft tab.

    Accepts either the full recon envelope (`raw["data"]`) or just the inner
    content-managed-page payload.
    """
    data = raw.get("data") if isinstance(raw, dict) and "data" in raw else raw
    if not isinstance(data, dict):
        return []

    markets = (data.get("attachments") or {}).get("markets") or {}
    tab_id = _find_draft_tab_id(data)
    if tab_id is None:
        # Fallback: if the tab structure is missing, scan all markets whose
        # name mentions "NFL Draft". Better than silently returning [].
        draft_market_ids = [
            mid for mid, m in markets.items()
            if "nfl draft" in (m.get("marketName") or "").lower()
        ]
    else:
        draft_market_ids = _collect_market_ids(data, tab_id)
        if not draft_market_ids:
            # Same fallback as above if the tab's cards didn't resolve.
            draft_market_ids = [
                mid for mid, m in markets.items()
                if "nfl draft" in (m.get("marketName") or "").lower()
            ]

    now = datetime.now()
    rows: List[OddsRow] = []

    for mid in draft_market_ids:
        m = markets.get(mid)
        if not m:
            continue
        m_name = (m.get("marketName") or "").strip()
        m_type = (m.get("marketType") or "").strip()
        if m.get("marketStatus") and m.get("marketStatus") != "OPEN":
            continue  # skip suspended / closed markets

        group = _MARKET_TYPE_GROUPS.get(m_type, f"prop_{m_type.lower()}")

        for runner in m.get("runners", []) or []:
            if runner.get("runnerStatus") and runner.get("runnerStatus") != "ACTIVE":
                continue
            name = (runner.get("runnerName") or "").strip()
            american = _american(runner)
            if not name or american is None:
                continue
            rows.append(OddsRow(
                book="fanduel",
                book_label=m_name,
                book_subject=name,
                american_odds=american,
                fetched_at=now,
                market_group=group,
            ))

    return rows


# ---------------------------------------------------------------------------
# Live fetch (delegates to recon_fd.py)
# ---------------------------------------------------------------------------

def fetch_raw() -> dict:
    """Call the recon script's REST fetch to refresh the FD envelope."""
    _THIS = Path(__file__).resolve()
    sys.path.insert(0, str(_THIS.parent.parent.parent))
    from nfl_draft.scrapers import recon_fd  # type: ignore[import-not-found]

    data, _url, meta = recon_fd.run_rest_phase()
    if data is None:
        raise RuntimeError("FD REST fetch returned no data; see recon_fd logs.")
    return {"book": "fanduel", "meta": meta, "data": data}


def fetch_draft_odds() -> List[OddsRow]:
    return parse_response(fetch_raw())


# ---------------------------------------------------------------------------
# MARKET_MAP helpers: regex helpers used by config/markets.py::_fd_entries()
# Kept close to the parser so any rename here stays consistent.
# ---------------------------------------------------------------------------

# Matches "2026 NFL Draft - Number 1 Overall Pick" OR "... : Number 1 Overall ..."
PICK_LABEL_RE = re.compile(r"Number\s+(\d+)\s+Overall\s+Pick", re.IGNORECASE)
# Matches "2026 NFL Draft - To Be a Top 5 Pick" / "To Be a 1st Round Pick"
TOP_N_LABEL_RE = re.compile(r"Top\s+(\d+)\s+Pick|1st Round|First Round", re.IGNORECASE)

# Matches "First Wide Receiver Drafted" / "First Offensive Lineman Drafted" etc.
FIRST_POS_LABEL_RE = re.compile(
    r"First\s+(.+?)\s+Drafted\s*$", re.IGNORECASE,
)
POSITION_MAP = {
    "quarterback": "QB", "qb": "QB",
    "running back": "RB", "rb": "RB",
    "wide receiver": "WR", "wr": "WR",
    "tight end": "TE", "te": "TE",
    "cornerback": "CB",
    "safety": "S",
    "linebacker": "LB",
    "defensive back": "DB",
    "defensive tackle": "DT",
    "defensive end": "DE",
    "edge": "EDGE",
    "offensive lineman": "OL", "offensive line": "OL", "ol": "OL",
}
