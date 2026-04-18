"""Wagerzon NFL Draft markets adapter.

WZ's NewScheduleHelper endpoint returns an envelope with
`result.listLeagues` - a list of buckets, each bucket a list of leagues.
Each league has a `Description` (e.g. 'NFL DRAFT 2026 - NUMBER 1 OVERALL
PICK') and `Games`.

Three main shapes surface across WZ's 32+ draft leagues:

  A) Single game, many GameLines (Picks 1-2): one market, N runners.
         game.htm     = market name
         GameLines[*] = runners via .tmname + .odds

  B) Many games, each with 1 GameLine (Picks 3-11, Position-Selected):
         game.htm     = market name (same across games in this league)
         GameLines[0].tmname = player name, .odds = American int str

  C) Many games, each with 1 over/under-style GameLine (Top N, Player
     Draft Position, Totals): the player is embedded in game.htm
     ('(DRAFTED TOP 5) ARVELL REESE'); GameLine carries odds via .odds.

We classify by league Description + htm parsing, then emit one OddsRow per
discrete (market, runner). Lines with empty/zero odds are skipped (WZ
seeds every market with an '*ALL BETS ACTION*' placeholder line that has
blank odds).

Envelope shape (recon_wz.py::run_rest_phase):
    {
      "book": "wagerzon",
      "meta": {"draft_league_ids": [...], "lg_param": "1270,2867,..."},
      "data": {
        "result": {
          "listLeagues": [
            [{"Description":"...", "IdLeague": 1270, "Games":[...]}],
            ...
          ]
        }
      }
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
# Helpers
# ---------------------------------------------------------------------------

def _american(raw) -> Optional[int]:
    """WZ ships odds like '-110' / '20000' / '' as strings. '' = no line."""
    if raw is None:
        return None
    s = str(raw).strip()
    if not s or s == "+":
        return None
    try:
        iv = int(s)
    except ValueError:
        return None
    if iv == 0:
        return None
    return iv


# "NFL DRAFT 2026 - NUMBER 1 OVERALL PICK" -> pick_number=1
_PICK_DESC_RE = re.compile(
    r"NUMBER\s+(\d+)\s+OVERALL\s+PICK", re.IGNORECASE,
)
# "NFL DRAFT 2026 - NUMBER 11-15 OVERALL PICK" -> range
_PICK_RANGE_DESC_RE = re.compile(
    r"NUMBER\s+(\d+)-(\d+)\s+OVERALL\s+PICK", re.IGNORECASE,
)
# "NFL DRAFT 2026 - 1ST CORNERBACK SELECTED" -> nth=1, pos='CORNERBACK'
_NTH_POS_HTM_RE = re.compile(
    r"^NFL\s+DRAFT\s+\d+\s*-\s*(\d+)(?:ST|ND|RD|TH)\s+(.+?)\s+SELECTED",
    re.IGNORECASE,
)
# "(DRAFTED TOP 5) ARVELL REESE" / "(DRAFTED TOP 10) ..."
_TOPN_HTM_RE = re.compile(
    r"\(DRAFTED\s+TOP\s+(\d+)\)\s+(.+)", re.IGNORECASE,
)
# "(DRAFTED IN R1) DILLON THIENEMAN"
_R1_HTM_RE = re.compile(r"\(DRAFTED\s+IN\s+R(\d+)\)\s+(.+)", re.IGNORECASE)
# Specific per-pick htm in the 11-15 league: "NUMBER 11 OVERALL PICK"
_PICK_GAME_HTM_RE = re.compile(
    r"NUMBER\s+(\d+)\s+OVERALL\s+PICK", re.IGNORECASE,
)

# WZ writes position names in all-caps (CORNERBACK, OL, WIDE RECEIVER, ...).
POSITION_MAP = {
    "cornerback": "CB", "cb": "CB",
    "wide receiver": "WR", "wr": "WR",
    "ol": "OL", "offensive lineman": "OL", "offensive line": "OL",
    "quarterback": "QB", "qb": "QB",
    "safety": "S",
    "running back": "RB", "rb": "RB",
    "tight end": "TE", "te": "TE",
    "linebacker": "LB", "lb": "LB",
    "defensive back": "DB", "db": "DB",
    "defensive tackle": "DT",
    "defensive end": "DE",
    "edge": "EDGE",
}


def _as_list(x) -> list:
    """Normalize WZ's 'might be dict or list' field shapes."""
    if x is None:
        return []
    if isinstance(x, list):
        return x
    return [x]


# ---------------------------------------------------------------------------
# Per-shape emitters
# ---------------------------------------------------------------------------

def _emit_pick_single_game(lg: dict, game: dict, now: datetime) -> list[OddsRow]:
    """Shape A: one game, many GameLines (Picks 1-2 on WZ).

    league Description encodes the pick number; each GameLine has
    `tmname` + `odds`.
    """
    desc = (lg.get("Description") or "").strip()
    m = _PICK_DESC_RE.search(desc)
    if not m:
        return []
    pick = int(m.group(1))
    book_label = desc  # use league Description as stable label
    rows: list[OddsRow] = []
    for ln in _as_list(game.get("GameLines")):
        name = (ln.get("tmname") or "").strip()
        if not name or name.startswith("*ALL BETS ACTION*"):
            continue
        american = _american(ln.get("odds"))
        if american is None:
            continue
        rows.append(OddsRow(
            book="wagerzon", book_label=book_label, book_subject=name,
            american_odds=american, fetched_at=now,
            market_group="pick_outright",
        ))
    return rows


def _emit_multi_game_pick(lg: dict, now: datetime) -> list[OddsRow]:
    """Shape B variant: one league, many games (one runner per game).

    The 11-15 league lumps picks 11/12/13/14/15 into one bucket, so pick
    number comes from each game's htm, not the league Description.
    """
    rows: list[OddsRow] = []
    for g in _as_list(lg.get("Games")):
        htm = (g.get("htm") or "").strip()
        m = _PICK_GAME_HTM_RE.search(htm)
        if not m:
            continue
        pick = int(m.group(1))
        book_label = f"NFL DRAFT 2026 - NUMBER {pick} OVERALL PICK"
        for ln in _as_list(g.get("GameLines")):
            name = (ln.get("tmname") or "").strip()
            if not name or name.startswith("*ALL BETS ACTION*"):
                continue
            american = _american(ln.get("odds"))
            if american is None:
                continue
            rows.append(OddsRow(
                book="wagerzon", book_label=book_label, book_subject=name,
                american_odds=american, fetched_at=now,
                market_group="pick_outright",
            ))
    return rows


def _emit_first_at_position(lg: dict, now: datetime) -> list[OddsRow]:
    """Shape B: league '2579' splits 8 first/second-at-position markets,
    each as its own `htm`."""
    rows: list[OddsRow] = []
    for g in _as_list(lg.get("Games")):
        htm = (g.get("htm") or "").strip()
        m = _NTH_POS_HTM_RE.match(htm)
        if not m:
            continue
        nth = int(m.group(1))
        group = "first_at_position" if nth == 1 else f"nth_at_position_{nth}"
        for ln in _as_list(g.get("GameLines")):
            name = (ln.get("tmname") or "").strip()
            if not name or name.startswith("*ALL BETS ACTION*"):
                continue
            american = _american(ln.get("odds"))
            if american is None:
                continue
            rows.append(OddsRow(
                book="wagerzon", book_label=htm, book_subject=name,
                american_odds=american, fetched_at=now,
                market_group=group,
            ))
    return rows


def _emit_top_n(lg: dict, now: datetime) -> list[OddsRow]:
    """Shape C: (DRAFTED TOP N) PLAYER; odds in GameLines[0].odds."""
    rows: list[OddsRow] = []
    desc = (lg.get("Description") or "").strip()
    for g in _as_list(lg.get("Games")):
        htm = (g.get("htm") or "").strip()
        m = _TOPN_HTM_RE.search(htm)
        if not m:
            m = _R1_HTM_RE.search(htm)
            if not m:
                continue
            # (DRAFTED IN R1) -> round 1 means top 32
            high = 32
            player = m.group(2).strip()
        else:
            high = int(m.group(1))
            player = m.group(2).strip()
        lines = _as_list(g.get("GameLines"))
        if not lines:
            continue
        american = _american(lines[0].get("odds"))
        if american is None:
            continue
        rows.append(OddsRow(
            book="wagerzon", book_label=desc, book_subject=player,
            american_odds=american, fetched_at=now,
            market_group=f"top_{high}_range",
        ))
    return rows


def _emit_team_position_of_first_pick(lg: dict, now: datetime) -> list[OddsRow]:
    """Shape B: '<TEAM> POSITION OF 1ST DRAFTED PLAYER' - which position
    does team X draft first. Prop type."""
    rows: list[OddsRow] = []
    for g in _as_list(lg.get("Games")):
        htm = (g.get("htm") or "").strip()
        for ln in _as_list(g.get("GameLines")):
            name = (ln.get("tmname") or "").strip()
            if not name or name.startswith("*ALL BETS ACTION*"):
                continue
            american = _american(ln.get("odds"))
            if american is None:
                continue
            rows.append(OddsRow(
                book="wagerzon", book_label=htm, book_subject=name,
                american_odds=american, fetched_at=now,
                market_group="prop_team_position_of_first_pick",
            ))
    return rows


def _emit_overunder_prop(lg: dict, now: datetime) -> list[OddsRow]:
    """Shape C (variant): player draft position / totals markets are
    binary over/under lines. Both sides emit as two rows under the same
    label with '<Over|Under> N.N' subject."""
    rows: list[OddsRow] = []
    for g in _as_list(lg.get("Games")):
        htm = (g.get("htm") or "").strip()
        gdesc = (g.get("gdesc") or "").strip()
        # Use gdesc when present (richer, e.g. 'NFL Draft 2026 - Caleb Downs Draft Position'),
        # else fall back to htm.
        label = gdesc or htm
        for ln in _as_list(g.get("GameLines")):
            ov_odds = _american(ln.get("ovoddst"))
            un_odds = _american(ln.get("unoddst"))
            ov_num = ln.get("ovt")   # e.g. "-9.5" (over number)
            un_num = ln.get("unt")   # e.g. "9.5"
            if ov_odds is not None and un_num:
                rows.append(OddsRow(
                    book="wagerzon", book_label=label,
                    book_subject=f"Over {un_num}",
                    american_odds=ov_odds, fetched_at=now,
                    market_group="prop_overunder",
                ))
            if un_odds is not None and un_num:
                rows.append(OddsRow(
                    book="wagerzon", book_label=label,
                    book_subject=f"Under {un_num}",
                    american_odds=un_odds, fetched_at=now,
                    market_group="prop_overunder",
                ))
    return rows


# ---------------------------------------------------------------------------
# Classifier
# ---------------------------------------------------------------------------

def _dispatch_league(lg: dict, now: datetime) -> list[OddsRow]:
    """Look at a league's Description + Games shape and route to the right
    emitter. Returns [] for buckets we don't handle yet."""
    desc = (lg.get("Description") or "").upper()
    games = _as_list(lg.get("Games"))

    # Empty leagues: skip silently.
    if not games:
        return []

    # 11-15 range league: multiple games, each a different pick number.
    if _PICK_RANGE_DESC_RE.search(desc):
        return _emit_multi_game_pick(lg, now)

    # Single-game pick (Picks 1-2 shape)
    if _PICK_DESC_RE.search(desc) and len(games) == 1 and len(
            _as_list(games[0].get("GameLines"))) > 1:
        return _emit_pick_single_game(lg, games[0], now)

    # Multi-game pick (Picks 3-10 shape: same pick number but per-runner games)
    if _PICK_DESC_RE.search(desc):
        rows: list[OddsRow] = []
        m = _PICK_DESC_RE.search(desc)
        pick = int(m.group(1))
        book_label = (lg.get("Description") or "").strip()
        for g in games:
            for ln in _as_list(g.get("GameLines")):
                name = (ln.get("tmname") or "").strip()
                if not name or name.startswith("*ALL BETS ACTION*"):
                    continue
                american = _american(ln.get("odds"))
                if american is None:
                    continue
                rows.append(OddsRow(
                    book="wagerzon", book_label=book_label,
                    book_subject=name, american_odds=american,
                    fetched_at=now, market_group="pick_outright",
                ))
        return rows

    # Position-selected league
    if "PLAYER SELECTED BY POSITION" in desc:
        return _emit_first_at_position(lg, now)

    # Top-N / Round-1 leagues
    if "TOP 5 PLAYERS DRAFTED" in desc or "TOP 10 PLAYERS DRAFTED" in desc \
            or "DRAFTED IN THE 1ST ROUND" in desc:
        return _emit_top_n(lg, now)

    # Teams' position-of-first-pick
    if "TEAMS POSITION OF 1ST DRAFTED PLAYER" in desc:
        return _emit_team_position_of_first_pick(lg, now)

    # Over/Under style: player draft position, totals
    if "PLAYERS DRAFT POSITION" in desc or "NFL DRAFT TOTALS" in desc:
        return _emit_overunder_prop(lg, now)

    # Fallback: treat as simple binary (each game is one runner).
    rows: list[OddsRow] = []
    for g in games:
        htm = (g.get("htm") or "").strip()
        for ln in _as_list(g.get("GameLines")):
            name = (ln.get("tmname") or "").strip()
            if not name or name.startswith("*ALL BETS ACTION*"):
                continue
            american = _american(ln.get("odds"))
            if american is None:
                continue
            rows.append(OddsRow(
                book="wagerzon", book_label=htm, book_subject=name,
                american_odds=american, fetched_at=now,
                market_group=f"prop_{(desc or 'wz').lower().replace(' ', '_')[:40]}",
            ))
    return rows


# ---------------------------------------------------------------------------
# Public parser
# ---------------------------------------------------------------------------

def parse_response(raw: dict) -> List[OddsRow]:
    """Walk WZ's listLeagues envelope and emit OddsRow per runner.

    Accepts either the full recon envelope or just the inner
    `result`-bearing payload.
    """
    data = raw.get("data") if isinstance(raw, dict) and "data" in raw else raw
    if not isinstance(data, dict):
        return []

    result = data.get("result") or {}
    ll = result.get("listLeagues") or []

    now = datetime.now()
    rows: List[OddsRow] = []
    for bucket in ll:
        # bucket may itself be a list (usual) or a single league dict.
        for lg in _as_list(bucket):
            rows.extend(_dispatch_league(lg, now))
    return rows


# ---------------------------------------------------------------------------
# Live fetch (delegates to recon_wz.py)
# ---------------------------------------------------------------------------

def fetch_raw() -> dict:
    """Refresh WZ via recon script's REST path (authenticated + multi-lg)."""
    _THIS = Path(__file__).resolve()
    sys.path.insert(0, str(_THIS.parent.parent.parent))
    from nfl_draft.scrapers import recon_wz  # type: ignore[import-not-found]

    data, _url, meta = recon_wz.run_rest_phase()
    if data is None:
        raise RuntimeError("WZ REST fetch returned no data; see recon_wz logs.")
    return {"book": "wagerzon", "meta": meta, "data": data}


def fetch_draft_odds() -> List[OddsRow]:
    return parse_response(fetch_raw())


# ---------------------------------------------------------------------------
# MARKET_MAP helpers
# ---------------------------------------------------------------------------

PICK_DESC_RE = _PICK_DESC_RE
PICK_RANGE_DESC_RE = _PICK_RANGE_DESC_RE
NTH_POS_HTM_RE = _NTH_POS_HTM_RE
PICK_GAME_HTM_RE = _PICK_GAME_HTM_RE
