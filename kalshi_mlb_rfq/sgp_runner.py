"""SGP scrape orchestration for the Kalshi MLB RFQ bot.

Owns the cycle that:
  1. Enumerates Kalshi MVE markets per open MLB game
  2. Writes target lines (one row per game x (spread, total)) to bot DB
  3. Spawns the 4 scraper subprocesses with MLB_SGP_DB_PATH redirect
  4. Reads back priced SGP odds into the bot's _SGP_ODDS_CACHE

This module is invoked from main_loop on the SGP cadence tick.
"""
from __future__ import annotations
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import duckdb

from kalshi_mlb_rfq import auth_client
from mlb_sgp._shared import TargetLine


def should_scrape(last_fetch_time: datetime | None,
                   now: datetime,
                   min_interval_sec: int) -> bool:
    """True if we should scrape this tick. Guards against tight cadences
    after crash-recovery restarts that hit an already-fresh DB.

    Both `last_fetch_time` and `now` are normalized to UTC when naive."""
    if last_fetch_time is None:
        return True
    if last_fetch_time.tzinfo is None:
        last_fetch_time = last_fetch_time.replace(tzinfo=timezone.utc)
    if now.tzinfo is None:
        now = now.replace(tzinfo=timezone.utc)
    age = (now - last_fetch_time).total_seconds()
    return age > min_interval_sec


def latest_sgp_fetch_time(bot_market_db: str) -> datetime | None:
    """Read MAX(fetch_time) from mlb_sgp_odds in bot_market_db.
    Returns None for missing DB / missing table / empty table."""
    if not Path(bot_market_db).exists():
        return None
    con = duckdb.connect(bot_market_db, read_only=True)
    try:
        tables = {t[0] for t in con.execute("SHOW TABLES").fetchall()}
        if "mlb_sgp_odds" not in tables:
            return None
        row = con.execute("SELECT MAX(fetch_time) FROM mlb_sgp_odds").fetchone()
        return row[0] if row and row[0] is not None else None
    finally:
        con.close()


# 3-letter Kalshi team code -> Odds-API canonical team name.
# Same mapping main.py uses (lifted to module level for testability).
_MLB_CODE_TO_TEAM = {
    "ARI": "Arizona Diamondbacks", "ATL": "Atlanta Braves", "BAL": "Baltimore Orioles",
    "BOS": "Boston Red Sox", "CHC": "Chicago Cubs", "CWS": "Chicago White Sox",
    "CIN": "Cincinnati Reds", "CLE": "Cleveland Guardians", "COL": "Colorado Rockies",
    "DET": "Detroit Tigers", "HOU": "Houston Astros", "KC": "Kansas City Royals",
    "LAA": "Los Angeles Angels", "LAD": "Los Angeles Dodgers", "MIA": "Miami Marlins",
    "MIL": "Milwaukee Brewers", "MIN": "Minnesota Twins", "NYM": "New York Mets",
    "NYY": "New York Yankees", "OAK": "Athletics", "ATH": "Athletics",
    "AZ": "Arizona Diamondbacks", "PHI": "Philadelphia Phillies",
    "PIT": "Pittsburgh Pirates", "SD": "San Diego Padres", "SF": "San Francisco Giants",
    "SEA": "Seattle Mariners", "STL": "St. Louis Cardinals", "TB": "Tampa Bay Rays",
    "TEX": "Texas Rangers", "TOR": "Toronto Blue Jays",
    "WAS": "Washington Nationals", "WSH": "Washington Nationals",
}


def _parse_event_suffix(suffix: str) -> tuple[str | None, str | None]:
    """Split a KXMLB* event suffix into (away_code, home_code).
    Matches main.py::_parse_event_suffix exactly."""
    if len(suffix) < 11 + 4:
        return None, None
    team_block = suffix[11:]
    for home_len in (3, 2):
        if len(team_block) <= home_len:
            continue
        home = team_block[-home_len:]
        away = team_block[:-home_len]
        if home in _MLB_CODE_TO_TEAM and away in _MLB_CODE_TO_TEAM:
            return away, home
    return None, None


def _fetch_kalshi_mlb_events() -> list[dict]:
    status, body, _ = auth_client.api(
        "GET", "/events?series_ticker=KXMLBGAME&status=open&limit=50")
    if status != 200 or not isinstance(body, dict):
        return []
    return body.get("events", [])


def _fetch_kalshi_spread_lines(suffix: str) -> list[tuple[float, str]]:
    status, body, _ = auth_client.api(
        "GET", f"/markets?event_ticker=KXMLBSPREAD-{suffix}&limit=50")
    if status != 200 or not isinstance(body, dict):
        return []
    out = []
    seen = set()
    for m in body.get("markets", []):
        ticker = m.get("ticker", "")
        prefix = f"KXMLBSPREAD-{suffix}-"
        if not ticker.startswith(prefix):
            continue
        spread_part = ticker[len(prefix):]
        digits = "".join(c for c in spread_part if c.isdigit())
        team_chars = "".join(c for c in spread_part if not c.isdigit())
        if not digits or not team_chars:
            continue
        n = int(digits)
        line = -(n - 0.5)
        key = round(line, 1)
        if key in seen:
            continue
        seen.add(key)
        out.append((line, "home"))
    return out


def _fetch_kalshi_total_lines(suffix: str) -> list[float]:
    status, body, _ = auth_client.api(
        "GET", f"/markets?event_ticker=KXMLBTOTAL-{suffix}&limit=50")
    if status != 200 or not isinstance(body, dict):
        return []
    out = []
    seen = set()
    for m in body.get("markets", []):
        ticker = m.get("ticker", "")
        try:
            n = int(ticker.rsplit("-", 1)[-1])
            line = n - 0.5
            key = round(line, 1)
            if key in seen:
                continue
            seen.add(key)
            out.append(line)
        except ValueError:
            continue
    return out


def _load_schedule(schedule_db_path: str) -> dict[tuple, dict]:
    """Read mlb_odds_temp from mlb.duckdb. Returns dict keyed by (home, away) tuple
    -> {game_id, home_team, away_team, commence_time}."""
    from datetime import datetime as _dt
    if not Path(schedule_db_path).exists():
        return {}
    con = duckdb.connect(schedule_db_path, read_only=True)
    try:
        tables = {t[0] for t in con.execute("SHOW TABLES").fetchall()}
        if "mlb_odds_temp" not in tables:
            return {}
        rows = con.execute(
            "SELECT id, home_team, away_team, commence_time FROM mlb_odds_temp"
        ).fetchall()
    finally:
        con.close()
    out = {}
    for game_id, home, away, ct_str in rows:
        normalized = ct_str.replace("Z", "+00:00") if ct_str and ct_str.endswith("Z") else ct_str
        try:
            ct = _dt.fromisoformat(normalized) if normalized else None
        except Exception:
            ct = None
        out[(home, away)] = {"game_id": game_id, "home_team": home,
                              "away_team": away, "commence_time": ct}
    return out


def enumerate_kalshi_targets(schedule_db_path: str) -> list[TargetLine]:
    """Enumerate all open Kalshi MVE (spread, total) tuples per MLB game.
    Returns a TargetLine per (game x spread x total) combination, FG only."""
    events = _fetch_kalshi_mlb_events()
    if not events:
        return []
    schedule = _load_schedule(schedule_db_path)
    targets: list[TargetLine] = []
    for ev in events:
        event_ticker = ev.get("event_ticker", "")
        if not event_ticker.startswith("KXMLBGAME-"):
            continue
        suffix = event_ticker.replace("KXMLBGAME-", "")
        away_code, home_code = _parse_event_suffix(suffix)
        if away_code is None or home_code is None:
            continue
        home_team = _MLB_CODE_TO_TEAM.get(home_code)
        away_team = _MLB_CODE_TO_TEAM.get(away_code)
        sched = schedule.get((home_team, away_team))
        if not sched:
            continue
        spreads = _fetch_kalshi_spread_lines(suffix)
        totals = _fetch_kalshi_total_lines(suffix)
        if not spreads or not totals:
            continue
        for spread, _who in spreads:
            for total in totals:
                targets.append(TargetLine(
                    game_id=sched["game_id"],
                    home_team=home_team, away_team=away_team,
                    commence_time=sched["commence_time"],
                    period="FG", spread=spread, total=total,
                ))
    return targets
