"""Shared SGP scrape orchestration for Kalshi MLB bots.

Enumerates Kalshi MVE markets, writes target lines to the caller's market DB,
spawns scraper subprocesses, and returns priced rows — callers wire results
into their own cache.
"""
from __future__ import annotations
import os
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Mapping

import duckdb

from kalshi_common import auth_client
from kalshi_common.leg_types import _MLB_CODE_TO_TEAM, _parse_event_suffix
from mlb_sgp._shared import TargetLine
from kalshi_common.sgp_service import SGPService  # noqa: F401  (re-export)


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


def _fetch_schedule_from_odds_api() -> dict[tuple, dict]:
    """Fetch today's MLB events from the Odds API.

    Returns dict keyed by (home_team, away_team) tuple →
    {game_id, home_team, away_team, commence_time}.

    Decouples bot from the dashboard's mlb.duckdb (R-locked) for schedule.
    ODDS_API_KEY is read from env first, falling back to ~/.Renviron (R's
    env file) per the run.py pattern (commit 75e02e9). Returns {} if the
    key is missing or the Odds API call fails — caller drops the cycle.
    """
    from datetime import datetime as _dt
    import json
    import urllib.request

    key = os.environ.get("ODDS_API_KEY")
    if not key:
        # Fallback: parse ~/.Renviron (R-side env file). Python doesn't read
        # it automatically the way Rscript does.
        renviron = Path.home() / ".Renviron"
        if renviron.exists():
            for raw_line in renviron.read_text().splitlines():
                line = raw_line.strip()
                if line.startswith("ODDS_API_KEY"):
                    _, _, val = line.partition("=")
                    key = val.strip().strip('"').strip("'")
                    break
    if not key:
        print("  _fetch_schedule: ODDS_API_KEY not set (env or ~/.Renviron)",
              flush=True)
        return {}

    url = f"https://api.the-odds-api.com/v4/sports/baseball_mlb/events?apiKey={key}"
    try:
        with urllib.request.urlopen(url, timeout=10) as resp:
            events = json.loads(resp.read().decode())
    except Exception as e:
        print(f"  _fetch_schedule: Odds API call failed — {e}", flush=True)
        return {}

    out: dict[tuple, dict] = {}
    for e in events:
        game_id = e.get("id")
        home = e.get("home_team")
        away = e.get("away_team")
        ct_str = e.get("commence_time")
        if not (game_id and home and away):
            continue
        normalized = (ct_str.replace("Z", "+00:00")
                      if ct_str and ct_str.endswith("Z") else ct_str)
        try:
            ct = _dt.fromisoformat(normalized) if normalized else None
        except Exception:
            ct = None
        out[(home, away)] = {
            "game_id": game_id, "home_team": home, "away_team": away,
            "commence_time": ct,
        }
    return out


def enumerate_kalshi_targets() -> list[TargetLine]:
    """Enumerate all open Kalshi MVE (spread, total) tuples per MLB game.
    Returns a TargetLine per (game x spread x total) combination, FG only.

    Schedule is fetched directly from the Odds API; we no longer read
    Answer Keys/mlb.duckdb (which is R-write-locked during pipeline runs).
    """
    events = _fetch_kalshi_mlb_events()
    if not events:
        print("  enumerate: 0 kalshi events", flush=True)
        return []
    schedule = _fetch_schedule_from_odds_api()
    print(f"  enumerate: {len(events)} kalshi events, {len(schedule)} schedule rows",
          flush=True)
    targets: list[TargetLine] = []
    matched_games = 0
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
        matched_games += 1
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
    print(f"  enumerate: matched_games={matched_games} → {len(targets)} target lines",
          flush=True)
    return targets


def write_target_lines(target_lines: list[TargetLine], db_path: str):
    """Atomic DELETE+INSERT of mlb_target_lines in bot market DB.
    Creates the table if missing."""
    con = duckdb.connect(db_path)
    try:
        con.execute("""
            CREATE TABLE IF NOT EXISTS mlb_target_lines (
                game_id        VARCHAR,
                home_team      VARCHAR,
                away_team      VARCHAR,
                commence_time  TIMESTAMP,
                period         VARCHAR,
                spread         DOUBLE,
                total          DOUBLE,
                written_at     TIMESTAMP
            )
        """)
        con.execute("BEGIN TRANSACTION")
        con.execute("DELETE FROM mlb_target_lines")
        if target_lines:
            # DuckDB TIMESTAMP columns are naive. Inserting a tz-aware datetime
            # converts it to LOCAL wall-clock and strips the tz (Python duckdb
            # driver behavior). That silently corrupts UTC hour-bucket matching
            # downstream (e.g. legacy match_events in scraper_*_sgp.py uses
            # _utc_bucket which reads .hour off the stored naive value).
            # Convert to UTC and strip tz before insert so the stored value
            # is wall-clock UTC.
            def _to_utc_naive(dt):
                if dt is None:
                    return None
                if dt.tzinfo is None:
                    return dt
                return dt.astimezone(timezone.utc).replace(tzinfo=None)

            now = _to_utc_naive(datetime.now(timezone.utc))
            values = []
            for t in target_lines:
                ct = _to_utc_naive(t.commence_time)
                values.extend([t.game_id, t.home_team, t.away_team,
                                ct, t.period, t.spread, t.total, now])
            placeholders = ",".join(["(?, ?, ?, ?, ?, ?, ?, ?)"] * len(target_lines))
            con.execute(f"INSERT INTO mlb_target_lines VALUES {placeholders}", values)
        con.execute("COMMIT")
    except Exception:
        con.execute("ROLLBACK")
        raise
    finally:
        con.close()


def run_scrapers(
    scraper_dir: str,
    scraper_names: list[str],
    venv_python: str,
    timeout_sec: int,
    env: Mapping[str, str] | None = None,
) -> dict[str, int]:
    """Spawn each scraper as a subprocess in parallel. Returns
    {scraper_name: return_code}.

    - All scrapers share `env` (defaulting to os.environ) and a global
      deadline of `timeout_sec` after which any still-running scraper
      is killed.
    - Hardened against handle leaks: subprocess Popen failures close
      log handles before re-raising.
    """
    if env is None:
        env_dict = dict(os.environ)
    else:
        env_dict = dict(os.environ)
        env_dict.update(env)

    log_handles: list = []
    procs: dict[str, subprocess.Popen] = {}
    try:
        for name in scraper_names:
            log_path = Path(scraper_dir) / f"{name}.runner.log"
            handle = open(log_path, "w")
            log_handles.append(handle)
            try:
                p = subprocess.Popen(
                    [venv_python, name],
                    cwd=scraper_dir,
                    env=env_dict,
                    stdout=handle,
                    stderr=subprocess.STDOUT,
                )
                procs[name] = p
            except Exception:
                handle.close()
                raise
        deadline = time.time() + timeout_sec
        rcs: dict[str, int] = {}
        for name, p in procs.items():
            remaining = deadline - time.time()
            if remaining <= 0:
                p.kill()
                rcs[name] = -1
                continue
            try:
                rcs[name] = p.wait(timeout=remaining)
            except subprocess.TimeoutExpired:
                p.kill()
                try:
                    p.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    pass
                rcs[name] = -1
        return rcs
    finally:
        for h in log_handles:
            try:
                h.close()
            except Exception:
                pass


def read_priced_rows(bot_market_db: str, max_age_sec: int):
    """Read mlb_sgp_odds rows fresher than max_age_sec from the bot DB.
    Returns a pandas DataFrame; empty if no fresh data."""
    import pandas as pd
    if not Path(bot_market_db).exists():
        return pd.DataFrame()
    con = duckdb.connect(bot_market_db, read_only=True)
    try:
        tables = {t[0] for t in con.execute("SHOW TABLES").fetchall()}
        if "mlb_sgp_odds" not in tables:
            return pd.DataFrame()
        df = con.execute(
            "SELECT game_id, combo, period, bookmaker, sgp_decimal, sgp_american, "
            "fetch_time, source, spread_line, total_line "
            "FROM mlb_sgp_odds WHERE fetch_time > NOW() - INTERVAL (CAST(? AS BIGINT)) SECOND",
            [max_age_sec],
        ).fetchdf()
        return df
    finally:
        con.close()


SCRAPER_NAMES = [
    "scraper_draftkings_sgp.py",
    "scraper_fanduel_sgp.py",
    "scraper_prophetx_sgp.py",
    "scraper_novig_sgp.py",
]


def sgp_cycle(
    bot_market_db: str,
    scraper_dir: str,
    venv_python: str,
    timeout_sec: int,
) -> dict[str, int]:
    """One full SGP scrape tick (atomic, serial):
      1. Enumerate Kalshi MVE → list[TargetLine] (schedule from Odds API)
      2. Write to mlb_target_lines in bot_market_db
      3. Spawn the 4 scrapers with MLB_SGP_DB_PATH=bot_market_db, MLB_SGP_PERIODS=FG

    Returns {scraper_name: return_code}.
    """
    targets = enumerate_kalshi_targets()
    write_target_lines(targets, db_path=bot_market_db)
    rcs = run_scrapers(
        scraper_dir=scraper_dir,
        scraper_names=SCRAPER_NAMES,
        venv_python=venv_python,
        timeout_sec=timeout_sec,
        env={
            "MLB_SGP_DB_PATH": bot_market_db,
            "MLB_SGP_PERIODS": "FG",
        },
    )
    return rcs
