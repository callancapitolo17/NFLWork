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


def _odds_api_key() -> str | None:
    """ODDS_API_KEY from env, falling back to ~/.Renviron (R's env file) per the
    run.py pattern (commit 75e02e9). Python doesn't read ~/.Renviron the way
    Rscript does, so we parse it ourselves."""
    key = os.environ.get("ODDS_API_KEY")
    if key:
        return key
    renviron = Path.home() / ".Renviron"
    if renviron.exists():
        for raw_line in renviron.read_text().splitlines():
            line = raw_line.strip()
            if line.startswith("ODDS_API_KEY"):
                _, _, val = line.partition("=")
                return val.strip().strip('"').strip("'")
    return None


def _fetch_schedule_from_odds_api() -> dict[tuple, dict]:
    """Fetch today's MLB events from the Odds API.

    Returns dict keyed by (home_team, away_team) tuple →
    {game_id, home_team, away_team, commence_time}.

    Decouples bot from the dashboard's mlb.duckdb (R-locked) for schedule.
    Returns {} if the key is missing or the Odds API call fails — caller
    drops the cycle.
    """
    from datetime import datetime as _dt
    import json
    import urllib.request

    key = _odds_api_key()
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


def _parse_singles_payload(events: list[dict], fetch_time: datetime) -> list[dict]:
    """Parse an Odds API /odds payload (markets=h2h,spreads,totals) into flat
    mlb_singles_odds rows. Pure function — no network — so it is unit-testable.

    Each row: {game_id, home_team, away_team, market, line, outcome, decimal,
    bookmaker, fetch_time, source}. market is 'moneyline'|'spread'|'total';
    outcome is 'home'|'away' (ml/spread) or 'over'|'under' (total). `line` is
    None for moneyline and the HOME-perspective points for spread/total (both
    sides of a market store the same line so the pricer can pair them).

    A market is emitted only if BOTH of its two sides are present and resolvable
    — an incomplete pair can't be devigged, so it is dropped silently.
    """
    rows: list[dict] = []
    for ev in events:
        gid = ev.get("id")
        home = ev.get("home_team")
        away = ev.get("away_team")
        if not (gid and home and away):
            continue
        for bk in ev.get("bookmakers", []):
            book = bk.get("key")
            if not book:
                continue
            for m in bk.get("markets", []):
                mkey = m.get("key")
                outs = m.get("outcomes", []) or []
                if mkey == "h2h":
                    by_team = {o.get("name"): o for o in outs}
                    if home not in by_team or away not in by_team:
                        continue
                    for outcome, team in (("home", home), ("away", away)):
                        rows.append(dict(
                            game_id=gid, home_team=home, away_team=away,
                            market="moneyline", line=None, outcome=outcome,
                            decimal=float(by_team[team]["price"]),
                            bookmaker=book, fetch_time=fetch_time, source="odds_api"))
                elif mkey == "spreads":
                    by_team = {o.get("name"): o for o in outs}
                    if home not in by_team or away not in by_team:
                        continue
                    home_point = by_team[home].get("point")
                    if home_point is None:
                        continue
                    line = float(home_point)  # home-perspective handicap
                    for outcome, team in (("home", home), ("away", away)):
                        rows.append(dict(
                            game_id=gid, home_team=home, away_team=away,
                            market="spread", line=line, outcome=outcome,
                            decimal=float(by_team[team]["price"]),
                            bookmaker=book, fetch_time=fetch_time, source="odds_api"))
                elif mkey == "totals":
                    by_side = {str(o.get("name", "")).lower(): o for o in outs}
                    if "over" not in by_side or "under" not in by_side:
                        continue
                    point = by_side["over"].get("point")
                    if point is None:
                        continue
                    line = float(point)
                    for outcome in ("over", "under"):
                        rows.append(dict(
                            game_id=gid, home_team=home, away_team=away,
                            market="total", line=line, outcome=outcome,
                            decimal=float(by_side[outcome]["price"]),
                            bookmaker=book, fetch_time=fetch_time, source="odds_api"))
    return rows


def fetch_mlb_singles_from_odds_api() -> list[dict]:
    """Fetch per-book moneyline / spread / total singles for all MLB games from
    the Odds API in one request. Returns flat rows (see _parse_singles_payload).
    Returns [] if the key is missing or the call fails — callers keep prior rows.
    """
    import json
    import urllib.request

    key = _odds_api_key()
    if not key:
        return []
    url = (f"https://api.the-odds-api.com/v4/sports/baseball_mlb/odds?apiKey={key}"
           "&regions=us&markets=h2h,spreads,totals&oddsFormat=decimal")
    try:
        with urllib.request.urlopen(url, timeout=15) as resp:
            events = json.loads(resp.read().decode())
    except Exception as e:
        print(f"  fetch_singles: Odds API call failed — {e}", flush=True)
        return []
    return _parse_singles_payload(events, datetime.now(timezone.utc))


def write_singles_odds(rows: list[dict], db_path: str):
    """Atomic DELETE+INSERT of mlb_singles_odds in the bot market DB. Creates
    the table if missing. No-op DELETE-only when rows is empty would wipe a good
    prior snapshot, so an empty `rows` is treated as 'keep prior' (early return).
    """
    con = duckdb.connect(db_path)
    try:
        con.execute("""
            CREATE TABLE IF NOT EXISTS mlb_singles_odds (
                game_id     VARCHAR,
                home_team   VARCHAR,
                away_team   VARCHAR,
                market      VARCHAR,
                line        DOUBLE,
                outcome     VARCHAR,
                decimal     DOUBLE,
                bookmaker   VARCHAR,
                fetch_time  TIMESTAMP,
                source      VARCHAR
            )
        """)
        if not rows:
            return
        con.execute("BEGIN TRANSACTION")
        con.execute("DELETE FROM mlb_singles_odds")
        values = []
        for r in rows:
            ft = r["fetch_time"]
            if ft is not None and ft.tzinfo is not None:
                ft = ft.astimezone(timezone.utc).replace(tzinfo=None)
            values.extend([r["game_id"], r["home_team"], r["away_team"],
                           r["market"], r["line"], r["outcome"], r["decimal"],
                           r["bookmaker"], ft, r["source"]])
        placeholders = ",".join(["(?,?,?,?,?,?,?,?,?,?)"] * len(rows))
        con.execute(f"INSERT INTO mlb_singles_odds VALUES {placeholders}", values)
        con.execute("COMMIT")
    except Exception:
        try:
            con.execute("ROLLBACK")
        except Exception:
            pass
        raise
    finally:
        con.close()


# Odds API player-prop market key -> our normalized market_type. The Kalshi
# leg prefixes map to these: KXMLBHR->home_runs, KXMLBHIT->hits,
# KXMLBTB->total_bases, KXMLBKS->strikeouts, KXMLBHRR->hits_runs_rbis.
_PROP_MARKET_MAP = {
    "batter_home_runs": "home_runs",
    "batter_hits": "hits",
    "batter_total_bases": "total_bases",
    "pitcher_strikeouts": "strikeouts",
    "batter_hits_runs_rbis": "hits_runs_rbis",
}
_PLAYER_SUFFIXES = {"jr", "sr", "ii", "iii", "iv", "v"}


def _norm_player(name: str) -> str:
    """Normalize a player name so Kalshi market titles and Odds API descriptions
    map to the same key: strip accents, lowercase, drop punctuation and generation
    suffixes (Jr/Sr/II...), collapse whitespace. Kalshi drops accented letters in
    ticker codes inconsistently, so we match on the clean title name, not the code.
    """
    import unicodedata
    decomposed = unicodedata.normalize("NFKD", name)
    ascii_only = "".join(c for c in decomposed if not unicodedata.combining(c))
    # Drop punctuation outright (so "J.D." -> "jd", "O'Neill" -> "oneill") while
    # keeping spaces as token boundaries. Both Kalshi titles and Odds API
    # descriptions go through this same function, so the key matches as long as
    # the rule is consistent.
    cleaned = "".join(c for c in ascii_only if c.isalnum() or c.isspace())
    tokens = [t for t in cleaned.lower().split() if t not in _PLAYER_SUFFIXES]
    return " ".join(tokens)


def _parse_props_event(ev: dict, fetch_time: datetime) -> list[dict]:
    """Parse one Odds API event-odds payload (player-prop markets) into flat
    mlb_player_props rows. Pure function. Each row: {game_id, market_type, player,
    line, outcome ('over'/'under'), decimal, bookmaker, fetch_time}. Only complete
    over+under pairs (per player/line/book) are emitted — an unpaired side can't
    be devigged. `player` is normalized via _norm_player."""
    gid = ev.get("id")
    if not gid:
        return []
    rows: list[dict] = []
    for bk in ev.get("bookmakers", []):
        book = bk.get("key")
        if not book:
            continue
        for m in bk.get("markets", []):
            mtype = _PROP_MARKET_MAP.get(m.get("key"))
            if not mtype:
                continue
            # group outcomes by (player, line) so we pair over/under
            by_key: dict[tuple, dict] = {}
            for o in m.get("outcomes", []) or []:
                player = o.get("description")
                point = o.get("point")
                side = str(o.get("name", "")).lower()
                if not player or point is None or side not in ("over", "under"):
                    continue
                by_key.setdefault((player, float(point)), {})[side] = o
            for (player, line), sides in by_key.items():
                if "over" not in sides or "under" not in sides:
                    continue
                pn = _norm_player(player)
                for outcome in ("over", "under"):
                    rows.append(dict(
                        game_id=gid, market_type=mtype, player=pn, line=line,
                        outcome=outcome, decimal=float(sides[outcome]["price"]),
                        bookmaker=book, fetch_time=fetch_time))
    return rows


def fetch_mlb_player_props_from_odds_api(game_ids: list[str]) -> list[dict]:
    """Fetch player props for the given Odds API event ids (one request per
    event — the props endpoint is per-event). Returns flat rows. Best-effort: a
    failed event is skipped; a missing key returns []."""
    import json
    import urllib.request

    key = _odds_api_key()
    if not key or not game_ids:
        return []
    markets = ",".join(_PROP_MARKET_MAP.keys())
    now = datetime.now(timezone.utc)
    out: list[dict] = []
    # The props endpoint is per-event, so this is one call per game and runs
    # synchronously inside the cycle. Use a short per-call timeout so a single
    # hung game can't stall the loop (and therefore quoting) for long.
    for gid in game_ids:
        url = (f"https://api.the-odds-api.com/v4/sports/baseball_mlb/events/{gid}"
               f"/odds?apiKey={key}&regions=us&markets={markets}&oddsFormat=decimal")
        try:
            with urllib.request.urlopen(url, timeout=6) as resp:
                ev = json.loads(resp.read().decode())
        except Exception:
            continue
        out.extend(_parse_props_event(ev, now))
    return out


def write_player_props(rows: list[dict], db_path: str):
    """Atomic DELETE+INSERT of mlb_player_props. Empty rows -> keep prior (early
    return), same as write_singles_odds."""
    con = duckdb.connect(db_path)
    try:
        con.execute("""
            CREATE TABLE IF NOT EXISTS mlb_player_props (
                game_id     VARCHAR,
                market_type VARCHAR,
                player      VARCHAR,
                line        DOUBLE,
                outcome     VARCHAR,
                decimal     DOUBLE,
                bookmaker   VARCHAR,
                fetch_time  TIMESTAMP
            )
        """)
        if not rows:
            return
        con.execute("BEGIN TRANSACTION")
        con.execute("DELETE FROM mlb_player_props")
        values = []
        for r in rows:
            ft = r["fetch_time"]
            if ft is not None and ft.tzinfo is not None:
                ft = ft.astimezone(timezone.utc).replace(tzinfo=None)
            values.extend([r["game_id"], r["market_type"], r["player"], r["line"],
                           r["outcome"], r["decimal"], r["bookmaker"], ft])
        placeholders = ",".join(["(?,?,?,?,?,?,?,?)"] * len(rows))
        con.execute(f"INSERT INTO mlb_player_props VALUES {placeholders}", values)
        con.execute("COMMIT")
    except Exception:
        try:
            con.execute("ROLLBACK")
        except Exception:
            pass
        raise
    finally:
        con.close()


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
    "scraper_betmgm_sgp.py",
    "scraper_caesars_sgp.py",
]


# Book name -> orchestrator module (source-label owner). Used by the
# service path to clear exactly the labels each book writes.
_BOOK_MODULES = {
    "draftkings": "mlb_sgp.draftkings",
    "fanduel": "mlb_sgp.fanduel",
    "prophetx": "mlb_sgp.prophetx",
    "novig": "mlb_sgp.novig",
}


def sgp_cycle(
    bot_market_db: str,
    scraper_dir: str | None = None,
    venv_python: str | None = None,
    timeout_sec: int | None = None,
    service=None,
) -> dict[str, int]:
    """One full SGP scrape tick.

    With `service` (an SGPService): in-process path —
      1. Enumerate Kalshi MVE -> list[TargetLine]
      2. Write mlb_target_lines (tipoff gating + debugging read it)
      3. service.refresh(targets); per successful book, clear that
         book's source labels and upsert its fresh rows. Failed books
         (None) keep their old rows — same outcome as a crashed
         subprocess today. Returns {book: row_count, failed: -1}.

    Without `service`: legacy subprocess path, unchanged. Returns
    {scraper_name: return_code}. Kept as the rollback hatch
    (scraper_dir / venv_python / timeout_sec are required then).
    """
    targets = enumerate_kalshi_targets()
    write_target_lines(targets, db_path=bot_market_db)

    # Singles odds (moneyline / spread / total) from the Odds API — one request,
    # all books. Powers the general cross-game pricer (moneyline parlays etc.).
    # Best-effort: a failed fetch keeps the prior snapshot, never breaks the cycle.
    try:
        singles = fetch_mlb_singles_from_odds_api()
        write_singles_odds(singles, db_path=bot_market_db)
        if singles:
            print(f"  singles: {len(singles)} odds-api rows "
                  f"({len({r['game_id'] for r in singles})} games)", flush=True)
    except Exception as e:
        singles = []
        print(f"  singles: write failed — {e}", flush=True)

    # Player props (HR / hits / total bases / strikeouts / hits+runs+RBIs) —
    # one Odds API request PER game (the props endpoint is per-event). Reuse the
    # game ids the singles fetch surfaced. Best-effort, never breaks the cycle.
    try:
        prop_game_ids = sorted({r["game_id"] for r in singles}) if singles else []
        props = fetch_mlb_player_props_from_odds_api(prop_game_ids)
        write_player_props(props, db_path=bot_market_db)
        if props:
            print(f"  props: {len(props)} odds-api rows "
                  f"({len({(r['game_id'], r['player']) for r in props})} player-games)",
                  flush=True)
    except Exception as e:
        print(f"  props: write failed — {e}", flush=True)

    if service is None:
        if not (scraper_dir and venv_python and timeout_sec):
            raise ValueError(
                "sgp_cycle: legacy subprocess path needs scraper_dir, "
                "venv_python and timeout_sec")
        return run_scrapers(
            scraper_dir=scraper_dir,
            scraper_names=SCRAPER_NAMES,
            venv_python=venv_python,
            timeout_sec=timeout_sec,
            env={
                "MLB_SGP_DB_PATH": bot_market_db,
                "MLB_SGP_PERIODS": "FG",
            },
        )

    import importlib
    from mlb_sgp import db as sgp_db

    results = service.refresh(targets)
    counts: dict[str, int] = {}
    for book, rows in results.items():
        if rows is None:
            counts[book] = -1
            continue
        mod = importlib.import_module(_BOOK_MODULES[book])
        sgp_db.clear_source(mod.SOURCE_LABEL, db_path=bot_market_db)
        # Keep the getattr guard: DK/FD/NV emit interpolated rows under a
        # SOURCE_LABEL_FALLBACK that must also be cleared so stale ones
        # don't linger; PX defines the label but never emits them. Do not
        # hardcode per-book — getattr handles "has a fallback or not".
        fallback_label = getattr(mod, "SOURCE_LABEL_FALLBACK", None)
        if fallback_label:
            sgp_db.clear_source(fallback_label, db_path=bot_market_db)
        sgp_db.upsert_priced_rows(rows, db_path=bot_market_db)
        counts[book] = len(rows)
    return counts
