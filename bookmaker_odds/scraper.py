#!/usr/bin/env python3
"""
Bookmaker.eu Odds Scraper
Fetches odds from Bookmaker.eu's internal API via curl_cffi (Chrome TLS fingerprint).
No browser needed — uses saved cookies from recon_bookmaker.py.

Usage:
    python scraper.py cbb
    python scraper.py nba

If cookies expire, run recon_bookmaker.py to refresh them.
"""

import json
import os
import sys
import time
import duckdb
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional
from zoneinfo import ZoneInfo

from curl_cffi import requests as cffi_requests
from dotenv import load_dotenv

# Add Answer Keys to path for shared team name resolution
sys.path.insert(0, str(Path(__file__).parent.parent / "Answer Keys"))
from canonical_match import load_team_dict, load_canonical_games, resolve_team_names

# Bookmaker's API returns naive Pacific wall-clock in (gmdt, gmtm) — confirmed
# by tools/TZ_AUDIT_FINDINGS.md (modal Δ vs sharp UTC was +7.02h = PDT). Prior
# R-side code in Tools.R::.parse_wz_game_dt interpreted these strings as ET,
# silently shifting every BKM game by 3 hours downstream and dropping ~half of
# BKM's MLB rows via .drop_past_games(). We now parse as America/Los_Angeles
# and store the UTC instant as game_start_time TIMESTAMPTZ so downstream code
# never has to guess the TZ. Same fix shape as Hoop88 (commit 5d10d24).
_BKM_TZ = ZoneInfo("America/Los_Angeles")


def _bkm_game_start_time(gmdt: str, gmtm: str) -> Optional[datetime]:
    """Parse BKM's naive PT (gmdt, gmtm) pair into a UTC-aware datetime.

    gmdt: 8-digit YYYYMMDD string (e.g. "20260522")
    gmtm: HH:MM or HH:MM:SS string (e.g. "18:30:00")
    Returns: tz-aware UTC datetime, or None if either input is malformed.
    """
    if not gmdt or len(gmdt) != 8:
        return None
    if not gmtm or len(gmtm) < 5:
        return None
    try:
        y, m, d = int(gmdt[0:4]), int(gmdt[4:6]), int(gmdt[6:8])
        hh, mm = int(gmtm[0:2]), int(gmtm[3:5])
    except (ValueError, IndexError):
        return None
    naive_pt = datetime(y, m, d, hh, mm, tzinfo=_BKM_TZ)
    return naive_pt.astimezone(timezone.utc)


# Load .env from bet_logger directory
env_path = Path(__file__).parent.parent / "bet_logger" / ".env"
load_dotenv(env_path)

BOOKMAKER_USERNAME = os.getenv("BOOKMAKER_USERNAME")
BOOKMAKER_PASSWORD = os.getenv("BOOKMAKER_PASSWORD")

SITE_URL = "https://be.bookmaker.eu/en/sports/"
API_BASE = "https://be.bookmaker.eu/gateway/BetslipProxy.aspx"
COOKIE_PATH = Path(__file__).parent / ".bookmaker_cookies.json"
DB_PATH = Path(__file__).parent / "bookmaker.duckdb"

# Sport configurations — markets are described by their BKM leagueDescEn label
# (the human-readable name visible on bookmaker.eu) plus the sportId code BKM
# uses for that sport. League IDs are NOT hardcoded: they are resolved from
# BKM's GetRoutingInfo endpoint at scrape time. This prevents the failure mode
# where someone guesses a leagueId and the scraper silently fetches the wrong
# market — that bug bit us once already (league 503 was labeled "1st 3 Innings"
# but is actually BKM's "2ND HALVES" market; commit aecc669).
#
# To add a new BKM market: pick its leagueDescEn from a fresh GetRoutingInfo
# response (or recon_bookmaker_api.json) and add an entry here. If BKM later
# re-labels the league, scraping returns 0 games for that entry with a clear
# warning — a visible break, not silent corruption.
SPORT_CONFIGS = {
    "cbb": {
        "sport_key": "basketball_ncaab",
        "table_name": "cbb_odds",
        # BKM groups leagues by sportId (broad category) and region (specific
        # league within the category). E.g., BASKETBALL/sportId=NBA has both
        # region=NBA (the real NBA) and region=WNBA — so sportId alone is NOT
        # enough to disambiguate. Always pin both.
        "bkm_sport_id": "NCB",  # UNVERIFIED: CBB was out of season at recon
                                # capture time, so this sportId+region are
                                # inferred. First in-season run WARN-and-skips
                                # if wrong; verify against routing then.
        "bkm_region": "NCAAB",
        "markets": [
            {"league_pattern": "GAME LINES",  "market": "spreads",     "period": "fg"},
            {"league_pattern": "1ST HALVES",  "market": "spreads_h1",  "period": "Half1"},
            # "Extra Games" / non-conference / tournament games may live
            # under a separate region — verify the label in season
            # before adding here.
        ],
    },
    "nba": {
        "sport_key": "basketball_nba",
        "table_name": "nba_odds",
        "bkm_sport_id": "NBA",
        "bkm_region": "NBA",  # NOT "WNBA" — same sportId, different region.
        "markets": [
            {"league_pattern": "GAME LINES",  "market": "spreads",     "period": "fg"},
            {"league_pattern": "1ST HALVES",  "market": "spreads_h1",  "period": "Half1"},
        ],
    },
    "mlb": {
        "sport_key": "baseball_mlb",
        "table_name": "mlb_odds",
        "bkm_sport_id": "MLB",
        "bkm_region": "MLB",
        "markets": [
            {"league_pattern": "GAME LINES",     "market": "spreads",     "period": "fg"},
            {"league_pattern": "1ST 5 INNINGS",  "market": "spreads_f5",  "period": "F5"},
            {"league_pattern": "2ND HALVES",     "market": "spreads_h2",  "period": "H2"},
        ],
    },
}


# =============================================================================
# API CLIENT (curl_cffi)
# =============================================================================


def _make_schedule_body(league_id: str) -> dict:
    return {
        "o": {
            "BORequestData": {
                "BOParameters": {
                    "BORt": {},
                    "LeaguesIdList": league_id,
                    "LanguageId": "0",
                    "LineStyle": "E",
                    "ScheduleType": "american",
                    "LinkDeriv": "true",
                }
            }
        }
    }


def _save_cookies(session):
    """Save session cookies to disk for reuse."""
    cookies = dict(session.cookies)
    with open(COOKIE_PATH, "w") as f:
        json.dump(cookies, f)


def _load_cookies(session):
    """Load saved cookies into session."""
    if not COOKIE_PATH.exists():
        return False
    try:
        with open(COOKIE_PATH) as f:
            cookies = json.load(f)
        for name, value in cookies.items():
            session.cookies.set(name, value, domain=".bookmaker.eu")
        return bool(cookies)
    except (json.JSONDecodeError, OSError):
        return False


def _create_session() -> cffi_requests.Session:
    """Create a curl_cffi session that impersonates Chrome."""
    session = cffi_requests.Session(impersonate="chrome")
    _load_cookies(session)
    return session


def fetch_schedule(session, league_id: str) -> dict | None:
    """Fetch GetSchedule via curl_cffi. Returns None if blocked."""
    try:
        resp = session.post(
            f"{API_BASE}/GetSchedule",
            json=_make_schedule_body(league_id),
            timeout=15,
        )
        if resp.status_code != 200:
            return None
        data = resp.json()
        if "Schedule" not in data:
            return None
        return data
    except Exception:
        return None


def fetch_routing_info(session) -> dict | None:
    """Fetch BKM's league catalog (GetRoutingInfo).

    The response shape is roughly:
        {
          "valid": true,
          "routedSports": [
            {"sportDescEn": "BASEBALL", "routedLeagues": [
                {"sportId": "MLB", "leagueDescEn": "GAME LINES",    "leagueId": "5"},
                {"sportId": "MLB", "leagueDescEn": "1ST 5 INNINGS", "leagueId": "6"},
                {"sportId": "MLB", "leagueDescEn": "2ND HALVES",    "leagueId": "503"},
                ...
            ]},
            ...
          ]
        }

    Returns the parsed JSON, or None on HTTP error / blocked request.
    """
    body = {"o": {"BORequestData": {"BOParameters": {"BORt": {}, "LanguageId": "0"}}}}
    try:
        resp = session.post(f"{API_BASE}/GetRoutingInfo", json=body, timeout=15)
        if resp.status_code != 200:
            return None
        data = resp.json()
        if not data.get("valid") or not data.get("routedSports"):
            return None
        return data
    except Exception:
        return None


def resolve_leagues(routing_info: dict | None, sport_id: str, region: str,
                    markets: list[dict]) -> list[dict]:
    """Resolve each wanted market against BKM's live league catalog.

    For each entry in `markets`, find the routedLeague where
    `(sportId, region, leagueDescEn) == (sport_id, region, market["league_pattern"])`
    and attach its `leagueId`. Markets that don't match anything are dropped
    with a warning — better to see "no such market" loudly than to silently
    fetch the wrong data.

    Why all three keys: BKM's catalog has e.g. region=NBA + region=WNBA both
    sitting under sportId=NBA, so sportId alone is ambiguous. The triple
    (sportId, region, leagueDescEn) is the actual primary key.

    Returns a list of `parse_schedule`-ready league configs (same shape as
    the old hardcoded SPORT_CONFIGS[sport]["leagues"]: id/name/market/period).
    """
    if not routing_info:
        print("WARNING: BKM GetRoutingInfo unavailable — cannot resolve leagues.")
        return []

    # Flatten all routedLeagues across sport buckets and index by
    # (sportId, region, leagueDescEn). sportId + region live on the inner
    # league, not the outer sport entry.
    catalog: dict[tuple[str, str, str], dict] = {}
    # Also track which leagueDescEn values are available per (sportId, region)
    # scope, so when a pattern miss happens we can tell the operator what they
    # COULD have asked for. Cheap and turns "WARN" into actionable feedback.
    patterns_by_scope: dict[tuple[str, str], list[str]] = {}
    for outer in routing_info.get("routedSports", []):
        for lg in (outer.get("routedLeagues") or []):
            # Skip malformed entries: a league with no leagueId is unusable
            # downstream (fetch_schedule would fail) and indexing it would
            # let us match against a phantom entry. Better to ignore it.
            if not lg.get("leagueId"):
                continue
            sid = lg.get("sportId") or ""
            rgn = lg.get("region") or ""
            ldEn = lg.get("leagueDescEn") or ""
            catalog[(sid, rgn, ldEn)] = lg
            patterns_by_scope.setdefault((sid, rgn), []).append(ldEn)

    resolved: list[dict] = []
    for mkt in markets:
        key = (sport_id, region, mkt["league_pattern"])
        lg = catalog.get(key)
        if lg is None:
            available = patterns_by_scope.get((sport_id, region), [])
            if available:
                hint = (f"available patterns under "
                        f"(sportId='{sport_id}', region='{region}'): {available}")
            else:
                hint = (f"no BKM leagues under "
                        f"(sportId='{sport_id}', region='{region}') at all "
                        f"— sport_id/region likely wrong or off-season")
            print(f"  WARNING: no BKM league matching "
                  f"leagueDescEn='{mkt['league_pattern']}' — skipping. {hint}")
            continue
        resolved.append({
            "id": str(lg["leagueId"]),
            "name": mkt["league_pattern"],
            "market": mkt["market"],
            "period": mkt["period"],
        })
    return resolved


def login(session, username: str, password: str) -> bool:
    """Login via curl_cffi. Returns True on success."""
    body = {
        "o": {
            "BORequestData": {
                "BOParameters": {
                    "BORt": {},
                    "Player": username,
                    "Password": password,
                    "loginKey": "",
                }
            }
        }
    }
    try:
        resp = session.post(f"{API_BASE}/Login", json=body, timeout=15)
        return resp.status_code == 200
    except Exception:
        return False


def _has_games(data: dict | None) -> bool:
    """Check if GetSchedule response contains games."""
    if not data:
        return False
    return bool(
        data.get("Schedule", {})
        .get("Data", {})
        .get("Leagues", {})
        .get("League", [])
    )


def _session_looks_healthy(*, blocked: bool, login_ok: bool) -> bool:
    """Did we pass Cloudflare AND authenticate successfully?

    If both are true, an empty schedule means no games posted right now —
    not a broken session. Caller should save empty and exit cleanly rather
    than escalating to the interactive recon browser.
    """
    return (not blocked) and login_ok


def _can_launch_interactive_recon() -> bool:
    """Only allow the Playwright browser popup when a human is at the keyboard.

    recon_bookmaker.py has three blocking input() calls. When the scraper is
    running as a piped subprocess of run.py, stdin is not a TTY and those
    prompts would hang the whole MLB/CBB pipeline indefinitely.
    """
    try:
        return sys.stdin.isatty()
    except (AttributeError, ValueError):
        return False


def _recon_rate_limited(sentinel: Path, *, min_gap_sec: int = 3600) -> bool:
    """Return True if recon was attempted less than `min_gap_sec` ago.

    Prevents a misbehaving caller (e.g. a user hammering ./run.sh after a
    real CF 403) from spawning back-to-back Chrome windows.
    """
    if not sentinel.exists():
        return False
    age = time.time() - sentinel.stat().st_mtime
    return age < min_gap_sec


def refresh_cookies():
    """Run recon_bookmaker.py to get fresh Cloudflare cookies."""
    import subprocess
    recon_script = Path(__file__).parent / "recon_bookmaker.py"
    python = sys.executable
    print(f"\nRunning recon script to refresh cookies...")
    print(f"  {python} {recon_script}\n")
    subprocess.run([python, str(recon_script)], check=True)

    # Extract cookies from recon output into our format
    recon_cookies = Path(__file__).parent / "recon_bookmaker_cookies.json"
    if recon_cookies.exists():
        with open(recon_cookies) as f:
            pw_cookies = json.load(f)
        cookie_dict = {}
        for c in pw_cookies:
            if "bookmaker.eu" in c.get("domain", ""):
                cookie_dict[c["name"]] = c["value"]
        with open(COOKIE_PATH, "w") as f:
            json.dump(cookie_dict, f)
        print(f"Saved {len(cookie_dict)} cookies from recon.")


# =============================================================================
# JSON PARSING
# =============================================================================


def parse_schedule(schedule_data: dict, league_config: dict, sport_key: str,
                   team_dict: dict, canonical_games: list,
                   fetch_time: str) -> list[dict]:
    """Parse GetSchedule response into 18-column DuckDB records.

    Response structure: Schedule.Data.Leagues.League[].dateGroup[].game[]
    Each game has Derivatives.line[] where index "0" is the main line.
    """
    records = []
    market = league_config["market"]
    period = league_config["period"]

    leagues = (schedule_data
               .get("Schedule", {})
               .get("Data", {})
               .get("Leagues", {})
               .get("League", []))

    for league in leagues:
        for date_group in league.get("dateGroup", []):
            for game in date_group.get("game", []):
                # _parse_game now returns a list (main + alt-line indices).
                records.extend(_parse_game(
                    game, market, period, sport_key,
                    team_dict, canonical_games, fetch_time
                ))

    return records


def _parse_game(game: dict, market: str, period: str, sport_key: str,
                team_dict: dict, canonical_games: list,
                fetch_time: str) -> list[dict]:
    """Parse a single game's `Derivatives.line` list into DuckDB records.

    BKM packs the main and ALL alt lines into one `Derivatives.line` array
    keyed by `index` (0 = main, ± non-zero = alternates). Previously this
    extracted only index 0 and discarded everything else; we now emit one
    record per usable index so the alt run-line / alt-total ladder reaches
    the dashboard. The main index keeps `market` and the ML; non-zero
    indices use `market = market.replace("spreads", "alternate_spreads")`
    and carry only spread + total (no ML — that's per-game, not per-line).
    Downstream `get_bookmaker_odds()` derives an `alternate_totals` canonical
    row from the same alt record via its existing `gsub("spreads","totals")`.
    """
    if game.get("Stat") != "O":
        return []

    away_raw = (game.get("vtm") or "").strip()
    home_raw = (game.get("htm") or "").strip()
    if not away_raw or not home_raw:
        return []

    # Resolve team names via shared canonical_match
    if team_dict or canonical_games:
        away_team, home_team = resolve_team_names(
            away_raw, home_raw, team_dict, canonical_games
        )
    else:
        away_team, home_team = away_raw, home_raw

    # Validate resolved teams against canonical games to filter out
    # cross-sport contamination (e.g., UFC fighters in MLB league 206)
    if canonical_games:
        matched = any(
            (away_team == cg["away_team"] and home_team == cg["home_team"]) or
            (away_team == cg["home_team"] and home_team == cg["away_team"])
            for cg in canonical_games
        )
        if not matched:
            return []

    # Game start time — parse BKM's naive PT (gmdt, gmtm) into UTC instant.
    game_start_time = _bkm_game_start_time(game.get("gmdt", ""), game.get("gmtm", ""))

    game_id = str(game.get("idgm", ""))

    lines = game.get("Derivatives", {}).get("line", [])
    if not lines:
        return []

    # market label for alternate rows: "spreads" -> "alternate_spreads",
    # "spreads_f5" -> "alternate_spreads_f5", "spreads_h2" -> "alternate_spreads_h2".
    alt_market = market.replace("spreads", "alternate_spreads")

    records: list[dict] = []
    for line in lines:
        is_main = str(line.get("index")) == "0"
        rec_market = market if is_main else alt_market

        away_spread = away_spread_price = home_spread = home_spread_price = None
        if line.get("s_sp") == 1:
            try:
                away_spread = float(line["vsprdt"])
                home_spread = float(line["hsprdt"])
                away_spread_price = int(line["vsprdoddst"])
                home_spread_price = int(line["hsprdoddst"])
            except (KeyError, ValueError, TypeError):
                pass

        total = over_price = under_price = None
        if line.get("s_tot") == 1:
            try:
                total = float(line["ovt"])
                over_price = int(line["ovoddst"])
                under_price = int(line["unoddst"])
            except (KeyError, ValueError, TypeError):
                pass

        # Moneyline is per-game (not per-line), so only the main index carries it.
        away_ml = home_ml = None
        if is_main and line.get("s_ml") == 1:
            try:
                away_ml = int(line["voddst"])
                home_ml = int(line["hoddst"])
            except (KeyError, ValueError, TypeError):
                pass

        # Skip lines with no usable odds. For the main line this preserves the
        # original behavior (game is dropped if main has nothing); for alts we
        # just drop that specific index.
        if away_spread is None and total is None and away_ml is None:
            continue

        records.append({
            "fetch_time": fetch_time,
            "sport_key": sport_key,
            "game_id": game_id,
            "game_start_time": game_start_time,
            "away_team": away_team,
            "home_team": home_team,
            "market": rec_market,
            "period": period,
            "away_spread": away_spread,
            "away_spread_price": away_spread_price,
            "home_spread": home_spread,
            "home_spread_price": home_spread_price,
            "total": total,
            "over_price": over_price,
            "under_price": under_price,
            "away_ml": away_ml,
            "home_ml": home_ml,
        })

    return records


# =============================================================================
# DATABASE
# =============================================================================


def init_database(sport: str):
    """Initialize DuckDB with the odds table for a sport."""
    config = SPORT_CONFIGS[sport]
    table_name = config["table_name"]

    conn = duckdb.connect(str(DB_PATH))
    try:
        # Self-healing migration: if the table predates the game_start_time
        # migration (old schema had game_date/game_time VARCHARs), drop it so
        # it's recreated with the new schema below. Ephemeral table — fully
        # repopulated every cycle, so dropping loses at most one stale snapshot.
        has_gst = conn.execute(
            "SELECT 1 FROM information_schema.columns "
            "WHERE table_name = ? AND column_name = 'game_start_time'",
            [table_name]
        ).fetchone()
        if has_gst is None:
            tbl_exists = conn.execute(
                "SELECT 1 FROM information_schema.tables WHERE table_name = ?",
                [table_name]
            ).fetchone()
            if tbl_exists is not None:
                print(f"[bkm] Migrating {table_name} to game_start_time schema "
                      f"(existing snapshot re-populated this run)")
                conn.execute(f"DROP TABLE {table_name}")

        conn.execute(f"""
            CREATE TABLE IF NOT EXISTS {table_name} (
                fetch_time TIMESTAMPTZ,
                sport_key VARCHAR,
                game_id VARCHAR,
                game_start_time TIMESTAMPTZ,
                away_team VARCHAR,
                home_team VARCHAR,
                market VARCHAR,
                period VARCHAR,
                away_spread FLOAT,
                away_spread_price INTEGER,
                home_spread FLOAT,
                home_spread_price INTEGER,
                total FLOAT,
                over_price INTEGER,
                under_price INTEGER,
                away_ml INTEGER,
                home_ml INTEGER
            )
        """)
    finally:
        conn.close()


def save_to_database(sport: str, odds_data: list):
    """Save scraped odds to DuckDB atomically (stage-and-swap)."""
    config = SPORT_CONFIGS[sport]
    table_name = config["table_name"]

    conn = duckdb.connect(str(DB_PATH))

    columns = [
        "fetch_time", "sport_key", "game_id", "game_start_time",
        "away_team", "home_team", "market", "period",
        "away_spread", "away_spread_price", "home_spread", "home_spread_price",
        "total", "over_price", "under_price", "away_ml", "home_ml"
    ]

    # Stage into a TEMP table cloned from the live schema, then atomically
    # swap. CREATE OR REPLACE TABLE ... AS SELECT is DuckDB's closest
    # equivalent of an atomic rename — readers see either the entire old
    # snapshot or the entire new one, never a half-written state. On an
    # empty scrape, leave the prior snapshot in place so consumers don't
    # see a momentarily empty table (matches BFA/DK/FD/WZ/Hoop88 pattern).
    if odds_data:
        conn.execute(
            f"CREATE OR REPLACE TEMP TABLE {table_name}_new AS "
            f"SELECT * FROM {table_name} LIMIT 0"
        )
        placeholders = ", ".join(["?" for _ in columns])
        tuples = [tuple(d.get(c) for c in columns) for d in odds_data]
        conn.executemany(
            f"INSERT INTO {table_name}_new VALUES ({placeholders})", tuples
        )
        conn.execute(
            f"CREATE OR REPLACE TABLE {table_name} AS SELECT * FROM {table_name}_new"
        )
        conn.execute(f"DROP TABLE IF EXISTS {table_name}_new")
    else:
        print(f"[bookmaker] empty scrape — leaving prior snapshot in place", flush=True)

    result = conn.execute(f"SELECT COUNT(*) FROM {table_name}").fetchone()
    print(f"Database now has {result[0]} total records in {table_name}")

    conn.close()


# =============================================================================
# MAIN
# =============================================================================


def scrape_bookmaker(sport: str):
    """Scrape Bookmaker.eu odds via curl_cffi and save to DuckDB.

    Uses saved cookies from recon_bookmaker.py for Cloudflare bypass.
    If cookies are expired, launches recon script to refresh them.
    """
    if sport not in SPORT_CONFIGS:
        raise ValueError(f"Unknown sport: {sport}. Available: {list(SPORT_CONFIGS.keys())}")

    config = SPORT_CONFIGS[sport]
    init_database(sport)

    # Aware UTC datetime — DuckDB stores it directly as TIMESTAMPTZ. A
    # naive strftime string would be interpreted as local-TZ on insert
    # and silently shift (same bug DK/FD hit pre-migration).
    fetch_time = datetime.now(timezone.utc)

    # Load team name resolution
    team_dict = load_team_dict(sport)
    canonical_games = load_canonical_games(sport)

    session = _create_session()

    # Hit the site to refresh Cloudflare cookies
    try:
        resp = session.get(SITE_URL, timeout=15)
        blocked = resp.status_code == 403
    except Exception:
        blocked = True

    # Session-health + league-discovery in one call: GetRoutingInfo returns
    # the full league catalog AND verifies the session is valid (it requires
    # cookies + a working CF bypass, same as GetSchedule).
    routing_info = None if blocked else fetch_routing_info(session)

    login_ok = True  # Default: assume prior-session cookies are still valid

    # If routing failed, try login before anything drastic
    if routing_info is None:
        login_ok = False
        if not blocked and BOOKMAKER_USERNAME and BOOKMAKER_PASSWORD:
            print("GetRoutingInfo returned nothing, logging in...")
            login_ok = login(session, BOOKMAKER_USERNAME, BOOKMAKER_PASSWORD)
            if login_ok:
                routing_info = fetch_routing_info(session)

    # If still no routing info, run the recon dance (same logic as before).
    if routing_info is None:
        sentinel = Path(__file__).parent / ".last_recon_attempt"
        if not _can_launch_interactive_recon():
            print("Cookies appear stale but stdin is not a TTY.")
            print("Skipping interactive recon. Run manually:")
            print(f"  cd {Path(__file__).parent} && ./venv/bin/python recon_bookmaker.py")
            save_to_database(sport, [])
            return []
        if _recon_rate_limited(sentinel):
            print("Recon was attempted recently (< 1h ago). Skipping to avoid popup loop.")
            save_to_database(sport, [])
            return []

        sentinel.touch()
        refresh_cookies()
        session = _create_session()
        resp = session.get(SITE_URL, timeout=15)
        if resp.status_code == 403:
            print("ERROR: Still blocked after recon. Check Cloudflare manually. Clearing stale data.")
            save_to_database(sport, [])
            return []
        if BOOKMAKER_USERNAME and BOOKMAKER_PASSWORD:
            login(session, BOOKMAKER_USERNAME, BOOKMAKER_PASSWORD)
        routing_info = fetch_routing_info(session)
        if routing_info is None:
            print("ERROR: GetRoutingInfo failed after recon + login. Clearing stale data.")
            save_to_database(sport, [])
            return []

    # Resolve wanted markets against BKM's live catalog. Any market that
    # doesn't match is logged and dropped — no silent mislabeling.
    resolved_leagues = resolve_leagues(
        routing_info,
        config["bkm_sport_id"],
        config["bkm_region"],
        config["markets"],
    )
    if not resolved_leagues:
        # Off-season, all markets renamed, or wrong sport_id. Session was
        # healthy enough to get routing info, so this is "BKM has nothing
        # for this sport right now" — save empty so stale rows clear.
        print(f"No {sport.upper()} markets resolved from BKM catalog — saving empty.")
        save_to_database(sport, [])
        return []

    # One-line summary so the operator can eyeball what got wired up before
    # the schedule fetches start. Format: "GAME LINES→5, 1ST 5 INNINGS→6, ..."
    routes = ", ".join(f"{lg['name']}→{lg['id']}" for lg in resolved_leagues)
    print(f"Resolved {len(resolved_leagues)} BKM {sport.upper()} market(s): {routes}")

    # Fetch each resolved league
    all_odds = []
    for league in resolved_leagues:
        print(f"Fetching {league['name']} (league {league['id']})...")
        data = fetch_schedule(session, league["id"])
        if data is None:
            print(f"  Failed to fetch league {league['id']}")
            continue

        records = parse_schedule(
            data, league, config["sport_key"],
            team_dict, canonical_games, fetch_time
        )
        all_odds.extend(records)
        print(f"  Parsed {len(records)} records")

    # Save cookies for next run
    _save_cookies(session)

    # Summary by market
    market_counts = {}
    for rec in all_odds:
        m = rec["market"]
        market_counts[m] = market_counts.get(m, 0) + 1
    for m, c in sorted(market_counts.items()):
        print(f"  {m}: {c} records")

    # Save to DuckDB (always save, even if empty, to clear stale data)
    save_to_database(sport, all_odds)

    print(f"\nScraped {len(all_odds)} total records for {sport.upper()}")
    return all_odds


if __name__ == "__main__":
    sport = sys.argv[1] if len(sys.argv) > 1 else "cbb"
    print(f"Starting Bookmaker.eu {sport.upper()} odds scraper...")
    scrape_bookmaker(sport)
