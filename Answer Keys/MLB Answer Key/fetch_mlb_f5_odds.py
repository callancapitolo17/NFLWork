#!/usr/bin/env python3
"""
Fetch historical MLB F5 odds from The Odds API for ROI backtesting.

2-step API approach (verified cost: 1 credit events list + 30 credits/game):
  1. GET /v4/historical/sports/baseball_mlb/events/?date=<snap>  → event IDs
  2. GET /v4/historical/sports/baseball_mlb/events/{id}/odds/    → F5 odds

Budget: keep 150,000 credits in reserve. At 30 credits/game, all 12,719
historical games (~381,000 total) fit within the budget.

Stores in pbp.duckdb/mlb_f5_odds_history (normalized: one row per
game_pk / bookmaker / market / outcome).

Usage:
    python fetch_mlb_f5_odds.py               # fetch all missing games
    python fetch_mlb_f5_odds.py --dry-run     # estimate cost only
    python fetch_mlb_f5_odds.py --status      # show progress by season
"""

import argparse
import os
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path

import duckdb
import requests

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
def _load_api_key() -> str:
    """Load Odds API key from env var or ~/.Renviron fallback."""
    key = os.environ.get("ODDS_API_KEY")
    if key:
        return key
    renviron = Path.home() / ".Renviron"
    if renviron.exists():
        for line in renviron.read_text().splitlines():
            if line.strip().startswith("ODDS_API_KEY"):
                return line.split("=", 1)[1].strip()
    raise RuntimeError("ODDS_API_KEY not found in env or ~/.Renviron")

API_KEY = _load_api_key()
BASE_URL = "https://api.the-odds-api.com/v4/historical/sports/baseball_mlb"
MARKETS = "h2h_1st_5_innings,totals_1st_5_innings,spreads_1st_5_innings"
REGIONS = "us"
RESERVE_CREDITS = 150_000
COST_PER_GAME = 30      # empirically: 3 F5 markets × us region = 30/game
COST_EVENTS_LIST = 1    # per events-list call
CALL_DELAY_SEC = 0.1    # throttle: ~10 calls/sec max (tested safe)
SNAPSHOT_OFFSET_MIN = 240  # 4 hours before game start (F5 archive has gaps at T-30min)

DB_PATH = Path(__file__).parent.parent / "pbp.duckdb"

# ---------------------------------------------------------------------------
# Database
# ---------------------------------------------------------------------------

CREATE_TABLE_SQL = """
CREATE TABLE IF NOT EXISTS mlb_f5_odds_history (
    game_pk       VARCHAR,
    snapshot_time TIMESTAMP,
    bookmaker     VARCHAR,
    market        VARCHAR,
    key           VARCHAR,    -- team name (h2h/spreads) or "Over"/"Under" (totals)
    price         INTEGER,    -- American odds
    point         FLOAT       -- spread or total line (NULL for h2h)
)
"""

CREATE_INDEX_SQL = """
CREATE UNIQUE INDEX IF NOT EXISTS idx_mlb_f5_odds_history
    ON mlb_f5_odds_history (game_pk, bookmaker, market, key)
"""


def init_db(conn):
    conn.execute(CREATE_TABLE_SQL)
    try:
        conn.execute(CREATE_INDEX_SQL)
    except duckdb.CatalogException:
        pass


def get_fetched_pks(conn) -> set:
    rows = conn.execute(
        "SELECT DISTINCT game_pk FROM mlb_f5_odds_history"
    ).fetchall()
    return {r[0] for r in rows}


def get_target_games(conn, fetched_pks: set) -> list[dict]:
    """All F5-eligible games not yet fetched, sorted most-recent first."""
    rows = conn.execute("""
        SELECT game_pk, game_date,
               COALESCE(game_start_time, CAST(game_date AS TIMESTAMP) + INTERVAL '18 hours') AS game_start_time,
               home_team, away_team,
               EXTRACT(YEAR FROM game_date)::INTEGER AS season
        FROM mlb_betting_pbp
        WHERE game_home_margin_inning_inning_5 IS NOT NULL
        ORDER BY game_date DESC, game_start_time ASC NULLS LAST
    """).fetchall()

    return [
        {
            "game_pk": str(game_pk),
            "game_date": game_date,
            "game_start_time": game_start_time,  # UTC datetime
            "home_team": home_team,
            "away_team": away_team,
            "season": season,
        }
        for (game_pk, game_date, game_start_time, home_team, away_team, season) in rows
        if str(game_pk) not in fetched_pks
    ]


def insert_records(conn, records: list[dict]):
    if not records:
        return
    conn.executemany("""
        INSERT OR IGNORE INTO mlb_f5_odds_history
            (game_pk, snapshot_time, bookmaker, market, key, price, point)
        VALUES (?, ?, ?, ?, ?, ?, ?)
    """, [
        (r["game_pk"], r["snapshot_time"], r["bookmaker"],
         r["market"], r["key"], r["price"], r["point"])
        for r in records
    ])


# ---------------------------------------------------------------------------
# API helpers
# ---------------------------------------------------------------------------

def _get(url: str, params: dict) -> tuple[dict | list | None, int, int]:
    """Make a GET request. Returns (data, credits_used, credits_remaining)."""
    try:
        resp = requests.get(url, params=params, timeout=15)
        used = int(resp.headers.get("x-requests-last", 0))
        remaining = int(resp.headers.get("x-requests-remaining", 0))
        if resp.status_code != 200:
            print(f"  HTTP {resp.status_code}: {resp.text[:200]}")
            return None, used, remaining
        return resp.json(), used, remaining
    except Exception as e:
        print(f"  Request error: {e}")
        return None, 0, 0


def fetch_events_at(snap_dt: datetime) -> tuple[list, int]:
    """
    Call /events/ endpoint at snap_dt.
    Returns (events_list, credits_remaining).
    events_list: [{id, home_team, away_team, commence_time}, ...]
    """
    iso = snap_dt.strftime("%Y-%m-%dT%H:%M:%SZ")
    data, _, remaining = _get(
        f"{BASE_URL}/events/",
        {"apiKey": API_KEY, "date": iso}
    )
    if data is None:
        return [], remaining
    events = data.get("data", []) if isinstance(data, dict) else []
    return events, remaining


def fetch_event_f5_odds(event_id: str, snap_dt: datetime) -> tuple[dict | None, int]:
    """
    Call /events/{id}/odds/ for F5 markets at snap_dt.
    Returns (event_data, credits_remaining).
    """
    iso = snap_dt.strftime("%Y-%m-%dT%H:%M:%SZ")
    data, _, remaining = _get(
        f"{BASE_URL}/events/{event_id}/odds/",
        {
            "apiKey": API_KEY,
            "regions": REGIONS,
            "markets": MARKETS,
            "date": iso,
            "oddsFormat": "american",
        }
    )
    if data is None:
        return None, remaining
    return data.get("data"), remaining


# ---------------------------------------------------------------------------
# Matching and parsing
# ---------------------------------------------------------------------------

def normalize(name: str) -> str:
    return name.lower().strip()


def build_team_index(events: list[dict]) -> dict:
    """Map (norm_away, norm_home) → event_id."""
    return {
        (normalize(e.get("away_team", "")), normalize(e.get("home_team", ""))): e["id"]
        for e in events
        if e.get("id") and e.get("home_team") and e.get("away_team")
    }


def parse_event_odds(event_data: dict, game_pk: str, snap_dt: datetime) -> list[dict]:
    """Flatten bookmaker/market/outcome structure into flat records."""
    records = []
    snap_str = snap_dt.strftime("%Y-%m-%d %H:%M:%S")
    for bookmaker in event_data.get("bookmakers", []):
        book_key = bookmaker.get("key", "")
        for market in bookmaker.get("markets", []):
            mkt_key = market.get("key", "")
            for outcome in market.get("outcomes", []):
                records.append({
                    "game_pk": game_pk,
                    "snapshot_time": snap_str,
                    "bookmaker": book_key,
                    "market": mkt_key,
                    "key": outcome.get("name", ""),
                    "price": outcome.get("price"),
                    "point": outcome.get("point"),
                })
    return records


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def check_current_balance() -> int:
    """Hit a cheap endpoint to get current credits remaining."""
    resp = requests.get(
        "https://api.the-odds-api.com/v4/sports/",
        params={"apiKey": API_KEY},
        timeout=10
    )
    return int(resp.headers.get("x-requests-remaining", 0))


def run(dry_run: bool = False, status_only: bool = False):
    conn = duckdb.connect(str(DB_PATH))
    init_db(conn)

    fetched_pks = get_fetched_pks(conn)
    target_games = get_target_games(conn, fetched_pks)

    # Count by season
    total_games = conn.execute("""
        SELECT COUNT(*) FROM mlb_betting_pbp
        WHERE game_home_margin_inning_inning_5 IS NOT NULL
    """).fetchone()[0]

    print("\n=== MLB F5 Historical Odds Acquisition ===")
    print(f"Total F5-eligible games: {total_games:,}")
    print(f"Already fetched:         {len(fetched_pks):,}")
    print(f"Remaining to fetch:      {len(target_games):,}")

    if status_only:
        season_rows = conn.execute("""
            SELECT EXTRACT(YEAR FROM game_date)::INTEGER AS season, COUNT(*) AS total
            FROM mlb_betting_pbp
            WHERE game_home_margin_inning_inning_5 IS NOT NULL
            GROUP BY 1 ORDER BY 1 DESC
        """).fetchall()
        fetched_rows = conn.execute("""
            SELECT EXTRACT(YEAR FROM CAST(snapshot_time AS DATE))::INTEGER AS season,
                   COUNT(DISTINCT game_pk)
            FROM mlb_f5_odds_history
            GROUP BY 1 ORDER BY 1 DESC
        """).fetchall() if fetched_pks else []

        fetched_map = {r[0]: r[1] for r in fetched_rows}
        print("\nSeason breakdown:")
        for season, total in season_rows:
            done = fetched_map.get(season, 0)
            pct = done / total * 100 if total else 0
            print(f"  {season}: {done:,}/{total:,} ({pct:.0f}%)")
        conn.close()
        return

    if not target_games:
        print("All games already fetched!")
        conn.close()
        return

    # Budget check
    current_remaining = check_current_balance()
    available = current_remaining - RESERVE_CREDITS
    # Estimate: 1 events-list call per unique snapshot + 30/game
    n_snaps = len({
        (g["game_start_time"] - timedelta(minutes=SNAPSHOT_OFFSET_MIN)).replace(second=0, microsecond=0)
        for g in target_games
    })
    estimated_cost = COST_EVENTS_LIST * n_snaps + COST_PER_GAME * len(target_games)
    max_affordable = max(0, (available - COST_EVENTS_LIST * n_snaps) // COST_PER_GAME)

    print(f"\nAPI credits remaining:  {current_remaining:,}")
    print(f"Reserve:                {RESERVE_CREDITS:,}")
    print(f"Available to spend:     {available:,}")
    print(f"Estimated total cost:   {estimated_cost:,}")
    print(f"Max games affordable:   {max_affordable:,}")
    games_to_process = min(len(target_games), max_affordable)
    print(f"Games to process:       {games_to_process:,}")

    if dry_run:
        print("\n[DRY RUN] First 10 games:")
        for g in target_games[:10]:
            snap = g["game_start_time"] - timedelta(minutes=SNAPSHOT_OFFSET_MIN)
            print(f"  {g['season']} {g['game_date']}  "
                  f"{g['away_team']} @ {g['home_team']}  "
                  f"snap={snap.strftime('%H:%M UTC')}")
        conn.close()
        return

    if games_to_process <= 0:
        print("Insufficient budget. Aborting.")
        conn.close()
        return

    games_queue = target_games[:games_to_process]

    # Group games by snapshot_dt (truncated to minute)
    # One events-list call per unique snapshot can serve multiple games
    from collections import defaultdict
    snap_groups: dict[datetime, list[dict]] = defaultdict(list)
    for g in games_queue:
        snap_dt = (g["game_start_time"] - timedelta(minutes=SNAPSHOT_OFFSET_MIN)
                   ).replace(second=0, microsecond=0)
        snap_groups[snap_dt].append(g)

    print(f"\nFetching {games_to_process:,} games across "
          f"{len(snap_groups):,} snapshot times...\n")

    remaining_credits = current_remaining
    fetched_count = 0
    no_odds_count = 0
    consecutive_422s = 0
    api_calls = 0
    pending = {g["game_pk"]: g for g in games_queue}  # global pending

    # Process newest-first: recent games have reliable F5 historical coverage
    # Stop early if we hit consecutive 422s (archive boundary reached)
    for snap_dt in sorted(snap_groups.keys(), reverse=True):
        if not pending:
            break
        if remaining_credits - RESERVE_CREDITS < COST_PER_GAME:
            print("Budget reserve reached. Stopping.")
            break

        # --- Step 1: get event IDs at this snapshot ---
        n_targeted = sum(1 for g in snap_groups[snap_dt] if g["game_pk"] in pending)
        if n_targeted == 0:
            continue  # all games at this snapshot already fetched
        print(f"\n[{snap_dt.strftime('%Y-%m-%d %H:%M UTC')}] "
              f"{n_targeted} game(s) | credits: {remaining_credits:,}", flush=True)
        events, ev_remaining = fetch_events_at(snap_dt)
        if ev_remaining > 0:
            remaining_credits = ev_remaining  # only trust non-zero values
        api_calls += 1
        time.sleep(CALL_DELAY_SEC)

        if not events:
            # No events at this snapshot — don't pop games, they may appear
            # at a different snapshot time. Just increment failure counter.
            consecutive_422s += 1
            continue

        consecutive_422s = 0  # events list succeeded
        team_index = build_team_index(events)

        # --- Step 2: per-event F5 odds for each pending game visible at snap ---
        for g in snap_groups[snap_dt]:
            game_pk = g["game_pk"]
            if game_pk not in pending:
                continue  # already captured by an earlier snapshot

            key = (normalize(g["away_team"]), normalize(g["home_team"]))
            event_id = team_index.get(key)

            if not event_id:
                # Try partial match (some books shorten team names)
                for (ev_away, ev_home), eid in team_index.items():
                    if (ev_away in normalize(g["away_team"]) or
                            normalize(g["away_team"]) in ev_away) and \
                       (ev_home in normalize(g["home_team"]) or
                            normalize(g["home_team"]) in ev_home):
                        event_id = eid
                        break

            if not event_id:
                pending.pop(game_pk, None)
                no_odds_count += 1
                continue

            if remaining_credits - RESERVE_CREDITS < COST_PER_GAME:
                print("Budget reserve reached mid-batch. Stopping.")
                break

            event_data, od_remaining = fetch_event_f5_odds(event_id, snap_dt)
            if od_remaining > 0:
                remaining_credits = od_remaining  # only trust non-zero values
            api_calls += 1
            time.sleep(CALL_DELAY_SEC)

            if event_data:
                records = parse_event_odds(event_data, game_pk, snap_dt)
                if records:
                    insert_records(conn, records)
                    conn.commit()
                    fetched_count += 1
                    consecutive_422s = 0
                    print(f"  [{fetched_count:>5}] {g['season']} {g['game_date']}  "
                          f"{g['away_team']} @ {g['home_team']}  "
                          f"{len(records)} records  "
                          f"({remaining_credits:,} credits left)")
                else:
                    no_odds_count += 1
                    consecutive_422s += 1
            else:
                consecutive_422s += 1

            pending.pop(game_pk, None)

        # Stop if we've hit a long streak of failures (archive boundary)
        if consecutive_422s >= 50:
            print(f"\nHit {consecutive_422s} consecutive failures at "
                  f"{snap_dt.strftime('%Y-%m-%d')} — archive boundary reached. Stopping.")
            break

    conn.close()

    print(f"\n=== Done ===")
    print(f"API calls:         {api_calls:,}")
    print(f"Games fetched:     {fetched_count:,}")
    print(f"Games w/no F5:     {no_odds_count:,}")
    print(f"Credits remaining: {remaining_credits:,}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true",
                        help="Show plan without calling API")
    parser.add_argument("--status", action="store_true",
                        help="Show fetch progress by season")
    args = parser.parse_args()
    run(dry_run=args.dry_run, status_only=args.status)
