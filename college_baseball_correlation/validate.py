#!/usr/bin/env python3
"""
College Baseball Historical Data Acquisition & Analysis

CLI tool to build the historical dataset needed for the college baseball
correlated parlay answer key.

Subcommands:
    --pull-scores 2024 2025    Scrape ESPN college baseball scores
    --pull-odds 2024 2025      Pull closing lines from Odds API
    --analyze                  Structural stats (mercy rate, extras, walk-off)

Data stored in college_baseball.duckdb:
    games  — ESPN game results (scores, regulation/extras, mercy)
    odds   — Odds API closing lines (ML, totals)

Usage:
    python validate.py --pull-scores 2024 2025
    python validate.py --pull-odds 2024 2025
    python validate.py --analyze
    python validate.py --pull-scores 2024 2025 --pull-odds 2024 2025 --analyze
"""

import argparse
import json
import math
import os
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path

import duckdb
import requests

DB_PATH = Path(__file__).parent / "college_baseball.duckdb"

# ESPN API endpoints
ESPN_SCOREBOARD = "https://site.api.espn.com/apis/site/v2/sports/baseball/college-baseball/scoreboard"

# Odds API
ODDS_API_BASE = "https://api.the-odds-api.com/v4"
ODDS_API_KEY = os.environ.get("ODDS_API_KEY", "")

# Fall back to ~/.Renviron if not in environment
if not ODDS_API_KEY:
    renviron = Path.home() / ".Renviron"
    if renviron.exists():
        for line in renviron.read_text().splitlines():
            if line.startswith("ODDS_API_KEY"):
                ODDS_API_KEY = line.split("=", 1)[1].strip().strip('"').strip("'")
                break


# =============================================================================
# DATABASE
# =============================================================================

def init_db():
    """Initialize DuckDB with games and odds tables."""
    con = duckdb.connect(str(DB_PATH))
    con.execute("""
        CREATE TABLE IF NOT EXISTS games (
            game_id VARCHAR PRIMARY KEY,
            date DATE,
            season INTEGER,
            home_team VARCHAR,
            away_team VARCHAR,
            home_score INTEGER,
            away_score INTEGER,
            total_runs INTEGER,
            home_margin INTEGER,
            home_won BOOLEAN,
            regulation_innings INTEGER,
            went_extras BOOLEAN,
            mercy_rule BOOLEAN,
            status VARCHAR
        )
    """)
    con.execute("""
        CREATE TABLE IF NOT EXISTS linescores (
            game_id VARCHAR,
            inning INTEGER,
            half VARCHAR,  -- 'top' or 'bottom'
            team VARCHAR,
            runs INTEGER,
            hits INTEGER,
            errors INTEGER,
            PRIMARY KEY (game_id, inning, half)
        )
    """)
    con.execute("""
        CREATE TABLE IF NOT EXISTS odds (
            game_id VARCHAR,
            commence_time TIMESTAMP,
            home_team VARCHAR,
            away_team VARCHAR,
            home_ml INTEGER,
            away_ml INTEGER,
            total_line FLOAT,
            over_odds INTEGER,
            under_odds INTEGER,
            bookmaker_key VARCHAR,
            sport_key VARCHAR,
            PRIMARY KEY (game_id, bookmaker_key)
        )
    """)
    con.close()


# =============================================================================
# ESPN SCORE SCRAPING
# =============================================================================

def pull_scores_for_date(session: requests.Session, date_str: str) -> list[dict]:
    """Pull all college baseball scores for a single date from ESPN."""
    params = {
        "dates": date_str,
        "limit": 500,
    }
    resp = session.get(ESPN_SCOREBOARD, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    games = []
    for event in data.get("events", []):
        competition = event.get("competitions", [{}])[0]
        competitors = competition.get("competitors", [])
        if len(competitors) != 2:
            continue

        status_detail = event.get("status", {})
        status_type = status_detail.get("type", {}).get("name", "")
        if status_type != "STATUS_FINAL":
            continue

        # Parse teams and scores
        home = away = None
        for comp in competitors:
            if comp.get("homeAway") == "home":
                home = comp
            else:
                away = comp

        if not home or not away:
            continue

        home_score = int(home.get("score", 0))
        away_score = int(away.get("score", 0))

        # Detect extras / mercy
        status_text = status_detail.get("type", {}).get("shortDetail", "")
        period = status_detail.get("period", 9)
        went_extras = period > 9
        # Mercy rule: game ended early with large margin (7+ after 7th, 10+ after 5th)
        margin = abs(home_score - away_score)
        mercy_rule = (period < 9 and margin >= 10) or (period == 7 and margin >= 7)

        game = {
            "game_id": str(event.get("id", "")),
            "date": date_str[:4] + "-" + date_str[4:6] + "-" + date_str[6:8],
            "season": int(date_str[:4]),
            "home_team": home.get("team", {}).get("displayName", ""),
            "away_team": away.get("team", {}).get("displayName", ""),
            "home_score": home_score,
            "away_score": away_score,
            "total_runs": home_score + away_score,
            "home_margin": home_score - away_score,
            "home_won": home_score > away_score,
            "regulation_innings": min(period, 9),
            "went_extras": went_extras,
            "mercy_rule": mercy_rule,
            "status": status_type,
        }
        games.append(game)

    return games


def pull_scores(seasons: list[int]):
    """Pull all college baseball scores for given seasons from ESPN."""
    init_db()
    session = requests.Session()
    session.headers.update({"User-Agent": "Mozilla/5.0"})

    total_games = 0
    for season in seasons:
        # College baseball: mid-February through late June
        start = datetime(season, 2, 14)
        end = datetime(season, 6, 30)

        current = start
        season_games = 0
        while current <= end:
            date_str = current.strftime("%Y%m%d")
            try:
                games = pull_scores_for_date(session, date_str)
                if games:
                    con = duckdb.connect(str(DB_PATH))
                    for g in games:
                        con.execute("""
                            INSERT OR REPLACE INTO games VALUES (
                                ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?
                            )
                        """, [
                            g["game_id"], g["date"], g["season"],
                            g["home_team"], g["away_team"],
                            g["home_score"], g["away_score"],
                            g["total_runs"], g["home_margin"], g["home_won"],
                            g["regulation_innings"], g["went_extras"],
                            g["mercy_rule"], g["status"],
                        ])
                    con.close()
                    season_games += len(games)

                    if season_games % 200 < len(games):
                        print(f"  {season} {date_str}: {season_games} games so far")
            except Exception as e:
                print(f"  Warning: {date_str} failed: {e}")

            current += timedelta(days=1)
            # Polite rate limiting
            time.sleep(0.15)

        total_games += season_games
        print(f"Season {season}: {season_games} games pulled")

    print(f"\nTotal: {total_games} games across {len(seasons)} seasons")

    con = duckdb.connect(str(DB_PATH))
    count = con.execute("SELECT COUNT(*) FROM games").fetchone()[0]
    con.close()
    print(f"Database now has {count} total games")


# =============================================================================
# ESPN LINESCORE PULLING
# =============================================================================

ESPN_SUMMARY = "https://site.api.espn.com/apis/site/v2/sports/baseball/college-baseball/summary"


def _fetch_linescore(session, game_id):
    """Fetch linescore for a single game. Returns list of row tuples or None."""
    try:
        resp = session.get(ESPN_SUMMARY, params={"event": game_id}, timeout=30)
        resp.raise_for_status()
        data = resp.json()

        header = data.get("header", {})
        comps = header.get("competitions", [])
        if not comps:
            return None

        competitors = comps[0].get("competitors", [])
        if len(competitors) != 2:
            return None

        rows = []
        for comp in competitors:
            team_name = comp.get("team", {}).get("displayName", "")
            is_home = comp.get("homeAway") == "home"
            linescores = comp.get("linescores", [])

            if not linescores:
                return None

            for inning_idx, ls in enumerate(linescores):
                inning = inning_idx + 1
                half = "bottom" if is_home else "top"
                runs = int(ls.get("displayValue", 0))
                hits = int(ls.get("hits", 0))
                errs = int(ls.get("errors", 0))
                rows.append((game_id, inning, half, team_name, runs, hits, errs))

        return rows if rows else None
    except Exception:
        return "error"


def pull_linescores(seasons: list[int]):
    """Pull inning-by-inning linescores from ESPN summary API for given seasons.

    Uses concurrent requests (8 workers) and batch DB inserts for speed.
    """
    from concurrent.futures import ThreadPoolExecutor, as_completed

    init_db()
    con = duckdb.connect(str(DB_PATH))

    # Get game IDs for requested seasons that don't already have linescores
    game_ids = [r[0] for r in con.execute("""
        SELECT g.game_id
        FROM games g
        WHERE g.season = ANY(?)
          AND g.game_id NOT IN (SELECT DISTINCT game_id FROM linescores)
        ORDER BY g.date
    """, [seasons]).fetchall()]

    con.close()
    print(f"Found {len(game_ids)} games needing linescores for seasons {seasons}")

    if not game_ids:
        print("All games already have linescores.")
        return

    session = requests.Session()
    session.headers.update({"User-Agent": "Mozilla/5.0"})

    success = 0
    skipped = 0
    errors = 0
    batch = []
    batch_size = 100

    def flush_batch(batch):
        if not batch:
            return
        con = duckdb.connect(str(DB_PATH))
        for row in batch:
            con.execute("INSERT OR REPLACE INTO linescores VALUES (?, ?, ?, ?, ?, ?, ?)", list(row))
        con.close()

    # Process in parallel with 8 workers
    with ThreadPoolExecutor(max_workers=8) as executor:
        futures = {executor.submit(_fetch_linescore, session, gid): gid for gid in game_ids}

        for i, future in enumerate(as_completed(futures)):
            result = future.result()

            if result is None:
                skipped += 1
            elif result == "error":
                errors += 1
            else:
                batch.extend(result)
                success += 1

            # Flush batch periodically
            if len(batch) >= batch_size * 17:  # ~17 rows per game
                flush_batch(batch)
                batch = []

            if (i + 1) % 500 == 0:
                flush_batch(batch)
                batch = []
                print(f"  Progress: {i+1}/{len(game_ids)} ({success} ok, {skipped} skipped, {errors} errors)")

    # Final flush
    flush_batch(batch)

    print(f"\nDone: {success} linescores pulled, {skipped} skipped, {errors} errors")

    con = duckdb.connect(str(DB_PATH))
    total = con.execute("SELECT COUNT(DISTINCT game_id) FROM linescores").fetchone()[0]
    con.close()
    print(f"Database now has linescores for {total} games")


# =============================================================================
# ODDS API PULLING
# =============================================================================

def pull_odds(seasons: list[int]):
    """Pull historical closing lines from Odds API."""
    if not ODDS_API_KEY:
        print("Error: ODDS_API_KEY environment variable not set")
        sys.exit(1)

    init_db()
    sport_key = "baseball_ncaa"

    total_records = 0
    for season in seasons:
        # College baseball season: Feb-Jun
        start = datetime(season, 2, 14)
        end = datetime(season, 6, 30)

        # Odds API historical endpoint: one day at a time
        current = start
        while current <= end:
            date_iso = current.strftime("%Y-%m-%dT12:00:00Z")
            try:
                resp = requests.get(
                    f"{ODDS_API_BASE}/sports/{sport_key}/odds-history",
                    params={
                        "apiKey": ODDS_API_KEY,
                        "regions": "us,us2,eu",
                        "markets": "h2h,totals",
                        "oddsFormat": "american",
                        "dateFormat": "iso",
                        "date": date_iso,
                    },
                    timeout=30,
                )

                if resp.status_code == 404:
                    # Sport may not exist for this date
                    current += timedelta(days=1)
                    continue
                elif resp.status_code == 422:
                    # No data for this date
                    current += timedelta(days=1)
                    continue
                resp.raise_for_status()

                data = resp.json().get("data", [])
                remaining = resp.headers.get("x-requests-remaining", "?")

                day_records = 0
                con = duckdb.connect(str(DB_PATH))

                for game in data:
                    game_id = game.get("id", "")
                    commence = game.get("commence_time", "")
                    home_team = game.get("home_team", "")
                    away_team = game.get("away_team", "")

                    for bookmaker in game.get("bookmakers", []):
                        bk_key = bookmaker.get("key", "")
                        home_ml = away_ml = None
                        total_line = over_odds = under_odds = None

                        for market in bookmaker.get("markets", []):
                            mk = market.get("key", "")
                            outcomes = {o["name"]: o for o in market.get("outcomes", [])}

                            if mk == "h2h":
                                home_ml = outcomes.get(home_team, {}).get("price")
                                away_ml = outcomes.get(away_team, {}).get("price")
                            elif mk == "totals":
                                over_data = outcomes.get("Over", {})
                                under_data = outcomes.get("Under", {})
                                total_line = over_data.get("point")
                                over_odds = over_data.get("price")
                                under_odds = under_data.get("price")

                        if home_ml is not None or total_line is not None:
                            con.execute("""
                                INSERT OR REPLACE INTO odds VALUES (
                                    ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?
                                )
                            """, [
                                game_id, commence, home_team, away_team,
                                home_ml, away_ml, total_line,
                                over_odds, under_odds, bk_key, sport_key,
                            ])
                            day_records += 1

                con.close()
                total_records += day_records

                if day_records > 0:
                    print(f"  {current.strftime('%Y-%m-%d')}: {day_records} odds records (API remaining: {remaining})")

            except requests.exceptions.HTTPError as e:
                if e.response.status_code == 429:
                    print("  Rate limited, waiting 60s...")
                    time.sleep(60)
                    continue
                print(f"  Warning: {current.strftime('%Y-%m-%d')} failed: {e}")
            except Exception as e:
                print(f"  Warning: {current.strftime('%Y-%m-%d')} failed: {e}")

            current += timedelta(days=1)
            time.sleep(0.5)  # Rate limit

    print(f"\nTotal: {total_records} odds records")

    con = duckdb.connect(str(DB_PATH))
    count = con.execute("SELECT COUNT(*) FROM odds").fetchone()[0]
    games_with_odds = con.execute("""
        SELECT COUNT(DISTINCT o.game_id)
        FROM odds o
        INNER JOIN games g ON g.home_team = o.home_team AND g.away_team = o.away_team
            AND g.date = CAST(o.commence_time AS DATE)
    """).fetchone()[0]
    con.close()
    print(f"Database now has {count} total odds records, {games_with_odds} games matched to scores")


# =============================================================================
# ANALYSIS
# =============================================================================

def normalize_team(name: str) -> str:
    """Normalize team name for fuzzy matching (ESPN vs Odds API differences)."""
    s = name.strip()
    # Common abbreviation differences
    replacements = [
        ("State ", "St "),
        ("App St ", "Appalachian St "),
        ("-Pine Bluff ", "-PB "),
        ("Little Rock", "Little Rock"),
        ("UConn Huskies", "Connecticut Huskies"),
        ("UMass ", "Massachusetts "),
        ("UNLV ", "UNLV "),
        ("UCF ", "UCF "),
        ("USC ", "Southern California "),
        ("LSU ", "LSU "),
        ("SMU ", "SMU "),
        ("TCU ", "TCU "),
        ("UTEP ", "UTEP "),
        ("UTSA ", "UT San Antonio "),
        ("UNC ", "North Carolina "),
        ("Ole Miss ", "Ole Miss "),
        ("Pitt ", "Pittsburgh "),
    ]
    return s


def build_team_mapping(con) -> dict:
    """Build ESPN -> Odds API team name mapping using fuzzy matching."""
    espn_teams = [r[0] for r in con.execute("SELECT DISTINCT home_team FROM games").fetchall()]
    odds_teams = [r[0] for r in con.execute("SELECT DISTINCT home_team FROM odds").fetchall()]

    mapping = {}

    # Pass 1: exact match
    odds_set = set(odds_teams)
    for et in espn_teams:
        if et in odds_set:
            mapping[et] = et

    # Pass 2: normalize common differences
    def make_key(name):
        k = name
        k = k.replace(" State ", " St ")
        k = k.replace("App St ", "Appalachian St ")
        # Hyphen vs space
        k = k.replace("-", " ")
        k = k.replace("SIU Edwardsville", "SIU-Edwardsville")
        return k.lower().strip()

    odds_by_key = {make_key(ot): ot for ot in odds_teams}
    for et in espn_teams:
        if et not in mapping:
            k = make_key(et)
            if k in odds_by_key:
                mapping[et] = odds_by_key[k]

    # Pass 3: hardcoded known mismatches (ESPN -> Odds API)
    hardcoded = {
        "Army Black Knights": "Army Knights",
        "Purdue Fort Wayne Mastodons": "Fort Wayne Mastodons",
        "Long Island University Sharks": "LIU Sharks",
        "Prairie View A&M Panthers": "Prairie View Panthers",
        "South Dakota State Jackrabbits": "South Dakota St Jackrabbits",
        "UT Arlington Mavericks": "UT-Arlington Mavericks",
        "Youngstown State Penguins": "Youngstown St Penguins",
        "Morehead State Eagles": "Morehead St Eagles",
    }
    for espn_name, odds_name in hardcoded.items():
        if espn_name in espn_teams and odds_name in odds_set and espn_name not in mapping:
            mapping[espn_name] = odds_name

    # Pass 4: match on last word (mascot) + first word similarity
    unmatched_espn = [et for et in espn_teams if et not in mapping]
    unmatched_odds = [ot for ot in odds_teams if ot not in mapping.values()]

    for et in unmatched_espn:
        et_parts = et.split()
        et_mascot = et_parts[-1].lower() if et_parts else ""
        for ot in unmatched_odds:
            ot_parts = ot.split()
            ot_mascot = ot_parts[-1].lower() if ot_parts else ""
            if et_mascot == ot_mascot and et_mascot:
                # Same mascot — check if location words overlap
                et_loc = set(w.lower() for w in et_parts[:-1])
                ot_loc = set(w.lower() for w in ot_parts[:-1])
                if et_loc & ot_loc:  # at least one location word in common
                    mapping[et] = ot
                    unmatched_odds.remove(ot)
                    break

    return mapping


def match_teams():
    """Normalize ESPN team names to match Odds API names in the games table."""
    con = duckdb.connect(str(DB_PATH))

    mapping = build_team_mapping(con)
    print(f"Built team mapping: {len(mapping)} ESPN -> Odds API pairs")

    # Show unmatched
    espn_teams = set(r[0] for r in con.execute(
        "SELECT DISTINCT home_team FROM games UNION SELECT DISTINCT away_team FROM games"
    ).fetchall())
    odds_teams = set(r[0] for r in con.execute(
        "SELECT DISTINCT home_team FROM odds UNION SELECT DISTINCT away_team FROM odds"
    ).fetchall())

    mapped_odds = set(mapping.values())
    unmatched_odds = odds_teams - mapped_odds
    if unmatched_odds:
        print(f"Unmatched Odds API teams ({len(unmatched_odds)}): {sorted(unmatched_odds)[:10]}...")

    # Apply mapping: update games table team names to match Odds API
    updated = 0
    for espn_name, odds_name in mapping.items():
        if espn_name != odds_name:
            n1 = con.execute(
                "UPDATE games SET home_team = ? WHERE home_team = ?",
                [odds_name, espn_name]
            ).fetchone()
            n2 = con.execute(
                "UPDATE games SET away_team = ? WHERE away_team = ?",
                [odds_name, espn_name]
            ).fetchone()
            updated += 1

    print(f"Updated {updated} team name mappings in games table")

    # Verify match count
    matched = con.execute("""
        SELECT COUNT(DISTINCT g.game_id)
        FROM games g
        INNER JOIN (
            SELECT home_team, away_team, CAST(commence_time AS DATE) as game_date
            FROM odds WHERE total_line IS NOT NULL AND home_ml IS NOT NULL
            GROUP BY 1,2,3
        ) o ON g.home_team = o.home_team AND g.away_team = o.away_team AND g.date = o.game_date
        WHERE NOT g.went_extras AND NOT g.mercy_rule
    """).fetchone()[0]
    print(f"Games with matched odds (regulation, no mercy): {matched}")

    con.close()


def analyze():
    """Run structural analysis on historical data."""
    con = duckdb.connect(str(DB_PATH))

    # Basic counts
    n_games = con.execute("SELECT COUNT(*) FROM games").fetchone()[0]
    n_odds = con.execute("SELECT COUNT(DISTINCT game_id) FROM odds").fetchone()[0]

    print(f"=== COLLEGE BASEBALL STRUCTURAL ANALYSIS ===")
    print(f"Total games: {n_games}")
    print(f"Games with odds: {n_odds}")

    if n_games == 0:
        con.close()
        print("No games found. Run --pull-scores first.")
        return

    # Structural stats
    stats = con.execute("""
        SELECT
            COUNT(*) as total,
            SUM(CASE WHEN went_extras THEN 1 ELSE 0 END) as extras,
            SUM(CASE WHEN mercy_rule THEN 1 ELSE 0 END) as mercy,
            AVG(total_runs) as avg_total,
            PERCENTILE_CONT(0.25) WITHIN GROUP (ORDER BY total_runs) as total_p25,
            PERCENTILE_CONT(0.50) WITHIN GROUP (ORDER BY total_runs) as total_p50,
            PERCENTILE_CONT(0.75) WITHIN GROUP (ORDER BY total_runs) as total_p75,
            AVG(CASE WHEN home_won THEN 1.0 ELSE 0.0 END) as home_win_pct,
            AVG(ABS(home_margin)) as avg_margin
        FROM games
        WHERE NOT mercy_rule AND NOT went_extras
    """).fetchone()

    total, extras, mercy, avg_total, p25, p50, p75, home_pct, avg_margin = stats
    print(f"\nRegulation games (no extras, no mercy): {total}")
    print(f"Extras rate: {extras}/{n_games} = {extras/n_games*100:.1f}%")
    print(f"Mercy rate: {mercy}/{n_games} = {mercy/n_games*100:.1f}%")
    print(f"Home win %: {home_pct*100:.1f}%")
    print(f"Avg total: {avg_total:.1f} (P25={p25:.0f}, P50={p50:.0f}, P75={p75:.0f})")
    print(f"Avg abs margin: {avg_margin:.1f}")

    # Correlation analysis: ML outcome vs Over/Under by total bucket
    print(f"\n=== ML + TOTAL CORRELATION ===")
    print("(Regulation games only, joined to odds for total_line)")

    # Match games to odds consensus total
    matched = con.execute("""
        WITH consensus AS (
            SELECT
                home_team, away_team,
                CAST(commence_time AS DATE) as game_date,
                MEDIAN(total_line) as total_line
            FROM odds
            WHERE total_line IS NOT NULL
            GROUP BY home_team, away_team, CAST(commence_time AS DATE)
        )
        SELECT
            g.game_id,
            g.home_team, g.away_team,
            g.home_score, g.away_score,
            g.total_runs, g.home_margin, g.home_won,
            c.total_line,
            CASE
                WHEN c.total_line < 12 THEN 'low (<12)'
                WHEN c.total_line <= 13.5 THEN 'medium (12-13.5)'
                ELSE 'high (>13.5)'
            END as total_bucket,
            CASE WHEN g.total_runs > c.total_line THEN 1 ELSE 0 END as went_over,
            CASE WHEN g.home_won THEN 1 ELSE 0 END as home_won_int
        FROM games g
        INNER JOIN consensus c
            ON g.home_team = c.home_team AND g.away_team = c.away_team
            AND g.date = c.game_date
        WHERE NOT g.went_extras AND NOT g.mercy_rule
            AND c.total_line IS NOT NULL
    """).fetchall()

    if not matched:
        print("No matched games found. Run --pull-odds first.")
        con.close()
        return

    print(f"Matched games: {len(matched)}")

    # Compute correlation factors by bucket
    # CF = P(away ML + over) / (P(away ML) * P(over))
    buckets = {}
    for row in matched:
        bucket = row[9]  # total_bucket
        if bucket not in buckets:
            buckets[bucket] = {"n": 0, "away_win": 0, "over": 0,
                               "away_over": 0, "away_under": 0,
                               "home_over": 0, "home_under": 0}
        b = buckets[bucket]
        b["n"] += 1
        away_won = not row[7]  # home_won
        went_over = row[10]

        if away_won:
            b["away_win"] += 1
        if went_over:
            b["over"] += 1
        if away_won and went_over:
            b["away_over"] += 1
        if away_won and not went_over:
            b["away_under"] += 1
        if not away_won and went_over:
            b["home_over"] += 1
        if not away_won and not went_over:
            b["home_under"] += 1

    print(f"\n{'Bucket':<20} {'N':>5} {'CF(Away+Over)':>14} {'CF(Away+Under)':>15} {'CF(Home+Over)':>14} {'CF(Home+Under)':>15}")
    print("-" * 90)
    for bucket_name in ["low (<12)", "medium (12-13.5)", "high (>13.5)"]:
        if bucket_name not in buckets:
            continue
        b = buckets[bucket_name]
        n = b["n"]
        p_away = b["away_win"] / n
        p_over = b["over"] / n
        p_home = 1 - p_away
        p_under = 1 - p_over

        # Correlation factors
        cf_away_over = (b["away_over"] / n) / (p_away * p_over) if p_away * p_over > 0 else 0
        cf_away_under = (b["away_under"] / n) / (p_away * p_under) if p_away * p_under > 0 else 0
        cf_home_over = (b["home_over"] / n) / (p_home * p_over) if p_home * p_over > 0 else 0
        cf_home_under = (b["home_under"] / n) / (p_home * p_under) if p_home * p_under > 0 else 0

        print(f"{bucket_name:<20} {n:>5} {cf_away_over:>14.3f} {cf_away_under:>15.3f} {cf_home_over:>14.3f} {cf_home_under:>15.3f}")

    con.close()


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="College Baseball Historical Data Acquisition & Analysis"
    )
    parser.add_argument("--pull-scores", nargs="+", type=int, metavar="YEAR",
                        help="Pull ESPN scores for given seasons")
    parser.add_argument("--pull-odds", nargs="+", type=int, metavar="YEAR",
                        help="Pull Odds API closing lines for given seasons")
    parser.add_argument("--pull-linescores", nargs="+", type=int, metavar="YEAR",
                        help="Pull ESPN inning-by-inning linescores for given seasons")
    parser.add_argument("--match-teams", action="store_true",
                        help="Normalize ESPN team names to match Odds API")
    parser.add_argument("--analyze", action="store_true",
                        help="Run structural analysis on historical data")

    args = parser.parse_args()

    if not any([args.pull_scores, args.pull_odds, args.pull_linescores, args.match_teams, args.analyze]):
        parser.print_help()
        return

    if args.pull_scores:
        print(f"Pulling ESPN scores for seasons: {args.pull_scores}")
        pull_scores(args.pull_scores)

    if args.pull_odds:
        print(f"\nPulling Odds API closing lines for seasons: {args.pull_odds}")
        pull_odds(args.pull_odds)

    if args.pull_linescores:
        print(f"\nPulling ESPN linescores for seasons: {args.pull_linescores}")
        pull_linescores(args.pull_linescores)

    if args.match_teams:
        print("\nMatching team names...")
        match_teams()

    if args.analyze:
        print()
        analyze()


if __name__ == "__main__":
    main()
