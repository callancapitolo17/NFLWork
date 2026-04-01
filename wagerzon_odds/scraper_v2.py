#!/usr/bin/env python3
"""
Wagerzon Odds Scraper v2 - REST API
Fetches odds via NewScheduleHelper.aspx JSON endpoint using requests.Session.
No browser required — auth via ASP.NET form POST, data via JSON API.
"""

import os
import re
import sys
import duckdb
import requests
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from dotenv import load_dotenv

from config import get_sport_config, WAGERZON_BASE_URL, WAGERZON_HELPER_URL
from team_mapping import normalize_team_name

# Add Answer Keys to path for shared team name resolution
sys.path.insert(0, str(Path(__file__).parent.parent / "Answer Keys"))
from canonical_match import load_team_dict, load_canonical_games, resolve_team_names

# Load environment
env_path = Path(__file__).parent.parent / "bet_logger" / ".env"
load_dotenv(env_path)

WAGERZON_USERNAME = os.getenv("WAGERZON_USERNAME")
WAGERZON_PASSWORD = os.getenv("WAGERZON_PASSWORD")

DB_PATH = Path(__file__).parent / "wagerzon.duckdb"


# =============================================================================
# HELPERS
# =============================================================================


def safe_float(val) -> Optional[float]:
    """Convert string to float, returning None for empty/invalid values."""
    if not val:
        return None
    try:
        return float(val)
    except (ValueError, TypeError):
        return None


def safe_int(val) -> Optional[int]:
    """Convert string to int, returning None for empty/invalid values."""
    if not val:
        return None
    try:
        return int(val)
    except (ValueError, TypeError):
        return None


# =============================================================================
# AUTH
# =============================================================================


def login(session: requests.Session):
    """Login to Wagerzon via ASP.NET form POST.

    1. GET the login page to capture __VIEWSTATE and other hidden fields
    2. POST credentials with the hidden fields
    3. Session cookie (ASP.NET_SessionId) maintains auth for subsequent requests
    """
    resp = session.get(WAGERZON_BASE_URL, timeout=15)
    resp.raise_for_status()

    # If already redirected to schedule, session is still valid
    if "NewSchedule" in resp.url:
        print("Already authenticated")
        return

    html = resp.text

    # Extract ASP.NET hidden fields from login form
    fields = {}
    for name in ["__VIEWSTATE", "__VIEWSTATEGENERATOR", "__EVENTVALIDATION",
                 "__EVENTTARGET", "__EVENTARGUMENT"]:
        match = re.search(rf'(?:name|id)="{name}"[^>]*value="([^"]*)"', html)
        if match:
            fields[name] = match.group(1)

    if "__VIEWSTATE" not in fields:
        raise RuntimeError("Could not find __VIEWSTATE on login page — page structure may have changed")

    fields["Account"] = WAGERZON_USERNAME
    fields["Password"] = WAGERZON_PASSWORD
    fields["BtnSubmit"] = ""

    resp = session.post(WAGERZON_BASE_URL, data=fields, timeout=15)
    resp.raise_for_status()
    print("Logged in successfully")


# =============================================================================
# API CLIENT
# =============================================================================


def fetch_odds_json(session: requests.Session, sport: str) -> dict:
    """Fetch odds JSON from NewScheduleHelper.aspx."""
    config = get_sport_config(sport)
    url = f"{WAGERZON_HELPER_URL}?WT=0&{config['url_params']}"

    resp = session.get(url, timeout=30, headers={
        "Accept": "application/json, text/plain, */*",
        "X-Requested-With": "XMLHttpRequest",
    })
    resp.raise_for_status()
    return resp.json()


# =============================================================================
# JSON PARSING
# =============================================================================


def parse_game_line(line: dict, game_id: str, period: str, market: str,
                    base: dict) -> Optional[dict]:
    """Parse a GameLine object (spread + total + ML) into a DuckDB record."""
    away_spread = safe_float(line.get("vsprdt"))
    away_spread_price = safe_int(line.get("vsprdoddst"))
    home_spread = safe_float(line.get("hsprdt"))
    home_spread_price = safe_int(line.get("hsprdoddst"))

    # unt is the positive total value; ovt is negative (same number)
    total = safe_float(line.get("unt"))
    over_price = safe_int(line.get("ovoddst"))
    under_price = safe_int(line.get("unoddst"))

    away_ml = safe_int(line.get("voddst"))
    home_ml = safe_int(line.get("hoddst"))

    if all(v is None for v in [away_spread, total, away_ml]):
        return None

    return {
        **base,
        "game_id": game_id,
        "market": market,
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
    }


def parse_team_total(line: dict, game_id: str, period: str, market: str,
                     base: dict) -> Optional[dict]:
    """Parse a team total GameLine (over/under only) into a DuckDB record."""
    total = safe_float(line.get("unt"))
    over_price = safe_int(line.get("ovoddst"))
    under_price = safe_int(line.get("unoddst"))

    if total is None:
        return None

    return {
        **base,
        "game_id": game_id,
        "market": market,
        "period": period,
        "away_spread": None,
        "away_spread_price": None,
        "home_spread": None,
        "home_spread_price": None,
        "total": total,
        "over_price": over_price,
        "under_price": under_price,
        "away_ml": None,
        "home_ml": None,
    }


def parse_moneyline_only(line: dict, game_id: str, period: str, market: str,
                         base: dict) -> Optional[dict]:
    """Parse a moneyline-only GameLine (no spread/total) into a DuckDB record.

    Used for prop markets like score-first, odd/even, score-in-1st-inning where
    the only data is two moneyline prices (away_ml / home_ml).
    """
    away_ml = safe_int(line.get("voddst"))
    home_ml = safe_int(line.get("hoddst"))

    if away_ml is None and home_ml is None:
        return None

    return {
        **base,
        "game_id": game_id,
        "market": market,
        "period": period,
        "away_spread": None,
        "away_spread_price": None,
        "home_spread": None,
        "home_spread_price": None,
        "total": None,
        "over_price": None,
        "under_price": None,
        "away_ml": away_ml,
        "home_ml": home_ml,
    }


def parse_odds(data: dict, sport: str) -> list[dict]:
    """Parse NewScheduleHelper JSON response into DuckDB records.

    Processes all leagues and extracts derivatives from each game's GameChilds:
      - idgmtyp 10: Full game (parent) — spread + total + ML
      - idgmtyp 15: First half — spread + total + ML
      - idgmtyp 19: Hits / H+R+E totals — total only
      - idgmtyp 25: Alt lines OR period totals (3 INN / 7 INN) — spread/total
      - idgmtyp 30: Pitcher props (outs, hits allowed, walks) — total only
      - idgmtyp 31: Odd/even total runs — ML only
      - idgmtyp 35: Team total (full game) — total only
      - idgmtyp 44: Score first / score first wins game — ML only
      - idgmtyp 47: Score in 1st inning (yes/no) — ML only
      - idgmtyp 66: Team total (1H) — total only
    """
    config = get_sport_config(sport)
    sport_key = config["sport_key"]
    fetch_time = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")

    team_dict = load_team_dict(sport) if sport != "nfl" else {}
    canonical_games = load_canonical_games(sport) if sport != "nfl" else []

    records = []

    leagues = data.get("result", {}).get("listLeagues", [[]])[0]
    if not leagues:
        print("No leagues in response")
        return records

    # Process all leagues — each may have a different period (FG, F5, H1)
    all_games = []
    for parent_league in leagues:
        league_desc = parent_league.get("Description", "")
        lg_games = parent_league.get("Games", [])

        # Detect period from league description
        desc_lower = league_desc.lower()
        if "1st 5" in desc_lower or "first 5" in desc_lower or "f5" in desc_lower:
            league_period = "F5"
        elif "1st half" in desc_lower or "first half" in desc_lower or "1h" in desc_lower:
            league_period = "Half1"
        else:
            league_period = "fg"

        if lg_games:
            print(f"Found {len(lg_games)} games in '{league_desc}' (period={league_period})")
        for g in lg_games:
            g["_league_period"] = league_period
            all_games.append(g)

    games = all_games

    for game in games:
        if not game.get("GameLines"):
            continue

        away_raw = game["vtm"]
        home_raw = game["htm"]

        # Format date: YYYYMMDD -> MM/DD
        gmdt = game.get("gmdt", "")
        game_date = f"{gmdt[4:6]}/{gmdt[6:8]}" if len(gmdt) == 8 else ""
        game_time = game.get("gmtm", "")[:5]  # "18:30:00" -> "18:30"

        # Resolve team names
        if team_dict or canonical_games:
            away_team, home_team = resolve_team_names(
                away_raw, home_raw, team_dict, canonical_games
            )
        else:
            away_team = normalize_team_name(away_raw, sport)
            home_team = normalize_team_name(home_raw, sport)

        away_rot = str(game["vnum"])
        home_rot = str(game["hnum"])
        game_id = f"{away_rot}-{home_rot}"

        base = {
            "fetch_time": fetch_time,
            "sport_key": sport_key,
            "game_date": game_date,
            "game_time": game_time,
            "away_team": away_team,
            "home_team": home_team,
        }

        # --- Full game line (parent) ---
        league_period = game.get("_league_period", "fg")
        line = game["GameLines"][0]
        rec = parse_game_line(line, game_id, league_period, "spreads", base)
        if rec:
            records.append(rec)
            print(f"  spreads: {away_team} @ {home_team} | "
                  f"{rec['away_spread']}/{rec['home_spread']} | {rec['total']}")

        # --- Child games (derivatives) ---
        alt_counter = 0
        for child in game.get("GameChilds", []):
            child_type = child.get("idgmtyp")
            child_lines = child.get("GameLines", [])
            if not child_lines:
                continue
            child_line = child_lines[0]
            child_vtm = child.get("vtm", "")
            child_vtm_upper = child_vtm.upper()

            # Helper: determine if a stripped child name belongs to the away
            # or home team.  Wagerzon child names drop the city abbreviation
            # (e.g. "GUARDIANS TEAM TOTAL" for parent "CLE GUARDIANS"), so we
            # check containment rather than exact equality.
            def _side(stripped_name: str) -> str:
                s = stripped_name.upper()
                if s in away_raw.upper():
                    return "away"
                return "home"

            if child_type == 15:
                # First half line (F5 in baseball) — spread + total + ML
                cid = f"{child['vnum']}-{child['hnum']}"
                rec = parse_game_line(child_line, cid, "h1", "spreads_h1", base)
                if rec:
                    records.append(rec)

            elif child_type == 19:
                # Hits totals or H+R+E totals
                # "GUARDIANS TOTAL HITS" → team hits total
                # "H+R+E (SF/SD)" → game-level hits+runs+errors total
                if "H+R+E" in child_vtm_upper:
                    rec = parse_team_total(
                        child_line, f"{game_id}-hre", "fg", "hre_total", base
                    )
                elif "TOTAL HITS" in child_vtm_upper:
                    # Wagerzon labels this with the away team name but it's
                    # the game total hits (one per game, ~13-16 range)
                    rec = parse_team_total(
                        child_line, f"{game_id}-totalhits", "fg",
                        "total_hits", base
                    )
                else:
                    rec = None
                if rec:
                    records.append(rec)

            elif child_type == 25:
                # Two sub-types share idgmtyp=25:
                #   1) Period totals: "3 INN TEAM" or "7 INN TEAM" — total only
                #   2) Alt lines: "TEAM ALT" — spread and/or total
                period_match = re.match(r"(\d+)\s+INN\s+", child_vtm_upper)
                if period_match:
                    # Period total (e.g. first 3 innings, first 7 innings)
                    innings = period_match.group(1)
                    rec = parse_team_total(
                        child_line, f"{game_id}-{innings}inn", f"f{innings}",
                        f"totals_f{innings}", base
                    )
                    if rec:
                        records.append(rec)
                else:
                    # Alt line — spread and total are paired in each entry
                    child_id = str(child.get("idgm", alt_counter))

                    # Alt spread
                    alt_away_spread = safe_float(child_line.get("vsprdt"))
                    if alt_away_spread is not None:
                        records.append({
                            **base,
                            "game_id": f"{game_id}-alts-{child_id}",
                            "market": "alternate_spreads_fg",
                            "period": "fg",
                            "away_spread": alt_away_spread,
                            "away_spread_price": safe_int(child_line.get("vsprdoddst")),
                            "home_spread": safe_float(child_line.get("hsprdt")),
                            "home_spread_price": safe_int(child_line.get("hsprdoddst")),
                            "total": None,
                            "over_price": None,
                            "under_price": None,
                            "away_ml": None,
                            "home_ml": None,
                        })

                    # Alt total
                    alt_total = safe_float(child_line.get("unt"))
                    if alt_total is not None:
                        records.append({
                            **base,
                            "game_id": f"{game_id}-altt-{child_id}",
                            "market": "alternate_totals_fg",
                            "period": "fg",
                            "away_spread": None,
                            "away_spread_price": None,
                            "home_spread": None,
                            "home_spread_price": None,
                            "total": alt_total,
                            "over_price": safe_int(child_line.get("ovoddst")),
                            "under_price": safe_int(child_line.get("unoddst")),
                            "away_ml": None,
                            "home_ml": None,
                        })

                alt_counter += 1

            elif child_type == 30:
                # Pitcher props — "L WEBB (SF) TOTAL OUTS", "L WEBB (SF) HITS ALLOWED"
                # Extract pitcher name and prop type from the child team name
                pitcher_match = re.match(
                    r"(.+?)\s*\([A-Z]+\)\s+(TOTAL OUTS|HITS ALLOWED|WALKS ALLOWED|STRIKEOUTS|EARNED RUNS)",
                    child_vtm_upper
                )
                if pitcher_match:
                    pitcher_name = pitcher_match.group(1).strip()
                    prop_type = pitcher_match.group(2).lower().replace(" ", "_")
                    # Use pitcher name in the game_id so each prop is unique
                    pitcher_slug = pitcher_name.lower().replace(" ", "_").replace(".", "")
                    rec = parse_team_total(
                        child_line,
                        f"{game_id}-pp-{pitcher_slug}-{prop_type}",
                        "fg",
                        f"pitcher_{prop_type}",
                        {**base, "away_team": pitcher_name, "home_team": f"{away_team} @ {home_team}"},
                    )
                    if rec:
                        records.append(rec)

            elif child_type == 31:
                # Odd/even total runs — ML only
                # "TOTAL RUNS ODD(CLE/LAD)" vs "TOTAL RUNS EVEN(CLE/LAD)"
                rec = parse_moneyline_only(
                    child_line, f"{game_id}-oddeven", "fg", "odd_even_runs", base
                )
                if rec:
                    records.append(rec)

            elif child_type == 35:
                # Full game team total
                tt_name = child_vtm_upper.replace(" TEAM TOTAL", "").strip()
                tt_side = _side(tt_name)
                rec = parse_team_total(
                    child_line, f"{game_id}-tt-{tt_side}", "fg",
                    f"team_totals_{tt_side}_fg", base
                )
                if rec:
                    records.append(rec)

            elif child_type == 44:
                # Two variants:
                #   "CLE GUARDIANS SC 1ST" — which team scores first (ML)
                #   "YES TM SCR 1ST WIN G(CLE/LAD)" — does the team that scores first win?
                if "SC 1ST" in child_vtm_upper and "WIN" not in child_vtm_upper:
                    rec = parse_moneyline_only(
                        child_line, f"{game_id}-scorefirst", "fg", "score_first", base
                    )
                else:
                    rec = parse_moneyline_only(
                        child_line, f"{game_id}-scorefirst-wins", "fg",
                        "score_first_wins_game", base
                    )
                if rec:
                    records.append(rec)

            elif child_type == 47:
                # Score in 1st inning — yes/no ML
                rec = parse_moneyline_only(
                    child_line, f"{game_id}-score1stinn", "fg",
                    "score_1st_inning", base
                )
                if rec:
                    records.append(rec)

            elif child_type == 66:
                # 1H team total
                tt_name = child_vtm_upper
                tt_name = tt_name.replace("1H ", "").replace(" TEAM TOTAL", "").strip()
                tt_side = _side(tt_name)
                rec = parse_team_total(
                    child_line, f"{game_id}-tt-{tt_side}-h1", "h1",
                    f"team_totals_{tt_side}_h1", base
                )
                if rec:
                    records.append(rec)

    # --- Race-to-10 leagues (separate from main games) ---
    # lg=1852 comes as a standalone league with idgmtyp=47 moneyline-only games.
    # Team names embed the prop: "DUKE GET 10PTS 1ST" or "DUKE GET 10 PTS 1ST"
    race_pattern = re.compile(r"(.+?)\s+GET\s+(\d+)\s*PTS\s+1ST", re.IGNORECASE)
    for league in leagues:
        desc = (league.get("Description", "") or "").upper()
        if "SCORE FIRST" not in desc:
            continue
        race_games = league.get("Games", [])
        print(f"Found {len(race_games)} race-to-X games in '{league.get('Description', '')}'")
        for game in race_games:
            if not game.get("GameLines"):
                continue
            line = game["GameLines"][0]
            away_odds = safe_int(line.get("voddst"))
            home_odds = safe_int(line.get("hoddst"))
            if away_odds is None or home_odds is None:
                continue

            # Strip prop suffix from team names and extract threshold
            away_raw = game["vtm"]
            home_raw = game["htm"]
            away_match = race_pattern.match(away_raw)
            home_match = race_pattern.match(home_raw)
            if not away_match or not home_match:
                continue
            away_clean = away_match.group(1).strip()
            home_clean = home_match.group(1).strip()
            threshold = away_match.group(2)

            # Resolve canonical names
            if team_dict or canonical_games:
                away_team, home_team = resolve_team_names(
                    away_clean, home_clean, team_dict, canonical_games
                )
            else:
                away_team = normalize_team_name(away_clean, sport)
                home_team = normalize_team_name(home_clean, sport)

            gmdt = game.get("gmdt", "")
            game_date = f"{gmdt[4:6]}/{gmdt[6:8]}" if len(gmdt) == 8 else ""
            game_time = game.get("gmtm", "")[:5]
            away_rot = str(game["vnum"])
            home_rot = str(game["hnum"])

            market_name = f"race_to_{threshold}"
            records.append({
                "fetch_time": fetch_time,
                "sport_key": sport_key,
                "game_id": f"{away_rot}-{home_rot}",
                "game_date": game_date,
                "game_time": game_time,
                "away_team": away_team,
                "home_team": home_team,
                "market": market_name,
                "period": "fg",
                "away_spread": None,
                "away_spread_price": None,
                "home_spread": None,
                "home_spread_price": None,
                "total": None,
                "over_price": None,
                "under_price": None,
                "away_ml": away_odds,
                "home_ml": home_odds,
            })
            print(f"  {market_name}: {away_team} @ {home_team} | {away_odds}/{home_odds}")

    return records


# =============================================================================
# DATABASE
# =============================================================================


def init_database(sport: str):
    """Initialize DuckDB with the odds table for a sport."""
    config = get_sport_config(sport)
    table_name = config["table_name"]

    conn = duckdb.connect(str(DB_PATH))
    conn.execute(f"""
        CREATE TABLE IF NOT EXISTS {table_name} (
            fetch_time TIMESTAMP,
            sport_key VARCHAR,
            game_id VARCHAR,
            game_date VARCHAR,
            game_time VARCHAR,
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
    conn.close()


def save_odds(odds_data: list[dict], sport: str):
    """Save odds to DuckDB (replaces previous scrape)."""
    if not odds_data:
        print("No odds data to save")
        return

    config = get_sport_config(sport)
    table_name = config["table_name"]

    conn = duckdb.connect(str(DB_PATH))

    columns = [
        "fetch_time", "sport_key", "game_id", "game_date", "game_time",
        "away_team", "home_team", "market", "period",
        "away_spread", "away_spread_price", "home_spread", "home_spread_price",
        "total", "over_price", "under_price", "away_ml", "home_ml"
    ]

    placeholders = ", ".join(["?" for _ in columns])

    # Clear old data before inserting fresh scrape
    conn.execute(f"DELETE FROM {table_name}")

    conn.executemany(f"""
        INSERT INTO {table_name} ({", ".join(columns)})
        VALUES ({placeholders})
    """, [
        tuple(d[col] for col in columns)
        for d in odds_data
    ])

    result = conn.execute(f"SELECT COUNT(*) FROM {table_name}").fetchone()
    print(f"\nDatabase now has {result[0]} total records in {table_name}")

    conn.close()


# =============================================================================
# MAIN
# =============================================================================


def scrape_wagerzon(sport: str, headless: bool = True):
    """Scrape Wagerzon odds via REST API and save to DuckDB.

    The headless parameter is kept for backward compatibility but is ignored —
    no browser is needed.
    """
    if not WAGERZON_USERNAME or not WAGERZON_PASSWORD:
        raise ValueError("WAGERZON_USERNAME and WAGERZON_PASSWORD must be set in bet_logger/.env")

    init_database(sport)

    session = requests.Session()

    # Step 1: Authenticate
    print("Logging in to Wagerzon...")
    login(session)

    # Step 2: Fetch odds JSON (main leagues)
    print(f"Fetching {sport.upper()} odds...")
    data = fetch_odds_json(session, sport)

    # Step 3: Parse JSON into records
    odds_data = parse_odds(data, sport)

    # Step 3b: Fetch prop leagues (race_to_10, etc.) — separate API calls
    config = get_sport_config(sport)
    for prop_name, prop_params in config.get("prop_params", {}).items():
        print(f"Fetching {sport.upper()} {prop_name} props...")
        prop_url = f"{WAGERZON_HELPER_URL}?WT=0&{prop_params}"
        try:
            resp = session.get(prop_url, timeout=30, headers={
                "Accept": "application/json, text/plain, */*",
                "X-Requested-With": "XMLHttpRequest",
            })
            resp.raise_for_status()
            prop_data = resp.json()
            prop_records = parse_odds(prop_data, sport)
            odds_data.extend(prop_records)
        except Exception as e:
            print(f"  Warning: failed to fetch {prop_name}: {e}")

    # Count by market type
    market_counts = {}
    for rec in odds_data:
        m = rec["market"]
        market_counts[m] = market_counts.get(m, 0) + 1
    for m, c in sorted(market_counts.items()):
        print(f"  {m}: {c} records")

    # Step 4: Save to DuckDB
    if odds_data:
        save_odds(odds_data, sport)
        print(f"\nSaved {len(odds_data)} total records for {sport.upper()}")
    else:
        print("No odds data scraped")

    return odds_data


if __name__ == "__main__":
    sport = sys.argv[1] if len(sys.argv) > 1 else "nfl"
    print(f"Starting Wagerzon {sport.upper()} odds scraper (REST API)...")
    scrape_wagerzon(sport)
