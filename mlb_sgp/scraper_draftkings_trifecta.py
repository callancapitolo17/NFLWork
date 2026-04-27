"""
DraftKings Trifecta SGP Scraper

Reads a JSON request file with one entry per (game_id, prop_type, side, legs),
resolves each leg into a DK selection ID via dk_leg_resolvers.LEG_RESOLVERS,
posts the combined SGP to DK's calculateBets endpoint, and writes one row
per request to mlb_trifecta_sgp_odds in the project's MLB DuckDB.

Usage:
    python3 scraper_draftkings_trifecta.py \\
        --input /tmp/trifecta_requests.json \\
        --db    "Answer Keys/mlb.duckdb"

Input JSON shape (one element per request):
    [
      {
        "game_id":   "0438fd942a9f1f07557cef989b1a4e4d",
        "home_team": "San Francisco Giants",
        "away_team": "Los Angeles Dodgers",
        "prop_type": "TRIPLE-PLAY",
        "side":      "home",
        "legs":      [
          {"type": "scores_first"},
          {"type": "wins_period", "period": "F5"},
          {"type": "wins_period", "period": "FG"}
        ]
      },
      ...
    ]

Exit codes:
    0  - success (some/all rows may have NULL sgp_decimal if DK couldn't price them)
    1  - usage error (bad CLI args or missing input file)
    2  - DK API failure (Akamai block, network, etc.)
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import duckdb

# Reuse the production DK SGP scraper's session + event-discovery functions
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from scraper_draftkings_sgp import init_session, fetch_dk_events  # noqa: E402
from dk_leg_resolvers import resolve_legs, find_market_by_name      # noqa: E402

# DK API URLs — same prefix the production scraper uses
SGP_EVENTS_URL_TPL = (
    "https://sportsbook-nash.draftkings.com/sites/US-SB/api/sportscontent/"
    "parlays/v1/sgp/events/{event_id}"
)
DK_CALCULATE_BETS_URL = "https://gaming-us-nj.draftkings.com/api/wager/v1/calculateBets"

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("scraper_draftkings_trifecta")


def map_odds_id_to_dk_event(events: list[dict], home_team: str, away_team: str) -> Optional[str]:
    """Find the DK event whose teams match the requested home/away. DK uses
    e.g. 'CHI Cubs' / 'SD Padres'; we use Odds API canonical names. We match
    by checking if the canonical team's last word is in DK's team string
    (e.g. 'San Francisco Giants' matches 'SF Giants' on 'Giants')."""
    home_last = home_team.split()[-1].lower()
    away_last = away_team.split()[-1].lower()
    for evt in events:
        dk_home = (evt.get('dk_home') or '').lower()
        dk_away = (evt.get('dk_away') or '').lower()
        if home_last in dk_home and away_last in dk_away:
            return str(evt['dk_event_id'])
    return None


def fetch_event_state(session, event_id: str) -> Optional[dict]:
    """Fetch the full SGP event payload (~1-6 MB)."""
    url = SGP_EVENTS_URL_TPL.format(event_id=event_id)
    try:
        r = session.get(url, timeout=30)
        r.raise_for_status()
        return r.json()
    except Exception as e:
        log.error("Failed to fetch event %s: %s", event_id, e)
        return None


def extract_team_names_from_event(event_state: dict) -> dict:
    """Walk the event payload to extract the DK team-name strings used in
    selection labels (e.g. 'CHI Cubs', 'SD Padres'). Picks them from the
    Moneyline market which always exists and always has both teams.

    DK Moneyline selection IDs end in _3 (home) or _1 (away), matching the
    convention used by the production DK SGP scraper.
    """
    market = find_market_by_name(event_state, 'Moneyline')
    if market is None:
        return {'home': '', 'away': ''}
    sels = market.get('selections', [])
    if len(sels) < 2:
        return {'home': '', 'away': ''}
    home_name = ''
    away_name = ''
    for sel in sels:
        sid = sel.get('id') or sel.get('selectionId') or ''
        nm = sel.get('name') or ''
        if sid.endswith('_3'):
            home_name = nm
        elif sid.endswith('_1'):
            away_name = nm
    return {'home': home_name, 'away': away_name}


def post_calculate_bets(session, sel_ids: list[str]) -> Optional[float]:
    """Post the combined SGP to DK's calculateBets endpoint. Returns the
    decimal odds, or None if DK rejected/errored.

    NOTE: This payload format was discovered in Task 1's combinability test.
    The DK calculateBets endpoint requires `selectionsForYourBet` (not the
    naive `bets/wagerOption/StraightSGP` wrapper). All legs share
    yourBetGroup=0 to mark them as one SGP.
    """
    body = {
        "selections": [],
        "selectionsForYourBet": [
            {"id": sid, "yourBetGroup": 0}
            for sid in sel_ids
        ],
        "selectionsForCombinator": [],
        "selectionsForProgressiveParlay": [],
        "oddsStyle": "american",
    }
    try:
        r = session.post(DK_CALCULATE_BETS_URL, json=body, timeout=30)
        if r.status_code != 200:
            log.warning("calculateBets returned %d: %s", r.status_code, r.text[:200])
            return None
        data = r.json()
        bets_arr = data.get('bets') or []
        if not bets_arr:
            log.warning("calculateBets returned no bets: %s", str(data)[:300])
            return None
        bet = bets_arr[0]
        # `trueOdds` is the canonical decimal odds field
        true_odds = bet.get('trueOdds')
        if true_odds is None:
            log.warning("No trueOdds in calculateBets response: %s", str(bet)[:300])
            return None
        return float(true_odds)
    except Exception as e:
        log.error("calculateBets POST failed: %s", e)
        return None


def decimal_to_american(dec: float) -> int:
    if dec >= 2.0:
        return int(round((dec - 1) * 100))
    return int(round(-100 / (dec - 1)))


def write_rows(con: duckdb.DuckDBPyConnection, rows: list[dict]) -> None:
    """Write all rows under one fetch_time so MAX(fetch_time) returns a
    consistent atomic snapshot.
    CREATE IF NOT EXISTS is idempotent — Plan #1 already created the table,
    but this protects against fresh checkouts and test fixtures."""
    if not rows:
        log.info("No rows to write.")
        return
    con.execute("""
        CREATE TABLE IF NOT EXISTS mlb_trifecta_sgp_odds (
            fetch_time     TIMESTAMP,
            game_id        VARCHAR,
            prop_type      VARCHAR,
            side           VARCHAR,
            legs_json      VARCHAR,
            selection_ids  VARCHAR,
            sgp_decimal    DOUBLE,
            sgp_american   INTEGER,
            source         VARCHAR
        )
    """)
    con.executemany(
        "INSERT INTO mlb_trifecta_sgp_odds VALUES (?,?,?,?,?,?,?,?,?)",
        [
            (
                r['fetch_time'], r['game_id'], r['prop_type'], r['side'],
                r['legs_json'], r['selection_ids'],
                r['sgp_decimal'], r['sgp_american'], r['source'],
            )
            for r in rows
        ],
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input', required=True, help='Path to trifecta requests JSON')
    ap.add_argument('--db', required=True, help='Path to mlb.duckdb')
    args = ap.parse_args()

    input_path = Path(args.input)
    db_path = Path(args.db)
    if not input_path.exists():
        log.error("Input file not found: %s", input_path)
        return 1

    requests_in = json.loads(input_path.read_text())
    log.info("Loaded %d trifecta requests", len(requests_in))

    session = init_session()
    events = fetch_dk_events(session)
    log.info("DK eventgroup returned %d MLB events", len(events))

    # Cache event-state fetches and team-name extractions per DK event id
    event_state_cache: dict[str, dict] = {}
    team_names_cache: dict[str, dict] = {}

    fetch_time = datetime.now(timezone.utc).replace(microsecond=0)
    rows = []

    for req in requests_in:
        game_id   = req['game_id']
        home_team = req['home_team']
        away_team = req['away_team']
        prop_type = req['prop_type']
        side      = req['side']
        legs      = req['legs']

        dk_event_id = map_odds_id_to_dk_event(events, home_team, away_team)
        if dk_event_id is None:
            log.warning("No DK event for %s vs %s — skipping (game_id=%s prop=%s side=%s)",
                        home_team, away_team, game_id, prop_type, side)
            continue

        if dk_event_id not in event_state_cache:
            state = fetch_event_state(session, dk_event_id)
            if state is None:
                log.warning("Could not fetch DK event %s — skipping its requests", dk_event_id)
                event_state_cache[dk_event_id] = None
                continue
            event_state_cache[dk_event_id] = state
            team_names_cache[dk_event_id] = extract_team_names_from_event(state)

        if event_state_cache[dk_event_id] is None:
            continue

        team_names = team_names_cache[dk_event_id]
        sel_ids = resolve_legs(legs, side, event_state_cache[dk_event_id], team_names)
        if sel_ids is None:
            log.info("Could not resolve all legs for game=%s prop=%s side=%s — writing NULL row",
                     game_id, prop_type, side)
            rows.append({
                'fetch_time': fetch_time, 'game_id': game_id,
                'prop_type': prop_type, 'side': side,
                'legs_json': json.dumps(legs), 'selection_ids': '',
                'sgp_decimal': None, 'sgp_american': None,
                'source': 'draftkings_direct',
            })
            continue

        sgp_decimal = post_calculate_bets(session, sel_ids)
        sgp_american = decimal_to_american(sgp_decimal) if sgp_decimal else None
        rows.append({
            'fetch_time': fetch_time, 'game_id': game_id,
            'prop_type': prop_type, 'side': side,
            'legs_json': json.dumps(legs), 'selection_ids': ','.join(sel_ids),
            'sgp_decimal': sgp_decimal, 'sgp_american': sgp_american,
            'source': 'draftkings_direct',
        })
        log.info("game=%s prop=%s side=%s sgp_decimal=%s",
                 game_id[:8], prop_type, side, sgp_decimal)

    log.info("Writing %d rows to %s", len(rows), db_path)
    con = duckdb.connect(str(db_path))
    try:
        write_rows(con, rows)
    finally:
        con.close()
    log.info("Done. fetch_time=%s", fetch_time.isoformat())
    return 0


if __name__ == '__main__':
    sys.exit(main())
