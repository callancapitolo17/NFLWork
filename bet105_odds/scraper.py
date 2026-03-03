#!/usr/bin/env python3
"""
Bet105.ag Odds Scraper
Fetches odds from LinePros platform via WebSocket (Socket.IO v4).
No browser needed — connects directly to the data feed at pandora.ganchrow.com.

Usage:
    python scraper.py cbb
    python scraper.py nba

Architecture:
    Bet105 is a LinePros white-label sportsbook. Odds are delivered via
    Socket.IO WebSocket from pandora.ganchrow.com with gzip-compressed
    binary payloads in decimal odds format. This scraper implements the
    Socket.IO v4 / Engine.IO v4 protocol at the raw WebSocket level for
    full control over binary message handling.
"""

import gzip
import json
import os
import sys
import time
import duckdb
from datetime import datetime, timezone
from pathlib import Path

import websocket
from dotenv import load_dotenv

# Add Answer Keys to path for shared team name resolution
sys.path.insert(0, str(Path(__file__).parent.parent / "Answer Keys"))
from canonical_match import load_team_dict, load_canonical_games, resolve_team_names

# Load credentials from shared .env
env_path = Path(__file__).parent.parent / "bet_logger" / ".env"
load_dotenv(env_path)

# =============================================================================
# CONFIGURATION
# =============================================================================

WS_URL = "wss://pandora.ganchrow.com/socket.io/?EIO=4&transport=websocket"
PARTNER_ID = os.getenv("BET105_PARTNER_ID", "111")
PREMATCH_KEY = os.getenv("BET105_PREMATCH_KEY")
PREMATCH_USER_ID = int(os.getenv("BET105_USER_ID", "0"))
PREMATCH_GROUP_ID = int(os.getenv("BET105_GROUP_ID", "0"))

DB_PATH = Path(__file__).parent / "bet105.duckdb"

SPORT_CONFIGS = {
    "cbb": {
        "sport_key": "basketball_ncaab",
        "sport_id": "102",
        "table_name": "cbb_odds",
        "leagues": {"7": "NCAA Division I", "878": "NCAA Extra"},
    },
    "nba": {
        "sport_key": "basketball_nba",
        "sport_id": "2",
        "table_name": "nba_odds",
        "leagues": {"6": "NBA"},
    },
}


# =============================================================================
# ODDS CONVERSION
# =============================================================================


def decimal_to_american(dec: float) -> int | None:
    """Convert decimal odds to American format."""
    if dec is None or dec <= 1.0:
        return None
    if dec >= 2.0:
        return round((dec - 1) * 100)
    else:
        return round(-100 / (dec - 1))


# =============================================================================
# SOCKET.IO v4 PROTOCOL (raw websocket)
# =============================================================================
#
# Engine.IO packet types: 0=open, 2=ping, 3=pong, 4=message
# Socket.IO packet types (inside EIO type 4):
#   0=CONNECT, 2=EVENT, 5=BINARY_EVENT
#
# Binary protocol:
#   Text frame: "451-[\"room.name\",{\"_placeholder\":true,\"num\":0}]"
#   Binary frame: <gzip compressed data>
#
# The 451- prefix means: EIO message(4) + SIO binary-event(5) + 1 attachment + separator(-)


class Bet105Scraper:
    """Low-level Socket.IO client for LinePros prematch odds feed."""

    def __init__(self, sport: str):
        self.sport = sport
        self.config = SPORT_CONFIGS[sport]
        self.ws = None
        self.events = {}        # eventId → {away_raw, home_raw, start_time, ...}
        self.coefficients = {}  # eventId → {fg: {...}, Half1: {...}}
        self.pending_binary_room = None  # Room name for next binary frame
        self.event_data_done = False
        self.coeff_subscribed = set()
        self.coeff_received = set()

    def _send_sio_event(self, event: str, data):
        """Send a Socket.IO event (type 42 = EIO message + SIO event)."""
        payload = json.dumps([event, data])
        self.ws.send(f"42{payload}")

    def _send_sio_subscribe(self, rooms):
        """Send a subscribe event."""
        self._send_sio_event("subscribe", rooms)

    def _on_text_message(self, msg: str):
        """Handle a text WebSocket frame."""
        if not msg:
            return

        # Engine.IO open packet
        if msg[0] == "0":
            config = json.loads(msg[1:])
            self.ping_interval = config.get("pingInterval", 25000) / 1000
            # Send Socket.IO connect
            self.ws.send("40")
            return

        # Engine.IO ping → pong
        if msg[0] == "2":
            self.ws.send("3")
            return

        # Engine.IO message (type 4)
        if msg[0] != "4":
            return

        inner = msg[1:]

        # Socket.IO CONNECT response
        if inner.startswith("0"):
            self._on_sio_connected()
            return

        # Socket.IO EVENT (42)
        if inner.startswith("2"):
            try:
                arr = json.loads(inner[1:])
                self._on_sio_event(arr[0], arr[1:])
            except (json.JSONDecodeError, IndexError) as e:
                print(f"  Warning: Failed to parse SIO event: {e}")
            return

        # Socket.IO BINARY_EVENT (51-)
        if inner.startswith("5"):
            # Parse: "51-[\"room.name\",{\"_placeholder\":true,\"num\":0}]"
            dash_idx = inner.index("-")
            try:
                arr = json.loads(inner[dash_idx + 1:])
                if isinstance(arr, list) and len(arr) >= 1:
                    self.pending_binary_room = arr[0]
            except (json.JSONDecodeError, ValueError) as e:
                print(f"  Warning: Failed to parse binary event header: {e}")
            return

    def _on_binary_message(self, data: bytes):
        """Handle a binary WebSocket frame (gzip-compressed payload)."""
        room = self.pending_binary_room
        self.pending_binary_room = None
        if not room:
            return

        try:
            decompressed = gzip.decompress(data).decode("utf-8")
            payload = json.loads(decompressed)
        except Exception as e:
            print(f"  Warning: Failed to decompress binary from {room}: {e}")
            return

        inner = payload.get("payload", {})
        is_diff = payload.get("isDiff", False)

        if ".eventData" in room and not is_diff:
            self._handle_event_data(inner)
        elif ".eventCoefficients." in room and not is_diff:
            event_id = room.rsplit(".", 1)[-1]
            self._handle_coefficients(event_id, inner)

    def _on_sio_connected(self):
        """Socket.IO connection established — start subscription chain."""
        print("  Connected to LinePros feed")

        # Set metadata
        self._send_sio_event("setSocketMetadata", {
            "partnerId": PARTNER_ID,
            "flavor": "prematch",
        })

        # Subscribe to system events
        self._send_sio_event("subscribeSystemEvents", {
            "partnerId": PARTNER_ID,
            "groupId": PREMATCH_GROUP_ID,
            "userId": PREMATCH_USER_ID,
        })

        # Subscribe to event data (game list)
        room = f"prematch.main.{PREMATCH_KEY}.eventData"
        self._send_sio_subscribe([room])

    def _on_sio_event(self, event_name: str, args: list):
        """Handle a Socket.IO text event."""
        # "subscribed" confirmations — trigger coefficient subscriptions
        if event_name == "subscribed" and self.event_data_done:
            # Check if this confirms an eventCoefficients subscription
            if args and isinstance(args[0], list):
                for room in args[0]:
                    if ".eventCoefficients." in str(room):
                        pass  # Subscription confirmed, data comes via binary

    def _handle_event_data(self, payload: dict):
        """Parse eventData to extract game list for target sport/leagues."""
        sport_id = self.config["sport_id"]
        target_leagues = set(self.config["leagues"].keys())

        s = payload.get("s", {})
        sport_data = s.get(sport_id, {})

        for country_id, leagues in sport_data.items():
            for league_id, events in leagues.items():
                if league_id not in target_leagues:
                    continue
                for event_id, event_arr in events.items():
                    if not isinstance(event_arr, list) or len(event_arr) < 3:
                        continue

                    away = event_arr[0]
                    home = event_arr[1]
                    start_time = event_arr[2]

                    self.events[event_id] = {
                        "away_raw": away[0].strip() if isinstance(away, list) and away else "",
                        "home_raw": home[0].strip() if isinstance(home, list) and home else "",
                        "start_time": start_time,
                        "league_id": league_id,
                    }

        print(f"  Found {len(self.events)} {self.sport.upper()} events")
        if len(self.events) == 0:
            print(f"  WARNING: 0 events found for {self.sport.upper()}. "
                  f"If games should exist, the PREMATCH_KEY may have rotated. "
                  f"Check BET105_PREMATCH_KEY in .env.")
        self.event_data_done = True

        # Subscribe to coefficients for each event
        for event_id in self.events:
            room = f"prematch.main.{PREMATCH_KEY}.eventCoefficients.{event_id}"
            self._send_sio_subscribe([room])
            self.coeff_subscribed.add(event_id)

    def _handle_coefficients(self, event_id: str, payload: dict):
        """Parse eventCoefficients to extract spread/total/ML odds."""
        self.coeff_received.add(event_id)

        coeffs = payload.get("c", {})
        if not coeffs or isinstance(coeffs, list):
            return

        parsed = {}
        for period_key, period_label in [("m", "fg"), ("h1", "Half1")]:
            period_data = coeffs.get(period_key, {})
            if not period_data:
                continue

            record = {}

            # Type 3 = Moneyline
            ml_data = period_data.get("3", {})
            ml_odds = ml_data.get("o", {})
            if isinstance(ml_odds, dict) and "1" in ml_odds and "2" in ml_odds:
                record["away_ml"] = decimal_to_american(ml_odds["1"])
                record["home_ml"] = decimal_to_american(ml_odds["2"])

            # Type 5 = Total
            total_data = period_data.get("5", {})
            total_odds = total_data.get("o", {})
            main_total = total_data.get("r")
            if main_total is not None:
                key = str(main_total)
                # Try int key if float key doesn't exist
                if key not in total_odds and "." not in key:
                    key = str(float(main_total))
                if key in total_odds:
                    arr = total_odds[key]
                    if isinstance(arr, list) and len(arr) >= 2:
                        record["total"] = float(main_total)
                        record["over_price"] = decimal_to_american(arr[0])
                        record["under_price"] = decimal_to_american(arr[1])

            # Type 6 = Spread
            spread_data = period_data.get("6", {})
            spread_odds = spread_data.get("o", {})
            main_spread = spread_data.get("r")
            if main_spread is not None:
                key = str(main_spread)
                if key not in spread_odds:
                    key = str(float(main_spread))
                if key in spread_odds:
                    arr = spread_odds[key]
                    if isinstance(arr, list) and len(arr) >= 2:
                        record["away_spread"] = float(main_spread)
                        record["away_spread_price"] = decimal_to_american(arr[0])
                        record["home_spread"] = -float(main_spread)
                        record["home_spread_price"] = decimal_to_american(arr[1])

            if record:
                parsed[period_label] = record

        if parsed:
            self.coefficients[event_id] = parsed

        # Progress
        total = len(self.coeff_subscribed)
        done = len(self.coeff_received)
        if done % 10 == 0 or done == total:
            print(f"  Coefficients: {done}/{total}", end="\r", flush=True)

    def scrape(self, timeout: int = 30) -> list[dict]:
        """Connect, subscribe, collect odds, disconnect. Returns DuckDB records."""
        print(f"Connecting to LinePros feed...")

        self.ws = websocket.WebSocket()
        try:
            self.ws.connect(WS_URL, timeout=10)
        except Exception as e:
            print(f"  Failed to connect: {e}")
            return []

        self.ws.settimeout(1.0)  # Non-blocking reads with 1s timeout

        deadline = time.time() + timeout

        while time.time() < deadline:
            try:
                opcode, data = self.ws.recv_data()
            except websocket.WebSocketTimeoutException:
                # Check if we're done
                if self.event_data_done and self.coeff_received >= self.coeff_subscribed:
                    break
                continue
            except websocket.WebSocketConnectionClosedException:
                print("  Connection closed by server")
                break

            if opcode == websocket.ABNF.OPCODE_TEXT:
                self._on_text_message(data.decode("utf-8"))
            elif opcode == websocket.ABNF.OPCODE_BINARY:
                self._on_binary_message(data)

            # Check completion
            if (self.event_data_done
                    and self.coeff_subscribed
                    and self.coeff_received >= self.coeff_subscribed):
                break

        try:
            self.ws.close()
        except Exception:
            pass

        if not self.event_data_done:
            print("  WARNING: Never received eventData. Connection may have been rejected. "
                  "Check BET105_USER_ID and BET105_GROUP_ID in .env.")

        print(f"\n  Received coefficients for {len(self.coefficients)} events "
              f"({len(self.coeff_received)}/{len(self.coeff_subscribed)} responses)")

        return self._build_records()

    def _build_records(self) -> list[dict]:
        """Convert collected data into 18-column DuckDB records."""
        fetch_time = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")
        sport_key = self.config["sport_key"]

        team_dict = load_team_dict(self.sport)
        canonical_games = load_canonical_games(self.sport)

        records = []
        for event_id, coeffs in self.coefficients.items():
            event = self.events.get(event_id)
            if not event:
                continue

            away_raw = event["away_raw"]
            home_raw = event["home_raw"]
            if not away_raw or not home_raw:
                continue

            # Resolve team names
            if team_dict or canonical_games:
                away_team, home_team = resolve_team_names(
                    away_raw, home_raw, team_dict, canonical_games
                )
            else:
                away_team, home_team = away_raw, home_raw

            # Game date/time from unix timestamp
            start_ts = event.get("start_time")
            if start_ts:
                dt = datetime.fromtimestamp(start_ts, tz=timezone.utc)
                game_date = dt.strftime("%m/%d")
                game_time = dt.strftime("%H:%M")
            else:
                game_date = game_time = ""

            for period_label, odds in coeffs.items():
                market = "spreads" if period_label == "fg" else "spreads_h1"

                record = {
                    "fetch_time": fetch_time,
                    "sport_key": sport_key,
                    "game_id": str(event_id),
                    "game_date": game_date,
                    "game_time": game_time,
                    "away_team": away_team,
                    "home_team": home_team,
                    "market": market,
                    "period": period_label,
                    "away_spread": odds.get("away_spread"),
                    "away_spread_price": odds.get("away_spread_price"),
                    "home_spread": odds.get("home_spread"),
                    "home_spread_price": odds.get("home_spread_price"),
                    "total": odds.get("total"),
                    "over_price": odds.get("over_price"),
                    "under_price": odds.get("under_price"),
                    "away_ml": odds.get("away_ml"),
                    "home_ml": odds.get("home_ml"),
                }
                records.append(record)

        return records


# =============================================================================
# DATABASE
# =============================================================================


def init_database(sport: str):
    """Initialize DuckDB with the odds table."""
    config = SPORT_CONFIGS[sport]
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


def save_to_database(sport: str, odds_data: list):
    """Save scraped odds to DuckDB (replaces previous scrape)."""
    config = SPORT_CONFIGS[sport]
    table_name = config["table_name"]

    conn = duckdb.connect(str(DB_PATH))

    columns = [
        "fetch_time", "sport_key", "game_id", "game_date", "game_time",
        "away_team", "home_team", "market", "period",
        "away_spread", "away_spread_price", "home_spread", "home_spread_price",
        "total", "over_price", "under_price", "away_ml", "home_ml"
    ]

    placeholders = ", ".join(["?" for _ in columns])

    conn.execute(f"DELETE FROM {table_name}")
    conn.executemany(f"""
        INSERT INTO {table_name} ({", ".join(columns)})
        VALUES ({placeholders})
    """, [
        tuple(d[col] for col in columns)
        for d in odds_data
    ])

    result = conn.execute(f"SELECT COUNT(*) FROM {table_name}").fetchone()
    print(f"Database now has {result[0]} total records in {table_name}")

    conn.close()


# =============================================================================
# MAIN
# =============================================================================


def scrape_bet105(sport: str):
    """Scrape Bet105 odds via LinePros WebSocket feed."""
    if sport not in SPORT_CONFIGS:
        raise ValueError(f"Unknown sport: {sport}. Available: {list(SPORT_CONFIGS.keys())}")

    if not PREMATCH_KEY or not PREMATCH_USER_ID or not PREMATCH_GROUP_ID:
        raise ValueError(
            "BET105_PREMATCH_KEY, BET105_USER_ID, and BET105_GROUP_ID must be set in bet_logger/.env"
        )

    init_database(sport)

    scraper = Bet105Scraper(sport)
    records = scraper.scrape(timeout=30)

    if not records:
        print(f"No records scraped for {sport.upper()}")
        return []

    # Summary by market
    market_counts = {}
    for rec in records:
        m = rec["market"]
        market_counts[m] = market_counts.get(m, 0) + 1
    for m, c in sorted(market_counts.items()):
        print(f"  {m}: {c} records")

    save_to_database(sport, records)
    print(f"\nScraped {len(records)} total records for {sport.upper()}")
    return records


if __name__ == "__main__":
    sport = sys.argv[1] if len(sys.argv) > 1 else "cbb"
    print(f"Starting Bet105 {sport.upper()} odds scraper...")
    scrape_bet105(sport)
