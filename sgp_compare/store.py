"""DuckDB storage for every parlay the page prices.

File: $SGP_COMPARE_DB, default sgp_compare/sgp_compare.duckdb (gitignored).
Table `parlay_quotes`: APPEND only, one row per venue per Price click, so the
history can later show how each venue prices correlation by market pair.
"""
from __future__ import annotations

import json
import os
import threading
from datetime import datetime, timezone
from pathlib import Path

import duckdb

from sgp_compare.legs import Game
from sgp_compare.pricing import VenueResult

DEFAULT_DB_PATH = Path(__file__).resolve().parent / "sgp_compare.duckdb"

CREATE_PARLAY_QUOTES = """
CREATE TABLE IF NOT EXISTS parlay_quotes (
    quoted_at                TIMESTAMPTZ NOT NULL,
    venue                    VARCHAR NOT NULL,
    status                   VARCHAR NOT NULL,
    detail                   VARCHAR,
    venue_game_id            VARCHAR NOT NULL,
    home_team                VARCHAR NOT NULL,
    away_team                VARCHAR NOT NULL,
    game_start_time          TIMESTAMPTZ NOT NULL,
    n_legs                   INTEGER NOT NULL,
    legs_json                VARCHAR NOT NULL,
    best_decimal_odds        DOUBLE,
    best_max_stake           DOUBLE,
    tiers_json               VARCHAR,
    partition_json           VARCHAR,
    parlay_raw_chance        DOUBLE,
    independent_raw_chance   DOUBLE,
    correlation_raw          DOUBLE,
    parlay_fair_chance       DOUBLE,
    independent_fair_chance  DOUBLE,
    correlation_fair         DOUBLE,
    partition_overround      DOUBLE,
    venue_response_json      VARCHAR
)
"""

INSERT_PARLAY_QUOTE = """
INSERT INTO parlay_quotes (quoted_at, venue, status, detail, venue_game_id, home_team, away_team,
    game_start_time, n_legs, legs_json, best_decimal_odds, best_max_stake, tiers_json, partition_json,
    parlay_raw_chance, independent_raw_chance, correlation_raw, parlay_fair_chance,
    independent_fair_chance, correlation_fair, partition_overround, venue_response_json)
VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
"""


class Store:
    def __init__(self, db_path: Path | str | None = None) -> None:
        self.db_path = Path(db_path or os.environ.get("SGP_COMPARE_DB") or DEFAULT_DB_PATH)
        self._lock = threading.Lock()
        with self._connect() as con:
            con.execute(CREATE_PARLAY_QUOTES)

    def _connect(self) -> duckdb.DuckDBPyConnection:
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        return duckdb.connect(str(self.db_path))

    def append_quote(self, game: Game, leg_count: int, result: VenueResult) -> None:
        """APPEND one row to parlay_quotes."""
        payload = result.to_json()
        report = result.report
        best = result.quote.best if result.quote else None
        row = [
            datetime.now(timezone.utc), result.venue, result.status, result.detail,
            game.venue_game_id, game.home_team, game.away_team, game.game_start_time,
            leg_count, json.dumps(payload["legs"]),
            best.decimal_odds if best else None, best.max_stake if best else None,
            json.dumps(payload["tiers"]), json.dumps(result.partition_decimal_odds),
            report.parlay_raw_chance if report else None, report.independent_raw_chance if report else None,
            report.correlation_raw if report else None, report.parlay_fair_chance if report else None,
            report.independent_fair_chance if report else None, report.correlation_fair if report else None,
            report.partition_overround if report else None,
            json.dumps(result.quote.raw) if result.quote else None,
        ]
        with self._lock, self._connect() as con:
            con.execute(INSERT_PARLAY_QUOTE, row)
