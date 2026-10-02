"""DuckDB state for the NFL fecta pricer (config.STATE_DB_PATH).

Tables:
    fecta_quotes   APPEND — one row per (refresh, special, book): the book's
                   partition fair or the reason it has none. The history that
                   lets us check later whether these fairs were right.
    placed_fectas  APPEND — one row per placement ATTEMPT, whatever Wagerzon said.
    settings       UPSERT — bankroll and Kelly fraction.
Each call opens and closes its own read-write connection: the app's threads
share one process, and DuckDB refuses a read-only connection to a file this
process already holds read-write.
"""
from __future__ import annotations

from datetime import datetime
from pathlib import Path

import duckdb

from nfl_specials import config

SCHEMA = [
    """CREATE TABLE IF NOT EXISTS fecta_quotes (
        quoted_at       TIMESTAMPTZ,
        wz_game_id      BIGINT,
        rotation        INTEGER,
        description     VARCHAR,
        team            VARCHAR,
        prop_type       VARCHAR,
        home_team       VARCHAR,
        away_team       VARCHAR,
        game_start_time TIMESTAMPTZ,
        wz_american     INTEGER,
        book            VARCHAR,
        fair_prob       DOUBLE,
        sgp_decimal     DOUBLE,
        overround       DOUBLE,
        n_cells         INTEGER,
        reason          VARCHAR
    )""",
    """CREATE TABLE IF NOT EXISTS placed_fectas (
        placed_at       TIMESTAMPTZ,
        account         VARCHAR,
        wz_game_id      BIGINT,
        rotation        INTEGER,
        description     VARCHAR,
        game_start_time TIMESTAMPTZ,
        wz_american     INTEGER,
        risk            DOUBLE,
        fair_prob       DOUBLE,
        ev              DOUBLE,
        status          VARCHAR,
        ticket_number   VARCHAR,
        error           VARCHAR
    )""",
    """CREATE TABLE IF NOT EXISTS settings (
        key   VARCHAR PRIMARY KEY,
        value DOUBLE
    )""",
]

QUOTE_COLUMNS = ("quoted_at", "wz_game_id", "rotation", "description", "team", "prop_type",
                 "home_team", "away_team", "game_start_time", "wz_american", "book",
                 "fair_prob", "sgp_decimal", "overround", "n_cells", "reason")
PLACEMENT_COLUMNS = ("placed_at", "account", "wz_game_id", "rotation", "description",
                     "game_start_time", "wz_american", "risk", "fair_prob", "ev",
                     "status", "ticket_number", "error")


class Store:
    def __init__(self, db_path: Path = config.STATE_DB_PATH) -> None:
        self.db_path = db_path
        con = duckdb.connect(str(db_path))
        try:
            for statement in SCHEMA:
                con.execute(statement)
        finally:
            con.close()

    def _insert(self, table: str, columns: tuple[str, ...], rows: list[dict]) -> None:
        if not rows:
            return
        placeholders = ", ".join("?" for _ in columns)
        con = duckdb.connect(str(self.db_path))
        try:
            con.executemany(f"INSERT INTO {table} ({', '.join(columns)}) VALUES ({placeholders})",
                            [tuple(row[c] for c in columns) for row in rows])
        finally:
            con.close()

    def append_quotes(self, rows: list[dict]) -> None:
        self._insert("fecta_quotes", QUOTE_COLUMNS, rows)

    def record_placement(self, row: dict) -> None:
        self._insert("placed_fectas", PLACEMENT_COLUMNS, [row])

    def placements_since(self, since: datetime) -> list[dict]:
        con = duckdb.connect(str(self.db_path))
        try:
            cursor = con.execute(
                f"SELECT {', '.join(PLACEMENT_COLUMNS)} FROM placed_fectas "
                "WHERE placed_at >= ? ORDER BY placed_at DESC", [since])
            return [dict(zip(PLACEMENT_COLUMNS, row)) for row in cursor.fetchall()]
        finally:
            con.close()

    def settings(self) -> dict[str, float]:
        values = {"bankroll": config.DEFAULT_BANKROLL, "kelly_fraction": config.DEFAULT_KELLY_FRACTION}
        con = duckdb.connect(str(self.db_path))
        try:
            for key, value in con.execute("SELECT key, value FROM settings").fetchall():
                values[key] = value
        finally:
            con.close()
        return values

    def save_settings(self, bankroll: float, kelly_fraction: float) -> None:
        con = duckdb.connect(str(self.db_path))
        try:
            con.executemany(
                "INSERT INTO settings (key, value) VALUES (?, ?) "
                "ON CONFLICT (key) DO UPDATE SET value = excluded.value",
                [("bankroll", bankroll), ("kelly_fraction", kelly_fraction)])
        finally:
            con.close()
