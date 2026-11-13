"""DuckDB state for the NFL fecta pricer (config.STATE_DB_PATH).

Tables:
    fecta_quotes   APPEND — one row per (refresh, special, book): a partition
                   fair (the special's, or a superfecta's trifecta part's) or
                   DraftKings' scores-first share, or the reason there is none;
                   priced_at is when the book's prices were read (earlier than
                   quoted_at when the refresh reused them, pricing.PriceMemory).
                   The history that lets us check later whether fairs were right.
    placed_fectas  APPEND — one row per placement ATTEMPT, whatever Wagerzon said;
                   free_play TRUE when it risked the free play balance, not cash.
    settings       UPSERT — bankroll, Kelly fraction and use_free_play (1 = on).
    half_parlay_quotes   APPEND: one row per (refresh, 1H game, book, corner):
                   a book's devigged corner fair, or Wagerzon's parlay price
                   (book 'Wagerzon', wz_decimal), or the reason there is none.
    placed_half_parlays  APPEND: one row per 1H Place past the page's own
                   checks, whatever happened (placed, refused, or may be placed).
Each call opens and closes its own read-write connection: the app's threads
share one process, and DuckDB refuses a read-only connection to a file this
process already holds read-write. One lock holds those calls to one at a
time: two threads attaching the file at once fail with "Unique file handle
conflict" (2026-11-13, a fecta pricing lane appending while the page read
placements killed that refresh's fecta board).
"""
from __future__ import annotations

import threading
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path
from typing import Iterator

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
        quote_kind      VARCHAR,     -- 'full' | 'trifecta_part' | 'scores_first_share'
        fair_prob       DOUBLE,      -- fair of what quote_kind names
        sgp_decimal     DOUBLE,
        overround       DOUBLE,
        n_cells         INTEGER,     -- SGP calls made
        sf_share        DOUBLE,      -- P(scores first | trifecta part), DK rows and the 'fixed' DK-down rows
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
    # Added 2026-11-12; rows placed before it read NULL (cash).
    "ALTER TABLE placed_fectas ADD COLUMN IF NOT EXISTS free_play BOOLEAN",
    """CREATE TABLE IF NOT EXISTS half_parlay_quotes (
        quoted_at       TIMESTAMPTZ,
        league          VARCHAR,     -- 'NFL' | 'CFB'
        wz_game_id      BIGINT,      -- Wagerzon idgm of the 1H line
        away_team       VARCHAR,
        home_team       VARCHAR,
        game_start_time TIMESTAMPTZ,
        away_spread     DOUBLE,
        total           DOUBLE,
        corner          VARCHAR,     -- 'away_over' | 'away_under' | 'home_over' | 'home_under'
        book            VARCHAR,     -- a sportsbook, or 'Wagerzon' for its parlay price
        fair_prob       DOUBLE,      -- whole-number line: the chance BOTH legs win
        overround       DOUBLE,
        wz_decimal      DOUBLE,
        reason          VARCHAR,
        push_return     DOUBLE       -- whole-number line: Wagerzon's expected $ back per $1 from pushes
    )""",
    # A file made before whole-number lines were priced lacks push_return.
    "ALTER TABLE half_parlay_quotes ADD COLUMN IF NOT EXISTS push_return DOUBLE",
    # Added 2026-11-13 with PriceMemory: when the book's prices behind a row
    # were read; earlier than quoted_at when a refresh reused them. NULL on
    # rows written before it (priced at quoted_at).
    "ALTER TABLE fecta_quotes ADD COLUMN IF NOT EXISTS priced_at TIMESTAMPTZ",
    "ALTER TABLE half_parlay_quotes ADD COLUMN IF NOT EXISTS priced_at TIMESTAMPTZ",
    """CREATE TABLE IF NOT EXISTS placed_half_parlays (
        placed_at       TIMESTAMPTZ,
        account         VARCHAR,
        league          VARCHAR,
        wz_game_id      BIGINT,
        corner          VARCHAR,
        away_team       VARCHAR,
        home_team       VARCHAR,
        game_start_time TIMESTAMPTZ,
        away_spread     DOUBLE,
        total           DOUBLE,
        wz_decimal      DOUBLE,
        wz_american     INTEGER,     -- the conservative whole odds it was sized at
        risk            DOUBLE,
        to_win          DOUBLE,
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
                 "home_team", "away_team", "game_start_time", "wz_american", "book", "quote_kind",
                 "fair_prob", "sgp_decimal", "overround", "n_cells", "sf_share", "reason", "priced_at")
HALF_QUOTE_COLUMNS = ("quoted_at", "league", "wz_game_id", "away_team", "home_team", "game_start_time",
                      "away_spread", "total", "corner", "book", "fair_prob", "overround", "wz_decimal", "reason",
                      "push_return", "priced_at")
HALF_PLACEMENT_COLUMNS = ("placed_at", "account", "league", "wz_game_id", "corner", "away_team", "home_team",
                          "game_start_time", "away_spread", "total", "wz_decimal", "wz_american", "risk",
                          "to_win", "fair_prob", "ev", "status", "ticket_number", "error")
PLACEMENT_COLUMNS = ("placed_at", "account", "wz_game_id", "rotation", "description",
                     "game_start_time", "wz_american", "risk", "fair_prob", "ev",
                     "status", "ticket_number", "error", "free_play")


class Store:
    def __init__(self, db_path: Path = config.STATE_DB_PATH) -> None:
        self.db_path = db_path
        self._lock = threading.Lock()
        with self._connect() as con:
            for statement in SCHEMA:
                con.execute(statement)

    @contextmanager
    def _connect(self) -> Iterator[duckdb.DuckDBPyConnection]:
        """A read-write connection, one at a time across the app's threads."""
        with self._lock:
            con = duckdb.connect(str(self.db_path))
            try:
                yield con
            finally:
                con.close()

    def _insert(self, table: str, columns: tuple[str, ...], rows: list[dict]) -> None:
        if not rows:
            return
        placeholders = ", ".join("?" for _ in columns)
        with self._connect() as con:
            con.executemany(f"INSERT INTO {table} ({', '.join(columns)}) VALUES ({placeholders})",
                            [tuple(row[c] for c in columns) for row in rows])

    def append_quotes(self, rows: list[dict]) -> None:
        self._insert("fecta_quotes", QUOTE_COLUMNS, rows)

    def record_placement(self, row: dict) -> None:
        self._insert("placed_fectas", PLACEMENT_COLUMNS, [row])

    def placements_since(self, since: datetime) -> list[dict]:
        with self._connect() as con:
            cursor = con.execute(
                f"SELECT {', '.join(PLACEMENT_COLUMNS)} FROM placed_fectas "
                "WHERE placed_at >= ? ORDER BY placed_at DESC", [since])
            return [dict(zip(PLACEMENT_COLUMNS, row)) for row in cursor.fetchall()]

    def last_placed(self, account: str, wz_game_id: int, since: datetime) -> dict | None:
        """The newest successful placement on this special by `account` since `since`."""
        with self._connect() as con:
            row = con.execute(
                f"SELECT {', '.join(PLACEMENT_COLUMNS)} FROM placed_fectas "
                "WHERE account = ? AND wz_game_id = ? AND status = 'placed' AND placed_at >= ? "
                "ORDER BY placed_at DESC LIMIT 1", [account, wz_game_id, since]).fetchone()
            return dict(zip(PLACEMENT_COLUMNS, row)) if row else None

    def placed_risk(self, account: str, wz_game_id: int, since: datetime) -> float:
        """Dollars `account` has placed on this special since `since` (Wagerzon took them)."""
        with self._connect() as con:
            total = con.execute(
                "SELECT COALESCE(SUM(risk), 0) FROM placed_fectas "
                "WHERE account = ? AND wz_game_id = ? AND status = 'placed' AND placed_at >= ?",
                [account, wz_game_id, since]).fetchone()[0]
            return float(total)

    def append_half_quotes(self, rows: list[dict]) -> None:
        self._insert("half_parlay_quotes", HALF_QUOTE_COLUMNS, rows)

    def record_half_placement(self, row: dict) -> None:
        self._insert("placed_half_parlays", HALF_PLACEMENT_COLUMNS, [row])

    def half_placements_since(self, since: datetime) -> list[dict]:
        with self._connect() as con:
            cursor = con.execute(
                f"SELECT {', '.join(HALF_PLACEMENT_COLUMNS)} FROM placed_half_parlays "
                "WHERE placed_at >= ? ORDER BY placed_at DESC", [since])
            return [dict(zip(HALF_PLACEMENT_COLUMNS, row)) for row in cursor.fetchall()]

    def settings(self) -> dict:
        """{bankroll, kelly_fraction, use_free_play (bool)}; defaults for keys never saved."""
        values = {"bankroll": config.DEFAULT_BANKROLL, "kelly_fraction": config.DEFAULT_KELLY_FRACTION,
                  "use_free_play": 0.0}
        with self._connect() as con:
            for key, value in con.execute("SELECT key, value FROM settings").fetchall():
                values[key] = value
        values["use_free_play"] = values["use_free_play"] == 1.0
        return values

    def save_settings(self, bankroll: float, kelly_fraction: float, use_free_play: bool) -> None:
        with self._connect() as con:
            con.executemany(
                "INSERT INTO settings (key, value) VALUES (?, ?) "
                "ON CONFLICT (key) DO UPDATE SET value = excluded.value",
                [("bankroll", bankroll), ("kelly_fraction", kelly_fraction),
                 ("use_free_play", 1.0 if use_free_play else 0.0)])
