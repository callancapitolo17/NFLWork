"""DuckDB store of Unabated's full-game fair ladders before kickoff (issue #264,
part of the alt-spread model epic #263), and of every edge's life to kickoff
(server/edgelog.js, 2026-11-07: whether an edge seen early is still there
later, scored against these ladders' closes).

Owns ladders.duckdb (config.LADDER_DB_PATH; /data/ladders.duckdb on the VM),
apart from bets.duckdb so the runner's ladder writes never wait on the lock
that /bets.json and Place share. This service is its only writer; nothing in
the service reads it back (analysis reads a copy).

Table:
  ladder_readings  one row per (reading, event_id, period_type_id,
               bet_type_id), POSTed by the server runner to /ladders.json
               (server/laddersnapshot.js) for every NFL / CFB / CBB game still
               to start:
                 reading 'opener'  the first row ever stored: INSERT, never
                                   overwritten (the runner re-sends it on
                                   every restart; the conflict keeps the first).
                 reading 'close'   UPSERT, kept only when its observed_at is
                                   newer than the stored one; observed_at is
                                   always before event_start (validated), so
                                   the row left at the start is the close.
                                   A close whose content (`content_hash`:
                                   mains, rungs, start, teams) equals the
                                   stored one moves only observed_at and
                                   received_at: a whole-row rewrite leaves
                                   dead space in the file (116 MB in the
                                   first week, 2026-11-13).
               side 0 = away / Over, side 1 = home / Under (Unabated's
               sides). `rungs_json` is [[side, points, fair_american, books,
               agreeing], ...]: Unabated's fair for that side at that number
               (raw, never devigged), how many board lines carried a fair
               there and how many agreed on it. A moneyline's points are null.
               Never pruned (a few thousand rows a week).
  edge_segments  one row per (line_key, segment_start), POSTed by the runner to
               /edge_segments.json: a run of league snapshots in which one NFL /
               CFB / CBB full-game board line showed Unabated's edge >= 1% at
               the same number, price (American and the book's own source
               price), fair and scheduled start. INSERT on its first row; a
               later row of the same segment only moves last_seen forward
               (UPSERT, kept when newer). last_seen is the last snapshot the
               line was seen as that edge, about 5 min plus one league refresh
               behind while the segment is open; a runner restart splits a
               segment in two (join back-to-back segments of a line that agree
               on all of those fields). Price
               and fair are Unabated's as shown (raw American, never devigged);
               liquidity is the resting amount at the segment's start.
               Never pruned.

The connection is shared by the HTTP handler threads behind one lock (DuckDB
connections are not thread-safe).
"""
import hashlib
import json
import logging
import os
import threading
from datetime import datetime
from pathlib import Path

import duckdb

READING_OPENER = "opener"
READING_CLOSE = "close"
READINGS = (READING_OPENER, READING_CLOSE)
# compact_file: a file this small is never worth a copy, and a copy is kept
# only when it is at most this share of the original.
COMPACT_MIN_BYTES = 32 * 1024 * 1024
COMPACT_KEEP_RATIO = 0.9

log = logging.getLogger(__name__)

_SCHEMA = """
CREATE TABLE IF NOT EXISTS ladder_readings (
    reading         VARCHAR NOT NULL,      -- 'opener' | 'close'
    event_id        BIGINT NOT NULL,       -- Unabated's event id
    period_type_id  INTEGER NOT NULL,      -- 1 = full game
    bet_type_id     INTEGER NOT NULL,      -- 1 moneyline, 2 spread, 3 total
    league_id       INTEGER NOT NULL,      -- 1 NFL, 2 CFB, 4 CBB
    event_name      VARCHAR,
    away_team_id    BIGINT,
    home_team_id    BIGINT,
    event_start     TIMESTAMPTZ NOT NULL,
    observed_at     TIMESTAMPTZ NOT NULL,  -- when the league snapshot it was read off landed
    main_points_0   DOUBLE,                -- Unabated's own main line, side 0
    main_fair_0     INTEGER,
    main_points_1   DOUBLE,                -- and side 1
    main_fair_1     INTEGER,
    rung_count      INTEGER NOT NULL,
    rungs_json      VARCHAR NOT NULL,
    received_at     TIMESTAMPTZ NOT NULL,
    content_hash    VARCHAR,               -- sha256 of what the row says besides its times (_content_hash)
    PRIMARY KEY (reading, event_id, period_type_id, bet_type_id)
)
"""

# Files made before content_hash: NULL there, so their first close rewrites once.
_ADD_CONTENT_HASH = "ALTER TABLE ladder_readings ADD COLUMN IF NOT EXISTS content_hash VARCHAR"

_COLUMNS = ("reading, event_id, period_type_id, bet_type_id, league_id, event_name, away_team_id, home_team_id, "
            "event_start, observed_at, main_points_0, main_fair_0, main_points_1, main_fair_1, rung_count, rungs_json, "
            "received_at, content_hash")


def _placeholders(columns: str) -> str:
    return ", ".join("?" for _ in columns.split(","))


_PLACEHOLDERS = _placeholders(_COLUMNS)

_INSERT_OPENER = f"""
INSERT INTO ladder_readings ({_COLUMNS}) VALUES ({_PLACEHOLDERS})
ON CONFLICT (reading, event_id, period_type_id, bet_type_id) DO NOTHING
"""

_UPSERT_CLOSE = f"""
INSERT INTO ladder_readings ({_COLUMNS}) VALUES ({_PLACEHOLDERS})
ON CONFLICT (reading, event_id, period_type_id, bet_type_id) DO UPDATE SET
    league_id = excluded.league_id, event_name = excluded.event_name, away_team_id = excluded.away_team_id,
    home_team_id = excluded.home_team_id, event_start = excluded.event_start, observed_at = excluded.observed_at,
    main_points_0 = excluded.main_points_0, main_fair_0 = excluded.main_fair_0,
    main_points_1 = excluded.main_points_1, main_fair_1 = excluded.main_fair_1,
    rung_count = excluded.rung_count, rungs_json = excluded.rungs_json, received_at = excluded.received_at,
    content_hash = excluded.content_hash
WHERE excluded.observed_at > ladder_readings.observed_at
"""

# An unchanged close: only its times move forward (no whole-row rewrite).
_TOUCH_CLOSE = """
UPDATE ladder_readings SET observed_at = ?, received_at = ?
WHERE reading = 'close' AND event_id = ? AND period_type_id = ? AND bet_type_id = ?
  AND content_hash = ? AND observed_at < ?
"""


_EDGE_SEGMENT_SCHEMA = """
CREATE TABLE IF NOT EXISTS edge_segments (
    line_key        VARCHAR NOT NULL,      -- the feed's board line key (book, market, side; an alt's number)
    segment_start   TIMESTAMPTZ NOT NULL,  -- the league snapshot it was first seen as this edge
    last_seen       TIMESTAMPTZ NOT NULL,  -- the last snapshot it was seen as this edge
    event_id        BIGINT NOT NULL,
    league_id       INTEGER NOT NULL,      -- 1 NFL, 2 CFB, 4 CBB
    period_type_id  INTEGER NOT NULL,      -- 1 = full game (ladder_readings' key)
    book_id         INTEGER NOT NULL,
    bet_type_id     INTEGER NOT NULL,      -- 1 moneyline, 2 spread, 3 total (full game)
    side_index      INTEGER NOT NULL,      -- 0 away / Over, 1 home / Under
    points          DOUBLE,                -- null on a moneyline
    is_alt          BOOLEAN NOT NULL,
    event_start     TIMESTAMPTZ NOT NULL,
    price_american  INTEGER NOT NULL,      -- the book's price as Unabated shows it
    source_price    DOUBLE,                -- the book's own price (e.g. an exchange's cents)
    source_format   INTEGER,
    fair_american   INTEGER NOT NULL,      -- Unabated's fair for this side at this number
    edge            DOUBLE NOT NULL,       -- Unabated's edge, a fraction (0.0156 = 1.56%)
    liquidity       DOUBLE,                -- resting at the start (exchanges); null elsewhere
    received_at     TIMESTAMPTZ NOT NULL,
    PRIMARY KEY (line_key, segment_start)
)
"""

_EDGE_SEGMENT_COLUMNS = ("line_key, segment_start, last_seen, event_id, league_id, period_type_id, book_id, bet_type_id, side_index, "
                         "points, is_alt, event_start, price_american, source_price, source_format, fair_american, edge, "
                         "liquidity, received_at")

_EDGE_SEGMENT_ROW_PLACEHOLDERS = f"({_placeholders(_EDGE_SEGMENT_COLUMNS)})"


def _upsert_edge_segments_sql(row_count: int) -> str:
    """One INSERT of `row_count` rows: a new segment is inserted, a known one
    only moves last_seen (and received_at) forward."""
    return f"""
INSERT INTO edge_segments ({_EDGE_SEGMENT_COLUMNS})
VALUES {", ".join([_EDGE_SEGMENT_ROW_PLACEHOLDERS] * row_count)}
ON CONFLICT (line_key, segment_start) DO UPDATE SET last_seen = excluded.last_seen, received_at = excluded.received_at
WHERE excluded.last_seen > edge_segments.last_seen
"""


def _rungs_json(row: dict) -> str:
    return json.dumps(row["rungs"], separators=(",", ":"))


def _content_hash(row: dict) -> str:
    """sha256 of everything a ladder row says besides its reading and times."""
    content = [row["leagueId"], row["eventName"], row["awayTeamId"], row["homeTeamId"], row["eventStart"].isoformat(),
               row["mainPoints0"], row["mainFair0"], row["mainPoints1"], row["mainFair1"], row["rungs"]]
    return hashlib.sha256(json.dumps(content, separators=(",", ":")).encode()).hexdigest()


def _table_counts(con: duckdb.DuckDBPyConnection, database: str) -> dict[str, int]:
    tables = [name for (name,) in con.execute(
        "SELECT table_name FROM duckdb_tables() WHERE database_name = ? ORDER BY table_name", [database]).fetchall()]
    return {table: con.execute(f'SELECT count(*) FROM "{database}"."{table}"').fetchone()[0] for table in tables}


def compact_file(db_path: Path) -> None:
    """Rewrite ladders.duckdb into a fresh file when that frees enough of it.

    DuckDB never gives back the space a rewritten row leaves behind, so a
    file that took many rewrites is mostly dead space. Run before the store
    opens the file (the service's start, so every deploy): a file under
    COMPACT_MIN_BYTES, or whose blocks in use fill more than
    COMPACT_KEEP_RATIO of it, is left alone; otherwise it is copied table by table
    (COPY FROM DATABASE keeps the schema and keys) to `<name>.compact`, the
    row counts are checked equal, and the copy replaces the file only when it
    is at most COMPACT_KEEP_RATIO of its size. Raises on any failure, leaving
    the original file in place (the caller logs and opens it as it is)."""
    if not db_path.exists():
        return
    with duckdb.connect(str(db_path)) as con:
        con.execute("CHECKPOINT")
        [(used_blocks, block_size)] = con.execute(
            "SELECT used_blocks, block_size FROM pragma_database_size()").fetchall()
    before = db_path.stat().st_size
    if before < COMPACT_MIN_BYTES:
        return
    # Most starts find a tight file: skip the copy (seconds on a big file)
    # unless the blocks in use would fit the size check. File size, not
    # total_blocks: DuckDB counts free blocks it has already cut off the end.
    used_bytes = used_blocks * block_size
    if used_bytes > before * COMPACT_KEEP_RATIO:
        log.info("ladder store %s not compacted: %d of its %d bytes in use", db_path, used_bytes, before)
        return
    compact_path = db_path.with_name(db_path.name + ".compact")
    for leftover in (compact_path, compact_path.with_name(compact_path.name + ".wal")):
        leftover.unlink(missing_ok=True)
    with duckdb.connect() as con:
        con.execute("ATTACH ? AS old_db (READ_ONLY)", [str(db_path)])
        con.execute("ATTACH ? AS new_db", [str(compact_path)])
        con.execute("COPY FROM DATABASE old_db TO new_db")
        old_counts = _table_counts(con, "old_db")
        new_counts = _table_counts(con, "new_db")
        con.execute("CHECKPOINT new_db")
        con.execute("DETACH new_db")
        con.execute("DETACH old_db")
    if old_counts != new_counts:
        compact_path.unlink(missing_ok=True)
        raise RuntimeError(f"compacting {db_path}: expected row counts {old_counts} in the copy, got {new_counts}")
    after = compact_path.stat().st_size
    if after > before * COMPACT_KEEP_RATIO:
        compact_path.unlink(missing_ok=True)
        log.info("ladder store %s not compacted: a copy is %d bytes of %d", db_path, after, before)
        return
    os.replace(compact_path, db_path)
    log.info("ladder store %s compacted from %d to %d bytes (%s)", db_path, before, after, new_counts)


class LadderStore:
    def __init__(self, db_path: str | Path):
        self._lock = threading.Lock()
        self._con = duckdb.connect(str(db_path))
        self._con.execute(_SCHEMA)
        self._con.execute(_ADD_CONTENT_HASH)
        self._con.execute(_EDGE_SEGMENT_SCHEMA)

    def close(self) -> None:
        with self._lock:
            self._con.close()

    def save_readings(self, rows: list[dict], received_at: datetime) -> int:
        """INSERT each opener (the first one stored stays) and write each close
        kept only when newer: an unchanged one moves only its observed_at and
        received_at, a changed one is UPSERTed whole. `rows` are
        service.validate_ladder_rows' output (event_start / observed_at
        already datetimes). Returns how many rows were written."""
        written = 0
        with self._lock:
            for row in rows:
                content_hash = _content_hash(row)
                if row["reading"] == READING_CLOSE:
                    [touched] = self._con.execute(_TOUCH_CLOSE, [
                        row["observedAt"], received_at, row["eventId"], row["periodTypeId"], row["betTypeId"],
                        content_hash, row["observedAt"],
                    ]).fetchone()
                    if touched:
                        written += touched
                        continue
                statement = _INSERT_OPENER if row["reading"] == READING_OPENER else _UPSERT_CLOSE
                [count] = self._con.execute(statement, [
                    row["reading"], row["eventId"], row["periodTypeId"], row["betTypeId"], row["leagueId"],
                    row["eventName"], row["awayTeamId"], row["homeTeamId"], row["eventStart"], row["observedAt"],
                    row["mainPoints0"], row["mainFair0"], row["mainPoints1"], row["mainFair1"], len(row["rungs"]),
                    _rungs_json(row), received_at, content_hash,
                ]).fetchone()
                written += count
        return written

    def save_edge_segments(self, rows: list[dict], received_at: datetime) -> int:
        """INSERT each new segment and move a known one's last_seen forward
        (never back), in one statement: a statement per row costs ~3 ms, and
        the runner sends up to 500 rows a request. `rows` are
        service.validate_edge_segment_rows' output (times already datetimes);
        two rows of one segment keep the later last_seen (one statement cannot
        update a row twice). Returns how many rows were written."""
        newest: dict[tuple, dict] = {}
        for row in rows:
            key = (row["lineKey"], row["segmentStart"])
            if key not in newest or row["lastSeen"] > newest[key]["lastSeen"]:
                newest[key] = row
        if not newest:
            return 0
        params = []
        for row in newest.values():
            params.extend([
                row["lineKey"], row["segmentStart"], row["lastSeen"], row["eventId"], row["leagueId"],
                row["periodTypeId"], row["bookId"], row["betTypeId"], row["sideIndex"], row["points"], row["isAlt"],
                row["eventStart"], row["price"], row["sourcePrice"], row["sourceFormat"], row["fairAmerican"], row["edge"],
                row["liquidity"], received_at,
            ])
        with self._lock:
            [written] = self._con.execute(_upsert_edge_segments_sql(len(newest)), params).fetchone()
        return written
