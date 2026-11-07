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
import json
import threading
from datetime import datetime
from pathlib import Path

import duckdb

READING_OPENER = "opener"
READING_CLOSE = "close"
READINGS = (READING_OPENER, READING_CLOSE)

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
    PRIMARY KEY (reading, event_id, period_type_id, bet_type_id)
)
"""

_COLUMNS = ("reading, event_id, period_type_id, bet_type_id, league_id, event_name, away_team_id, home_team_id, "
            "event_start, observed_at, main_points_0, main_fair_0, main_points_1, main_fair_1, rung_count, rungs_json, "
            "received_at")
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
    rung_count = excluded.rung_count, rungs_json = excluded.rungs_json, received_at = excluded.received_at
WHERE excluded.observed_at > ladder_readings.observed_at
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


class LadderStore:
    def __init__(self, db_path: str | Path):
        self._lock = threading.Lock()
        self._con = duckdb.connect(str(db_path))
        self._con.execute(_SCHEMA)
        self._con.execute(_EDGE_SEGMENT_SCHEMA)

    def close(self) -> None:
        with self._lock:
            self._con.close()

    def save_readings(self, rows: list[dict], received_at: datetime) -> int:
        """INSERT each opener (the first one stored stays) and UPSERT each close
        (kept only when newer). `rows` are service.validate_ladder_rows' output
        (event_start / observed_at already datetimes). Returns how many rows
        were written."""
        written = 0
        with self._lock:
            for row in rows:
                statement = _INSERT_OPENER if row["reading"] == READING_OPENER else _UPSERT_CLOSE
                [count] = self._con.execute(statement, [
                    row["reading"], row["eventId"], row["periodTypeId"], row["betTypeId"], row["leagueId"],
                    row["eventName"], row["awayTeamId"], row["homeTeamId"], row["eventStart"], row["observedAt"],
                    row["mainPoints0"], row["mainFair0"], row["mainPoints1"], row["mainFair1"], len(row["rungs"]),
                    json.dumps(row["rungs"], separators=(",", ":")), received_at,
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
