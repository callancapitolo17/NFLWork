"""DuckDB store of Unabated's full-game fair ladders before kickoff (issue #264,
part of the alt-spread model epic #263).

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
_PLACEHOLDERS = ", ".join("?" for _ in _COLUMNS.split(","))

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


class LadderStore:
    def __init__(self, db_path: str | Path):
        self._lock = threading.Lock()
        self._con = duckdb.connect(str(db_path))
        self._con.execute(_SCHEMA)

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
