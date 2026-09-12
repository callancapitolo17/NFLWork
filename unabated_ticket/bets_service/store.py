"""DuckDB persistence for the bets service.

Owns unabated_ticket/bets_service/bets.duckdb (this service is its only
writer). Tables:
  bets         one row per normalised record, UPSERT on `id` (the venue-native
               bet id) — every record ever seen, never pruned: the future CLV
               work needs the full history. `record` holds the JSON the panel
               receives; the other columns exist for the retention window and
               ad-hoc queries.
  source_runs  APPEND, one row per source poll whether it succeeded or not —
               the per-source freshness the panel shows.
A failed poll writes a source_runs row and touches nothing in `bets`, so a
dark source keeps serving its previous records.

The DuckDB connection is shared between the poll thread and the HTTP handler
threads behind one lock (DuckDB connections are not thread-safe).
"""
import json
import logging
import threading
from datetime import datetime, timedelta, timezone
from pathlib import Path

import duckdb

log = logging.getLogger(__name__)

_SCHEMA = """
CREATE TABLE IF NOT EXISTS bets (
    id             VARCHAR PRIMARY KEY,
    source         VARCHAR,
    venue          VARCHAR,
    league         VARCHAR,
    status         VARCHAR,
    placed_at      TIMESTAMPTZ,
    closed_at      TIMESTAMPTZ,
    event_date     VARCHAR,
    record         VARCHAR,      -- JSON, the normalised record as served
    first_seen_at  TIMESTAMPTZ,
    last_seen_at   TIMESTAMPTZ
);
CREATE TABLE IF NOT EXISTS source_runs (
    source       VARCHAR,
    started_at   TIMESTAMPTZ,
    finished_at  TIMESTAMPTZ,
    ok           BOOLEAN,
    error        VARCHAR,
    n_records    INTEGER
);
"""

_UPSERT_BET = """
INSERT INTO bets (id, source, venue, league, status, placed_at, closed_at, event_date,
                  record, first_seen_at, last_seen_at)
VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
ON CONFLICT (id) DO UPDATE SET
    source = excluded.source, venue = excluded.venue, league = excluded.league,
    status = excluded.status, placed_at = excluded.placed_at, closed_at = excluded.closed_at,
    event_date = excluded.event_date, record = excluded.record,
    last_seen_at = excluded.last_seen_at
"""

# Open bets always; settled/closed ones within the window; a non-open bet with
# no closedAt is kept (dropping it would be silent) — same rule as bets.js.
_SELECT_WINDOW = """
SELECT record FROM bets
WHERE status = 'open' OR closed_at IS NULL OR closed_at >= ?
ORDER BY placed_at DESC, id
"""

# Timestamps come back as epoch seconds: reading a TIMESTAMPTZ into Python
# needs pytz, which the service's venv does not carry.
_SELECT_LATEST_RUNS = """
SELECT source, epoch(started_at), ok, error, n_records
FROM source_runs
QUALIFY row_number() OVER (PARTITION BY source ORDER BY started_at DESC) = 1
"""

_SELECT_LATEST_OK_RUNS = """
SELECT source, epoch(started_at), n_records
FROM source_runs
WHERE ok
QUALIFY row_number() OVER (PARTITION BY source ORDER BY started_at DESC) = 1
"""


def _parse_iso(value: object) -> datetime | None:
    if not isinstance(value, str):
        return None
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        return None
    return parsed if parsed.tzinfo else parsed.replace(tzinfo=timezone.utc)


def _epoch_to_iso(epoch_seconds: float | None) -> str | None:
    if epoch_seconds is None:
        return None
    return datetime.fromtimestamp(epoch_seconds, timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


class BetsStore:
    def __init__(self, db_path: str | Path):
        self._lock = threading.Lock()
        self._con = duckdb.connect(str(db_path))
        self._con.execute(_SCHEMA)

    def close(self) -> None:
        with self._lock:
            self._con.close()

    def upsert_bets(self, records: list[dict], seen_at: datetime) -> None:
        rows = [(
            record["id"], record.get("source"), record.get("venue"), record.get("league"),
            record.get("status"), _parse_iso(record.get("placedAt")), _parse_iso(record.get("closedAt")),
            record.get("eventDate"), json.dumps(record, separators=(",", ":")), seen_at, seen_at,
        ) for record in records]
        with self._lock:
            self._con.executemany(_UPSERT_BET, rows)

    def log_source_run(self, source: str, started_at: datetime, finished_at: datetime,
                       ok: bool, error: str | None, n_records: int) -> None:
        with self._lock:
            self._con.execute(
                "INSERT INTO source_runs (source, started_at, finished_at, ok, error, n_records) "
                "VALUES (?, ?, ?, ?, ?, ?)",
                [source, started_at, finished_at, ok, error, n_records])

    def load_bets(self, days: int, now: datetime) -> list[dict]:
        """Records in the retention window: open, or closed within `days` of `now`."""
        cutoff = now - timedelta(days=days)
        with self._lock:
            rows = self._con.execute(_SELECT_WINDOW, [cutoff]).fetchall()
        return [json.loads(row[0]) for row in rows]

    def source_status(self) -> dict[str, dict]:
        """{source: {fetchedAt, ok, error, count}} — fetchedAt/count from the
        latest SUCCESSFUL run (freshness), ok/error from the latest run of any kind."""
        with self._lock:
            latest = self._con.execute(_SELECT_LATEST_RUNS).fetchall()
            latest_ok = {row[0]: row for row in self._con.execute(_SELECT_LATEST_OK_RUNS).fetchall()}
        status = {}
        for source, started_at, ok, error, _n_records in latest:
            ok_row = latest_ok.get(source)
            status[source] = {
                "fetchedAt": _epoch_to_iso(ok_row[1]) if ok_row else None,
                "ok": bool(ok),
                "error": error,
                "count": ok_row[2] if ok_row else 0,
            }
        return status
