"""DuckDB persistence for the bets service.

Owns unabated_ticket/bets_service/bets.duckdb (this service is its only
writer). Tables:
  bets         one row per normalised record, UPSERT on `id` (the venue-native
               bet id) — every record ever seen, never pruned: the future CLV
               work needs the full history. `record` holds the JSON the panel
               receives; the other columns exist for the retention window and
               ad-hoc queries. A row whose content is unchanged is not
               rewritten (see upsert_bets; `content_hash` is the compare key),
               so `last_seen_at` and the record's own `sourceFetchedAt` are
               the last poll that CHANGED the record, not the last poll that
               saw it — no column records the latter.
  source_runs  APPEND, one row per source poll whether it succeeded or not —
               the per-source freshness the panel shows. Pruned to
               `source_runs_retention_days` (see _prune_source_runs_locked),
               always keeping each source's latest run and latest OK run.
  team_crosswalk  one row per (venue, league, venue team) -> Unabated team id,
               learned by the panel from a bet the board joined by its venue
               id (#118 step 4; the panel POSTs them, see service.py). INSERT
               only: a held key is never rewritten — the same id is a no-op,
               a different id is a conflict the caller is told about. Served
               with every /bets.json; DELETE /crosswalk.json empties it.
               The one exception is a manual attach (pin_bet): Cal is the
               authority, so its rows REPLACE a held key and carry
               `pinned_bet_id`; unpin_bet deletes exactly those rows.
  bet_pins     one row per bet Cal attached to a board event by hand (the
               panel's Attach control), UPSERT on `bet_id`. Never pruned: a
               few a week, and the history of manual matches.
A failed poll writes a source_runs row and touches nothing in `bets`, so a
dark source keeps serving its previous records.

The DuckDB connection is shared between the poll thread and the HTTP handler
threads behind one lock (DuckDB connections are not thread-safe).
"""
import hashlib
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
-- #125: sha256 of the record content the UPSERT compares (see _content_hash).
-- Rows from before the column carry NULL and are rewritten once.
ALTER TABLE bets ADD COLUMN IF NOT EXISTS content_hash VARCHAR;
CREATE TABLE IF NOT EXISTS source_runs (
    source       VARCHAR,
    started_at   TIMESTAMPTZ,
    finished_at  TIMESTAMPTZ,
    ok           BOOLEAN,
    error        VARCHAR,
    n_records    INTEGER
);
CREATE TABLE IF NOT EXISTS team_crosswalk (
    venue               VARCHAR NOT NULL,
    league              VARCHAR NOT NULL,
    venue_team_key      VARCHAR NOT NULL,   -- Novig team id; Kalshi event-title name
    venue_team_name     VARCHAR,
    unabated_team_id    VARCHAR NOT NULL,   -- the id the board's lines carry; key = "<league>:<id>"
    unabated_team_name  VARCHAR,
    learned_from        VARCHAR,            -- "<bet id> on board event <event id>"
    learned_at          TIMESTAMPTZ NOT NULL,
    PRIMARY KEY (venue, league, venue_team_key)
);
-- The bet whose manual attach wrote the row; NULL on rows learned from an id join.
ALTER TABLE team_crosswalk ADD COLUMN IF NOT EXISTS pinned_bet_id VARCHAR;
CREATE TABLE IF NOT EXISTS bet_pins (
    bet_id          VARCHAR PRIMARY KEY,   -- bets.id, the venue-native bet id
    venue           VARCHAR NOT NULL,
    league          VARCHAR NOT NULL,
    event_id        VARCHAR NOT NULL,      -- Unabated's eventId, as a string like the team ids
    event_start     TIMESTAMPTZ,
    away_team_id    VARCHAR,
    home_team_id    VARCHAR,
    away_team_name  VARCHAR,
    home_team_name  VARCHAR,
    pinned_at       TIMESTAMPTZ NOT NULL
);
"""

_INSERT_CROSSWALK = """
INSERT INTO team_crosswalk (venue, league, venue_team_key, venue_team_name, unabated_team_id,
                            unabated_team_name, learned_from, learned_at)
VALUES (?, ?, ?, ?, ?, ?, ?, ?)
"""

_SELECT_CROSSWALK_KEYS = "SELECT venue, league, venue_team_key, unabated_team_id FROM team_crosswalk"

_SELECT_CROSSWALK = """
SELECT venue, league, venue_team_key, venue_team_name, unabated_team_id, unabated_team_name,
       learned_from, epoch(learned_at), pinned_bet_id
FROM team_crosswalk
ORDER BY learned_at DESC, venue, league, venue_team_key
"""

# A manual attach replaces whatever is held for the key (see the module docstring).
_REPLACE_CROSSWALK = """
INSERT INTO team_crosswalk (venue, league, venue_team_key, venue_team_name, unabated_team_id,
                            unabated_team_name, learned_from, learned_at, pinned_bet_id)
VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
ON CONFLICT (venue, league, venue_team_key) DO UPDATE SET
    venue_team_name = excluded.venue_team_name, unabated_team_id = excluded.unabated_team_id,
    unabated_team_name = excluded.unabated_team_name, learned_from = excluded.learned_from,
    learned_at = excluded.learned_at, pinned_bet_id = excluded.pinned_bet_id
"""

_UPSERT_PIN = """
INSERT INTO bet_pins (bet_id, venue, league, event_id, event_start, away_team_id, home_team_id,
                      away_team_name, home_team_name, pinned_at)
VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
ON CONFLICT (bet_id) DO UPDATE SET
    venue = excluded.venue, league = excluded.league, event_id = excluded.event_id,
    event_start = excluded.event_start, away_team_id = excluded.away_team_id,
    home_team_id = excluded.home_team_id, away_team_name = excluded.away_team_name,
    home_team_name = excluded.home_team_name, pinned_at = excluded.pinned_at
"""

_SELECT_PINS = """
SELECT bet_id, venue, league, event_id, epoch(event_start), away_team_id, home_team_id,
       away_team_name, home_team_name, epoch(pinned_at)
FROM bet_pins
ORDER BY pinned_at DESC, bet_id
"""

# One statement per UPSERT_ROWS_PER_STATEMENT rows ({values} is that many
# 12-placeholder tuples): DuckDB runs executemany one statement per row, which
# measured 440 ms for a 528-row poll against 20-28 ms for multi-row VALUES.
_UPSERT_BETS = """
INSERT INTO bets (id, source, venue, league, status, placed_at, closed_at, event_date,
                  record, content_hash, first_seen_at, last_seen_at)
VALUES {values}
ON CONFLICT (id) DO UPDATE SET
    source = excluded.source, venue = excluded.venue, league = excluded.league,
    status = excluded.status, placed_at = excluded.placed_at, closed_at = excluded.closed_at,
    event_date = excluded.event_date, record = excluded.record, content_hash = excluded.content_hash,
    last_seen_at = excluded.last_seen_at
WHERE bets.content_hash IS DISTINCT FROM excluded.content_hash
"""

_UPSERT_ROW_PLACEHOLDERS = "(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"
UPSERT_ROWS_PER_STATEMENT = 200

# Keeps every source's latest run AND latest successful run whatever their
# age: source_status() reads exactly those two, and a source failing for
# longer than the window must keep showing its last good read, not blank out.
_PRUNE_SOURCE_RUNS = """
DELETE FROM source_runs
WHERE started_at < ?
  AND rowid NOT IN (SELECT rowid FROM source_runs
                    QUALIFY row_number() OVER (PARTITION BY source, ok ORDER BY started_at DESC) = 1)
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


# `sourceFetchedAt` is the poll's own clock — every source restamps every
# record with it on every poll (sources/kalshi.py:182) — so an unchanged-row
# check that compared it would never skip anything. Excluding it means the
# STORED stamp is the last poll that changed the record. That is the right
# value to serve: the panel breaks merge ties on it (bets.js dedupeByNativeId,
# newest wins), and "content current as of its last change" is the honest
# claim — restamping it with the latest poll would also mark records the
# source no longer returns (Novig/BetOnline pull a rolling window; `bets`
# serves open rows forever) as fresh, and a stale copy would then beat a
# newer content-script read for good.
VOLATILE_RECORD_FIELDS = ("sourceFetchedAt",)

# source_runs is pruned at most this often: source_status() only ever reads the
# latest run per source, so a DELETE on every poll would scan the table under
# the store lock for nothing.
SOURCE_RUNS_PRUNE_INTERVAL_SEC = 3600


def _content_hash(record: dict) -> str:
    """sha256 of the record with the volatile fields removed, keys sorted —
    the UPSERT's compare key: an exact text compare (no Python `==`, where
    1 == 1.0 and True == 1). The hash covers the record only — a change to
    how the derived columns (status, closed_at, ...) are extracted does not
    rewrite existing rows, and a build from before the column leaves the
    hash stale while rewriting `record`; after either, force a one-time
    rewrite with `UPDATE bets SET content_hash = NULL`."""
    content = {key: value for key, value in record.items() if key not in VOLATILE_RECORD_FIELDS}
    return hashlib.sha256(json.dumps(content, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _last_record_per_id(records: list[dict]) -> list[dict]:
    """The records with any duplicate id collapsed to its LAST occurrence.
    No source emits one today; if one did, a multi-row INSERT ... ON CONFLICT
    silently keeps the FIRST row for a key repeated within the statement
    (verified on DuckDB 1.4.4) — the old per-row executemany applied them in
    order, so last-wins is the behaviour kept, and the duplicate is logged."""
    by_id: dict[str, dict] = {}
    duplicates: set[str] = set()
    for record in records:
        if record["id"] in by_id:
            duplicates.add(record["id"])
        by_id[record["id"]] = record
    if duplicates:
        log.warning("upsert_bets: %d record id(s) repeated in one poll, last occurrence kept: %s",
                    len(duplicates), sorted(duplicates))
    return list(by_id.values())


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
    def __init__(self, db_path: str | Path, source_runs_retention_days: int):
        """`source_runs_retention_days` 0 disables the prune (the repo's
        convention for switching a guard off); negative is a config error."""
        if source_runs_retention_days < 0:
            raise ValueError(
                f"source_runs_retention_days must be 0 (no prune) or more, got {source_runs_retention_days}")
        self._lock = threading.Lock()
        self._source_runs_retention_days = source_runs_retention_days
        self._last_source_runs_prune_at: datetime | None = None
        self._con = duckdb.connect(str(db_path))
        self._con.execute(_SCHEMA)

    def close(self) -> None:
        with self._lock:
            self._con.close()

    def upsert_bets(self, records: list[dict], seen_at: datetime) -> int:
        """UPSERT the records; a row whose content_hash is unchanged is left
        alone by the UPSERT's WHERE (measured: zero WAL bytes over 20
        identical polls of 500 rows, against 1.3 MB unconditionally). Returns
        how many rows were inserted or updated. A source re-sends every record
        it has ever seen on every poll (Kalshi: ~111 records every 60 s, of
        which ~109 are unchanged), and a DuckDB update is copy-on-write —
        rewriting an identical row appends a new version and turns the old
        one into garbage, which is what grew the WAL to 9.8 MB against a
        4.7 MB file holding 500 rows (#125). `last_seen_at` is written only
        on those rows — it is the last poll that CHANGED the record."""
        rows = [(
            record["id"], record.get("source"), record.get("venue"), record.get("league"),
            record.get("status"), _parse_iso(record.get("placedAt")), _parse_iso(record.get("closedAt")),
            record.get("eventDate"), json.dumps(record, separators=(",", ":")), _content_hash(record),
            seen_at, seen_at,
        ) for record in _last_record_per_id(records)]
        written = 0
        with self._lock:
            for start in range(0, len(rows), UPSERT_ROWS_PER_STATEMENT):
                chunk = rows[start:start + UPSERT_ROWS_PER_STATEMENT]
                statement = _UPSERT_BETS.format(values=", ".join([_UPSERT_ROW_PLACEHOLDERS] * len(chunk)))
                # DuckDB's INSERT result is the rows it inserted or updated —
                # the ones the WHERE let through (verified on 1.4.4).
                [written_in_chunk] = self._con.execute(statement, [value for row in chunk for value in row]).fetchone()
                written += written_in_chunk
        return written

    def log_source_run(self, source: str, started_at: datetime, finished_at: datetime,
                       ok: bool, error: str | None, n_records: int) -> None:
        with self._lock:
            self._con.execute(
                "INSERT INTO source_runs (source, started_at, finished_at, ok, error, n_records) "
                "VALUES (?, ?, ?, ?, ?, ?)",
                [source, started_at, finished_at, ok, error, n_records])
            self._prune_source_runs_locked(finished_at)

    def _prune_source_runs_locked(self, now: datetime) -> int:
        """DELETE source_runs rows older than the retention window — never a
        source's latest run of either outcome, so its latest run and latest
        OK run always survive — at most once per
        SOURCE_RUNS_PRUNE_INTERVAL_SEC; returns how many went. The store lock
        is already taken. Nothing reads a run but source_status(), which
        full-scans this table twice under the lock on every /bets.json and
        /health, so an unbounded table (~2,000 rows/day measured) slows every
        poll and every panel refresh. Never raises: the run that triggered it
        is already logged, and an error here must not re-enter
        run_source_once's failure path and log that poll a second time as
        failed. The throttle advances first, so a persistent error is
        reported hourly, not every poll."""
        if self._source_runs_retention_days == 0:
            return 0
        if (self._last_source_runs_prune_at is not None
                and (now - self._last_source_runs_prune_at).total_seconds() < SOURCE_RUNS_PRUNE_INTERVAL_SEC):
            return 0
        self._last_source_runs_prune_at = now
        cutoff = now - timedelta(days=self._source_runs_retention_days)
        try:
            [deleted] = self._con.execute(_PRUNE_SOURCE_RUNS, [cutoff]).fetchone()
        except Exception:  # noqa: BLE001 — a failed prune is logged, never a failed poll
            log.exception("source_runs: prune of rows older than %s failed; retrying in %ds",
                          cutoff.isoformat(), SOURCE_RUNS_PRUNE_INTERVAL_SEC)
            return 0
        if deleted:
            log.info("source_runs: pruned %d row(s) older than %s", deleted, cutoff.isoformat())
        return deleted

    def load_bets(self, days: int, now: datetime) -> list[dict]:
        """Records in the retention window: open, or closed within `days` of `now`."""
        cutoff = now - timedelta(days=days)
        with self._lock:
            rows = self._con.execute(_SELECT_WINDOW, [cutoff]).fetchall()
        return [json.loads(row[0]) for row in rows]

    def learn_crosswalk(self, rows: list[dict], learned_at: datetime) -> dict:
        """INSERT the crosswalk rows not yet held (validated shape: see
        service.validate_crosswalk_rows). A held key is never rewritten: the
        same Unabated id is a no-op, a different one is returned in
        `conflicts` with what is held. Returns {learned, conflicts}."""
        learned = 0
        conflicts: list[dict] = []
        with self._lock:
            held = {(venue, league, key): team_id
                    for venue, league, key, team_id in self._con.execute(_SELECT_CROSSWALK_KEYS).fetchall()}
            for row in rows:
                key = (row["venue"], row["league"], row["venueTeamKey"])
                if key in held:
                    if held[key] != row["unabatedTeamId"]:
                        conflicts.append({"venue": row["venue"], "league": row["league"],
                                          "venueTeamKey": row["venueTeamKey"],
                                          "held": held[key], "proposed": row["unabatedTeamId"]})
                    continue
                self._con.execute(_INSERT_CROSSWALK, [
                    row["venue"], row["league"], row["venueTeamKey"], row.get("venueTeamName"),
                    row["unabatedTeamId"], row.get("unabatedTeamName"), row.get("learnedFrom"), learned_at])
                held[key] = row["unabatedTeamId"]
                learned += 1
        if conflicts:
            log.warning("crosswalk: %d row(s) refused, a different Unabated id is held: %s", len(conflicts), conflicts)
        return {"learned": learned, "conflicts": conflicts}

    def load_crosswalk(self) -> list[dict]:
        """Every crosswalk row, newest first, in the panel's camelCase shape."""
        with self._lock:
            rows = self._con.execute(_SELECT_CROSSWALK).fetchall()
        return [{
            "venue": venue, "league": league, "venueTeamKey": venue_team_key, "venueTeamName": venue_team_name,
            "unabatedTeamId": unabated_team_id, "unabatedTeamName": unabated_team_name,
            "learnedFrom": learned_from, "learnedAt": _epoch_to_iso(learned_at), "pinnedBetId": pinned_bet_id,
        } for venue, league, venue_team_key, venue_team_name, unabated_team_id, unabated_team_name,
              learned_from, learned_at, pinned_bet_id in rows]

    def clear_crosswalk(self) -> int:
        """DELETE every crosswalk row; returns how many were held."""
        with self._lock:
            count = self._con.execute("SELECT count(*) FROM team_crosswalk").fetchone()[0]
            self._con.execute("DELETE FROM team_crosswalk")
        log.info("crosswalk: cleared %d row(s)", count)
        return count

    def bet_venue(self, bet_id: str) -> str | None:
        """The venue of the stored bet with this id, or None when `bets` holds
        no such record (a pin must name a real bet, of its own venue)."""
        with self._lock:
            row = self._con.execute("SELECT venue FROM bets WHERE id = ?", [bet_id]).fetchone()
        return row[0] if row else None

    def pin_bet(self, pin: dict, crosswalk_rows: list[dict], pinned_at: datetime) -> None:
        """Cal's manual attach, in one transaction: UPSERT the pin (validated
        shape: see service.validate_pin_request), drop the crosswalk rows an
        earlier attach of the same bet taught, then write this attach's rows,
        REPLACING a held key (a manual lesson overrides an automatic one)."""
        with self._lock:
            self._con.execute("BEGIN TRANSACTION")
            try:
                self._con.execute(_UPSERT_PIN, [
                    pin["betId"], pin["venue"], pin["league"], pin["eventId"], _parse_iso(pin.get("eventStart")),
                    pin.get("awayTeamId"), pin.get("homeTeamId"), pin.get("awayTeamName"), pin.get("homeTeamName"),
                    pinned_at])
                self._con.execute("DELETE FROM team_crosswalk WHERE pinned_bet_id = ?", [pin["betId"]])
                for row in crosswalk_rows:
                    self._con.execute(_REPLACE_CROSSWALK, [
                        row["venue"], row["league"], row["venueTeamKey"], row.get("venueTeamName"),
                        row["unabatedTeamId"], row.get("unabatedTeamName"), row.get("learnedFrom"), pinned_at,
                        pin["betId"]])
                self._con.execute("COMMIT")
            except Exception:
                self._con.execute("ROLLBACK")
                raise
        log.info("pin: bet %s -> %s event %s, %d crosswalk row(s) taught",
                 pin["betId"], pin["league"], pin["eventId"], len(crosswalk_rows))

    def unpin_bet(self, bet_id: str) -> dict:
        """Undo an attach: DELETE the pin and the crosswalk rows it taught.
        A row it replaced is not restored. Returns {removedPin, removedRows}."""
        with self._lock:
            self._con.execute("BEGIN TRANSACTION")
            try:
                [removed_pin] = self._con.execute("DELETE FROM bet_pins WHERE bet_id = ?", [bet_id]).fetchone()
                [removed_rows] = self._con.execute(
                    "DELETE FROM team_crosswalk WHERE pinned_bet_id = ?", [bet_id]).fetchone()
                self._con.execute("COMMIT")
            except Exception:
                self._con.execute("ROLLBACK")
                raise
        log.info("pin: bet %s unpinned, %d crosswalk row(s) removed", bet_id, removed_rows)
        return {"removedPin": removed_pin > 0, "removedRows": removed_rows}

    def load_pins(self) -> list[dict]:
        """Every pin, newest first, in the panel's camelCase shape."""
        with self._lock:
            rows = self._con.execute(_SELECT_PINS).fetchall()
        return [{
            "betId": bet_id, "venue": venue, "league": league, "eventId": event_id,
            "eventStart": _epoch_to_iso(event_start), "awayTeamId": away_team_id, "homeTeamId": home_team_id,
            "awayTeamName": away_team_name, "homeTeamName": home_team_name, "pinnedAt": _epoch_to_iso(pinned_at),
        } for bet_id, venue, league, event_id, event_start, away_team_id, home_team_id,
              away_team_name, home_team_name, pinned_at in rows]

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
