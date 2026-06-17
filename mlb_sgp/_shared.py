"""Shared types and utilities for the MLB SGP library.

Module-private (`_` prefix) but stable API consumed by per-book modules
(draftkings.py, fanduel.py, prophetx.py, novig.py) and by callers
(dashboard scraper shims, kalshi_mlb_rfq/sgp_runner.py).
"""
from __future__ import annotations
import logging
import threading
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Literal

import duckdb

logger = logging.getLogger(__name__)

Period = Literal["FG", "F5"]


@dataclass(frozen=True)
class TargetLine:
    """One (game, period, spread, total) tuple the bot/dashboard wants priced."""
    game_id: str
    home_team: str
    away_team: str
    commence_time: datetime
    period: Period
    spread: float   # signed, home-perspective (negative = home favored)
    total: float


@dataclass(frozen=True)
class PricedRow:
    """One book's SGP price for a (game, period, spread, total, combo) tuple."""
    game_id: str
    combo: str
    period: Period
    spread_line: float
    total_line: float
    bookmaker: str
    source: str
    sgp_decimal: float
    sgp_american: int
    fetch_time: datetime


def decimal_to_american(dec: float) -> int:
    """Convert decimal odds to American format. Favorites are negative."""
    if dec >= 2.0:
        return int(round((dec - 1.0) * 100))
    return int(round(-100.0 / (dec - 1.0)))


def american_to_decimal(am: int) -> float:
    """Inverse of decimal_to_american."""
    if am > 0:
        return 1.0 + am / 100.0
    return 1.0 + 100.0 / abs(am)


def _utc_bucket(ts: datetime | str) -> str:
    """Extract a UTC "YYYY-MM-DDTHH" bucket string from a timestamp.

    Used as a match key when correlating events across data sources at
    date+hour granularity (avoids spurious matches across doubleheaders).
    Accepts both datetime objects and ISO-8601 strings.

    Empty/None input returns "" to match existing scraper behavior
    (callers pass `lines.get("commence_time", "")` and rely on "" to
    trigger their team-only fallback matcher).
    """
    if not ts:
        return ""
    if isinstance(ts, str):
        # Normalize Z suffix to +00:00 for fromisoformat
        normalized = ts.replace("Z", "+00:00") if ts.endswith("Z") else ts
        ts = datetime.fromisoformat(normalized)
    if ts.tzinfo is None:
        ts = ts.replace(tzinfo=timezone.utc)
    else:
        ts = ts.astimezone(timezone.utc)
    return ts.strftime("%Y-%m-%dT%H")


def load_target_lines(db_path: str) -> list[TargetLine]:
    """Load target lines from a DuckDB file.

    Prefers `mlb_target_lines` (bot-written, multi-row per game) when
    present. Falls back to legacy `mlb_parlay_lines` (dashboard-written,
    one row per game with FG+F5 columns) and emits one TargetLine per
    period per game (FG always; F5 only when F5 columns are non-NULL).

    Returns [] for missing DB or missing tables — callers decide what
    to do with empty input.

    Note: Returns all rows without filtering by `written_at` or recency.
    If the caller needs fresh data, it is the caller's responsibility to
    filter at the SQL level or post-process the returned list.
    """
    if not Path(db_path).exists():
        return []
    # Retry on lock conflict — parallel SGP scrapers may be writing to the same
    # bot market DB and briefly hold an exclusive lock during upsert_priced_rows.
    # Without retry, a read_only open during that window raises IOException.
    import time as _time
    last_err = None
    for attempt in range(10):
        try:
            con = duckdb.connect(db_path, read_only=True)
            break
        except duckdb.IOException as e:
            msg = str(e).lower()
            if "lock" not in msg and "in use" not in msg:
                raise
            last_err = e
            if attempt == 9:
                raise
            _time.sleep(min(0.1 * (2 ** attempt), 1.5))
    try:
        tables = {t[0] for t in con.execute("SHOW TABLES").fetchall()}
        if "mlb_target_lines" in tables:
            return _load_from_target_lines(con)
        if "mlb_parlay_lines" in tables:
            return _load_from_parlay_lines(con)
        return []
    finally:
        con.close()


def _load_from_target_lines(con: duckdb.DuckDBPyConnection) -> list[TargetLine]:
    rows = con.execute("""
        SELECT game_id, home_team, away_team, commence_time, period, spread, total
        FROM mlb_target_lines
        ORDER BY game_id, period, spread, total
    """).fetchall()
    return [
        TargetLine(
            game_id=r[0], home_team=r[1], away_team=r[2],
            commence_time=r[3], period=r[4], spread=r[5], total=r[6],
        ) for r in rows
    ]


def _load_from_parlay_lines(con: duckdb.DuckDBPyConnection) -> list[TargetLine]:
    rows = con.execute("""
        SELECT game_id, home_team, away_team, commence_time,
               fg_spread, fg_total, f5_spread, f5_total
        FROM mlb_parlay_lines
        ORDER BY game_id
    """).fetchall()
    out: list[TargetLine] = []
    for r in rows:
        game_id, home, away, ct_raw, fg_s, fg_t, f5_s, f5_t = r
        # commence_time in legacy table may be VARCHAR or TIMESTAMP.
        # Parse VARCHAR defensively — skip rows with malformed timestamps.
        if isinstance(ct_raw, str):
            try:
                normalized = ct_raw.replace("Z", "+00:00") if ct_raw.endswith("Z") else ct_raw
                ct = datetime.fromisoformat(normalized)
            except (ValueError, AttributeError):
                # Malformed/empty timestamp — skip this game's rows entirely.
                logger.warning(
                    "load_target_lines: skipping game_id=%s with malformed commence_time=%r",
                    game_id, ct_raw,
                )
                continue
        else:
            ct = ct_raw
        if fg_s is not None and fg_t is not None:
            out.append(TargetLine(
                game_id=game_id, home_team=home, away_team=away,
                commence_time=ct, period="FG", spread=fg_s, total=fg_t,
            ))
        if f5_s is not None and f5_t is not None:
            out.append(TargetLine(
                game_id=game_id, home_team=home, away_team=away,
                commence_time=ct, period="F5", spread=f5_s, total=f5_t,
            ))
    return out


class TTLCache:
    """Tiny per-key TTL cache for structure fetches (event lists,
    selection-id dictionaries). NEVER cache prices with this.

    Thread-safe for the lock around store access; concurrent misses on
    the same key may both fetch (last write wins) — acceptable for
    idempotent GETs, and in practice each book's structure fetches run
    single-threaded (the hoisting phase of price_sgps).
    """

    def __init__(self, ttl_sec: float, now_fn=time.monotonic):
        self.ttl_sec = ttl_sec
        self._now = now_fn
        self._lock = threading.Lock()
        self._store: dict = {}

    def get_or_fetch(self, key, fetch_fn):
        with self._lock:
            ent = self._store.get(key)
            if ent is not None and (self._now() - ent[0]) < self.ttl_sec:
                return ent[1]
        val = fetch_fn()
        with self._lock:
            self._store[key] = (self._now(), val)
        return val

    def clear(self):
        with self._lock:
            self._store.clear()
