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

# Moneyline × total is a 4-cell grid exactly like spread × total: the four cells
# partition the outcome space (exactly one of home/away wins, exactly one of
# over/under), so they sum to 1 before vig and the existing devig_book n-way
# probit devig works on them unchanged. Stored with total_line set and
# spread_line = NULL (there is no spread leg) — db.upsert_priced_rows dedups
# NULL-safely (IS NOT DISTINCT FROM), so NULL is the correct, leak-free marker.
ML_TOTAL_COMBO_NAMES = (
    "Home ML + Over",
    "Home ML + Under",
    "Away ML + Over",
    "Away ML + Under",
)


class BookTransportError(Exception):
    """A book's HTTP transport failed — the book is DOWN, not empty.

    Raised by the per-book clients on a non-200 (other than a per-event 404),
    a connection error/timeout, malformed JSON, or a failed auth gate. It is
    deliberately NOT caught by the orchestrators' ``price_sgps``: it must
    propagate so callers can tell "the book has no markets today" (``[]``)
    apart from "the book is unreachable".

    Every caller MUST, on seeing this:
      1. preserve the book's previously-written ``mlb_sgp_odds`` rows
         (do NOT run ``db.clear_source``), and
      2. count a failure (``SGPService`` strikes the book and reinits its
         client after ``MAX_FAILURES_BEFORE_REINIT``).

    ``stage`` is one of:
      "auth"      — session-level credential/token gate (MGM accessid,
                    CZR AWS-WAF token)
      "events"    — the book's event/fixture listing call
      "structure" — a per-event market/selection structure fetch
      "price"     — the combo pricing call(s); see ``PriceCallTally``
    """

    def __init__(self, book: str, stage: str, status_code: int | None = None,
                 cause: BaseException | None = None, detail: str = ""):
        self.book = book
        self.stage = stage
        self.status_code = status_code
        self.cause = cause
        self.detail = detail
        parts = [f"{book} transport failed at stage={stage}"]
        if status_code is not None:
            parts.append(f"status={status_code}")
        if cause is not None:
            parts.append(f"cause={type(cause).__name__}: {cause}")
        if detail:
            parts.append(detail)
        super().__init__(" | ".join(parts))


def check_response(book: str, stage: str, resp, allow_404: bool = False) -> bool:
    """Raise ``BookTransportError`` unless ``resp`` is a usable 200.

    Returns False (instead of raising) for a 404 when ``allow_404`` — a
    per-event structure fetch for a game the book has removed (postponed,
    delisted) is a legitimate skip, not a dead book. Without that carve-out
    one stale game_id would abort a whole book's cycle.
    """
    status = getattr(resp, "status_code", 200)
    if status == 200:
        return True
    if allow_404 and status == 404:
        return False
    raise BookTransportError(book, stage, status_code=status,
                             detail=(getattr(resp, "text", "") or "")[:200])


def json_or_raise(book: str, stage: str, resp):
    """``resp.json()`` with malformed bodies converted to BookTransportError.

    A WAF challenge page and a Cloudflare interstitial both come back 200 with
    HTML — decoding them is the only way to notice.
    """
    try:
        return resp.json()
    except Exception as e:
        raise BookTransportError(book, stage, cause=e,
                                 detail="response body is not JSON") from e


class PriceCallTally:
    """Per-cycle tally of one book's combo-price calls.

    Motivation: a book's event listing and market structure can be perfectly
    healthy while every PRICE call is blocked — DK, for instance, reads from
    ``sportsbook-nash.draftkings.com`` but prices on
    ``gaming-us-nj.draftkings.com``. Individual declined combos are normal
    (partial row drops), but a whole cycle of price calls returning nothing is
    a transport verdict, not an empty slate.

    Lives on the persistent book client, so ``price_sgps`` takes a
    ``snapshot()`` at the top and calls ``verdict(snapshot)`` before returning
    — scoping the judgement to THIS cycle. Thread-safe: every orchestrator
    prices combos through a thread pool.
    """

    # Two games' worth of combos. Below this a legitimately thin slate (one
    # game the book declines to build) would look like a dead book.
    MIN_ATTEMPTS_FOR_VERDICT = 8

    def __init__(self, book: str):
        self.book = book
        self._lock = threading.Lock()
        self._attempted = 0
        self._succeeded = 0

    def record(self, ok: bool) -> None:
        with self._lock:
            self._attempted += 1
            if ok:
                self._succeeded += 1

    def snapshot(self) -> tuple[int, int]:
        with self._lock:
            return self._attempted, self._succeeded

    def wrap(self, price_fn, ok=bool):
        """Return ``price_fn`` with each call's outcome recorded.

        Used for DK/FD/NV, whose price call is a legacy module-level function
        threaded through the orchestrator rather than a client method. ``ok``
        maps the return value to success/failure — override it for helpers
        that return a tuple (Novig's ``submit_parlay`` returns
        ``(priced, auth_failed)``, whose plain truthiness is always True).
        """
        def tallied(*args, **kwargs):
            result = price_fn(*args, **kwargs)
            self.record(bool(ok(result)))
            return result
        return tallied

    def verdict(self, since: tuple[int, int]) -> None:
        """Raise if every price call since ``since`` failed."""
        attempted, succeeded = self.snapshot()
        n_attempted = attempted - since[0]
        n_succeeded = succeeded - since[1]
        if n_attempted >= self.MIN_ATTEMPTS_FOR_VERDICT and n_succeeded == 0:
            raise BookTransportError(
                self.book, "price",
                detail=f"all {n_attempted} price calls this cycle failed")


def price_tally_for(client, book: str) -> PriceCallTally:
    """The client's own ``PriceCallTally``, or a throwaway one.

    Orchestrator tests inject ``MagicMock`` clients, whose ``.price_calls``
    is another mock — wrapping a price function through it would silently
    replace the function. A throwaway tally keeps those callers working:
    nothing records into it, so ``verdict()`` is a no-op.
    """
    tally = getattr(client, "price_calls", None)
    return tally if isinstance(tally, PriceCallTally) else PriceCallTally(book)


class PriceCallTallyMixin:
    """Gives a book client a lazily-created ``.price_calls`` tally.

    Lazy (rather than set in ``__init__``) because the client tests build
    instances via ``Client.__new__`` to skip the real HTTP session.
    """

    BOOK: str = ""

    @property
    def price_calls(self) -> PriceCallTally:
        tally = self.__dict__.get("_price_calls")
        if tally is None:
            tally = PriceCallTally(self.BOOK)
            self.__dict__["_price_calls"] = tally
        return tally


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
    """One book's SGP price for a (game, period, spread, total, combo) tuple.

    spread_line / total_line are None for a combo that has no leg of that kind
    (e.g. a moneyline×total combo has total_line set and spread_line=None).
    """
    game_id: str
    combo: str
    period: Period
    spread_line: float | None
    total_line: float | None
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


# ------------------------------------------------------------------ #
# Phase 2 on-demand pricing types (kalshi_mlb_mm on-demand engine).  #
# ------------------------------------------------------------------ #

@dataclass(frozen=True)
class ResolvedLeg:
    """One canonical leg resolved at one book.

    ref / opposite_ref are book-specific selection descriptors (opaque to
    callers). opposite_ref=None means the book one-sides this line — that
    routes the book to Route B (correlation transfer). single/opposite
    decimals are the leg's two-sided single odds where the book's structure
    carries them (PX/NV/MGM/CZR, FD after odds capture); None at DK.
    """
    ref: object
    opposite_ref: object | None = None
    single_decimal: float | None = None
    opposite_decimal: float | None = None


@dataclass(frozen=True)
class OnDemandBookResult:
    """One book's on-demand fair for one same-game leg set."""
    book: str
    fair: float
    route: str              # "partition" | "transfer"
    n_cells_priced: int
    latency_sec: float
    # transfer_fair - partition_fair, computed only where both routes came
    # free (Route B vig-cancellation live measurement); None otherwise.
    route_gap: float | None = None


@dataclass(frozen=True)
class GameRef:
    """Game identity handed to per-book event matching."""
    game_id: str
    home_team: str
    away_team: str
    commence_time: object   # datetime | None
