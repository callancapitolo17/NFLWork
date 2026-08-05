"""Shared types and utilities for the MLB SGP library.

Module-private (`_` prefix) but stable API consumed by per-book modules
(draftkings.py, fanduel.py, prophetx.py, novig.py) and by callers
(dashboard scraper shims, kalshi_mlb_rfq/sgp_runner.py).
"""
from __future__ import annotations
import logging
import random
import threading
import time
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Literal

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


# --------------------------------------------------------------------------- #
# Bounded retry / backoff (issue #34)                                          #
# --------------------------------------------------------------------------- #

# A retryable status is one where the SAME request may well succeed if we ask
# again: the book is overloaded (5xx) or throttling us (429). Everything else
# non-200 — 401/403 auth, 404 delisted, 422 non-combinable — is a real verdict
# and retrying it only burns the caller's latency budget.
RETRYABLE_STATUS_CODES = frozenset({429})
RETRYABLE_STATUS_FLOOR = 500

# curl_cffi's RequestException and requests' RequestException BOTH subclass
# OSError (curl_cffi 0.14: RequestException -> CurlError -> OSError), as do
# socket.timeout and the builtin ConnectionError family. One base therefore
# covers DNS failure, connect/read timeout, TLS error and connection reset
# across every transport this package uses.
_RETRYABLE_EXCEPTION_TYPES = (OSError,)

# Belt-and-braces for a future transport whose exceptions do NOT subclass
# OSError: silently NOT retrying is the exact failure mode this ticket exists
# to remove, so fall back to matching the class name.
_RETRYABLE_EXCEPTION_NAME_HINTS = (
    "timeout", "connection", "ssl", "proxy", "dns", "curl", "socket",
)


@dataclass(frozen=True)
class RetryProfile:
    """How hard to retry one wire call, and how much delay that may add.

    Two named profiles exist (below). Both are bounded in ATTEMPTS and in
    cumulative SLEEP, because the on-demand pricing path runs inside an RFQ
    quote budget — an unbounded backoff would trade a lost book for a lost
    quote, which is the worse failure.
    """
    name: str
    max_attempts: int
    base_delay_sec: float
    backoff_multiplier: float
    max_total_delay_sec: float      # hard ceiling on cumulative sleep


# Sweep / structure-warming: the caller has no deadline worth protecting, and
# losing a book costs a whole refresh cycle. Worst case ~5.6s of added sleep
# (1.5s then 3.0s, each +/-25% jitter), hard-capped at 10s.
RETRY_BACKGROUND = RetryProfile(
    name="background", max_attempts=3, base_delay_sec=1.5,
    backoff_multiplier=2.0, max_total_delay_sec=10.0,
)

# On-demand pricing (inside the RFQ quote budget, #50 targets p95 <= 8s warm)
# and EVERY price call regardless of path — price calls are the high-fan-out
# surface (DK runs targets x combos concurrently), so a 3-attempt backoff
# there would stall a degraded sweep for minutes. One fast retry only:
# worst case ~0.63s of added sleep, hard-capped at 1.0s.
RETRY_LIVE = RetryProfile(
    name="live", max_attempts=2, base_delay_sec=0.5,
    backoff_multiplier=2.0, max_total_delay_sec=1.0,
)

# Indirected so tests can neutralize/observe backoff without burning wall
# clock (see mlb_sgp/tests/conftest.py). Resolved at CALL time, not bound as
# a default argument, so monkeypatching the module attribute works.
_sleep = time.sleep
_jitter = random.random


def is_retryable_exception(exc: BaseException) -> bool:
    """True when ``exc`` is a transport blip worth one more attempt."""
    if isinstance(exc, _RETRYABLE_EXCEPTION_TYPES):
        return True
    name = type(exc).__name__.lower()
    return any(hint in name for hint in _RETRYABLE_EXCEPTION_NAME_HINTS)


def is_retryable_status(status: int) -> bool:
    """True for 429 and 5xx. Explicitly False for every 4xx incl. 404."""
    return status in RETRYABLE_STATUS_CODES or status >= RETRYABLE_STATUS_FLOOR


def _retry_after_seconds(resp) -> float | None:
    """``Retry-After`` in seconds, or None when absent/unparseable.

    Only the delta-seconds form is honored. The HTTP-date form is legal but
    no SGP book emits it, and mis-parsing a date into a huge sleep is a worse
    outcome than falling back to our own backoff.
    """
    headers = getattr(resp, "headers", None) or {}
    try:
        raw = headers.get("Retry-After")
    except AttributeError:
        return None
    if raw is None:
        return None
    try:
        seconds = float(str(raw).strip())
    except (TypeError, ValueError):
        return None
    return seconds if seconds > 0 else None


def request_with_retry(call: Callable[[], object], *, profile: RetryProfile,
                       book: str, stage: str, detail: str = "",
                       body_ok: Callable[[object], bool] | None = None,
                       sleep_fn: Callable[[float], None] | None = None,
                       jitter_fn: Callable[[], float] | None = None):
    """Run ``call()`` (a zero-arg wire request) with bounded retries.

    Returns the LAST response. It deliberately does NOT judge that response:
    status classification stays with ``check_response`` / ``json_or_raise``
    so ``BookTransportError`` is constructed in exactly one place and the
    per-event 404 carve-out keeps working. A 404 is never retried.

    Raises ``BookTransportError`` only when every attempt raised — with the
    last exception as ``cause``, matching the pre-retry behavior of the
    ``except Exception: raise BookTransportError(...)`` blocks it replaces.
    A NON-retryable exception (a ``TypeError`` from our own bug) still raises
    on attempt 1; it is simply never retried.

    ``body_ok`` marks a 200 whose BODY is unusable as retryable — Caesars'
    AWS-WAF challenge comes back 200 with an HTML interstitial, and asking
    again is how that book recovers.

    Worst-case added latency, measured by
    ``test_retry_backoff.py::test_worst_case_added_delay_stays_within_profile``:
    RETRY_BACKGROUND <= 5.63s (hard cap 10s), RETRY_LIVE <= 0.63s (hard cap
    1.0s). Both are per WIRE CALL, not per pricing operation.
    """
    sleep = sleep_fn if sleep_fn is not None else _sleep
    jitter = jitter_fn if jitter_fn is not None else _jitter

    slept_total = 0.0
    attempts_made = 0
    last_exc: BaseException | None = None
    last_status: int | None = None
    last_resp = None

    def _detail() -> str:
        """Caller's detail plus how hard we actually tried — the attempt
        count is what tells a bot.log reader "transient blip" from "book is
        gone"."""
        suffix = (f"failed after {attempts_made} attempt(s) "
                  f"[{profile.name} profile]")
        return f"{detail} | {suffix}" if detail else suffix

    # A price-stage retry is routine and HIGH VOLUME (one per combo across a
    # thread pool), so it logs at DEBUG — hundreds of WARNINGs per cycle would
    # drown bot.log, and PriceCallTally's all-failed verdict is the aggregate
    # signal that actually matters. A session-level retry (auth/events/
    # structure) is rare and consequential, so it stays at WARNING.
    retry_log = logger.debug if stage == "price" else logger.warning

    # max(1, ...) so a misconfigured profile degrades to a single attempt
    # rather than returning None and handing callers a non-response.
    for attempt in range(max(1, profile.max_attempts)):
        retry_after = None
        attempts_made += 1
        try:
            resp = call()
        except Exception as e:                      # noqa: BLE001 — reraised below
            if not is_retryable_exception(e):
                raise BookTransportError(book, stage, status_code=last_status,
                                         cause=e, detail=_detail()) from e
            last_exc = e
        else:
            status = getattr(resp, "status_code", 200)
            last_resp, last_status, last_exc = resp, status, None
            body_unusable = (status == 200 and body_ok is not None
                             and not body_ok(resp))
            if not (is_retryable_status(status) or body_unusable):
                return resp                          # 200, or a real verdict
            retry_after = _retry_after_seconds(resp) if status == 429 else None

        # Out of attempts: hand the caller the last response so
        # check_response raises, or raise the last transport exception.
        if attempt == max(1, profile.max_attempts) - 1:
            break

        delay = _next_delay(profile, attempt, retry_after, jitter, slept_total)
        if delay is None:
            break                                    # sleep budget exhausted
        retry_log(
            "%s: retrying stage=%s attempt %d/%d (status=%s exc=%s) in "
            "%.2fs [%s profile]",
            book, stage, attempt + 1, profile.max_attempts, last_status,
            type(last_exc).__name__ if last_exc else None, delay, profile.name)
        sleep(delay)
        slept_total += delay

    if last_exc is not None:
        raise BookTransportError(book, stage, status_code=last_status,
                                 cause=last_exc, detail=_detail()) from last_exc
    return last_resp


def _next_delay(profile: RetryProfile, attempt: int, retry_after: float | None,
                jitter: Callable[[], float], slept_total: float) -> float | None:
    """Seconds to sleep before the next attempt, or None to stop retrying.

    ``Retry-After`` wins when it FITS the remaining budget; a book asking for
    30s when the profile allows 1s is not "sane", so we fall back to our own
    backoff rather than blowing a quote deadline honoring it.
    """
    remaining = profile.max_total_delay_sec - slept_total
    if remaining <= 0:
        return None
    if retry_after is not None and retry_after <= remaining:
        return retry_after
    # +/-25% jitter around the exponential step, so a book that 5xx'd every
    # concurrent combo at once doesn't get the whole pool back in lockstep.
    step = profile.base_delay_sec * (profile.backoff_multiplier ** attempt)
    delay = step * (0.75 + 0.5 * jitter())
    return min(delay, remaining)


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
        # Last HTTP status a price call declined with, and the attempt number
        # it came from. The attempt number is what scopes it to ONE cycle:
        # without it, a previous cycle's 403 would be reported against a later
        # cycle whose declines carried no status at all.
        self._last_status: int | None = None
        self._last_status_at = 0

    def record(self, ok: bool) -> None:
        with self._lock:
            self._attempted += 1
            if ok:
                self._succeeded += 1

    def note_status(self, status: int) -> None:
        """Record the status of the decline currently in flight.

        Separate from ``record`` because DK's ``calculate_sgp`` swallows a
        non-200 into ``None`` — by the time ``wrap`` tallies the result the
        status is gone. Passed in as ``calculate_sgp(..., on_decline=...)`` so
        the value reaches the verdict rather than being logged and dropped
        (issue #39: the missing ``:403`` is what made a blocked price HOST look
        like an ordinary run of unpriceable combos).

        Stamps the in-flight attempt (``_attempted + 1``): this fires inside
        the price call, before ``wrap`` increments the counter.
        """
        with self._lock:
            self._last_status = status
            self._last_status_at = self._attempted + 1

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
        """Raise if every price call since ``since`` failed.

        Carries the last decline status so ``sgp_fetch_health.error_class``
        reads ``BookTransportError:price:403`` — stage AND status, which is
        what distinguishes a blocked price host from a slate of combos the
        book simply won't build.
        """
        attempted, succeeded = self.snapshot()
        n_attempted = attempted - since[0]
        n_succeeded = succeeded - since[1]
        if n_attempted >= self.MIN_ATTEMPTS_FOR_VERDICT and n_succeeded == 0:
            with self._lock:
                # Only report a status recorded during THIS cycle.
                status = (self._last_status
                          if self._last_status_at > since[0] else None)
            raise BookTransportError(
                self.book, "price", status_code=status,
                detail=f"all {n_attempted} price calls this cycle failed")


def accepts_on_decline(price_fn) -> bool:
    """True if ``price_fn`` takes the ``on_decline`` status callback.

    The per-book price function is re-imported per call so tests and the
    dashboard can monkeypatch it, and a stand-in written before issue #39 does
    not accept the kwarg. Binding it unconditionally would make every such
    substitute raise TypeError on every combo — i.e. a silent zero-row cycle.
    Mocks and builtins have no introspectable signature; those simply do not
    get the callback.
    """
    import inspect
    try:
        return "on_decline" in inspect.signature(price_fn).parameters
    except (TypeError, ValueError):
        return False


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
            # setdefault is atomic under the GIL, so concurrent price calls
            # on a pooled client all end up sharing ONE tally.
            tally = self.__dict__.setdefault("_price_calls",
                                             PriceCallTally(self.BOOK))
        return tally


# --------------------------------------------------------------------------- #
# Per-fetch drop counters + parse tripwires (issue #35)                        #
# --------------------------------------------------------------------------- #

# One tick per LOGICAL fetch, never per wire attempt: these counters sit
# outside request_with_retry, so a combo that needed two attempts still counts
# once. The two are complementary — PriceCallTally answers "is the book's
# price endpoint reachable?", these answer "does our parser still understand
# what the book sends back?".
COUNTER_NAMES = (
    "events_seen",         # events the book listed at all
    "events_matched",      # of those, ones matched to a game on our slate
    "structure_fetches",   # WIRE events/structure fetches (TTL-cache misses);
                           # 0 = fully cache-served — #50 derives cold/warm
    "targets_attempted",   # UNITS DIFFER BY PATH — see the warning below
    "legs_attempted",      # candidate lines/legs handed to structure lookup
    "legs_resolved",       # of those, ones the book actually offers
    "prices_returned",     # UNITS DIFFER BY PATH — see the warning below
    "sanity_drops",        # combos dropped by SANITY_MULT_RATIO
    "parse_failures",      # non-transport broad-except catches, by named stage
    "transport_errors",    # BookTransportError seen during this fetch
)

# !! UNITS WARNING for #38 (per-book fetch-health history) and anything else
# that aggregates these ACROSS the two paths:
#
#   sweep      targets_attempted = TARGET LINES priced (each fans out to ~4
#                                  combo price calls)
#              prices_returned   = PricedRows produced (includes rows derived
#                                  by the integer-line fallback, which are not
#                                  1:1 with price calls)
#   on-demand  targets_attempted = PRICE CALLS issued (2^N partition cells,
#                                  plus any Route B singles)
#              prices_returned   = those calls that returned a decimal
#
# Each path's own tripwire is correct in its own units ("we priced targets and
# got no rows" / "we issued calls and got no prices"), but a ratio computed
# across paths is NOT meaningful: the sweep's is ~1:4 by construction and
# on-demand's is 1:1. Store the path alongside the counts and compare within a
# path, never across. Unifying the units would mean counting price calls in the
# sweep too, which loses "rows produced" — the number that actually says
# whether the sweep wrote data — so the split is deliberate, not an oversight.


@dataclass(frozen=True)
class FetchCountersSnapshot:
    """Immutable per-fetch counter record. Stable shape — #38 persists it."""
    book: str
    path: str                       # "sweep" | "on_demand" | "warming"
    events_seen: int
    events_matched: int
    structure_fetches: int
    targets_attempted: int
    legs_attempted: int
    legs_resolved: int
    prices_returned: int
    sanity_drops: int
    parse_failures: int
    transport_errors: int
    by_stage: tuple                 # ((stage, count), ...) sorted by stage


class FetchCounters:
    """What one book's one fetch actually did — the parse-side counterpart
    to ``PriceCallTally``.

    Motivation: the FanDuel alt-total paren regression produced ZERO rows for
    weeks with no error anywhere. Transport was healthy, the book was healthy;
    our parser had simply stopped recognizing the format. Nothing in the
    pipeline could tell that apart from "FanDuel has no markets today".

    Both the sweep (``price_sgps``) and the on-demand path
    (``SGPService.price_on_demand``) fill the SAME counter shape, so a reader
    — or #38's health table — reads one vocabulary for both.

    Thread-safe: sweep pricing fans targets and combos across thread pools.
    """

    def __init__(self, book: str, path: str = "sweep"):
        self.book = book
        self.path = path
        self._lock = threading.Lock()
        self._counts = dict.fromkeys(COUNTER_NAMES, 0)
        self._by_stage: dict[str, int] = {}
        self._emitted = False

    def bump(self, name: str, n: int = 1) -> None:
        """Add ``n`` to a counter. Unknown names raise — a typo'd counter
        that silently vanished would recreate the exact bug this fixes."""
        if name not in self._counts:
            raise AttributeError(
                f"unknown counter {name!r}; expected one of {COUNTER_NAMES}")
        with self._lock:
            self._counts[name] += n

    def record_exception(self, stage: str, logger_: logging.Logger,
                         exc: BaseException | None = None) -> int:
        """Count one broad-``except`` catch at ``stage``, in the RIGHT bucket.

        Use this — not ``record_parse_failure`` — wherever the caught
        exception could have come off the wire. The book clients differ: CZR
        and MGM swallow a price-stage ``BookTransportError`` internally and
        return None, but DK/FD/NV/PX price through module-level helpers whose
        ``request_with_retry(stage="price")`` RAISES it, and that exception
        lands in the orchestrators' broad excepts. Filing it under
        ``parse_failures`` would point a fixer at the parser while the real
        story is a dead endpoint — the exact misdiagnosis this epic exists to
        remove — and it would also leave the ``prices_empty`` tripwire armed
        when #33 already owns that failure.
        """
        if isinstance(exc, BookTransportError):
            self.bump("transport_errors")
            return self._log_at_stage(stage, logger_, exc)
        return self.record_parse_failure(stage, logger_, exc)

    def record_parse_failure(self, stage: str, logger_: logging.Logger,
                             exc: BaseException | None = None) -> int:
        """Count one PARSE failure at ``stage`` and log it.

        Only for exceptions that cannot be transport (pure structure/parse
        code). Anywhere the wire is reachable, call ``record_exception``.

        Returns the running count for this stage.
        """
        with self._lock:
            self._counts["parse_failures"] += 1
        return self._log_at_stage(stage, logger_, exc)

    def _log_at_stage(self, stage: str, logger_: logging.Logger,
                      exc: BaseException | None) -> int:
        """Tally the per-stage breakdown and log at the volume-safe level.

        Volume policy (same reasoning as #34's price-stage DEBUG decision):
        the FIRST failure of a stage in this fetch is a WARNING — that is the
        signal a human needs — and every repeat drops to DEBUG. A wholesale
        break would otherwise emit one WARNING per target per book per cycle
        and drown bot.log; the aggregate rides in the summary line and trips
        the tripwire anyway.
        """
        with self._lock:
            count = self._by_stage.get(stage, 0) + 1
            self._by_stage[stage] = count
        log_at = logger_.warning if count == 1 else logger_.debug
        log_at("%s: caught %s at stage=%s (#%d this fetch): %s",
               self.book, type(exc).__name__ if exc else "error", stage,
               count, exc)
        return count

    def snapshot(self) -> FetchCountersSnapshot:
        with self._lock:
            counts = dict(self._counts)
            by_stage = tuple(sorted(self._by_stage.items()))
        return FetchCountersSnapshot(book=self.book, path=self.path,
                                     by_stage=by_stage, **counts)

    def tripwires(self) -> list[str]:
        """Names of the "we went silent" conditions that currently hold.

        Each compares a stage's OUTPUT against its INPUT, so a genuinely
        empty slate (no input) never fires.
        """
        with self._lock:
            c = dict(self._counts)
        fired = []
        # The book listed games but none mapped onto our slate: a team-name
        # or canonical-match regression.
        if c["events_seen"] > 0 and c["events_matched"] == 0:
            fired.append("events_unmatched")
        # THE FanDuel case: candidate lines existed, the parser recognized none.
        if c["legs_attempted"] > 0 and c["legs_resolved"] == 0:
            fired.append("legs_unresolved")
        # Priced targets, got nothing back, and transport never complained —
        # a transport failure is #33's story and is already loud there.
        if (c["targets_attempted"] > 0 and c["prices_returned"] == 0
                and c["transport_errors"] == 0):
            fired.append("prices_empty")
        return fired

    def counter_fields(self) -> str:
        """The counters as ``name=value`` pairs — shared by both log lines."""
        snap = self.snapshot()
        parts = [f"book={snap.book}", f"path={snap.path}"]
        parts += [f"{name}={getattr(snap, name)}" for name in COUNTER_NAMES]
        if snap.by_stage:
            parts.append("stages=" + ",".join(f"{s}:{n}"
                                              for s, n in snap.by_stage))
        return " ".join(parts)

    def summary_line(self) -> str:
        return "sgp fetch " + self.counter_fields()

    def emit(self, logger_: logging.Logger, level: int = logging.INFO) -> None:
        """Log the summary line plus one WARNING if any tripwire fired.

        Idempotent — ``price_sgps`` emits from a ``finally``, and a nested
        caller emitting again must not double-log.
        """
        with self._lock:
            if self._emitted:
                return
            self._emitted = True
        # isEnabledFor FIRST: the on-demand path emits at DEBUG on every RFQ
        # for every book, and building the summary string only to have the
        # handler drop it is pure overhead on the quote path (#50's budget).
        if logger_.isEnabledFor(level):
            logger_.log(level, "%s", self.summary_line())
        fired = self.tripwires()
        if fired:
            logger_.warning(
                "PARSE TRIPWIRE book=%s path=%s tripped=%s — %s",
                self.book, self.path, ",".join(fired), self.counter_fields())


@contextmanager
def fetch_counters(book: str, path: str, logger_: logging.Logger,
                   level: int = logging.INFO, *,
                   counters: "FetchCounters | None" = None):
    """Scope one book's one fetch, guaranteeing the summary is emitted.

    Wraps ``price_sgps`` bodies (which have several ``return []``
    short-circuits) and ``price_on_demand``. A ``BookTransportError`` on the
    way out is tallied first — that suppresses the ``prices_empty`` tripwire,
    which would otherwise point a fixer at the parser when the real story is
    a dead endpoint.

    ``counters`` (issue #38) lets a CALLER own the object instead of having
    one created here. ``SGPService`` uses that to read a sweep's counters
    after ``price_sgps`` returns — including on the exception paths, where
    there is no return value to ride back on — and persist them to
    ``sgp_fetch_health``. Passing a counter whose ``book``/``path`` disagree
    with the arguments would silently mislabel a health row, so it raises.
    """
    if counters is None:
        counters = FetchCounters(book, path)
    elif (counters.book, counters.path) != (book, path):
        raise ValueError(
            f"fetch_counters: caller passed counters for "
            f"({counters.book!r}, {counters.path!r}) but this fetch is "
            f"({book!r}, {path!r})")
    try:
        yield counters
    except BookTransportError:
        counters.bump("transport_errors")
        raise
    finally:
        counters.emit(logger_, level=level)


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

    def get_or_fetch(self, key, fetch_fn, *, miss_cb=None):
        """Cached value for ``key``, fetching (and storing) on a miss.

        ``miss_cb`` (keyword-only) fires exactly when ``fetch_fn`` is about
        to run — i.e. this call pays the wire cost. #50 uses it to tag
        on-demand health rows cold vs. warm; passing None keeps the
        pre-#50 behavior for every existing call site.
        """
        with self._lock:
            ent = self._store.get(key)
            if ent is not None and (self._now() - ent[0]) < self.ttl_sec:
                return ent[1]
        if miss_cb is not None:
            miss_cb()
        val = fetch_fn()
        with self._lock:
            self._store[key] = (self._now(), val)
        return val

    def evict_older_than(self, max_age_sec: float) -> int:
        """Drop entries older than ``max_age_sec``; returns how many.

        #50's warming pass calls this with its own cadence before re-walking
        the slate: get_or_fetch never refreshes on a hit, so without the
        evict an entry would silently age to the TTL and expire BETWEEN
        passes — handing the cold fetch to whatever RFQ arrives first.
        """
        now = self._now()
        with self._lock:
            stale = [k for k, (t, _) in self._store.items()
                     if (now - t) > max_age_sec]
            for k in stale:
                del self._store[k]
        return len(stale)

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
    # Per-fetch drop counters for THIS on-demand call (issue #35). Optional so
    # existing constructions keep working; #38 persists it to the health table.
    counters: "FetchCountersSnapshot | None" = None


@dataclass(frozen=True)
class GameRef:
    """Game identity handed to per-book event matching."""
    game_id: str
    home_team: str
    away_team: str
    commence_time: object   # datetime | None
