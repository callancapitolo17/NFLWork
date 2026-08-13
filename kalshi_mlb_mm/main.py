"""Kalshi MLB MM (maker) bot — REST-polling daemon.

Composes the prior modules into three timed loops:
  - discovery: poll open RFQs, scope/price/quote in-scope spread×total combos
  - confirm:   last-look gate on accepted quotes before confirming the fill
  - risk:      kill-switch, book-freshness auto-pull, tipoff-cancel, drift sweep

Fair value is now pure book-consensus median — the model was removed in the v1
hardening pass (`fairs.py` docstring explains why). Book fairs come from live
on-demand fetches triggered by the RFQ being priced (#54; the background sweep
was removed entirely in #81 — the maker's only background book traffic is #50's
structure warming) and are gated by `_consensus_filter` (issue #20: z-space
dispersion threshold — quote only when the books' fairs agree within
SIGMA_Z_MAX in probit space; a loudly dissenting book declines the combo
rather than being outvoted).
"""
import argparse
import json
import logging
import os
import signal
import statistics
import threading
import time
import uuid
from datetime import datetime, timedelta, timezone

import duckdb
from scipy.stats import norm

from kalshi_common import auth_client, legset, sgp_runner
from kalshi_common.ev_calc import maker_fee_per_contract
from kalshi_common.leg_types import (
    _parse_event_suffix,
    _MLB_CODE_TO_TEAM,
)
from kalshi_mlb_mm import (config, db, notify, pricing, research, risk, scope,
                           router, settlement, singles)
from kalshi_mlb_mm.rfq_source import RestRFQSource
from kalshi_mlb_mm.ws_rfq_source import _parse_created_ts
from kalshi_mlb_mm.quote_gateway import RestQuoteGateway
# #23: the constituent-singles helpers moved to kalshi_mlb_mm/singles.py so the
# quote-time sanity gate, the risk sweep, and #17's confirm veto all share ONE
# implementation of "read a Kalshi single and compare it". Imported under their
# original private names so existing call sites and test monkeypatches
# (main._leg_market_prices, main._singles_moved) keep working unchanged.
from kalshi_mlb_mm.singles import leg_market_prices as _leg_market_prices
from kalshi_mlb_mm.singles import singles_moved as _singles_moved

log = logging.getLogger("kalshi_mlb_mm")

_running = threading.Event()
_running.set()
_PREV_BOOK_FAIR = {}     # combo_market_ticker -> last blended book fair (circuit breaker)
_CONSTITUENT_POLL_CURSOR = 0   # #23: rotates the budgeted constituent poll (see _rotated_poll_order)
_SCOPE_CACHE = {}        # market_ticker -> (in_scope, game_id, legs)
# B5 fix (issue #18): bound the scope cache. Overflow evicts the OLDEST
# quarter, never clear(): the WS mirror serves the exchange-wide open-RFQ
# set (4k+ tickers on 2026-08-10), so this is a working-set limit — a full
# clear would re-fetch the entire universe at ~0.5s/ticker every pass.
_SCOPE_CACHE_MAX = 5000
_DISCOVERY_CURSOR = 0    # rotating start offset for the time-boxed discovery pass


def _scope_cache_put(ticker: str, in_scope: bool, game_id, legs) -> None:
    """Insert a scope verdict, evicting the OLDEST quarter on overflow
    (dict = insertion order) — never clear(): with 4k+ open exchange-wide
    tickers a full clear forces a full re-resolve of the universe every
    pass. max(1, ...) keeps the cap a hard bound even for tiny test-sized
    caps where MAX // 4 == 0."""
    if len(_SCOPE_CACHE) >= _SCOPE_CACHE_MAX:
        evict_n = max(1, _SCOPE_CACHE_MAX // 4)
        for stale_ticker in list(_SCOPE_CACHE)[:evict_n]:
            del _SCOPE_CACHE[stale_ticker]
    _SCOPE_CACHE[ticker] = (in_scope, game_id, legs)


_VOID_HALT_ACTIVE = False  # N12: track prior void-rate halt state for edge-triggered notify

# Phase 2 on-demand pricing: engine constructed in main_loop (None in tests /
# before startup — every use is None-guarded, failing safe to "don't quote").
_ENGINE = None
_OD_RESULT_EMITTED = {}  # leg_set_hash -> landed_at of last on_demand_result emit
# Confirm-window budget for the confirm-tick live re-fetch: the ~30s Kalshi
# confirm window minus the poll gap and a buffer for the confirm API call.
CONFIRM_REFETCH_BUDGET_SEC = 20.0


def _signal_handler(_s, _f):
    _running.clear()


def _configure_auth():
    auth_client.configure(config.KALSHI_API_KEY_ID, config.KALSHI_PRIVATE_KEY_PATH,
                          config.KALSHI_BASE_URL, config.PROJECT_ROOT)


def _target_line_tick():
    """Refresh `mlb_target_lines` from Kalshi MVE enumeration + the Odds API
    schedule (issue #81 — the sweep-free half of the old sgp cycle; ZERO
    book requests). Game resolution, tipoff gating and #50's warming all
    read the table this writes. O-1 heir: new games / corrected schedule
    rows appear here now, so the game-metadata caches drop in `finally` —
    even when Kalshi enumeration fails mid-tick (the old sweep refresh
    invalidated before its early returns for the same reason)."""
    try:
        sgp_runner.target_line_cycle(bot_market_db=str(config.MARKET_DB),
                                     both_teams=True)
    finally:
        _invalidate_game_caches()


def _coverage_summary_tick(*, sgp_service, rfq_source=None):
    """Drain the service's on-demand coverage tally into one
    `on_demand_coverage` research event (issue #81 — the heir to the old
    sweep's per-pass book counts: with no sweep there is no periodic book
    count, so this is how a silently-dead book stays visible in research,
    not just #37 alerts).
    An idle window emits nothing — zero flights is the poll's normal state
    on an empty slate, not a coverage fact.

    Also emits `rfq_ingestion_summary` (2026-08-11): the WS ingestion filter
    drops the non-MLB firehose with NO per-RFQ record (that record volume is
    what it exists to prevent), so this periodic aggregate is where non-MLB
    flow stays stored — cumulative drop count + current mirror size, one row
    per COVERAGE_SUMMARY_SEC. Analysts diff consecutive rows for rates."""
    if rfq_source is not None and hasattr(rfq_source, "ingestion_dropped"):
        research.emit("rfq_ingestion_summary",
                      payload=dict(
                          non_mlb_dropped_total=rfq_source.ingestion_dropped,
                          drops_by_reason=dict(
                              getattr(rfq_source, "ingestion_drops", {}) or {}),
                          kept_total=getattr(rfq_source, "ingestion_kept", None),
                          mirror_size=len(rfq_source.poll()),
                          window_sec=config.COVERAGE_SUMMARY_SEC))
    snapshot = sgp_service.coverage.snapshot_and_reset()
    if not snapshot:
        return
    research.emit("on_demand_coverage",
                  payload=dict(books=snapshot,
                               window_sec=config.COVERAGE_SUMMARY_SEC))


def _consensus_filter(book_fairs: dict[str, float]) -> dict[str, float]:
    """Z-space dispersion gate (issue #20) — regression-test oracle.

    Algorithm (kept as an INDEPENDENT implementation of router.consensus so
    the two can cross-check each other in tests):
      1. If fewer than MIN_AGREEING_BOOKS supplied → return {} (no quote).
         A single book has sigma == 0 by construction, so the count check
         must refuse it before dispersion is considered.
      2. Transform each book's devigged fair via norm.ppf (probit z-space).
      3. If the sample stddev (ddof=1) of the z-values exceeds SIGMA_Z_MAX
         → return {} (books genuinely disagree; the dissenter may be the
         informed one, so decline instead of outvoting it).
      4. Otherwise return ALL supplied books — there is no outlier removal;
         the caller medians the full set to get the fair.
    """
    if len(book_fairs) < max(config.MIN_AGREEING_BOOKS, 2):
        return {}
    zs = [float(norm.ppf(min(max(f, 1e-6), 1.0 - 1e-6)))
          for f in book_fairs.values()]
    if statistics.stdev(zs) > config.SIGMA_Z_MAX:
        return {}
    return dict(book_fairs)


# O-1 fix (issue #18): _commence_time and _resolve_game_for_legs used to open
# a fresh read-only DuckDB connection per RFQ per 2s tick (lock churn against
# the scraper's writes and the monitor's readers). Game→game_id resolution and
# commence times are static per slate, so cache them in memory and invalidate
# on every SGP refresh (each scrape cycle + startup warm-up). Lookup FAILURES
# (None) are never cached — a transient DB lock must not stick until the next
# refresh (the risk sweep fail-safe cancels on None).
_COMMENCE_CACHE: dict = {}       # game_id -> commence_time
_RESOLVE_CACHE: dict = {}        # event_ticker -> game_id
_GAME_REF_CACHE: dict = {}       # game_id -> GameRef (issue #89: was an
                                 # uncached MARKET_DB open per pending
                                 # sub-combo per tick)


def _invalidate_game_caches():
    _COMMENCE_CACHE.clear()
    _RESOLVE_CACHE.clear()
    _GAME_REF_CACHE.clear()


def _commence_time(game_id):
    if game_id in _COMMENCE_CACHE:
        return _COMMENCE_CACHE[game_id]
    ct = _commence_time_uncached(game_id)
    if ct is not None:
        _COMMENCE_CACHE[game_id] = ct
    return ct


def _commence_time_uncached(game_id):
    # read from mlb_target_lines (written by sgp_runner) for tipoff gating
    if not config.MARKET_DB.exists():
        return None
    try:
        con = duckdb.connect(str(config.MARKET_DB), read_only=True)
    except duckdb.IOException:
        return None
    try:
        row = con.execute(
            "SELECT commence_time FROM mlb_target_lines WHERE game_id=? LIMIT 1",
            [game_id]).fetchone()
        if not row or row[0] is None:
            return None
        ct = row[0]
        # 2026-08-11 bug: mlb_target_lines stores commence_time as a NAIVE
        # TIMESTAMP whose instants are UTC (sgp_runner writes Odds API
        # commence_time), but risk._now_matching treats naive datetimes as
        # LOCAL (a contract for the R pipeline's naive-local columns). The
        # mismatch made every tipoff look ~7h later than reality — the gate
        # never fired and in-play games generated live-flight traffic all
        # day. Attach UTC here so the comparison is aware-vs-aware.
        # Durable fix (column -> TIMESTAMPTZ in sgp_runner) is tracked
        # separately; it touches the taker's copy of the table too.
        if ct.tzinfo is None:
            ct = ct.replace(tzinfo=timezone.utc)
        return ct
    except duckdb.CatalogException:
        return None
    finally:
        con.close()


def _today_fills():
    # N8: treat unreconciled fills conservatively — count as the per-fill cap
    # so a fill whose confirm-response shape is still unverified during the 30s
    # reconcile window can't under-count and let extra quotes through the cap gates.
    start = datetime.now(timezone.utc).replace(hour=0, minute=0, second=0, microsecond=0)
    with db.connect(read_only=True) as con:
        rows = con.execute(
            "SELECT game_id, "
            "CASE WHEN reconciled THEN price * contracts ELSE ? END AS exposure "
            "FROM fills WHERE filled_at >= ?",
            [config.max_fill_exposure_usd(), start]).fetchall()
    return [{"game_id": g, "price": exp} for g, exp in rows]


def _today_fills_by_game():
    """Today's exposure fanned out to one row per (fill × game), for the PER-GAME
    cap. A cross-game combo counts its FULL exposure against EVERY game it touches
    (correlated risk: a bust in any game loses the whole combo), via the
    `fill_games` map. `_today_fills` stays one-row-per-combo so the DAILY cap and
    P&L are never double-counted. Same reconciled/unreconciled exposure basis as
    `_today_fills`, computed from `fills` so reconciliation still applies.
    """
    start = datetime.now(timezone.utc).replace(hour=0, minute=0, second=0, microsecond=0)
    with db.connect(read_only=True) as con:
        rows = con.execute(
            "SELECT fg.game_id, "
            "CASE WHEN f.reconciled THEN f.price * f.contracts ELSE ? END AS exposure "
            "FROM fills f JOIN fill_games fg ON f.fill_id = fg.fill_id "
            "WHERE f.filled_at >= ?",
            [config.max_fill_exposure_usd(), start]).fetchall()
    return [{"game_id": g, "price": exp} for g, exp in rows]


class _PassSnapshot:
    """One discovery pass's read-side gate state, loaded in ONE connection.

    Issue #89: the 2026-08-10 incident fix batched the seen-lookup and the
    decision writes, but the IN-SCOPE gate path still paid 4-5 separate
    `db.connect()` opens per RFQ (~17ms each on the live state DB, measured
    ~61ms/RFQ) — ~10s of pure connection churn per pass at ~200 open RFQs.
    That made the effective RFQ re-visit cadence the PASS duration (8-12s),
    not DISCOVERY_SEC (2s), so the #55 quorum landings sat unread and the
    median RFQ-created -> quote-submitted latency was 18.5s on RFQs whose
    median lifetime is 10s. Every per-RFQ gate read now resolves from these
    in-memory maps.

    Read-your-writes: quotes submitted / replaced / cancelled DURING the
    pass are folded in via record_submitted / remove_quote, preserving the
    N7 + R-1 rule that quotes submitted earlier in the same pass count
    against the caps immediately. Fills, cooldowns and creator counts are
    written only by OTHER loop arms (confirm/settlement), which never run
    mid-pass — the loop is single-threaded — so a pass-start snapshot of
    those is exact, not approximate.
    """

    def __init__(self, *, open_quotes: dict, quote_games: dict,
                 cooldowns: dict, last_fills: dict, combo_exposure: dict,
                 creator_fill_counts: dict):
        self._open_quotes = open_quotes            # qid -> row dict
        self._quote_games = quote_games            # qid -> [game_id, ...]
        self._cooldowns = cooldowns                # ticker -> cooled_until
        self._last_fills = last_fills              # ticker -> (filled_at, fair)
        self._combo_exposure = combo_exposure      # ticker -> filled $ exposure
        self._creator_fill_counts = creator_fill_counts  # creator_id -> n

    @classmethod
    def load(cls, now_utc):
        fallback_exposure = config.max_fill_exposure_usd()
        creator_cutoff = now_utc - timedelta(
            hours=config.PER_CREATOR_WINDOW_HOURS)
        with db.connect(read_only=True) as con:
            open_rows = con.execute(
                "SELECT quote_id, rfq_id, combo_market_ticker, game_id, "
                "worst_exposure_usd, yes_bid, no_bid "
                "FROM live_quotes WHERE status='open'").fetchall()
            game_rows = con.execute(
                "SELECT qg.quote_id, qg.game_id FROM quote_games qg "
                "JOIN live_quotes lq ON lq.quote_id = qg.quote_id "
                "WHERE lq.status='open'").fetchall()
            cooldown_rows = con.execute(
                "SELECT combo_market_ticker, cooled_until "
                "FROM combo_cooldown").fetchall()
            cooled_tickers = [t for t, _ in cooldown_rows]
            fill_rows = []
            if cooled_tickers:
                # Newest-first so dict construction below keeps the LATEST
                # fill per ticker (matches the old ORDER BY ... LIMIT 1).
                fill_rows = con.execute(
                    "SELECT combo_market_ticker, filled_at, "
                    "COALESCE(fair_at_confirm, blended_fair_at_quote) "
                    "FROM fills "
                    "WHERE combo_market_ticker IN (SELECT unnest(?)) "
                    "ORDER BY filled_at DESC", [cooled_tickers]).fetchall()
            exposure_rows = con.execute(
                "SELECT combo_market_ticker, "
                "SUM(CASE WHEN reconciled THEN price * contracts ELSE ? END) "
                "FROM fills GROUP BY combo_market_ticker",
                [fallback_exposure]).fetchall()
            creator_rows = con.execute(
                "SELECT s.creator_id, COUNT(*) FROM fills f "
                "JOIN seen_rfqs s ON f.rfq_id = s.rfq_id "
                "WHERE f.filled_at >= ? GROUP BY s.creator_id",
                [creator_cutoff]).fetchall()
        open_quotes = {}
        for qid, rid, ticker, game_id, worst, yb, nb in open_rows:
            # NULL worst_exposure_usd (pre-migration rows) counts at the
            # per-fill cap — same conservative convention as N8.
            open_quotes[qid] = dict(
                rfq_id=rid, ticker=ticker, game_id=game_id,
                worst_exposure_usd=float(worst if worst is not None
                                         else fallback_exposure),
                yes_bid=yb, no_bid=nb)
        quote_games: dict = {}
        for qid, gid in game_rows:
            quote_games.setdefault(qid, []).append(gid)
        last_fills = {}
        for ticker, filled_at, fair in fill_rows:
            last_fills.setdefault(ticker, (filled_at, fair))
        return cls(open_quotes=open_quotes, quote_games=quote_games,
                   cooldowns=dict(cooldown_rows), last_fills=last_fills,
                   combo_exposure={t: float(e or 0) for t, e in exposure_rows},
                   creator_fill_counts={c: int(n or 0)
                                        for c, n in creator_rows})

    # -- per-creator halt (H4) ------------------------------------------ #

    def creator_fill_count(self, creator_id: str) -> int:
        return self._creator_fill_counts.get(creator_id, 0)

    def creator_halt_active(self, creator_id: str) -> bool:
        """H4 (per-creator): refuse a counterparty that has already farmed
        us for >= PER_CREATOR_FILL_HALT fills in the recent window."""
        if not creator_id:
            return False
        return self.creator_fill_count(creator_id) >= config.PER_CREATOR_FILL_HALT

    # -- open-quote exposure (R-1, issue #22) --------------------------- #

    def _open_rows(self, exclude_rfq_id):
        # `!=` mirrors the old SQL IS DISTINCT FROM: exclude_rfq_id=None
        # excludes nothing (every real rfq_id differs from None).
        return [(qid, row) for qid, row in self._open_quotes.items()
                if row["rfq_id"] != exclude_rfq_id]

    def open_exposure_rows(self, exclude_rfq_id: str | None = None):
        """Worst-case exposure of OPEN quotes, one row per quote — for the
        DAILY cap. The cap gates must see money that could fill in one burst
        (a counterparty sweeping resting quotes), not just money that already
        filled. Shape matches `_today_fills` rows so callers can concatenate
        the two lists straight into `risk.daily_cap_ok`.

        `exclude_rfq_id`: the RFQ currently being (re-)quoted. Its existing
        open quote is SUPERSEDED by the new one (hysteresis keeps it, or
        Kalshi auto-cancels it when the replacement lands — issue #16), so
        counting it would self-block every replace near the caps."""
        return [{"game_id": row["game_id"], "price": row["worst_exposure_usd"]}
                for _, row in self._open_rows(exclude_rfq_id)]

    def open_exposure_by_game(self, exclude_rfq_id: str | None = None):
        """Open-quote worst-case exposure fanned out to one row per
        (quote × game) via `quote_games` — for the PER-GAME cap. A cross-game
        quote counts its FULL worst case against EVERY game it touches; a
        quote with no quote_games rows (pre-migration) falls back to its
        primary game_id so it is never invisible to the per-game cap."""
        rows = []
        for qid, row in self._open_rows(exclude_rfq_id):
            for gid in self._quote_games.get(qid) or [row["game_id"]]:
                rows.append({"game_id": gid,
                             "price": row["worst_exposure_usd"]})
        return rows

    # -- per-combo / per-RFQ open-quote reads --------------------------- #

    def open_count_for_ticker(self, ticker: str) -> int:
        return sum(1 for row in self._open_quotes.values()
                   if row["ticker"] == ticker)

    def open_quote_ids_for_ticker(self, ticker: str) -> list:
        return [qid for qid, row in self._open_quotes.items()
                if row["ticker"] == ticker]

    def existing_open_quote(self, rfq_id: str):
        """(quote_id, yes_bid, no_bid) of this RFQ's resting quote, or None."""
        for qid, row in self._open_quotes.items():
            if row["rfq_id"] == rfq_id:
                return (qid, row["yes_bid"], row["no_bid"])
        return None

    # -- cooldown (H9/P-7) + per-combo fill exposure (H8/N7/N8) --------- #

    def has_cooldown_row(self, ticker: str) -> bool:
        return ticker in self._cooldowns

    def cooldown_active(self, ticker: str, now_utc) -> bool:
        cooled_until = self._cooldowns.get(ticker)
        return cooled_until is not None and cooled_until > now_utc

    def last_fill(self, ticker: str):
        """(filled_at, fair_at_confirm|blended_fair_at_quote) of the latest
        fill on this combo, or None. Loaded only for cooled tickers — the
        only callers are behind has_cooldown_row."""
        return self._last_fills.get(ticker)

    def combo_fill_exposure(self, ticker: str) -> float:
        return self._combo_exposure.get(ticker, 0.0)

    # -- in-pass maintenance (read-your-writes for the cap gates) ------- #

    def record_submitted(self, quote_id: str, rfq_id: str, ticker: str,
                         primary_game_id, game_ids: list,
                         worst_exposure_usd: float,
                         yes_bid: float, no_bid: float) -> None:
        self._open_quotes[quote_id] = dict(
            rfq_id=rfq_id, ticker=ticker, game_id=primary_game_id,
            worst_exposure_usd=float(worst_exposure_usd),
            yes_bid=yes_bid, no_bid=no_bid)
        self._quote_games[quote_id] = list(game_ids)

    def remove_quote(self, quote_id: str) -> None:
        """A quote left 'open' mid-pass (replaced / cancelled / pulled)."""
        self._open_quotes.pop(quote_id, None)
        self._quote_games.pop(quote_id, None)


def _open_quote_exposure_rows(exclude_rfq_id: str | None = None):
    """One-off wrapper over `_PassSnapshot.open_exposure_rows` (the pass
    itself loads one snapshot and reuses it — issue #89)."""
    snapshot = _PassSnapshot.load(datetime.now(timezone.utc))
    return snapshot.open_exposure_rows(exclude_rfq_id)


def _open_quote_exposure_by_game(exclude_rfq_id: str | None = None):
    """One-off wrapper over `_PassSnapshot.open_exposure_by_game` (see
    `_open_quote_exposure_rows`)."""
    snapshot = _PassSnapshot.load(datetime.now(timezone.utc))
    return snapshot.open_exposure_by_game(exclude_rfq_id)


def _post_fill_live_refresh_landed(by_game: dict, filled_at) -> bool:
    """True iff EVERY game leg set of the combo has a COMPLETED on-demand
    flight with >=1 book fair that landed AFTER `filled_at` (P-7 issue #21,
    re-pointed by #57). The COMBO_COOLDOWN_SEC floor (60s) alone would let a
    re-quote be priced off the exact book snapshot that produced the
    picked-off fill; the old check waited for a full sweep pass, which at
    the slow cadence would stretch every cooldown to minutes — the gate now
    keys on the TARGETED fetch `_ensure_post_fill_fetches` fires instead.

    The engine stamps landings on the monotonic clock; the age converts to a
    wall instant here for the comparison against `fills.filled_at`
    (TIMESTAMPTZ — aware; a naive value is treated as session-local,
    mirroring `risk._now_matching`). Fail-CLOSED: no engine, no completed
    priced flight, a pre-fill landing, or any comparison error keeps the
    combo cooled — a stale re-quote is worse than a skipped one.
    """
    try:
        if filled_at is None or _ENGINE is None:
            return False
        now_utc = datetime.now(timezone.utc)
        filled_cmp = (filled_at if filled_at.tzinfo is not None
                      else filled_at.astimezone())
        for gl in by_game.values():
            age = _ENGINE.completed_fetch_age_sec(legset.leg_set_hash(gl))
            if age is None:
                return False
            landed_wall = now_utc - timedelta(seconds=age)
            if landed_wall <= filled_cmp:
                return False
        return True
    except Exception:
        log.warning("[post_fill_refresh_check_failed] — keeping combo cooled")
        return False


def _ensure_post_fill_fetches(rfq_id, ticker, by_game: dict) -> None:
    """Queue a TARGETED on-demand fetch for every game leg set of a filled
    combo (#57). Fired eagerly by the confirm tick the moment a fill is
    recorded (so the awaiting-refresh gate usually clears right at the 60s
    floor) and lazily by the cooldown gate whenever the landing is missing
    (process restart / engine store prune self-heal). Reuses the live-feed
    machinery wholesale — per-book pacing gates, job budget, dropped-book
    accounting, #38 health rows. Fail-safe: never raises into its tick."""
    try:
        if _ENGINE is None:
            return
        for gl in by_game.values():
            gid = _resolve_game_for_legs(gl)
            gref = _game_ref(gid) if gid is not None else None
            if gref is None:
                continue
            h = legset.leg_set_hash(gl)
            if _ENGINE.landed_empty(h):
                # The books just answered "nothing" for this leg set — while
                # that zero-book landing is fresh, re-feeding would hammer a
                # dead slate every tick. Retry resumes once it ages out
                # (~QUOTE_FRESH_SEC), mirroring the pricing path's re-feed
                # rule after live_fetch_timeout.
                continue
            if _ENGINE.ensure_fetch(h, gref, gl):
                research.emit("on_demand_requested", rfq_id=rfq_id,
                              ticker=ticker,
                              payload=dict(leg_set_hash=h, game_id=gid,
                                           n_legs=len(gl),
                                           trigger="post_fill"))
    except Exception:
        log.warning("[post_fill_fetch_failed] ticker=%s", ticker)


def _resolve_game_for_legs(game_legs: list) -> str | None:
    """Cached wrapper (O-1) — keyed by the game code; failures (None)
    are not cached; invalidated on SGP refresh."""
    if not game_legs:
        return None
    game_code = game_legs[0].game_id
    if game_code in _RESOLVE_CACHE:
        return _RESOLVE_CACHE[game_code]
    gid = _resolve_game_for_legs_uncached(game_legs)
    if gid is not None:
        _RESOLVE_CACHE[game_code] = gid
    return gid


def _resolve_game_for_legs_uncached(game_legs: list) -> str | None:
    """Resolve ONE game's CanonicalLegs to the mlb_target_lines game_id, or None.

    game_legs[0].game_id is the family-independent GAME CODE shared by all legs
    of that game (issue #71 — it used to be the whole event ticker, which
    differs per market family). Parses it with _parse_event_suffix +
    _MLB_CODE_TO_TEAM and looks up mlb_target_lines. Fail-safe: returns None on
    any error — never raises.
    """
    try:
        if not game_legs:
            return None
        # Already the trailing code; tolerate a stray prefix rather than
        # returning None, so a future ticker shape can't silently unprice a game.
        suffix = game_legs[0].game_id.rsplit("-", 1)[-1]
        if not suffix:
            return None
        away_code, home_code = _parse_event_suffix(suffix)
        if not away_code or not home_code:
            return None
        home_name = _MLB_CODE_TO_TEAM.get(home_code)
        away_name = _MLB_CODE_TO_TEAM.get(away_code)
        if not home_name or not away_name or not config.MARKET_DB.exists():
            return None
        try:
            con = duckdb.connect(str(config.MARKET_DB), read_only=True)
        except (duckdb.IOException, duckdb.CatalogException):
            return None
        try:
            row = con.execute(
                "SELECT game_id FROM mlb_target_lines "
                "WHERE home_team=? AND away_team=? LIMIT 1",
                [home_name, away_name]).fetchone()
            return row[0] if row else None
        except (duckdb.IOException, duckdb.CatalogException):
            return None
        finally:
            con.close()
    except Exception:
        return None


def _game_ref(game_id):
    """Cached wrapper (O-1 pattern, issue #89) — the pending-RFQ feed block
    calls this every tick per pending sub-combo, and the uncached MARKET_DB
    open added up at ~200 open RFQs. Failures (None) are not cached;
    invalidated with the other game-metadata caches on target-line refresh."""
    if game_id in _GAME_REF_CACHE:
        return _GAME_REF_CACHE[game_id]
    ref = _game_ref_uncached(game_id)
    if ref is not None:
        _GAME_REF_CACHE[game_id] = ref
    return ref


def _game_ref_uncached(game_id):
    """GameRef (id + team names + commence) from mlb_target_lines, or None.

    The per-book on-demand fetch needs team names to match the book's own
    event list. Fail-safe: None on any error — the RFQ then just stays
    on_demand_pending (never raises into the tick)."""
    try:
        from mlb_sgp._shared import GameRef
        if not game_id or not config.MARKET_DB.exists():
            return None
        try:
            con = duckdb.connect(str(config.MARKET_DB), read_only=True)
        except (duckdb.IOException, duckdb.CatalogException):
            return None
        try:
            row = con.execute(
                "SELECT home_team, away_team, commence_time "
                "FROM mlb_target_lines WHERE game_id=? LIMIT 1",
                [game_id]).fetchone()
        except (duckdb.IOException, duckdb.CatalogException):
            return None
        finally:
            con.close()
        if not row:
            return None
        return GameRef(game_id=game_id, home_team=row[0], away_team=row[1],
                       commence_time=row[2])
    except Exception:
        return None


def _out_of_scope_reason(legs, canon) -> str:
    """Granular out-of-scope classification (measurement blind-spot fix)."""
    if not legs:
        return "out_of_scope_unparseable"
    if canon is None:
        return "out_of_scope_non_mlb"      # some leg untypeable (other sport etc.)
    if len(canon) < 2:
        return "out_of_scope_lone_single"
    # F5-winner TEAM legs are in scope since #86 (re-encoded as +-0.5 run
    # lines at parse). Only the TIE market still parses as (ml, F5) — no
    # book leg can express it. Labeled distinctly so the declined demand
    # stays measurable.
    if any(l.market_type == "ml" and l.period == "F5" for l in canon):
        return "out_of_scope_f5_tie_leg"
    return "out_of_scope"


def _maybe_emit_on_demand_result(rfq_id, ticker, hash_, legs=()):
    """Emit on_demand_result once per landing (not per tick it gets read)."""
    try:
        if _ENGINE is None:
            return
        landed = _ENGINE.landed_at(hash_)
        if landed is None or _OD_RESULT_EMITTED.get(hash_) == landed:
            return
        _OD_RESULT_EMITTED[hash_] = landed
        # bounded: prune oldest entries (insertion-ordered dict) — days-long
        # uptime with real flow must not leak (adversarial review #9)
        while len(_OD_RESULT_EMITTED) > 512:
            _OD_RESULT_EMITTED.pop(next(iter(_OD_RESULT_EMITTED)))
        res = _ENGINE.lookup_results(hash_) or {}
        # #55: budget-dropped books ride along (getattr — duck-typed test
        # engines may not implement the accessor). latency_sec per book is
        # the arrival-order/latency record the quorum design is judged by.
        dropped = getattr(_ENGINE, "dropped_books", lambda h: ())(hash_)
        research.emit("on_demand_result", rfq_id=rfq_id, ticker=ticker,
                      payload=dict(
                          leg_set_hash=hash_,
                          # #85: which periods the flight's legs span
                          # ("FG"/"F5") — additive key, report.py readers
                          # select fields by name and ignore it.
                          periods=sorted({l.period for l in legs}),
                          dropped_books=list(dropped or ()),
                          books={b: dict(fair=r.fair, route=r.route,
                                         n_cells=r.n_cells_priced,
                                         latency_sec=r.latency_sec)
                                 for b, r in res.items()}))
    except Exception:                      # research must never break the tick
        pass


def _live_games_detail(by_game):
    """Per-game live-fetch trace for quote_priced (issue #54 acceptance: a
    live-mode quote's fair must be provably backed by a post-RFQ fetch).
    {hash: {age_sec, books: {book: {fair, route, latency_sec}}}} for every
    game with a stored engine result; None when there is none (cached mode).
    Fail-safe: research decoration only, never raises into the tick."""
    try:
        if _ENGINE is None:
            return None
        out = {}
        for gl in by_game.values():
            h = legset.leg_set_hash(gl)
            res = _ENGINE.lookup_results(h)
            if not res:
                continue
            out[h] = dict(
                age_sec=_ENGINE.result_age_sec(h),
                # #85: periods the game's legs span ("FG"/"F5") — additive
                # key, safe for report.py's live_games readers.
                periods=sorted({l.period for l in gl}),
                books={b: dict(fair=r.fair, route=r.route,
                               latency_sec=r.latency_sec)
                       for b, r in res.items()})
        return out or None
    except Exception:
        return None


def _on_demand_fill_info(canon):
    """Per-on-demand-game route/consensus detail for the fill research event."""
    try:
        if _ENGINE is None or not canon:
            return None
        out = {}
        for gl in legset.partition_by_game(canon).values():
            if legset.classify_subcombo(gl) != "on_demand":
                continue
            h = legset.leg_set_hash(gl)
            res = _ENGINE.lookup_results(h)
            if not res:
                out[h] = None
                continue
            fairs = {b: r.fair for b, r in res.items()}
            det = router.consensus_detail(fairs, config.MIN_AGREEING_BOOKS,
                                          config.SIGMA_Z_MAX)
            out[h] = dict(
                books={b: dict(fair=r.fair, route=r.route) for b, r in res.items()},
                consensus_books=det[1] if det else [])
        return out or None
    except Exception:
        return None


def _quote_age_sec(submitted_at) -> float | None:
    """Seconds between quote submission and now, for confirm_singles_check.
    O-2 convention: submitted_at is TIMESTAMPTZ (aware) post-migration; a
    naive value (pre-migration row) is naive-LOCAL, compared accordingly.
    Fail-safe: any surprise shape returns None rather than raising."""
    try:
        if submitted_at is None:
            return None
        if submitted_at.tzinfo is not None:
            return max(0.0, (datetime.now(timezone.utc) - submitted_at).total_seconds())
        return max(0.0, (datetime.now() - submitted_at).total_seconds())
    except Exception:
        return None


def _fill_game_ids(legs, primary_game_id):
    """All game_ids a filled combo touches, for per-game exposure attribution.

    Parses `legs`, partitions by game, resolves each. Fail-safe: ALWAYS includes
    `primary_game_id`, and the whole derivation is wrapped so ANY unexpected error
    degrades to primary-only rather than raising. A recorded fill must never be
    lost from the per-game cap, and this function must never throw into the fill
    path — so callers compute it BEFORE the DB write block (see `_confirm_tick`),
    guaranteeing `fills` is never written without matching `fill_games` rows.
    """
    ids = set()
    try:
        canon = legset.parse_legs(legs) if legs else None
        if canon:
            for gl in legset.partition_by_game(canon).values():
                g = _resolve_game_for_legs(gl)
                if g:
                    ids.add(g)
    except Exception:  # fail-safe: degrade to primary-only, never raise
        ids = set()
    if primary_game_id:
        ids.add(primary_game_id)
    return sorted(ids)


def _priceable_in_phase1(canon: list) -> bool:
    """True iff every game's sub-combo is priceable via Phase-1 grid routes.

    Phase 1 = 2-leg spread×total / ml×total grid cells OR cross-game combos
    (≥2 legs across games). Lone single-leg RFQs are explicitly excluded: a
    correlated-combo maker quoting a straight single is unintended.
    on_demand (3+ legs same game) and unpriceable (0 legs) are not yet supported.
    """
    if len(canon) < 2:
        return False
    by_game = legset.partition_by_game(canon)
    if not by_game:
        return False
    return all(
        legset.classify_subcombo(gl) in {"single", "grid_spread_total", "grid_ml_total"}
        for gl in by_game.values()
    )


def _priceable(canon: list) -> bool:
    """Phase 2 scope: every game's sub-combo has SOME pricing route.

    Adds "on_demand" to the Phase 1 routes (live book queries via the
    OnDemandEngine). The RFQ-level >=2-legs rule stays — lone singles remain
    excluded (unchanged Phase 1 decision). _priceable_in_phase1 above is kept
    verbatim as the Phase 1 regression oracle (test_router_integration.py).
    """
    if len(canon) < 2:
        return False
    by_game = legset.partition_by_game(canon)
    if not by_game:
        return False
    return all(
        legset.classify_subcombo(gl) in
        {"single", "grid_spread_total", "grid_ml_total", "on_demand"}
        for gl in by_game.values()
    )


def _get_position_contracts(market_ticker: str, timeout: int = 5) -> int | None:
    """Authoritative current position via /portfolio/positions. Returns a SIGNED
    int (positive = long YES, negative = long NO). None on API failure or timeout.
    Short default timeout so a stuck positions call can't block the confirm loop.
    Mirrors kalshi_mlb_rfq.rfq_client.get_position_contracts."""
    try:
        status, body, _ = auth_client.api(
            "GET", f"/portfolio/positions?ticker={market_ticker}&limit=10",
            timeout=timeout)
    except Exception:
        return None
    if status != 200 or not isinstance(body, dict):
        return None
    for p in body.get("market_positions") or []:
        if p.get("ticker") == market_ticker:
            fp = p.get("position_fp")
            if fp is not None:
                return int(float(fp))
            return int(p.get("position", 0))
    return 0


def _rfq_requested_contracts(rfq: dict) -> float:
    """Requested size in contracts, from the live RFQ shape.

    Kalshi returns `contracts_fp` as a fixed-point STRING (e.g. "21.00") —
    verified from live rfq_received payloads 2026-06-08; the legacy int
    `contracts` field is kept as a fallback. Returns 0.0 when absent or zero,
    which means the RFQ is dollar-denominated (`target_cost_dollars`) and the
    size gate must run post-pricing instead (see _worst_fill_exposure_usd)."""
    fp = rfq.get("contracts_fp")
    if fp is not None:
        try:
            return float(fp)
        except (TypeError, ValueError):
            pass
    try:
        return float(rfq.get("contracts") or 0)
    except (TypeError, ValueError):
        return 0.0


def _accept_fill_contracts(q: dict | None) -> float | None:
    """Filled size from the quote-status response, or None when unparseable.

    Kalshi live payloads denominate sizes as fixed-point STRINGS
    (`contracts_fp`) — mirror _rfq_requested_contracts: parse that first,
    fall back to the legacy int `contracts`. None means "size unknown";
    callers must record conservatively (0 + reconciled=FALSE), never 1."""
    if not isinstance(q, dict):
        return None
    fp = q.get("contracts_fp")
    if fp is not None:
        try:
            parsed = float(fp)
            if parsed > 0:
                return parsed
        except (TypeError, ValueError):
            pass
    legacy = q.get("contracts")
    if legacy is not None:
        try:
            parsed = float(legacy)
            if parsed > 0:
                return parsed
        except (TypeError, ValueError):
            pass
    return None


def _worst_fill_exposure_usd(rfq: dict, q) -> float:
    """Worst-case DOLLARS at risk if this RFQ fills at quote q, maxed over both
    sides the creator could take. Our cost basis on the side we'd hold IS our
    max loss (a binary settles to 0/1). Mirrors the confirm-path side logic:
      accepted_side=='no'  -> we hold YES at yes_ask = 1 - no_bid
      accepted_side=='yes' -> we hold NO  at no_ask  = 1 - yes_bid
    """
    yes_ask = 1.0 - q.no_bid   # our cost if we end up holding YES
    no_ask = 1.0 - q.yes_bid   # our cost if we end up holding NO
    n = _rfq_requested_contracts(rfq)
    if n > 0:
        # contract-denominated: we hold n of whichever side; worst = pricier side
        return n * max(yes_ask, no_ask)
    target = float(rfq.get("target_cost_dollars") or 0)
    if target > 0:
        # dollar-denominated: creator spends `target` on their side at their ask;
        # contracts = target / creator_ask; we hold the opposite at our_ask.
        # creator buys YES (pays yes_ask) -> we hold NO  : (target/yes_ask)*no_ask
        # creator buys NO  (pays no_ask)  -> we hold YES : (target/no_ask)*yes_ask
        if yes_ask <= 0 or no_ask <= 0:
            return float("inf")  # degenerate quote -> always skip
        return max((target / yes_ask) * no_ask, (target / no_ask) * yes_ask)
    return 0.0


_DECISION_INSERT_SQL = (
    "INSERT INTO quote_decisions (decision_id, rfq_id, quote_id, "
    "combo_market_ticker, game_id, decision, reason, model_fair, book_fair, "
    "blended_fair, yes_bid, no_bid, observed_at) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)")


def _log_decision(decision, *, rfq_id=None, quote_id=None, ticker=None, game_id=None,
                  reason=None, model=None, book=None, blended=None, yb=None, nb=None,
                  buffer: list | None = None):
    """Record a quote decision. With buffer=None (default) the row is written
    immediately on its own connection. The discovery pass instead passes its
    per-pass buffer: on the ~1GB state DB EVERY connection close pays a full
    CHECKPOINT (~0.5s, sampled 2026-08-10), so thousands of per-RFQ writes
    must collapse into one flush per pass or the tick starves the loop."""
    row = [str(uuid.uuid4()), rfq_id, quote_id, ticker, game_id, decision, reason,
           model, book, blended, yb, nb, datetime.now(timezone.utc)]
    if buffer is not None:
        buffer.append(row)
    else:
        with db.connect() as con:
            con.execute(_DECISION_INSERT_SQL, row)
    # Event 3: decision — research mirror of every quote_decisions row.
    research.emit("decision", rfq_id=rfq_id, quote_id=quote_id, ticker=ticker,
                  payload=dict(decision=decision, game_id=game_id, reason=reason,
                               model=model, book=book, blended=blended, yb=yb, nb=nb))


def _void_rate_halt_triggered() -> bool:
    """H4 (global): if voided/(voided+confirmed) > threshold over the recent
    window, halt new quotes. A high void rate means we're either getting picked
    off (book moves we cancel on) or chronically can't re-price — both are
    signals to step back instead of farming volume."""
    cutoff = datetime.now(timezone.utc) - timedelta(hours=config.VOID_RATE_WINDOW_HOURS)
    with db.connect(read_only=True) as con:
        row = con.execute(
            "SELECT "
            "  SUM(CASE WHEN decision LIKE 'voided_%' THEN 1 ELSE 0 END), "
            "  SUM(CASE WHEN decision='confirmed' THEN 1 ELSE 0 END) "
            "FROM quote_decisions WHERE observed_at >= ?",
            [cutoff]).fetchone()
    voided = int(row[0] or 0)
    confirmed = int(row[1] or 0)
    denom = voided + confirmed
    if denom == 0:
        return False
    return (voided / denom) > config.VOID_RATE_HALT_THRESHOLD


def _pull_quorum_quote(gateway, rid, ticker, game_id, snapshot) -> bool:
    """#55 straggler pull: a late book blew the #20 dispersion gate while a
    quote from an earlier (thinner) quorum rests on this RFQ. The resting
    price was formed from LESS information than we now have, and what we
    have says the books disagree — the stale side of that disagreement is
    us, so the quote comes down. Returns True iff a resting quote was
    pulled (the caller then skips the ordinary gate-decline log)."""
    existing = snapshot.existing_open_quote(rid)
    if existing is None:
        return False
    qid = existing[0]
    snapshot.remove_quote(qid)
    try:
        gateway.cancel(qid)
    except Exception:
        pass
    with db.connect() as con:
        con.execute(
            "UPDATE live_quotes SET status='cancelled', closed_at=? WHERE quote_id=?",
            [datetime.now(timezone.utc), qid])
    _log_decision("pulled", rfq_id=rid, quote_id=qid, ticker=ticker,
                  game_id=game_id, reason="quorum_dispersion_bust")
    research.emit("quote_pulled", rfq_id=rid, ticker=ticker,
                  payload=dict(quote_id=qid, game_id=game_id,
                               reason="quorum_dispersion_bust"))
    return True


def _discovery_tick(source, gateway, dry_run):
    global _VOID_HALT_ACTIVE
    # FIX M2: kill-switch — stop quoting immediately if the kill file exists
    if config.KILL_FILE.exists():
        return
    # No freshness top gate (#57/#81): every RFQ's data need is enforced
    # per-RFQ by the live-feed block below (on_demand_pending /
    # live_fetch_timeout / live_too_few_books) — there is no background
    # sweep, so nothing market-wide exists to gate on.
    # H4 + N12: global void-rate halt — step back if we're voiding too often.
    # N12: edge-triggered notification so operator gets a single ping on state
    # transitions (False→True = halt; True→False = resume) rather than silence.
    void_halt_now = _void_rate_halt_triggered()
    if void_halt_now and not _VOID_HALT_ACTIVE:
        # Transition False → True: halt started.
        cutoff = datetime.now(timezone.utc) - timedelta(hours=config.VOID_RATE_WINDOW_HOURS)
        with db.connect(read_only=True) as con:
            row = con.execute(
                "SELECT "
                "  SUM(CASE WHEN decision LIKE 'voided_%' THEN 1 ELSE 0 END), "
                "  SUM(CASE WHEN decision='confirmed' THEN 1 ELSE 0 END) "
                "FROM quote_decisions WHERE observed_at >= ?",
                [cutoff]).fetchone()
        voided = int(row[0] or 0)
        confirmed = int(row[1] or 0)
        denom = voided + confirmed
        rate = voided / denom if denom else 0.0
        notify.halt("void_rate",
                    detail=f"void_rate={rate:.0%} in last {config.VOID_RATE_WINDOW_HOURS}h")
        # Event 7: halt_event (fire transition)
        research.emit("halt_event",
                      payload=dict(halt_type="void_rate", transition="fire",
                                   current_rate=rate,
                                   window_hours=config.VOID_RATE_WINDOW_HOURS))
    elif not void_halt_now and _VOID_HALT_ACTIVE:
        # Transition True → False: halt lifted.
        notify.resume("void_rate", detail="void_rate fell back below threshold")
        # Event 7: halt_event (clear transition)
        research.emit("halt_event",
                      payload=dict(halt_type="void_rate", transition="clear",
                                   current_rate=0.0,
                                   window_hours=config.VOID_RATE_WINDOW_HOURS))
    _VOID_HALT_ACTIVE = void_halt_now
    if void_halt_now:
        _log_decision("halted_high_void_rate")
        return
    rfqs = source.poll()
    now_utc = datetime.now(timezone.utc)
    # Issue #89: ONE read connection loads every per-RFQ gate input (open
    # quotes, cooldowns, per-combo/creator fill aggregates). The per-RFQ
    # db.connect() opens this replaces were ~61ms/RFQ — the pass ran 8-12s
    # at ~200 open RFQs, which IS the quote-latency regression (#55's
    # quorum landings sat unread until the next slow pass).
    snapshot = _PassSnapshot.load(now_utc)
    open_count = len(snapshot.open_exposure_rows())
    # FIX I3: load today's fills once before the loop for cap checks
    fills_today = _today_fills()
    # Per-game exposure fanned out across every game each combo touches (cross-game
    # combos count against BOTH games). Daily cap keeps using fills_today (one row
    # per combo) so it is never double-counted; only the per-game cap uses this.
    game_exposure_rows = _today_fills_by_game()
    scope_fetches_this_tick = 0
    scope_fetches_deferred = 0
    # Per-pass write buffers + one bulk seen-lookup (2026-08-10 incident,
    # round 2): on the ~1GB state DB every connection close pays a full
    # CHECKPOINT (~0.5s — confirmed by stack sampling), so the old
    # per-RFQ seen-SELECT + per-decision INSERT pattern cost the pass
    # O(thousands) of checkpoints and starved the loop even after scope
    # checks went wire-free. One read + one write connection per pass.
    decision_rows: list[list] = []
    seen_insert_rows: list[list] = []

    def _decide(decision, **kw):
        _log_decision(decision, buffer=decision_rows, **kw)

    pass_rfq_ids = [r.get("id") for r in rfqs if r.get("id")]
    seen_ids: set[str] = set()
    if pass_rfq_ids:
        with db.connect(read_only=True) as con:
            seen_ids = {row[0] for row in con.execute(
                "SELECT rfq_id FROM seen_rfqs WHERE rfq_id IN (SELECT unnest(?))",
                [pass_rfq_ids]).fetchall()}
    # Time-boxed, rotating pass (2026-08-10 round 3): a Sunday-peak mirror
    # held 9k+ open RFQs — no per-RFQ cost is low enough that an unbounded
    # pass can't starve the loop. The lap gets a wall-clock budget with
    # guaranteed progress (>=1 RFQ), and starts where the last lap left off
    # so tail RFQs are not systematically starved.
    global _DISCOVERY_CURSOR
    start_at = _DISCOVERY_CURSOR % len(rfqs) if rfqs else 0
    if start_at:
        rfqs = rfqs[start_at:] + rfqs[:start_at]
    pass_started = time.monotonic()
    pass_deadline = pass_started + config.DISCOVERY_PASS_BUDGET_SEC
    iterated = 0
    timebox_deferred = 0
    stale_skipped = 0
    tipoff_skipped = 0
    horizon_skipped = 0
    try:
        for rfq in rfqs:
            if iterated and time.monotonic() > pass_deadline:
                timebox_deferred = len(rfqs) - iterated
                break
            iterated += 1
            # A poll() snapshot can be the entire exchange-wide open-RFQ set
            # (thousands of tickers) — without this check a SIGTERM only takes
            # effect after the full pass, which is the 2026-08-09/10 shutdown
            # wedge and loop starvation. Bail promptly; nothing here is
            # half-done (each RFQ is handled atomically).
            if not _running.is_set():
                break
            if open_count >= config.MAX_OPEN_QUOTES:
                break
            rid = rfq.get("id")
            ticker = rfq.get("market_ticker")
            if not rid or not ticker:
                continue
            # Age gate BEFORE any other work: an RFQ resting past
            # MAX_RFQ_AGE_SEC un-quoted is near expiry or already picked
            # over — quoting it late is adverse selection, and classifying
            # it is wasted lap budget. No decision row (it would re-bloat
            # the state DB on every restart backlog); the periodic pass
            # summary carries the count. Unparseable/missing created_ts
            # passes through (fail-open).
            rfq_created = _parse_created_ts(rfq.get("created_ts"))
            if rfq_created is not None and (
                    (now_utc - rfq_created).total_seconds() > config.MAX_RFQ_AGE_SEC):
                stale_skipped += 1
                continue
            # Kalshi RFQ poll responses tag the creator as creator_user_id; fall
            # back to creator_id for forward-compatibility if the API name changes.
            creator_id = rfq.get("creator_user_id") or rfq.get("creator_id") or ""
            # scope (cache the market lookup verdict). Resolved BEFORE the
            # seen_rfqs read so a budget-deferred ticker costs zero DB opens —
            # during a 4k-backlog drain the deferred majority is re-visited
            # every tick, and a per-RFQ connect there is ~400k pointless opens.
            if ticker in _SCOPE_CACHE:
                in_scope, game_id, legs = _SCOPE_CACHE[ticker]
                canon = legset.parse_legs(legs) if legs else None
            elif (rfq_legs := scope.decode_legs(rfq)) is not None:
                # The RFQ payload itself carries mve_selected_legs (verified
                # against WS rfq_created frames and the REST gap-fill,
                # 2026-08-10) in the same {event_ticker, market_ticker, side}
                # shape as GET /markets/{ticker} — so scope resolves with ZERO
                # wire calls. The budgeted get_market below survives only as a
                # fallback for payloads missing the field.
                legs = rfq_legs
                canon = legset.parse_legs(legs)
                in_scope = bool(canon and _priceable(canon))
                game_id = None
                _scope_cache_put(ticker, in_scope, game_id, legs)
            else:
                if scope_fetches_this_tick >= config.SCOPE_FETCH_BUDGET_PER_TICK:
                    # Budget exhausted: leave this ticker unresolved for a later
                    # tick (the mirror is level-triggered, so it will be offered
                    # again). This bounds a tick's REST work so a giant backlog
                    # drains across many breathing loop iterations instead of
                    # one starving mega-pass.
                    scope_fetches_deferred += 1
                    continue
                market = source.get_market(ticker)
                scope_fetches_this_tick += 1
                if market is None:
                    # B5 fix (issue #18): a transient fetch failure (HTTP blip)
                    # must not be cached as in_scope=False — that permanently
                    # blacklisted the combo until restart — nor recorded as a
                    # decided RFQ. Skip this tick; a later poll retries.
                    log.warning("[market_fetch_failed] ticker=%s — will retry "
                                "next tick", ticker)
                    continue
                legs = scope.decode_legs(market)
                canon = legset.parse_legs(legs) if legs else None
                in_scope = bool(canon and _priceable(canon))
                game_id = None
                _scope_cache_put(ticker, in_scope, game_id, legs)
            seen = rid in seen_ids
            if not in_scope:
                if not seen:
                    # Phase 2 instrumentation: granular reason + the legs
                    # themselves, so out-of-scope flow is finally measurable
                    # (previously legs_json was NULL and every reason was the
                    # single opaque "out_of_scope").
                    oos_reason = _out_of_scope_reason(legs, canon)
                    _decide("skipped", rfq_id=rid, ticker=ticker, reason=oos_reason)
                    seen_insert_rows.append(
                        [rid, ticker, False, None,
                         json.dumps(legs) if legs else None, now_utc,
                         oos_reason, creator_id])
                    seen_ids.add(rid)
                continue
            # Event 1: rfq_received — in-scope candidate, first time we act on it.
            if not seen:
                research.emit("rfq_received", rfq_id=rid, ticker=ticker,
                              payload=dict(rfq_keys=list(rfq.keys()), rfq_raw=rfq))
            # H4 (per-creator): refuse to quote a counterparty that's already
            # farmed us this window. Check BEFORE size gate so we don't waste cycles.
            if snapshot.creator_halt_active(creator_id):
                _decide("skipped", rfq_id=rid, ticker=ticker, reason="skipped_creator_halt")
                # Event 10: creator_halt_skip
                research.emit("creator_halt_skip", rfq_id=rid, ticker=ticker,
                              payload=dict(creator_id=creator_id,
                                           fill_count_window=snapshot.creator_fill_count(creator_id),
                                           window_hours=config.PER_CREATOR_WINDOW_HOURS))
                continue
            # Size gate. Kalshi RFQs are denominated EITHER in contracts
            # (contracts_fp, fixed-point string) OR in dollars (target_cost_dollars,
            # with contracts_fp == "0.00"). Live data 2026-06-08: ~85% of in-scope
            # flow was dollar-denominated, so both paths matter. The per-fill dollar
            # cap runs post-pricing (needs the quote — see _worst_fill_exposure_usd).
            req_contracts = _rfq_requested_contracts(rfq)
            target_cost = float(rfq.get("target_cost_dollars") or 0)
            if req_contracts <= 0 and target_cost <= 0:
                _decide("skipped", rfq_id=rid, ticker=ticker, reason="size_unknown")
                continue
            # Multi-game resolution: canon must be parseable and all games resolvable.
            if canon is None:
                _decide("skipped", rfq_id=rid, ticker=ticker, reason="no_game")
                continue
            by_game = legset.partition_by_game(canon)
            game_ids_list = [_resolve_game_for_legs(gl) for gl in by_game.values()]
            if any(gid is None for gid in game_ids_list):
                _decide("skipped", rfq_id=rid, ticker=ticker, reason="no_game")
                continue
            game_id = game_ids_list[0]  # primary id stored in schema (single column)
            # FIX I3: exposure cap gates (wired — were defined but never called).
            # R-1 (issue #22): both gates ALSO count open quotes' worst-case
            # exposure — a burst sweep of resting quotes could otherwise fill to
            # many multiples of the caps before the first fill registers. Reads
            # the pass snapshot, which record_submitted keeps current, so quotes
            # submitted earlier in this same loop count immediately (N7).
            open_quote_rows = snapshot.open_exposure_rows(exclude_rfq_id=rid)
            if not risk.daily_cap_ok(fills_today + open_quote_rows,
                                     config.daily_exposure_cap_usd()):
                _decide("skipped", rfq_id=rid, ticker=ticker, game_id=game_id,
                              reason="daily_cap")
                continue
            # Per-game cap: all games must pass (cross-game combos checked per-game).
            # game_exposure_rows attributes each combo's full stake to EVERY game it
            # touches, so a game can't slip past its cap by only ever being a 2nd leg;
            # open_quotes_by_game does the same for still-resting quotes via quote_games.
            open_quotes_by_game = snapshot.open_exposure_by_game(exclude_rfq_id=rid)
            if any(not risk.per_game_cap_ok(gid, game_exposure_rows + open_quotes_by_game,
                                            config.BANKROLL, config.MAX_GAME_EXPOSURE_PCT)
                   for gid in game_ids_list):
                _decide("skipped", rfq_id=rid, ticker=ticker, game_id=game_id,
                              reason="per_game_cap")
                continue
            # H9 + P-7 (issue #21, re-pointed by #57): per-combo cooldown, two
            # layers after a fill:
            # (1) hard floor — COMBO_COOLDOWN_SEC past the fill (unchanged);
            # (2) data-aware — even after the floor, stay cooled until EVERY game
            #     this combo touches has a completed post-fill TARGETED fetch
            #     (fired at fill time by the confirm tick, re-fired lazily here).
            #     Without it a re-quote could be priced before the books were
            #     ever re-asked after the picked-off fill.
            # cooled_until is TIMESTAMPTZ (O-2), so the aware bind compares as a
            # true instant. last_fill_fair feeds the same-price block below (only
            # set once the refresh gate has passed).
            last_fill_fair = None
            if snapshot.has_cooldown_row(ticker):
                if snapshot.cooldown_active(ticker, now_utc):
                    _decide("skipped", rfq_id=rid, ticker=ticker, game_id=game_id,
                                  reason="in_cooldown")
                    continue
                last_fill = snapshot.last_fill(ticker)
                if last_fill is not None:
                    if not _post_fill_live_refresh_landed(by_game, last_fill[0]):
                        _ensure_post_fill_fetches(rid, ticker, by_game)
                        _decide("skipped", rfq_id=rid, ticker=ticker,
                                      game_id=game_id,
                                      reason="in_cooldown_awaiting_refresh")
                        continue
                    last_fill_fair = last_fill[1]
            # Tipoff gate: earliest commence across all games. Counter, not a
            # decision row — in-play RFQ flow runs ~100/min at evening peak
            # and per-pass rows at that rate is the DB bloat the 08-11 prune
            # removed (pass-summary log + rfq_ingestion_summary carry it).
            commence_times = [_commence_time(gid) for gid in game_ids_list]
            earliest_ct = min((ct for ct in commence_times if ct is not None), default=None)
            if not risk.tipoff_ok(earliest_ct, config.TIPOFF_CANCEL_MIN):
                tipoff_skipped += 1
                continue
            # Flight horizon: books price SGP combos only near game time —
            # a combo whose farthest game is hours out comes back
            # too_few_books after 6 wasted fetches. Skip (counter only)
            # until every game is inside the horizon. 0 disables.
            if config.FLIGHT_HORIZON_HOURS > 0:
                latest_ct = max((ct for ct in commence_times if ct is not None),
                                default=None)
                if latest_ct is not None and (
                        (latest_ct - now_utc).total_seconds()
                        > config.FLIGHT_HORIZON_HOURS * 3600):
                    horizon_skipped += 1
                    continue
            # Live pricing (#54): EVERY in-scope sub-combo — grids, cross-game
            # singles, novel shapes — rides a live feed keyed off this very poll.
            # A fresh (<=QUOTE_FRESH_SEC) result prices this tick; otherwise
            # ensure a fetch is queued and skip. The tick re-enters every open
            # RFQ every 2s, so this single rule IS the feed: re-fetch fires each
            # time the result ages out, and stops the moment the RFQ leaves the
            # poll. The sweep cache is never in the quote path: a fetch that
            # lands with ZERO books is a DECLINE (live_fetch_timeout), and while
            # that empty landing is fresh we deliberately do not re-feed, so a
            # dead slate retries every ~QUOTE_FRESH_SEC. Placed after the
            # caps/cooldown/tipoff gates so gated RFQs never generate book
            # traffic.
            od_pending = False
            od_timed_out = False
            for gl, od_gid in zip(by_game.values(), game_ids_list):
                if _ENGINE is None:
                    od_pending = True          # engine absent -> fail-safe skip
                    continue
                od_hash = legset.leg_set_hash(gl)
                if _ENGINE.lookup(od_hash) is not None:
                    _maybe_emit_on_demand_result(rid, ticker, od_hash, gl)
                    continue
                if _ENGINE.landed_empty(od_hash):
                    od_timed_out = True
                    continue
                gref = _game_ref(od_gid)
                if gref is not None and _ENGINE.ensure_fetch(od_hash, gref, gl):
                    research.emit("on_demand_requested", rfq_id=rid, ticker=ticker,
                                  payload=dict(leg_set_hash=od_hash, game_id=od_gid,
                                               n_legs=len(gl)))
                od_pending = True
            if od_pending:
                _decide("skipped", rfq_id=rid, ticker=ticker, game_id=game_id,
                              reason="on_demand_pending")
                continue
            if od_timed_out:
                # Only claimed once nothing is still pending: every sub-combo
                # either priced or landed empty — "we asked live and the books
                # didn't answer" is a decline, not a wait.
                _decide("skipped", rfq_id=rid, ticker=ticker, game_id=game_id,
                              reason="live_fetch_timeout")
                continue
            # Fair value: router prices via legset/grid across all games; fresh
            # on-demand fairs are injected via the engine's pure lookup. The
            # detail carries sigma_pts (cross-book dispersion) and n_games for
            # the uncertainty-scaled margin (#19); gate declines get their own
            # skip reasons (#20) so the report can tell "not enough books" from
            # "books disagree" — everything else stays "no_fair".
            od_lookup = _ENGINE.lookup if _ENGINE is not None else None
            fair_detail, gate_reason = router.combo_fair_detail(
                legs, None, _resolve_game_for_legs,
                config.MIN_AGREEING_BOOKS, config.SIGMA_Z_MAX,
                on_demand_fairs=od_lookup, live_routing=True)
            blended = fair_detail.fair if fair_detail is not None else None
            book_med = blended  # single consensus fair; book_med == blended
            if blended is None or not (config.MIN_FAIR_PROB <= blended <= config.MAX_FAIR_PROB):
                skip_reason = (gate_reason if fair_detail is None and gate_reason
                               in ("too_few_books", "consensus_dispersion")
                               else "no_fair")
                if skip_reason == "too_few_books":
                    # "The live fetch answered too thin" — named so the monitor
                    # separates a thin live slate from the pre-#54 cache era.
                    skip_reason = "live_too_few_books"
                # #55: a dispersion bust doesn't just block a NEW quote — it
                # pulls the RESTING one this RFQ may already have (a straggler
                # book just contradicted the thinner quorum that priced it).
                # too_few_books deliberately does NOT pull: a set that merely
                # thinned was not contradicted, and the risk sweep's drift and
                # constituent-jump breakers still cover the resting quote.
                if (skip_reason == "consensus_dispersion"
                        and _pull_quorum_quote(gateway, rid, ticker, game_id,
                                               snapshot)):
                    continue
                _decide("skipped", rfq_id=rid, ticker=ticker, game_id=game_id,
                              reason=skip_reason, model=None, book=book_med,
                              blended=blended)
                continue
            # circuit breaker bookkeeping
            prev = _PREV_BOOK_FAIR.get(ticker)
            _PREV_BOOK_FAIR[ticker] = book_med
            if (prev is not None and book_med is not None
                    and risk.book_move_triggered(prev, book_med, config.BOOK_MOVE_CB_THRESHOLD)):
                # H1: circuit breaker now actually CANCELS open quotes on this combo
                # (previously it only blocked re-quotes — resting quotes were exposed).
                opens = snapshot.open_quote_ids_for_ticker(ticker)
                for open_qid in opens:
                    try:
                        gateway.cancel(open_qid)
                    except Exception:
                        pass
                    with db.connect() as con:
                        con.execute(
                            "UPDATE live_quotes SET status='cancelled', closed_at=? WHERE quote_id=?",
                            [datetime.now(timezone.utc), open_qid])
                    snapshot.remove_quote(open_qid)
                _decide("circuit_breaker", rfq_id=rid, ticker=ticker, game_id=game_id,
                              reason=f"book_move_pulled_{len(opens)}_quotes")
                # Event 9: circuit_breaker
                research.emit("circuit_breaker", rfq_id=rid, ticker=ticker,
                              payload=dict(prev_book_med=prev, cur_book_med=book_med,
                                           threshold=config.BOOK_MOVE_CB_THRESHOLD,
                                           opens_cancelled=len(opens),
                                           game_id=game_id))
                continue
            # P-7 same-price block (issue #21): books have refreshed since the
            # last fill on this combo (refresh gate above), but if consensus fair
            # hasn't actually moved we would re-post the exact price that just got
            # picked off. Placed after the circuit breaker so a big move still
            # cancels resting quotes first.
            if last_fill_fair is not None and risk.same_price_blocked(
                    blended, last_fill_fair, config.QUOTE_HYSTERESIS):
                _decide("skipped", rfq_id=rid, ticker=ticker, game_id=game_id,
                              reason="same_price_block", book=book_med, blended=blended)
                continue
            # #19: uncertainty-scaled margin — cross-book dispersion widens the
            # cushion, game count compounds the ROI divisor. #55: a quote whose
            # thinnest per-game consensus is exactly 2 books rides a 2-sample
            # sigma, so it carries QUORUM_MARGIN_ADDON as one extra floor term;
            # the add-on drops out on its own when a straggler thickens the set
            # and the quote refines.
            quorum_addon_pts = (config.QUORUM_MARGIN_ADDON
                                if fair_detail.min_n_books == 2 else 0.0)
            q = pricing.quote(blended, config.TARGET_ROI,
                              sigma_pts=fair_detail.sigma_pts,
                              n_games=fair_detail.n_games,
                              min_margin_pts=config.MIN_MARGIN_PTS,
                              k_sigma=config.K_SIGMA,
                              quorum_addon_pts=quorum_addon_pts)
            if q is None:
                _decide("skipped", rfq_id=rid, ticker=ticker, game_id=game_id,
                              reason="unpriceable")
                continue
            # Per-fill dollar cap (replaces the old contract cap). Worst-case over
            # both sides the creator could take; our held-side cost basis = max loss.
            fill_exposure = _worst_fill_exposure_usd(rfq, q)
            if fill_exposure > config.max_fill_exposure_usd():
                _decide("skipped", rfq_id=rid, ticker=ticker, game_id=game_id,
                              reason="size_gate")
                continue
            # H8 + N7 + N8: per-combo exposure cap — block multi-RFQ concentration
            # on one combo. N7: also count outstanding open live_quotes' worst-case
            # exposure so a burst of RFQs on the same combo can't all be accepted
            # before any fills register. N8: treat unreconciled fills conservatively
            # (counted at the per-fill cap) for the same reason.
            combo_exp = snapshot.combo_fill_exposure(ticker)
            inflight_count = snapshot.open_count_for_ticker(ticker)
            worst_inflight = float(inflight_count) * config.max_fill_exposure_usd()
            if combo_exp + worst_inflight + fill_exposure > config.MAX_COMBO_EXPOSURE_USD:
                _decide("skipped", rfq_id=rid, ticker=ticker, game_id=game_id,
                              reason="per_combo_cap")
                continue
            # FIX I1 + M1: dedup + hysteresis — skip re-quote if price unchanged,
            # else mark the old row 'replaced' so we never leave ghost-open rows.
            # #55: computed BEFORE the quote_priced emit so the refine action
            # (initial / hold / replace) rides in the research payload.
            existing = snapshot.existing_open_quote(rid)
            replaced_qid = None
            refine_action = "initial"
            if existing:
                eqid, eyb, enb = existing
                if (abs(q.yes_bid - eyb) < config.QUOTE_HYSTERESIS
                        and abs(q.no_bid - enb) < config.QUOTE_HYSTERESIS):
                    refine_action = "hold"  # price essentially unchanged
                else:
                    # price moved beyond hysteresis → replace. B2 fix (issue #16):
                    # submit FIRST, mark the old row 'replaced' only on success.
                    # Kalshi auto-cancels the prior quote only when a new one
                    # lands, so marking 'replaced' before a submit that then fails
                    # would leave a live exchange quote invisible to the confirm
                    # loop and risk sweep. On failure the old row stays 'open'
                    # (still tracked); the next tick retries the replace.
                    refine_action = "replace"
                    replaced_qid = eqid
            # Event 2: quote_priced — we have a valid quote from the pricer.
            # #19 item 5: margin components ride in the firehose payload (NOT in
            # quote_decisions — the #14 report's column semantics stay intact).
            # #55: quorum_size + quorum_addon_pts + refine_action make quorum
            # quoting and every refinement measurable from research data.
            research.emit("quote_priced", rfq_id=rid, ticker=ticker,
                          payload=dict(leg_set_hash=legset.leg_set_hash(canon),
                                       n_games=len(by_game),
                                       blended_fair=blended, yes_bid=q.yes_bid, no_bid=q.no_bid,
                                       game_id=game_id,
                                       sigma_pts=q.sigma_pts, floor_pts=q.floor_pts,
                                       quorum_size=fair_detail.min_n_books,
                                       quorum_addon_pts=q.quorum_addon_pts,
                                       refine_action=refine_action,
                                       resting_quote_id=existing[0] if existing else None,
                                       roi_pts_yes=q.roi_pts_yes, roi_pts_no=q.roi_pts_no,
                                       margin_pts_yes=q.margin_pts_yes,
                                       margin_pts_no=q.margin_pts_no,
                                       live_games=_live_games_detail(by_game)))
            if dry_run:
                _decide("dry_run_quote", rfq_id=rid, ticker=ticker, game_id=game_id,
                              model=None, book=book_med, blended=blended, yb=q.yes_bid, nb=q.no_bid)
                continue
            if refine_action == "hold":
                continue  # leave the resting quote in place
            # #17 singles baseline: snapshot the raw Kalshi odds of every leg at
            # quote time — the confirm-window veto compares a fresh read against
            # this and voids on any tick. No baseline ⇒ don't quote: a quote
            # without one is doomed to void on accept, and chronic non-confirms
            # are abusive behavior Kalshi can throttle.
            leg_snapshot = _leg_market_prices(legs)
            if leg_snapshot is None:
                _decide("skipped", rfq_id=rid, ticker=ticker, game_id=game_id,
                              reason="no_leg_snapshot")
                continue
            # #23 item 3: correlation sanity against Kalshi's live singles. The leg
            # snapshot we just fetched for #17's veto IS the marginal anchor, so
            # this costs ZERO extra API calls. Independent of the #20 gate: that
            # one asks whether the books agree with EACH OTHER, this asks whether
            # their consensus is consistent with the real-time single-leg prices —
            # tightly-agreeing books can still be jointly wrong. Degenerate books
            # (yes_ask=100, empty, crossed) yield no marginals: we log the miss and
            # quote on book consensus alone, exactly as before this ticket, rather
            # than declining on missing information.
            marginals = singles.marginals_for_legs(leg_snapshot, legs)
            sanity = singles.corr_sanity(blended, marginals,
                                         config.CORR_PREMIUM_MIN,
                                         config.CORR_PREMIUM_MAX) if marginals else None
            # Item 5: premium + marginals per quote — the tuning dataset for the
            # band and the calibration dataset for a Phase-3 correlation model.
            research.emit("corr_sanity_check", rfq_id=rid, ticker=ticker,
                          payload=dict(game_id=game_id, combo_fair=blended,
                                       n_legs=len(legs), marginals=marginals,
                                       baseline_independent=(
                                           sanity.baseline_independent if sanity else None),
                                       premium=sanity.premium if sanity else None,
                                       frechet_lo=sanity.frechet_lo if sanity else None,
                                       frechet_hi=sanity.frechet_hi if sanity else None,
                                       reason=sanity.reason if sanity else None))
            gated = sanity is not None and (
                (sanity.reason == "frechet" and config.CORR_SANITY_FRECHET_ENABLED)
                or (sanity.reason == "premium" and config.CORR_SANITY_PREMIUM_ENABLED))
            if gated:
                _decide("skipped", rfq_id=rid, ticker=ticker, game_id=game_id,
                              reason=f"corr_sanity_{sanity.reason}",
                              book=book_med, blended=blended)
                continue
            qid = gateway.submit_quote(rid, q.yes_bid, q.no_bid)
            if qid:
                with db.connect() as con:
                    # One explicit transaction: DuckDB autocommits per statement,
                    # and a crash between marking the old row 'replaced' and
                    # inserting the new one would leave the freshly-submitted
                    # exchange quote untracked.
                    con.execute("BEGIN TRANSACTION")
                    if replaced_qid is not None:
                        con.execute(
                            "UPDATE live_quotes SET status='replaced', closed_at=? WHERE quote_id=?",
                            [datetime.now(timezone.utc), replaced_qid])
                    con.execute(
                        "INSERT INTO live_quotes (quote_id, rfq_id, combo_market_ticker, "
                        "game_id, yes_bid, no_bid, model_fair, book_fair, blended_fair, "
                        "status, submitted_at, closed_at, worst_exposure_usd, leg_prices_json) "
                        "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                        [qid, rid, ticker, game_id, q.yes_bid, q.no_bid, None, book_med,
                         blended, "open", datetime.now(timezone.utc), None, fill_exposure,
                         json.dumps(leg_snapshot)])
                    # R-1 (issue #22): game attribution for the OPEN-quote exposure
                    # gates — one row per game the combo touches, same rule as
                    # fill_games. Written in this transaction so an open quote can
                    # never exist without its game map.
                    for quote_game_id in sorted(set(game_ids_list)):
                        con.execute(
                            "INSERT OR REPLACE INTO quote_games (quote_id, game_id) "
                            "VALUES (?, ?)", [qid, quote_game_id])
                    con.execute(
                        "INSERT OR REPLACE INTO seen_rfqs "
                        "(rfq_id, market_ticker, in_scope, game_id, legs_json, "
                        "first_seen_at, last_decision, creator_id) "
                        "VALUES (?,?,?,?,?,?,?,?)",
                        [rid, ticker, True, game_id, json.dumps(legs),
                         datetime.now(timezone.utc), "quoted", creator_id])
                    con.execute("COMMIT")
                _decide("quoted", rfq_id=rid, quote_id=qid, ticker=ticker, game_id=game_id,
                              model=None, book=book_med, blended=blended, yb=q.yes_bid, nb=q.no_bid)
                open_count += 1
                # Fold the new quote into the pass snapshot so the cap gates
                # see it for the REST of this pass (N7/R-1 read-your-writes).
                if replaced_qid is not None:
                    snapshot.remove_quote(replaced_qid)
                snapshot.record_submitted(qid, rid, ticker, game_id,
                                          sorted(set(game_ids_list)),
                                          fill_exposure, q.yes_bid, q.no_bid)
    finally:
        # One write connection — one checkpoint — for the whole pass, and it
        # runs even when the pass breaks early (SIGTERM, MAX_OPEN_QUOTES).
        if decision_rows or seen_insert_rows:
            with db.connect() as con:
                if decision_rows:
                    con.executemany(_DECISION_INSERT_SQL, decision_rows)
                if seen_insert_rows:
                    con.executemany(
                        "INSERT OR REPLACE INTO seen_rfqs "
                        "(rfq_id, market_ticker, in_scope, game_id, legs_json, "
                        "first_seen_at, last_decision, creator_id) "
                        "VALUES (?,?,?,?,?,?,?,?)", seen_insert_rows)
        # Advance the rotation so the next lap resumes past this lap's work.
        _DISCOVERY_CURSOR = (start_at + iterated) % len(rfqs) if rfqs else 0
    # pass_sec is the #89 health metric: the RFQ re-visit cadence IS the
    # pass duration whenever it exceeds DISCOVERY_SEC, so a slow pass is a
    # quote-latency regression even when nothing was deferred. Log it with
    # the usual counters, and unconditionally when it breaches the cadence.
    pass_sec = time.monotonic() - pass_started
    if (scope_fetches_deferred or timebox_deferred or stale_skipped
            or tipoff_skipped or horizon_skipped
            or pass_sec > config.DISCOVERY_SEC):
        log.info("discovery pass: iterated=%d/%d pass_sec=%.2f stale_skipped=%d "
                 "tipoff_skipped=%d horizon_skipped=%d "
                 "timebox_deferred=%d scope_fetch_deferred=%d",
                 iterated, len(rfqs), pass_sec, stale_skipped, tipoff_skipped,
                 horizon_skipped, timebox_deferred, scope_fetches_deferred)


def _confirm_tick(gateway, dry_run):
    with db.connect(read_only=True) as con:
        # FIX I2: also select model_fair and book_fair so we can carry them into fills
        live = con.execute(
            "SELECT quote_id, rfq_id, combo_market_ticker, game_id, yes_bid, no_bid, "
            "blended_fair, model_fair, book_fair, submitted_at, leg_prices_json "
            "FROM live_quotes WHERE status='open'").fetchall()
    for (qid, rid, ticker, game_id, yb, nb, prev_fair, model_fair_at_q,
         book_fair_at_q, submitted_at, leg_prices_json) in live:
        status, body, _ = auth_client.api("GET", f"/communications/quotes/{qid}")
        q = body.get("quote") if isinstance(body, dict) else None
        st = (q or {}).get("status")
        if st == "accepted":
            # Event 4: accept_observed — BEFORE side-inference logic, to capture
            # the raw Kalshi response shape for first-launch field verification.
            research.emit("accept_observed", rfq_id=rid, quote_id=qid, ticker=ticker,
                          payload=dict(response_keys=list((q or {}).keys()),
                                       response_raw=q))
            # Which side do we end up HOLDING? The taker accepted one side of our
            # two-sided quote; we hold the opposite. If they took our yes_bid (they
            # sold us YES), we hold YES; the API marks the accepted_side as theirs.
            accepted_side = (q or {}).get("accepted_side")
            side_held = "no" if accepted_side == "yes" else "yes"
            price = (1 - nb) if side_held == "yes" else (1 - yb)  # ask we transact at
            fee = maker_fee_per_contract(price)
            with db.connect(read_only=True) as con:
                legs_row = con.execute(
                    "SELECT legs_json FROM seen_rfqs WHERE rfq_id=?", [rid]).fetchone()
            legs = json.loads(legs_row[0]) if legs_row and legs_row[0] else []
            # H2: last-look HARD-FAIL when we can't re-price.
            # Principle: "can't re-price ⇒ don't confirm." Previously a missing
            # legs row or empty book_fairs silently fell back to prev_fair, which
            # made the drift check a no-op for exactly the cases that matter most.
            if not legs:
                _log_decision("voided_no_legs", rfq_id=rid, quote_id=qid, ticker=ticker,
                              game_id=game_id)
                with db.connect() as con:
                    con.execute(
                        "UPDATE live_quotes SET status='voided', closed_at=? WHERE quote_id=?",
                        [datetime.now(timezone.utc), qid])
                continue
            # #17 (reworked): singles-move veto — the fast last look. Re-read
            # the raw Kalshi odds of every leg (each leg is its own singles
            # market, so this is <~1s regardless of combo shape) and refuse
            # the fill if ANY leg's bid or ask moved even one tick since the
            # quote-time snapshot. Rationale: over the seconds between quote
            # and accept the legs' correlation is structural and static — the
            # marginals carry the pickoff signal, so unmoved raw leg odds ⇒
            # unmoved combo fair. Missing baseline, unparseable snapshot, or
            # a failed fresh read all void ("can't verify ⇒ don't confirm").
            try:
                snapshot = json.loads(leg_prices_json) if leg_prices_json else None
            except (TypeError, ValueError):
                snapshot = None
            fresh_leg_prices = _leg_market_prices(legs)
            singles_moved = (snapshot is None or fresh_leg_prices is None
                             or _singles_moved(snapshot, fresh_leg_prices))
            # Firehose: quote age at accept + snapshot-vs-fresh per-leg odds.
            # Measures how often the veto saves us AND (via the deltas) what a
            # non-zero tolerance would have passed — the tuning dataset.
            research.emit("confirm_singles_check", rfq_id=rid, quote_id=qid, ticker=ticker,
                          payload=dict(quote_age_sec=_quote_age_sec(submitted_at),
                                       snapshot=snapshot, fresh=fresh_leg_prices,
                                       moved=singles_moved, prev_fair=prev_fair))
            if singles_moved:
                _log_decision("voided_singles_moved", rfq_id=rid, quote_id=qid,
                              ticker=ticker, game_id=game_id)
                with db.connect() as con:
                    con.execute(
                        "UPDATE live_quotes SET status='voided', closed_at=? WHERE quote_id=?",
                        [datetime.now(timezone.utc), qid])
                continue
            # Phase 2 live last look: for combos with on-demand games, trigger
            # a synchronous priority re-fetch budgeted against the confirm
            # window (on-demand shapes have no cache fairs at all — without a
            # live fetch they cannot be re-priced for the EV check below).
            canon_c = legset.parse_legs(legs)
            refetch_jobs = []
            # refetch_ok is the engine's own "every job landed AFTER this
            # refetch entered" guarantee. It must gate the confirm: the result
            # store serves lookups for QUOTE_FRESH_SEC, so after a failed
            # re-fetch a PRE-entry result (e.g. from another accept on the
            # same combo seconds ago) could still satisfy the lookup — relying
            # on "failed fetch ⇒ stale lookup" alone would confirm on a
            # previously-fetched number.
            refetch_ok = True
            if _ENGINE is not None and canon_c:
                for gl in legset.partition_by_game(canon_c).values():
                    # #54: EVERY sub-combo re-fetches — the quote's fair came
                    # from the engine, whose result has aged out by accept
                    # time; the cache must not stand in.
                    gid_c = _resolve_game_for_legs(gl)
                    gref = _game_ref(gid_c) if gid_c else None
                    if gref is None:
                        refetch_jobs = None  # unresolvable -> stale -> void path
                        refetch_ok = False
                        break
                    refetch_jobs.append((legset.leg_set_hash(gl), gref, gl))
                if refetch_jobs:
                    refetch_ok = _ENGINE.refetch_now(refetch_jobs,
                                                     CONFIRM_REFETCH_BUDGET_SEC)
            od_lookup = _ENGINE.lookup if _ENGINE is not None else None
            cur_fair = router.combo_fair(legs, None, _resolve_game_for_legs,
                                         config.MIN_AGREEING_BOOKS, config.SIGMA_Z_MAX,
                                         on_demand_fairs=od_lookup,
                                         live_routing=True)
            if cur_fair is None or not refetch_ok:
                _log_decision("voided_no_fresh_books", rfq_id=rid, quote_id=qid, ticker=ticker,
                              game_id=game_id)
                with db.connect() as con:
                    con.execute(
                        "UPDATE live_quotes SET status='voided', closed_at=? WHERE quote_id=?",
                        [datetime.now(timezone.utc), qid])
                continue
            if risk.last_look_ok(side_held, price, fee, cur_fair, prev_fair,
                                 config.FAIR_DRIFT_TOLERANCE):
                if not dry_run and gateway.confirm(qid):
                    # B4 fix (issue #18): sizes arrive as contracts_fp strings;
                    # an unknown size records 0 (NOT 1) — the N8 convention
                    # already counts unreconciled fills at the per-fill cap,
                    # and the reconcile sweep sets the true size from Kalshi
                    # positions.
                    contracts = _accept_fill_contracts(q)
                    if contracts is None:
                        log.warning(
                            "[fill_size_unparseable] quote_id=%s ticker=%s "
                            "response_keys=%s — recording 0 contracts pending "
                            "reconcile", qid, ticker, list((q or {}).keys()))
                        contracts = 0.0
                    # N5: record the fill IMMEDIATELY with EXPECTED side/size
                    # (no positions retry in the hot path — that lived here
                    # before and could spend up to 15s/fill, blowing the 30s
                    # confirm window during accept bursts). reconciled=FALSE
                    # marks the row for the background sweep to verify against
                    # Kalshi positions and correct side/size if needed.
                    now_ts = datetime.now(timezone.utc)
                    fill_id = str(uuid.uuid4())
                    # Derive game attribution BEFORE the write: _fill_game_ids is
                    # fail-safe (never raises), so the fills + fill_games inserts
                    # below stay consistent — a fills row is never committed
                    # without its matching fill_games rows (which would let the
                    # per-game cap fail open / under-count).
                    fill_game_ids = _fill_game_ids(legs, game_id)
                    with db.connect() as con:
                        con.execute(
                            "INSERT INTO fills (fill_id, quote_id, rfq_id, combo_market_ticker, "
                            "game_id, side_held, contracts, price, fee, model_fair_at_quote, "
                            "book_fair_at_quote, blended_fair_at_quote, fair_at_confirm, "
                            "realized_pnl, filled_at, reconciled) "
                            "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                            # FIX I2: carry model_fair_at_q / book_fair_at_q from the live_quotes row
                            [fill_id, qid, rid, ticker, game_id, side_held,
                             contracts, price, fee, model_fair_at_q, book_fair_at_q,
                             prev_fair, cur_fair, None, now_ts, False])
                        # Per-game exposure ledger: one row per game the combo
                        # touches (cross-game combos land 2+). Powers the per-game
                        # cap without double-counting the daily cap / P&L.
                        for g in fill_game_ids:
                            con.execute(
                                "INSERT OR REPLACE INTO fill_games (fill_id, game_id) "
                                "VALUES (?, ?)", [fill_id, g])
                        con.execute(
                            "UPDATE live_quotes SET status='filled', closed_at=? WHERE quote_id=?",
                            [now_ts, qid])
                        # H9: arm the per-combo cooldown — block fresh quotes
                        # on the same combo for COMBO_COOLDOWN_SEC.
                        con.execute(
                            "INSERT OR REPLACE INTO combo_cooldown "
                            "(combo_market_ticker, cooled_until) VALUES (?, ?)",
                            [ticker, now_ts + timedelta(seconds=config.COMBO_COOLDOWN_SEC)])
                    # Event 5: fill_recorded — after INSERT, before notify.
                    research.emit("fill_recorded", rfq_id=rid, quote_id=qid, ticker=ticker,
                                  payload=dict(side_held_expected=side_held,
                                               contracts_expected=contracts,
                                               price=price, fee=fee,
                                               prev_fair=prev_fair, cur_fair=cur_fair,
                                               model_fair_at_q=model_fair_at_q,
                                               book_fair_at_q=book_fair_at_q,
                                               # Phase 2: per-book route + which
                                               # books agreed (DK+NV-only fills
                                               # are a named risk metric).
                                               on_demand=_on_demand_fill_info(canon_c)))
                    notify.fill(ticker, side_held, contracts, price)
                    _log_decision("confirmed", rfq_id=rid, quote_id=qid, ticker=ticker,
                                  game_id=game_id)
                    # #57: fire the targeted post-fill fetch NOW — the
                    # cooldown's awaiting-refresh gate clears on this
                    # flight's landing, usually right at the 60s floor,
                    # with no sweep pass involved.
                    if canon_c:
                        _ensure_post_fill_fetches(
                            rid, ticker, legset.partition_by_game(canon_c))
            else:
                _log_decision("voided_last_look", rfq_id=rid, quote_id=qid, ticker=ticker,
                              game_id=game_id)
                with db.connect() as con:
                    con.execute(
                        "UPDATE live_quotes SET status='voided', closed_at=? WHERE quote_id=?",
                        [datetime.now(timezone.utc), qid])
        elif st in ("cancelled", "expired", "executed"):
            with db.connect() as con:
                con.execute("UPDATE live_quotes SET status=?, closed_at=? WHERE quote_id=?",
                            [st, datetime.now(timezone.utc), qid])
        elif status == 404 and isinstance(body, dict):
            # Kalshi DELETES the quote object when its RFQ expires: the GET
            # returns 404 {'error': {'code': 'not_found'}} (live-verified
            # 2026-08-12). Without this branch st is None, nothing matches,
            # and the row wedges as 'open' forever — counting against
            # MAX_OPEN_QUOTES and both exposure caps until the bot silently
            # stops quoting. Only an explicit 404 with a parsed JSON body is
            # terminal; transient failures (5xx, timeout, non-dict body,
            # connection errors) fall through and leave the quote open.
            log.info("[quote_gone] quote_id=%s ticker=%s — 404 on quote GET; "
                     "RFQ expired and Kalshi deleted the quote; closing as "
                     "expired", qid, ticker)
            with db.connect() as con:
                con.execute(
                    "UPDATE live_quotes SET status='expired', closed_at=? WHERE quote_id=?",
                    [datetime.now(timezone.utc), qid])


def _reconcile_sweep_tick():
    """N5: background reconciliation sweep. Picks up fills written with
    reconciled=FALSE and verifies side/size against Kalshi `/portfolio/positions`.

    For each unreconciled fill on ticker T:
      - prior_total = signed sum of OTHER fills on T (excluding this one)
      - actual      = current Kalshi position on T (signed)
      - delta       = actual - prior_total  → THIS fill's true side/size
    If delta disagrees with what we recorded, UPDATE the fill row. Either way,
    set reconciled=TRUE so we don't re-check it.

    N10: phantom-fill detection — if delta==0 but expected was non-zero, the
    confirm succeeded but no position materialized (e.g. cancel race, internal
    Kalshi mismatch). Mark contracts=0 + reconciled=TRUE + log [phantom_fill].
    Caps are no longer inflated by confirms that didn't actually fill.

    N11: max-age fallback — if positions API is persistently down, a fill older
    than MAX_RECONCILE_AGE_SEC gets marked reconciled=TRUE with its recorded
    values + [reconcile_max_age_fallback] warning. Prevents indefinite
    cap-blocking on Kalshi outages.

    Lifted out of the confirm hot path (was up to 15s/fill via retries on
    eventual consistency). Sweep runs every RECONCILE_SWEEP_SEC.
    """
    now = datetime.now(timezone.utc)
    with db.connect(read_only=True) as con:
        rows = con.execute(
            "SELECT fill_id, combo_market_ticker, side_held, contracts, filled_at "
            "FROM fills WHERE reconciled = FALSE").fetchall()
    if not rows:
        return
    for fill_id, ticker, recorded_side, recorded_ct, filled_at in rows:
        # N11: compute fill age before the API call so we can apply the max-age
        # fallback if the positions API is persistently unreachable.
        # O-2: filled_at is TIMESTAMPTZ (aware) post-migration — compare aware
        # UTC. Naive branch kept only for rows read before a restart applies
        # the migration; those are naive-LOCAL, so compare naive-local.
        if filled_at is not None:
            if filled_at.tzinfo is not None:
                fill_age = max(0.0, (now - filled_at).total_seconds())
            else:
                fill_age = max(0.0, (datetime.now() - filled_at).total_seconds())
        else:
            fill_age = 0.0
        actual = _get_position_contracts(ticker)
        if actual is None:
            if fill_age > config.MAX_RECONCILE_AGE_SEC:
                # N11: API persistently down — fall back to recorded values so
                # caps don't block new quotes indefinitely during an outage.
                with db.connect() as con:
                    con.execute("UPDATE fills SET reconciled=TRUE WHERE fill_id=?",
                                [fill_id])
                log.warning("[reconcile_max_age_fallback] fill_id=%s ticker=%s "
                            "age=%.0fs — positions API unreachable, marking "
                            "reconciled with recorded values", fill_id, ticker, fill_age)
                research.emit("reconcile_done", ticker=ticker,
                              payload=dict(actual_total=None, prior_total=None,
                                           delta=None, recorded_signed=None,
                                           outcome="max_age_fallback"))
            else:
                # API down — leave reconciled=FALSE; next sweep will retry.
                log.warning("[position_reconcile_unavailable] ticker=%s fill_id=%s "
                            "— will retry next sweep", ticker, fill_id)
            continue
        # prior_total = signed sum of all OTHER fills on this ticker.
        with db.connect(read_only=True) as con:
            prior_total = con.execute(
                "SELECT COALESCE(SUM(CASE WHEN side_held='yes' THEN contracts "
                "ELSE -contracts END), 0) FROM fills "
                "WHERE combo_market_ticker=? AND fill_id<>?",
                [ticker, fill_id]).fetchone()[0]
        delta = int(actual) - int(prior_total)
        recorded_signed = (int(recorded_ct) if recorded_side == "yes"
                           else -int(recorded_ct))
        # N10: phantom-fill — confirm succeeded but no position materialized.
        if delta == 0 and recorded_signed != 0:
            with db.connect() as con:
                con.execute(
                    "UPDATE fills SET contracts=0, reconciled=TRUE WHERE fill_id=?",
                    [fill_id])
            log.warning("[phantom_fill] fill_id=%s ticker=%s "
                        "expected=%d actual_delta=0 (kalshi_total=%d "
                        "prior=%d) — marking contracts=0",
                        fill_id, ticker, recorded_signed, actual, prior_total)
            research.emit("reconcile_done", ticker=ticker,
                          payload=dict(actual_total=int(actual),
                                       prior_total=int(prior_total),
                                       delta=delta, recorded_signed=recorded_signed,
                                       outcome="phantom"))
            continue
        if delta != recorded_signed:
            log.warning("[position_mismatch] ticker=%s fill_id=%s "
                        "recorded=%d actual_delta=%d "
                        "(kalshi_total=%d prior=%d) "
                        "— trusting Kalshi",
                        ticker, fill_id, recorded_signed, delta, actual, prior_total)
            if delta > 0:
                new_side, new_ct = "yes", delta
            elif delta < 0:
                new_side, new_ct = "no", abs(delta)
            else:
                # delta == 0 AND recorded_signed == 0: already correct, just mark done.
                new_side, new_ct = recorded_side, 0
            with db.connect() as con:
                con.execute(
                    "UPDATE fills SET side_held=?, contracts=?, reconciled=TRUE "
                    "WHERE fill_id=?",
                    [new_side, new_ct, fill_id])
            research.emit("reconcile_done", ticker=ticker,
                          payload=dict(actual_total=int(actual),
                                       prior_total=int(prior_total),
                                       delta=delta, recorded_signed=recorded_signed,
                                       outcome="mismatch_corrected"))
        else:
            with db.connect() as con:
                con.execute("UPDATE fills SET reconciled=TRUE WHERE fill_id=?",
                            [fill_id])
            research.emit("reconcile_done", ticker=ticker,
                          payload=dict(actual_total=int(actual),
                                       prior_total=int(prior_total),
                                       delta=delta, recorded_signed=recorded_signed,
                                       outcome="matched"))


def _quote_game_ids(legs_json: str | None) -> list[str] | None:
    """All game_ids a resting quote's combo touches, or None if ANY is
    unresolvable (missing/unparseable legs, unresolved game). B1 fix (issue
    #15): live_quotes.game_id holds only the PRIMARY game, so the sweep
    re-derives the full set from legs_json via the same legset path discovery
    uses. Callers must treat None as fail-safe → cancel, never fail open."""
    if not legs_json:
        return None
    try:
        canon = legset.parse_legs(json.loads(legs_json))
        if not canon:
            return None
        gids = [_resolve_game_for_legs(gl)
                for gl in legset.partition_by_game(canon).values()]
        if any(gid is None for gid in gids):
            return None
        return gids
    except Exception:
        return None


def _current_consensus_fair(legs_json: str | None) -> float | None:
    """Current book-consensus fair for a resting quote's combo, or None.

    Extracted so the drift check and the #23 constituent-jump breaker's
    "were the books quiet?" test read the SAME number. Live-routed (#57):
    a resting quote's RFQ is still in the poll, so the engine holds fairs at
    most ~QUOTE_FRESH_SEC old — sweep rows would be blind between the slow
    passes. Fail-safe: any parse or pricing error is None (no signal), never
    an exception into the sweep.
    """
    if not legs_json:
        return None
    try:
        legs = json.loads(legs_json)
        return router.combo_fair(legs, None, _resolve_game_for_legs,
                                 config.MIN_AGREEING_BOOKS,
                                 config.SIGMA_Z_MAX,
                                 on_demand_fairs=(_ENGINE.lookup
                                                  if _ENGINE is not None else None),
                                 live_routing=True)
    except Exception:
        return None


def _open_quote_constituent_tickers(live_rows) -> set:
    """Every DISTINCT constituent ticker across all open quotes.

    This set IS the sweep's API cost: one GET per element, per sweep. It is
    bounded by (open quotes × legs) and is usually far smaller, because quotes
    on the same game share legs. Zero open quotes → zero calls.
    """
    tickers = set()
    for row in live_rows:
        leg_prices_json = row[6]
        if not leg_prices_json:
            continue
        try:
            tickers.update(str(t) for t in json.loads(leg_prices_json))
        except Exception:
            continue
    return tickers


def _rotated_poll_order(tickers: set) -> list:
    """Constituent tickers in a deterministic order, rotated each sweep.

    The poll is wall-clock budgeted (see singles.fetch_market_prices), so on a
    large resting book some tickers go unpolled. Rotating the start point means
    the SAME tail is not starved every sweep — over a few sweeps every ticker
    gets looked at. Sorted first so the order is reproducible in tests.
    """
    global _CONSTITUENT_POLL_CURSOR
    ordered = sorted(tickers)
    if not ordered:
        return ordered
    start = _CONSTITUENT_POLL_CURSOR % len(ordered)
    return ordered[start:] + ordered[:start]


def _games_for_tickers(raw_legs: list, tickers: set) -> set:
    """game_ids of the legs whose market_ticker is in `tickers`.

    CanonicalLeg carries the game code but not the market_ticker, so the raw
    leg dicts map ticker → event_ticker, `legset.game_id_of` reduces that to
    the same family-independent code `partition_by_game` keys on (#71), and
    that yields the CanonicalLegs `_resolve_game_for_legs` consumes.
    """
    out = set()
    try:
        canon = legset.parse_legs(raw_legs)
        if not canon:
            return out
        by_game = legset.partition_by_game(canon)
        jumped_games = {legset.game_id_of(str(leg.get("event_ticker") or ""))
                        for leg in raw_legs
                        if str(leg.get("market_ticker") or "") in tickers}
        for game_code in jumped_games:
            game_legs = by_game.get(game_code)
            gid = _resolve_game_for_legs(game_legs) if game_legs else None
            if gid:
                out.add(gid)
    except Exception:
        return out
    return out


def _constituent_jumped_games(live_rows, current_prices) -> tuple[set, dict]:
    """(game_ids whose constituents jumped since quote time, per-quote detail).

    #23 items 1-2. Our books refresh every ~150-165s; the constituent Kalshi
    singles trade in real time, so a constituent that moved since we quoted is
    the fastest evidence that the market left our resting price behind.

    The placement baseline is `live_quotes.leg_prices_json` — the raw Kalshi
    odds #17 already snapshots at submission — so this needs no schema of its
    own. The jump is unconditional (#54): every quote was live-priced at
    placement, so any constituent move since then means the market moved
    after us.

    Fail-safe: any row we cannot evaluate contributes NO jump signal. It is
    still covered by the tipoff / staleness / book-drift gates, and #17's
    confirm veto still fails closed on any accept.
    """
    jumped_games, detail = set(), {}
    if not current_prices:
        return jumped_games, detail
    for (qid, _game_id, _ticker, book_fair_at_q, _rid,
         legs_json, leg_prices_json) in live_rows:
        if not legs_json or not leg_prices_json:
            continue
        try:
            baseline = json.loads(leg_prices_json)
            raw_legs = json.loads(legs_json)
        except Exception:
            continue
        moves = singles.jumped_tickers(baseline, current_prices,
                                       config.CONSTITUENT_JUMP_THRESHOLD)
        if not moves:
            continue
        detail[qid] = moves
        jumped_games.update(_games_for_tickers(raw_legs, set(moves)))
    return jumped_games, detail


def _risk_sweep_tick(gateway, *, book_health=None):
    if config.KILL_FILE.exists():
        notify.halt("kill_switch")
        return
    # "We can no longer trust our fair" auto-pull, re-keyed by #57: the old
    # rule read the sweep view's age, which under the slow sweep was empty
    # most of every cycle and would cancel the whole book constantly. Live-fetch
    # health (#38 failure streaks via BookHealthAlerter.consensus_dark) is
    # cadence-invariant: fewer than MIN_AGREEING_BOOKS books answering
    # on-demand fetches means resting quotes can no longer be re-verified,
    # so stop resting risk. Fail-safe: no alerter wired (tests, alerts
    # disabled) → no pull; the #23 breakers below still protect the book.
    books_unhealthy = book_health is not None and book_health.consensus_dark()
    with db.connect(read_only=True) as con:
        # Pull the per-quote fields we need for the drift-since-quote check too.
        live = con.execute(
            "SELECT lq.quote_id, lq.game_id, lq.combo_market_ticker, lq.book_fair, "
            "       lq.rfq_id, sr.legs_json, lq.leg_prices_json "
            "FROM live_quotes lq LEFT JOIN seen_rfqs sr ON lq.rfq_id = sr.rfq_id "
            "WHERE lq.status='open'").fetchall()
    # #23: real-time constituent poll. ONE GET per DISTINCT ticker across all
    # open quotes — skipped entirely when flat, and when the books are dark
    # (every quote is being cancelled anyway, so the calls would be wasted).
    jumped_games, jump_detail = set(), {}
    if live and not books_unhealthy:
        global _CONSTITUENT_POLL_CURSOR
        poll_order = _rotated_poll_order(_open_quote_constituent_tickers(live))
        current_prices = singles.fetch_market_prices(
            poll_order, budget_sec=config.CONSTITUENT_POLL_BUDGET_SEC)
        # Advance past what we actually got, so the next sweep starts where
        # this one ran out of budget (approximate: failed reads are not
        # counted, which just re-polls them sooner).
        _CONSTITUENT_POLL_CURSOR += max(1, len(current_prices))
        jumped_games, jump_detail = _constituent_jumped_games(live, current_prices)
        if jumped_games:
            research.emit("constituent_jump",
                          payload=dict(games=sorted(jumped_games),
                                       threshold=config.CONSTITUENT_JUMP_THRESHOLD,
                                       moves=jump_detail))
    for qid, game_id, ticker, book_fair_at_q, rid, legs_json, leg_prices_json in live:
        cancel = False
        cancel_reason = None
        if books_unhealthy:
            cancel, cancel_reason = True, "books_unhealthy"
        else:
            # B1 fix: tipoff-check the EARLIEST commence time across ALL the
            # combo's games, not just the primary game_id column. Any game we
            # can't resolve or clock we can't read → cancel (fail-safe).
            gids = _quote_game_ids(legs_json)
            if gids is None:
                cancel, cancel_reason = True, "unresolvable_game"
            else:
                commence_times = [_commence_time(gid) for gid in gids]
                if any(ct is None for ct in commence_times):
                    cancel, cancel_reason = True, "unresolvable_game"
                else:
                    earliest_ct = min(commence_times)
                    if not risk.tipoff_ok(earliest_ct, config.TIPOFF_CANCEL_MIN):
                        cancel, cancel_reason = True, "tipoff"
        if not cancel and jumped_games:
            # #23 item 2: a constituent single moved past the threshold since
            # this quote was placed. Evaluated BEFORE the book-drift check
            # because it is the more specific and much faster signal, so it
            # wins the logged reason when both would fire. Game-level: a quote
            # is pulled when ANY game it touches is compromised, even if its
            # own legs are not the ones that moved.
            gids = _quote_game_ids(legs_json)
            if gids and (set(gids) & jumped_games):
                cancel, cancel_reason = True, "constituent_jump"
        if not cancel:
            # H1 (part 2): per-open-quote drift-since-quote sweep.
            # Recompute current book consensus for this combo; if it has drifted
            # more than BOOK_MOVE_CB_THRESHOLD from book_fair-at-quote, cancel.
            # Catches gradual drift that the per-tick (last vs current) circuit
            # breaker misses (e.g., 1¢ moves over several ticks adding to 4¢).
            if book_fair_at_q is not None and legs_json:
                cur_med = _current_consensus_fair(legs_json)
                try:
                    if cur_med is not None and risk.book_move_triggered(
                            book_fair_at_q, cur_med, config.BOOK_MOVE_CB_THRESHOLD):
                        cancel, cancel_reason = True, "book_drift"
                except Exception:
                    pass
        if cancel:
            try:
                gateway.cancel(qid)
            except Exception:
                pass
            with db.connect() as con:
                con.execute(
                    "UPDATE live_quotes SET status='cancelled', closed_at=? WHERE quote_id=?",
                    [datetime.now(timezone.utc), qid])
            _log_decision("sweep_cancel", rfq_id=rid, quote_id=qid, ticker=ticker,
                          game_id=game_id, reason=cancel_reason)


def build_rfq_source():
    """WS discovery source (issue #56, WS-only by user decision 2026-08-07 —
    no mode switch; rollback is a git revert).

    The WS source owns its loud automatic REST fallback (watchdog trip →
    ERROR + notify + rest serving; recovery flips back), so the tick loop
    never needs to know how discovery is feeling. Must run after
    _configure_auth(): the WS handshake signs with the same injected
    credentials.
    """
    from kalshi_mlb_mm import ws_rfq_source
    source = ws_rfq_source.WebSocketRFQSource(rest_source=RestRFQSource())
    source.start()
    log.info("WS RFQ discovery started (REST is automatic fallback only)")
    return source


def main_loop(dry_run: bool):
    from kalshi_mlb_mm.log_setup import setup_logging
    setup_logging()
    research.init_research_db()
    _configure_auth()
    db.init_database()
    sid = db.start_session(pid=os.getpid(), dry_run=dry_run)
    research.set_session_id(str(sid))
    log.info("=== MM bot session %s dry_run=%s ===", sid, dry_run)
    source, gateway = build_rfq_source(), RestQuoteGateway()
    from kalshi_common.book_health import build_alerter
    from kalshi_common.sgp_service import SGPService
    # health_db_path (issue #38): per-book fetch health lands in the SAME
    # market DB (and write lock) the SGP rows go to. Buffered — see
    # flush_health() in the tick loop below.
    # book_health (issue #37): run-time alerting on consecutive failed
    # fetches. In-memory streaks fed by the same two call sites that write
    # health rows — never a query, which on duckdb 1.4.4 could cost a quote
    # (see kalshi_common/book_health.py). Rule B reuses MIN_AGREEING_BOOKS so
    # the alert fires at exactly the count the quoting gate cares about.
    # The alerter is shared: the service feeds it per-fetch outcomes, and the
    # risk sweep reads consensus_dark() as its "stop resting risk" signal
    # (#57). build_alerter returns None when alerting is disabled — the
    # sweep's book_health=None default then simply skips the pull.
    book_alerter = build_alerter(
        label="MM",
        enabled=config.BOOK_ALERT_ENABLED,
        streak_threshold=config.BOOK_ALERT_STREAK,
        min_healthy_books=config.MIN_AGREEING_BOOKS,
        paths=config.BOOK_ALERT_PATHS)
    # #81: no per_book_deadline_sec — that knob budgets service.refresh(),
    # the full-slate sweep the maker no longer runs. On-demand flights are
    # budgeted by on_demand_deadline_sec; warming gets its own explicit
    # wall budget at the warm_cycle call (STRUCTURE_WARM_BUDGET_SEC).
    sgp_service = SGPService(on_demand_deadline_sec=config.ON_DEMAND_DEADLINE_SEC,
                             health_db_path=str(config.MARKET_DB),
                             book_health=book_alerter)
    # Phase 2: on-demand pricing engine shares the service's persistent
    # clients + structure caches. Always on — no switch (user decision);
    # the bot-wide kill file remains the emergency stop.
    global _ENGINE
    from kalshi_mlb_mm.on_demand import OnDemandEngine
    _ENGINE = OnDemandEngine(sgp_service)
    # Synchronous warm-up: populate mlb_target_lines before the loop starts
    # (game resolution, tipoff gating and the first warming pass all read
    # it). #81: ZERO book requests here — the maker's book traffic is now
    # purely on-demand flights + #50's structure warming. both_teams=True
    # (issue #70) enumerates BOTH teams' margin lines (signed spread_line:
    # home margin negative, away margin positive), matching the taker.
    warmup_targets_ok = True
    try:
        _target_line_tick()
    except Exception as e:
        warmup_targets_ok = False
        log.warning("warmup target-line refresh failed: %s", e)
    # "warm" starts at 0.0 so the first loop iteration warms immediately:
    # every book's on-demand structure caches (and the CZR WAF token) are
    # cold until the first warming pass (issue #50). A FAILED warmup
    # target-line refresh retries on the first loop iteration too (then
    # falls back to the 300s cadence) — without it, a Kalshi blip at boot
    # leaves game resolution blind for a full period.
    last = {"disc": 0.0, "conf": 0.0, "risk": 0.0, "reconcile": 0.0,
            "settle": 0.0, "warm": 0.0, "coverage": time.time(),
            "targets": time.time() if warmup_targets_ok else 0.0}
    try:
        while _running.is_set():
            now = time.time()
            if now - last["disc"] >= config.DISCOVERY_SEC:
                try:
                    _discovery_tick(source, gateway, dry_run)
                except Exception as e:
                    log.error("disc err: %s", e)
                last["disc"] = now
            if now - last["conf"] >= config.CONFIRM_SEC:
                try:
                    _confirm_tick(gateway, dry_run)
                except Exception as e:
                    log.error("conf err: %s", e)
                last["conf"] = now
            if now - last["risk"] >= config.RISK_SWEEP_SEC:
                try:
                    _risk_sweep_tick(gateway, book_health=book_alerter)
                except Exception as e:
                    log.error("risk err: %s", e)
                last["risk"] = now
            if (not dry_run
                    and now - last["reconcile"] >= config.RECONCILE_SWEEP_SEC):
                try:
                    _reconcile_sweep_tick()
                except Exception as e:
                    log.error("reconcile err: %s", e)
                last["reconcile"] = now
            # Settlement sweep (issue #12): populate realized_pnl once markets
            # settle. Own slow cadence — only matters hours post-game. Gated
            # off dry-run like reconcile (both read /portfolio-side truth that
            # only exists for live sessions).
            if (not dry_run
                    and now - last["settle"] >= config.SETTLEMENT_SWEEP_SEC):
                try:
                    settlement.settlement_sweep_tick()
                except Exception as e:
                    log.error("settle err: %s", e)
                last["settle"] = now
            # Issue #81: target-line refresh — Kalshi + Odds API only, zero
            # book requests. New games appear here, so the game-metadata
            # caches drop inside the tick.
            if now - last["targets"] >= config.TARGET_LINE_REFRESH_SEC:
                try:
                    _target_line_tick()
                except Exception as e:
                    log.error("target-line err: %s", e)
                last["targets"] = now
            # Issue #81: periodic on-demand coverage summary — the research
            # record of which books are actually answering live fetches.
            if now - last["coverage"] >= config.COVERAGE_SUMMARY_SEC:
                try:
                    _coverage_summary_tick(sgp_service=sgp_service, rfq_source=source)
                except Exception as e:
                    log.error("coverage err: %s", e)
                last["coverage"] = now
            # Issue #50: structure-only warming — no pricing calls, keeps
            # every book's on-demand structure caches + CZR token hot so an
            # RFQ never pays the cold-start penalty. Runs inline on the tick
            # thread on its own FAST cadence: warming is what keeps live
            # fetches quick.
            if now - last["warm"] >= config.STRUCTURE_WARM_SEC:
                try:
                    warmed = sgp_runner.warm_cycle(
                        bot_market_db=str(config.MARKET_DB),
                        service=sgp_service,
                        cadence_sec=config.STRUCTURE_WARM_SEC,
                        wall_budget_sec=config.STRUCTURE_WARM_BUDGET_SEC)
                    log.info("structure warm: %s", warmed)
                except Exception as e:
                    log.error("warm err: %s", e)
                last["warm"] = now
            research.flush()  # per-tick batched drain
            sgp_service.flush_health()   # issue #38: drain fetch-health rows
            sgp_service.check_book_health()   # issue #37: book-death alerts
            time.sleep(0.25)   # short sleep → responsive SIGTERM
    finally:
        # Stop the WS reader first (no-op for RestRFQSource) so shutdown's
        # quote-cancel sweep doesn't race reconnect attempts and their logs.
        if hasattr(source, "stop"):
            source.stop()
        with db.connect(read_only=True) as con:
            live = [r[0] for r in con.execute(
                "SELECT quote_id FROM live_quotes WHERE status='open'").fetchall()]
        for qid in live:
            try:
                gateway.cancel(qid)
                with db.connect() as con:
                    con.execute(
                        "UPDATE live_quotes SET status='cancelled', closed_at=? WHERE quote_id=?",
                        [datetime.now(timezone.utc), qid])
            except Exception:
                pass
        # #81: drain the partial coverage window into the buffer, then the
        # buffer to disk — otherwise up to COVERAGE_SUMMARY_SEC of per-book
        # outcomes die with the process.
        try:
            _coverage_summary_tick(sgp_service=sgp_service, rfq_source=source)
        except Exception:
            pass
        research.flush()  # final drain before shutdown
        db.end_session(sid)
        log.info("=== shutdown complete ===")


def cli():
    p = argparse.ArgumentParser()
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args()
    signal.signal(signal.SIGINT, _signal_handler)
    signal.signal(signal.SIGTERM, _signal_handler)
    main_loop(dry_run=args.dry_run)


if __name__ == "__main__":
    cli()
