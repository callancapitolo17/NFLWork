"""ProphetX SGP orchestration library.

Pure function: `price_sgps(target_lines)` -> `list[PricedRow]`.

Consumes ``ProphetXClient`` (HTTP transport) and calls the existing
helper functions in ``scraper_prophetx_sgp.py`` (``match_events``,
``_find_market``, ``_pick_selection``, ``_verify_competitor_ids``) to
compose the four-combo SGP pricing per (game, period, spread_line,
total_line).

No DB I/O. Caller (scraper shim or bot sgp_runner) handles persistence.

Source-label tracking
---------------------
Rows produced by the main RFQ path are tagged ``prophetx_direct``.
The legacy scraper also produces ``prophetx_interpolated`` rows via
its integer-line fallback, but those are essentially never emitted in
practice (no PX-interpolated rows have ever appeared on main). To keep
this orchestrator small and focused on the live path, we expose the
fallback label as a module-level constant for parity with DK/FD but do
NOT invoke the integer-fallback path here — that path will continue to
live in the legacy scraper module for the dashboard shim's use.

Sanity filter (F5-Over defense)
-------------------------------
ProphetX's parlay pricer has a known systematic bug where F5-Over
combos occasionally return decimal odds 5-7x larger than the naive
independent leg-multiply would suggest. We apply the same defense the
legacy scraper uses: drop any combo whose ``parlay_decimal / (leg1_dec
* leg2_dec) > SANITY_MULT_RATIO``. Legitimate anti-correlated combos
reach ~1.15x; 1.5x leaves headroom while still catching the bug.

Helper-signature deviations from the original plan
--------------------------------------------------
The plan-spec didn't match the actual helper signatures shipped in
``scraper_prophetx_sgp.py``. Specifically:

* ``match_events(px_events, parlay_lines)`` consumes the legacy
  px_events dict shape (keys: ``px_event_id``, ``px_home``, ``px_away``,
  ``px_home_competitor_id``, ``px_away_competitor_id``, ``scheduled``)
  and the legacy parlay-lines dict shape. It returns matched dicts with
  ``px_event_id`` / ``px_home_competitor_id`` / ``px_away_competitor_id``
  / ``fg_spread_line`` / ``f5_spread_line`` / ``fg_total_line`` /
  ``f5_total_line``. We translate ``client.list_events()`` ``Event``
  dataclasses into that dict shape, then collapse per-period
  ``TargetLine`` rows into the parlay-lines dict shape.

* ``_pick_selection`` predicates receive raw selection dicts with keys
  ``id`` (outcomeId), ``lineID``, ``line``, ``competitorId``, ``name``
  (e.g. "Over 8.5"), and ``odds`` (single-leg American). We mirror the
  legacy scraper's predicates exactly: competitorId + signed line for
  spread legs, name-prefix + line for total legs.

* ``client.submit_parlay_rfq`` returns ``(chosen_offer, used_fallback)``
  where ``chosen_offer["odds"]`` is the **American** parlay odds, not
  decimal. We convert via ``american_to_decimal`` and round.

* No integer-line fallback path. PX's interpolated path has produced
  zero rows on main; we skip it here to keep the orchestrator small.
"""
from __future__ import annotations

import logging
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone

from mlb_sgp._shared import (BookTransportError, FetchCounters, PricedRow,
                             TargetLine, american_to_decimal,
                             decimal_to_american, fetch_counters,
                             price_tally_for)
from mlb_sgp.prophetx_client import ProphetXClient, SelectionLeg

logger = logging.getLogger(__name__)

BOOK_NAME = "prophetx"
SOURCE_LABEL = "prophetx_direct"
SOURCE_LABEL_FALLBACK = "prophetx_interpolated"

# F5-Over systematic-bug defense. See module docstring.
SANITY_MULT_RATIO = 1.5

# PX RFQs are a real market footprint, but width does NOT change RFQs
# *per cycle* (same targets priced) — only how bursty they are. PX is the
# cycle's long pole at low width; 6-wide gives comfortable ≤60s margin at
# production scale. Live ramp 2026-06-16 found ZERO backoff (429/403) up
# to 12-wide, so 6 is well within the safe range. Env-overridable; raise
# toward the (untriggered, >12) ceiling only if a future probe confirms.
PX_TARGET_PARALLELISM_DEFAULT = 6


def _resolve_parallelism(parallelism: int | None) -> int:
    if parallelism is not None:
        return parallelism
    return int(os.environ.get("MLB_SGP_PX_PARALLELISM",
                              str(PX_TARGET_PARALLELISM_DEFAULT)))

# Combo names — byte-identical to scraper_prophetx_sgp.py so the
# dashboard / kalshi_mlb_rfq leg lookups keep matching. The "F5 "
# prefix is applied to F5-period rows in the orchestrator.
_COMBO_NAMES = (
    "Home Spread + Over",
    "Home Spread + Under",
    "Away Spread + Over",
    "Away Spread + Under",
)


def _extract_offered_lines_px(
    markets: list[dict],
    period_key: str,
    home_id,
    find_market=None,
    spread_market_name: str | None = None,
    total_market_name: str | None = None,
) -> dict[str, set]:
    """Return ``{"spreads": set, "totals": set}`` for one period from PX
    raw event markets.

    PX stores each selection's ``line`` from the *outcome's* own
    perspective. To get a single home-perspective signed line per
    market, we restrict the spread iteration to selections whose
    ``competitorId == home_id``.

    Total selections are filtered to those whose name starts with
    "over" or "under" (defensive — some auxiliary selections sit in the
    same market).

    ``find_market`` is the market-lookup callable (signature
    ``(markets, name) -> dict | None``). In production we pass
    ``scraper_prophetx_sgp._find_market`` which honors PX's
    ``NAME_ALIASES`` for renamed markets — keeping that behavior in the
    Filter A path matters because a missed market would silently drop
    every target. ``spread_market_name`` / ``total_market_name`` are
    the PX market display names for this period (``"Run Line"`` /
    ``"Total Runs"`` for FG, ``"1st-5th Inning Spread"`` etc. for F5).
    All three default to ``None`` and are resolved from
    ``scraper_prophetx_sgp`` when omitted — that lazy lookup keeps the
    production call site unchanged while letting tests pass everything
    in directly (avoiding ``sys.path`` mucking that would shadow the
    project-root ``db`` module).

    Used by Filter A in ``price_sgps`` to drop targets PX doesn't offer
    before the pricing loop runs.
    """
    if find_market is None or spread_market_name is None or total_market_name is None:
        from scraper_prophetx_sgp import _find_market as _fm, MARKET_NAMES
        if find_market is None:
            find_market = _fm
        if spread_market_name is None:
            spread_market_name = MARKET_NAMES[period_key]["spread"]
        if total_market_name is None:
            total_market_name = MARKET_NAMES[period_key]["total"]

    spread_mkt = find_market(markets, spread_market_name)
    total_mkt = find_market(markets, total_market_name)
    spreads: set = set()
    totals: set = set()
    if spread_mkt:
        for sel in _iter_selections(spread_mkt):
            if sel.get("competitorId") != home_id:
                continue
            line_val = sel.get("line")
            if line_val is not None:
                spreads.add(float(line_val))
    if total_mkt:
        for sel in _iter_selections(total_mkt):
            nm = (sel.get("name") or "").lower()
            if not (nm.startswith("over") or nm.startswith("under")):
                continue
            line_val = sel.get("line")
            if line_val is not None:
                totals.add(float(line_val))
    return {"spreads": spreads, "totals": totals}


def price_sgps(
    target_lines: list[TargetLine],
    periods: tuple[str, ...] = ("FG",),
    client: ProphetXClient | None = None,
    verbose: bool = False,
    parallelism: int | None = None,
    *,
    counters: FetchCounters | None = None,
) -> list[PricedRow]:
    """Price every target line against the ProphetX RFQ pricer.

    Parameters
    ----------
    target_lines
        List of ``TargetLine`` records — one per (game, period, spread,
        total) tuple the caller wants priced.
    periods
        Which periods to price. ``("FG",)`` (bot default) or
        ``("FG", "F5")`` (dashboard). Targets in other periods are
        silently dropped.
    client
        Optional pre-built ``ProphetXClient``. If ``None``, a new one is
        created (which opens a real curl_cffi Chrome-TLS session). Tests
        should always pass a mock.
    verbose
        Forwarded to leg-level helpers for debug printing.
    parallelism
        Target-pool width. When ``None``, resolves to the
        ``MLB_SGP_PX_PARALLELISM`` env var, else
        ``PX_TARGET_PARALLELISM_DEFAULT``. Default kept at 2 because PX
        RFQs are a real market footprint — do not widen without a probe.

    Returns
    -------
    list[PricedRow]
        Up to four PricedRow per (game, period) — one per combo. Combos
        that fail the SANITY_MULT_RATIO check are silently dropped.

    This is a thin wrapper over ``_price_sgps`` so issue #35's per-fetch
    counter summary and parse tripwire are emitted on every exit path,
    early returns included.

    ``counters`` (issue #38, keyword-only) prices into a CALLER-owned
    ``FetchCounters`` instead of a fresh one, so ``SGPService`` can read
    this fetch's counts afterwards and persist them to
    ``sgp_fetch_health``. Omitting it keeps the pre-#38 behavior.
    """
    with fetch_counters(BOOK_NAME, "sweep", logger,
                        counters=counters) as counters:
        return _price_sgps(target_lines, periods, client, verbose,
                           parallelism, counters)


def _price_sgps(
    target_lines: list[TargetLine],
    periods: tuple[str, ...],
    client: ProphetXClient | None,
    verbose: bool,
    parallelism: int | None,
    counters: FetchCounters,
) -> list[PricedRow]:
    if not target_lines:
        return []

    # Filter to requested periods early; if nothing remains, bail before
    # touching the network / building a client.
    targets = [t for t in target_lines if t.period in periods]
    if not targets:
        return []

    # Build the client lazily so tests that filter everything out via
    # the `periods` argument don't need to mock any HTTP at all.
    if client is None:
        client = ProphetXClient(verbose=verbose)

    # Import legacy helpers lazily — keeps the period-filter early-exit
    # test free of any HTTP imports.
    from scraper_prophetx_sgp import (
        match_events,
        _find_market,
        _pick_selection,
        _verify_competitor_ids,
        MARKET_NAMES,
    )

    # Issue #33: the RFQ endpoint can be blocked while the read endpoints
    # still work. client.submit_parlay_rfq records into this tally; a cycle
    # where EVERY price call failed raises a "price"-stage transport error.
    price_calls = price_tally_for(client, BOOK_NAME)
    tally_start = price_calls.snapshot()

    # ----- Group target lines into the legacy parlay_lines dict shape ----- #
    # match_events expects: {game_id: {home_team, away_team, commence_time,
    # fg_spread_line, fg_total_line, f5_spread_line, f5_total_line}}.
    parlay_lines: dict[str, dict] = {}
    for t in targets:
        ent = parlay_lines.setdefault(t.game_id, {
            "home_team": t.home_team,
            "away_team": t.away_team,
            "commence_time": t.commence_time,
            "fg_spread_line": None,
            "fg_total_line": None,
            "f5_spread_line": None,
            "f5_total_line": None,
        })
        if t.period == "FG":
            ent["fg_spread_line"] = t.spread
            ent["fg_total_line"] = t.total
        elif t.period == "F5":
            ent["f5_spread_line"] = t.spread
            ent["f5_total_line"] = t.total

    # ----- Phase 1: list PX events, translate to legacy dict shape, match ----- #
    events = client.list_events()
    # match_events consumes the legacy px_events shape — translate Events.
    # Legacy keys: px_event_id, px_home, px_away, px_home_competitor_id,
    # px_away_competitor_id, scheduled.
    px_events = [
        {
            "px_event_id": e.event_id,
            "px_home": e.home_team,
            "px_away": e.away_team,
            "px_home_competitor_id": e.home_id,
            "px_away_competitor_id": e.away_id,
            "scheduled": e.start_time,
        }
        for e in events
    ]
    counters.bump("events_seen", len(px_events))
    matched = match_events(px_events, parlay_lines)
    counters.bump("events_matched", len(matched or []))
    if not matched:
        return []

    matched_by_gid = {m["game_id"]: m for m in matched}

    out: list[PricedRow] = []
    fetch_now = datetime.now(timezone.utc)

    # Cache markets per game_id so two TargetLines on the same game
    # (e.g. FG + F5) only hit the event-markets endpoint once.
    markets_cache: dict[str, list[dict]] = {}

    # ----- Phase 1.6: group + cache + filter (Filter A) ----- #
    # Fetch per-game markets up-front so we can know what lines PX
    # offers before iterating targets, then drop targets PX doesn't
    # carry. PX has no integer-line fallback path in this orchestrator,
    # so the filter is strict: spread and total must both be offered.
    targets_by_game: dict[str, list[TargetLine]] = {}
    for t in targets:
        targets_by_game.setdefault(t.game_id, []).append(t)

    filtered_targets: list[TargetLine] = []
    skip_games: set[str] = set()
    for game_id, game_targets in targets_by_game.items():
        game = matched_by_gid.get(game_id)
        if game is None:
            skip_games.add(game_id)
            continue
        if game_id not in markets_cache:
            market_objs = client.fetch_event_markets(game["px_event_id"])
            markets_cache[game_id] = [
                {
                    "id": m.market_id,
                    "name": m.name,
                    "marketLines": m.selections,
                    "outcomes": m.outcomes,
                    # Moneyline order book (top-level `selections`), for the ML
                    # standalone price used by the SANITY_MULT_RATIO check.
                    "selections": m.ml_selections,
                }
                for m in market_objs
            ]
        markets = markets_cache[game_id]

        home_id = game["px_home_competitor_id"]
        away_id = game["px_away_competitor_id"]
        if not _verify_competitor_ids(markets, home_id, away_id):
            skip_games.add(game_id)
            continue

        # Build offered (spread, total) sets per period from the cached
        # raw-dict markets via the testable _extract_offered_lines_px
        # helper. PX stores `line` on each selection from the outcome's
        # own perspective, so we restrict spread lines to the home_id
        # outcomes to get a home-perspective set.
        offered_per_period: dict[str, dict[str, set]] = {
            period_key: _extract_offered_lines_px(markets, period_key, home_id)
            for period_key in ("fg", "f5")
        }

        pre = 0
        post = 0
        for t in game_targets:
            pre += 1
            offered = offered_per_period.get(t.period.lower())
            if offered is None:
                continue
            if float(t.spread) not in offered["spreads"]:
                continue
            if float(t.total) not in offered["totals"]:
                continue
            filtered_targets.append(t)
            post += 1
        counters.bump("legs_attempted", pre)
        counters.bump("legs_resolved", post)
        logger.debug("%s: game %s: %d kalshi -> %d offered",
                     BOOK_NAME, game_id, pre, post)
    counters.bump("targets_attempted", len(filtered_targets))

    # ----- Phase 2: per target row, build 4 legs and price 4 combos ----- #
    # Layer target-level concurrency on top of the per-combo parallelism
    # inside _price_combos_parallel. 2 concurrent targets × 4 combos = 8
    # in-flight RFQs per book — conservative given PX's rate limit.
    def _price_one_target(t: TargetLine) -> list[PricedRow]:
        """Compute up to 4 PricedRows for a single (game, period, spread, total).

        Pure function over the surrounding closure (markets_cache,
        matched_by_gid, skip_games, client, fetch_now, verbose). Returns
        an empty list when the target can't be priced (any missing leg,
        any combo dropped by SANITY_MULT_RATIO).
        """
        if t.game_id in skip_games:
            return []
        game = matched_by_gid.get(t.game_id)
        if game is None:
            return []
        period_key = "fg" if t.period == "FG" else "f5"
        markets = markets_cache[t.game_id]
        home_id = game["px_home_competitor_id"]
        away_id = game["px_away_competitor_id"]

        spread_mkt = _find_market(markets, MARKET_NAMES[period_key]["spread"])
        total_mkt = _find_market(markets, MARKET_NAMES[period_key]["total"])
        if not spread_mkt or not total_mkt:
            return []

        # Pick the four legs at the target spread/total lines. PX stores
        # each outcome's line from its own perspective: home_spread = +1.5
        # means home outcome has line == +1.5.
        home_sel = _pick_selection(
            spread_mkt,
            lambda s, lv=t.spread, hid=home_id: (
                s.get("competitorId") == hid and _line_eq(s.get("line"), lv)
            ),
        )
        away_sel = _pick_selection(
            spread_mkt,
            lambda s, lv=-t.spread, aid=away_id: (
                s.get("competitorId") == aid and _line_eq(s.get("line"), lv)
            ),
        )
        over_sel = _pick_selection(
            total_mkt,
            lambda s, lv=t.total: (
                (s.get("name") or "").lower().startswith("over")
                and _line_eq(s.get("line"), lv)
            ),
        )
        under_sel = _pick_selection(
            total_mkt,
            lambda s, lv=t.total: (
                (s.get("name") or "").lower().startswith("under")
                and _line_eq(s.get("line"), lv)
            ),
        )
        if not (home_sel and away_sel and over_sel and under_sel):
            return []

        prefix = "" if t.period == "FG" else "F5 "
        combos = (
            ("Home Spread + Over",  home_sel, over_sel),
            ("Home Spread + Under", home_sel, under_sel),
            ("Away Spread + Over",  away_sel, over_sel),
            ("Away Spread + Under", away_sel, under_sel),
        )

        spread_mid = spread_mkt.get("id")
        total_mid = total_mkt.get("id")

        priced_by_combo = _price_combos_parallel(
            client, game["px_event_id"], spread_mid, total_mid,
            combos, counters=counters,
        )

        target_rows: list[PricedRow] = []
        for combo_name, _sp_sel, _tot_sel in combos:
            sgp_decimal = priced_by_combo.get(combo_name)
            if sgp_decimal is None:
                continue
            target_rows.append(PricedRow(
                game_id=t.game_id,
                combo=prefix + combo_name,
                period=t.period,
                spread_line=t.spread,
                total_line=t.total,
                bookmaker=BOOK_NAME,
                source=SOURCE_LABEL,
                sgp_decimal=round(sgp_decimal, 4),
                sgp_american=decimal_to_american(sgp_decimal),
                fetch_time=fetch_now,
            ))
        return target_rows

    with ThreadPoolExecutor(max_workers=max(1, _resolve_parallelism(parallelism))) as pool:
        futures = [pool.submit(_price_one_target, t) for t in filtered_targets]
        for f in as_completed(futures):
            try:
                out.extend(f.result())
            except Exception as e:
                counters.record_exception("price_target", logger, e)

    # ----- Phase 3: moneyline × total combos (FG) ----- #
    out.extend(_price_ml_total_for_games(
        client, markets_cache, matched_by_gid, filtered_targets,
        fetch_now, parallelism, verbose, counters=counters,
    ))

    counters.bump("prices_returned", len(out))
    price_calls.verdict(tally_start)   # raises if EVERY price call failed
    return out


def _price_ml_total_for_games(
    client, markets_cache, matched_by_gid, filtered_targets,
    fetch_now, parallelism, verbose, *,
    counters: FetchCounters | None = None,
) -> list[PricedRow]:
    """Price the 4 moneyline×total combos per (game, distinct FG total) on PX.

    PX nests the moneyline ORDER BOOK under the market's `selections` key (one
    best-first ladder per outcome), NOT under `marketLines` like spread/total.
    The top ladder entry carries the best standalone odds plus id/lineID/line —
    so the ML leg is a normal priced selection, and ML×total flows through the
    SAME _price_combos_parallel + SANITY_MULT_RATIO path as spread×total (full
    F5-Over-style defense, no special-casing). spread_line=None.
    """
    from scraper_prophetx_sgp import _find_market, _pick_selection, MARKET_NAMES

    totals_by_game: dict[str, set] = {}
    for t in filtered_targets:
        if t.period == "FG":
            totals_by_game.setdefault(t.game_id, set()).add(t.total)

    def _pick_ml_selection(ml_market: dict, competitor_id):
        """Best (first) order-book entry for a competitor from the ML book.

        The entry carries odds + id + lineID, so it serves as both the RFQ leg
        and the single-leg sanity price. A side with no resting orders (empty
        book) returns None → that game's ML is skipped (can't sanity-check)."""
        for ladder in ml_market.get("selections") or []:
            if (ladder and isinstance(ladder, list) and isinstance(ladder[0], dict)
                    and ladder[0].get("competitorId") == competitor_id):
                return ladder[0]
        return None

    def _price_one_game(game_id: str, totals: set) -> list[PricedRow]:
        rows: list[PricedRow] = []
        game = matched_by_gid.get(game_id)
        markets = markets_cache.get(game_id)
        if not game or not markets:
            return rows
        ml_mkt = _find_market(markets, MARKET_NAMES["fg"]["moneyline"])
        total_mkt = _find_market(markets, MARKET_NAMES["fg"]["total"])
        if not ml_mkt or not total_mkt:
            return rows
        home_sel = _pick_ml_selection(ml_mkt, game["px_home_competitor_id"])
        away_sel = _pick_ml_selection(ml_mkt, game["px_away_competitor_id"])
        if not (home_sel and away_sel):
            return rows
        ml_mid = ml_mkt.get("id")
        total_mid = total_mkt.get("id")
        event_id = game["px_event_id"]
        for total in totals:
            over_sel = _pick_selection(
                total_mkt,
                lambda s, lv=total: (
                    (s.get("name") or "").lower().startswith("over")
                    and _line_eq(s.get("line"), lv)))
            under_sel = _pick_selection(
                total_mkt,
                lambda s, lv=total: (
                    (s.get("name") or "").lower().startswith("under")
                    and _line_eq(s.get("line"), lv)))
            if not (over_sel and under_sel):
                continue
            ml_combos = (
                ("Home ML + Over", home_sel, over_sel),
                ("Home ML + Under", home_sel, under_sel),
                ("Away ML + Over", away_sel, over_sel),
                ("Away ML + Under", away_sel, under_sel),
            )
            priced = _price_combos_parallel(
                client, event_id, ml_mid, total_mid, ml_combos,
                counters=counters)
            for combo_name, _ml, _tot in ml_combos:
                dec = priced.get(combo_name)
                if dec is None:
                    continue
                rows.append(PricedRow(
                    game_id=game_id, combo=combo_name, period="FG",
                    spread_line=None, total_line=total,
                    bookmaker=BOOK_NAME, source=SOURCE_LABEL,
                    sgp_decimal=round(dec, 4),
                    sgp_american=decimal_to_american(dec),
                    fetch_time=fetch_now))
        return rows

    out: list[PricedRow] = []
    if not totals_by_game:
        return out
    with ThreadPoolExecutor(max_workers=max(1, _resolve_parallelism(parallelism))) as pool:
        futures = [pool.submit(_price_one_game, gid, totals)
                   for gid, totals in totals_by_game.items()]
        for f in as_completed(futures):
            try:
                out.extend(f.result())
            except Exception as e:
                counters.record_exception("price_ml_target", logger, e)
    return out


def _price_combos_parallel(
    client: ProphetXClient,
    px_event_id: str,
    spread_mid,
    total_mid,
    combos: tuple,
    *,
    counters: FetchCounters | None = None,
) -> dict:
    """Price the 4 combo flavors in parallel and return {combo_name: sgp_decimal}.

    Each combo submits a parlay RFQ via ``client.submit_parlay_rfq`` on
    the shared curl_cffi session. The SANITY_MULT_RATIO filter (F5-Over
    defense) runs *inside* the worker so a combo that fails sanity is
    omitted from the result dict — same behavior as the sequential path,
    just executed concurrently.

    ``spread_mid`` is the first leg's market id — a spread market for
    spread×total combos, or the moneyline market for ML×total combos (the
    ML pricer reuses this same path now that the ML leg carries odds). The
    ``combos`` tuples are ``(name, leg1_sel, total_sel)`` either way.

    ``counters`` is keyword-only and optional (direct-call tests).
    """
    counters = counters or FetchCounters(BOOK_NAME, "sweep")
    results: dict = {}

    def _price_one(combo_name, sp_sel, tot_sel):
        try:
            legs = [
                _selection_to_leg(px_event_id, spread_mid, sp_sel),
                _selection_to_leg(px_event_id, total_mid, tot_sel),
            ]
            offer, _used_fallback = client.submit_parlay_rfq(legs)
            if offer is None:
                return combo_name, None

            am_parlay = offer.get("odds")
            if am_parlay is None:
                return combo_name, None
            try:
                sgp_decimal = american_to_decimal(int(am_parlay))
            except (TypeError, ValueError):
                return combo_name, None

            # SANITY filter — block F5-Over systematic mispricing.
            leg1_dec = _safe_leg_decimal(sp_sel)
            leg2_dec = _safe_leg_decimal(tot_sel)
            if not _passes_sanity_mult_ratio(sgp_decimal, leg1_dec, leg2_dec):
                counters.bump("sanity_drops")
                # DEBUG, not WARNING: one line per dropped combo across a
                # thread pool. The per-fetch count rides in the summary.
                logger.debug("%s: SANITY-DROP %s parlay=%.2f naive=%.2f",
                             BOOK_NAME, combo_name, sgp_decimal,
                             leg1_dec * leg2_dec)
                return combo_name, None

            return combo_name, sgp_decimal
        except Exception as e:
            counters.record_exception("price_combo", logger, e)
            return combo_name, None

    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = [
            pool.submit(_price_one, combo_name, sp_sel, tot_sel)
            for combo_name, sp_sel, tot_sel in combos
        ]
        for fut in as_completed(futures):
            try:
                combo_name, sgp_decimal = fut.result()
            except Exception as e:
                counters.record_exception("price_combo_future", logger, e)
                continue
            if sgp_decimal is not None:
                results[combo_name] = sgp_decimal
    return results


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _passes_sanity_mult_ratio(
    parlay_decimal: float, leg1_dec: float, leg2_dec: float
) -> bool:
    """ProphetX F5-Over defensive filter.

    The PX parlay pricer occasionally returns decimals 5-7x larger than
    the naive independent leg-multiply would suggest. We block any
    combo whose parlay decimal exceeds ``SANITY_MULT_RATIO`` times the
    independent multiply. Legitimate anti-correlated combos reach ~1.15x;
    1.5x is enough headroom while still catching the bug.

    Returns False (drop the combo) for non-positive naive multiplies —
    defensive against missing single-leg odds upstream.
    """
    naive = leg1_dec * leg2_dec
    if naive <= 0:
        return False
    return parlay_decimal / naive <= SANITY_MULT_RATIO


def _iter_selections(market: dict):
    """Yield every concrete selection dict under a PX market.

    PX market lines come in two shapes (mirrors ``_pick_selection``):
    - ``marketLines[].selections`` — nested per-side lists (with Nones
      sprinkled in when only one side is offered),
    - ``marketLines[].outcomes`` — flat list,
    - top-level ``market.selections`` — used by moneyline.

    We walk all three shapes and yield every non-None selection dict, so
    callers can collect line values without caring which shape PX
    returned for this particular market.
    """
    for line_grp in market.get("marketLines") or []:
        for side in line_grp.get("selections") or []:
            if side is None:
                continue
            for sel in side:
                if sel is None:
                    continue
                yield sel
        for sel in line_grp.get("outcomes") or []:
            if sel is None:
                continue
            yield sel
    for side in market.get("selections") or []:
        if side is None:
            continue
        for sel in side:
            if sel is None:
                continue
            yield sel


def _line_eq(a, b, eps: float = 1e-6) -> bool:
    """Float-safe equality for line values (halves or quarters)."""
    if a is None or b is None:
        return False
    try:
        return abs(float(a) - float(b)) < eps
    except (TypeError, ValueError):
        return False


def _selection_to_leg(event_id: str, market_id, sel: dict) -> SelectionLeg:
    """Build a ``SelectionLeg`` from a raw PX selection dict.

    PX selection-dict keys: ``id`` (outcomeId), ``lineID``, ``line``.
    """
    return SelectionLeg(
        sport_event_id=event_id,
        market_id=str(market_id) if market_id is not None else "",
        outcome_id=str(sel.get("id", "")),
        line_id=str(sel.get("lineID", "")),
        line=sel.get("line", 0.0),
    )


def _safe_leg_decimal(sel: dict) -> float:
    """Extract a single-leg decimal price from a PX selection dict.

    Returns 0.0 when the selection has no ``odds`` field — the SANITY
    filter then treats the combo as suspect (naive=0 -> drops it).
    """
    am = sel.get("odds")
    if am is None:
        return 0.0
    try:
        return american_to_decimal(int(am))
    except (TypeError, ValueError):
        return 0.0


# ---------------------------------------------------------------------------
# Phase 2 on-demand pricing (kalshi_mlb_mm) — resolve_legs + price_selection_set
#
# Append-only additions: nothing above this banner (in particular
# ``price_sgps``) is touched. The ``_od_*`` helpers mirror the exact
# selection-picking mechanics the orchestrator uses (``_find_market`` /
# ``_pick_selection`` from scraper_prophetx_sgp, and the nested
# ``_pick_ml_selection`` in ``_price_ml_total_for_games``) WITHOUT importing
# the legacy scraper module — its top-level imports (curl_cffi,
# canonical_match, db) are heavyweight and shadow the project-root ``db``
# module in tests. The FG market names have no NAME_ALIASES entries, so an
# exact-name match is behavior-identical to ``_find_market`` here.
# ---------------------------------------------------------------------------

# (market_type, period) -> market display-name aliases, first hit wins (#85).
# FG names are scraper_prophetx_sgp.MARKET_NAMES["fg"] (no aliases exist).
# F5 names are that module's 2026-04-24 catalogue + its NAME_ALIASES — the
# 2026-08-12 recon board carried ZERO F5 markets at PX, so these are the best
# known names; an F5 leg with no matching market declines the book cleanly.
# No F5 moneyline name has ever been catalogued (and F5-winner legs are
# service-blocked until #86 anyway), so ("ml", "F5") is deliberately absent.
_OD_MARKET_NAMES = {
    ("spread", "FG"): ("Run Line",),
    ("total", "FG"): ("Total Runs",),
    ("ml", "FG"): ("Moneyline",),
    ("spread", "F5"): ("1st-5th Inning Spread", "1st-5th Inning Run Line",
                       "1st 5 Innings Spread", "F5 Run Line"),
    ("total", "F5"): ("1st-5th Inning Total Runs", "1st 5 Innings Total Runs",
                      "F5 Total Runs"),
}


def _od_normalize_markets(markets) -> list[dict]:
    """``Market`` dataclasses (fetch_event_markets output) or raw dicts ->
    the raw-dict shape all pickers consume. Mirrors the markets_cache
    construction in ``price_sgps`` (moneyline order book under top-level
    ``selections``, spread/total lines under ``marketLines``)."""
    out: list[dict] = []
    for m in markets or []:
        if isinstance(m, dict):
            out.append(m)
        else:
            out.append({
                "id": m.market_id,
                "name": m.name,
                "marketLines": m.selections,
                "outcomes": m.outcomes,
                "selections": m.ml_selections,
            })
    return out


def _od_find_market(markets: list[dict], names: tuple) -> dict | None:
    """Exact-name market lookup over an alias tuple, first hit wins —
    the scraper_prophetx_sgp.NAME_ALIASES pattern."""
    for name in names:
        for m in markets:
            if m.get("name") == name:
                return m
    return None


def _od_pick_selection(market: dict, predicate) -> dict | None:
    """First selection matching ``predicate`` — same walk order as
    scraper_prophetx_sgp._pick_selection (marketLines selections/outcomes,
    then top-level selections), via the shared ``_iter_selections``."""
    for sel in _iter_selections(market):
        if predicate(sel):
            return sel
    return None


def _od_best_ml_entry(ml_market: dict, competitor_id) -> dict | None:
    """Best (first) order-book entry for a competitor from the ML book.

    Mirrors ``_pick_ml_selection`` inside ``_price_ml_total_for_games``:
    PX nests the moneyline order book under top-level ``selections`` as one
    best-first ladder per outcome. A side with no resting orders returns
    None."""
    for ladder in ml_market.get("selections") or []:
        if (ladder and isinstance(ladder, list) and isinstance(ladder[0], dict)
                and ladder[0].get("competitorId") == competitor_id):
            return ladder[0]
    return None


def _od_leg_decimal(sel: dict) -> float | None:
    """Single-leg decimal odds for a ResolvedLeg, or None.

    Same american->decimal conversion as ``_safe_leg_decimal`` but maps the
    "unpriceable" cases (odds absent, zero, unparseable) to None instead of
    0.0 — ResolvedLeg's contract is ``single_decimal=None`` when the book's
    structure carries no usable odds. Zero is guarded up front because
    ``american_to_decimal(0)`` divides by zero."""
    am = sel.get("odds")
    if not am:
        return None
    try:
        return american_to_decimal(int(am))
    except (TypeError, ValueError, ZeroDivisionError):
        return None


def resolve_legs(structure, legs, home_team, away_team, *,
                 counters: FetchCounters | None = None):
    """Resolve CanonicalLegs to ProphetX SelectionLeg refs for one event.

    Parameters
    ----------
    structure : dict
        Per-event PX structure::

            {"event_id": str,          # PX sportEvent id (RFQ payload key)
             "markets":  list,         # ProphetXClient.fetch_event_markets()
                                       # output (Market dataclasses) or the
                                       # equivalent raw dicts
             "home_id":  int,          # PX competitor ids from the matched
             "away_id":  int}          # Event (list_events / match_events)

        The competitor ids and event id ride along because PX's leg pickers
        select spread/ML sides by ``competitorId`` (not team name) and every
        RFQ leg embeds the sportEvent id — the market list alone is not
        enough context.
    legs : list[CanonicalLeg]
        market_type in {"spread","total","ml"}; spread/total ``line`` is
        SIGNED home-perspective (negative = home favored); side in
        {"home","away"} (spread/ml) or {"over","under"} (total).
    home_team, away_team : str
        Unused for PX (side selection is competitor-id based); accepted for
        cross-book interface parity.

    Returns
    -------
    list[ResolvedLeg] | None
        One ResolvedLeg per input leg, order preserved. ``ref`` /
        ``opposite_ref`` are ``SelectionLeg``s ready for
        ``price_selection_set``. None when any leg's CHOSEN side is
        unresolvable; a missing OPPOSITE side yields ``opposite_ref=None``
        (routes the book to Route B). Never raises.
    """
    from mlb_sgp._shared import ResolvedLeg

    try:
        event_id = structure["event_id"]
        markets = _od_normalize_markets(structure["markets"])
        home_id = structure["home_id"]
        away_id = structure["away_id"]

        out: list[ResolvedLeg] = []
        for leg in legs:
            # Period dispatch — the ONE place period is applied at PX (#85):
            # an unknown (market_type, period) pair declines the whole book.
            mkt_names = _OD_MARKET_NAMES.get((leg.market_type, leg.period))
            if mkt_names is None:
                return None
            market = _od_find_market(markets, mkt_names)
            if market is None:
                return None

            if leg.market_type == "spread":
                if leg.side not in ("home", "away"):
                    return None
                # PX stores each outcome's line from its OWN perspective:
                # home outcome carries leg.line, away outcome -leg.line
                # (same predicates as price_sgps' home_sel / away_sel).
                home_pred = (lambda s, lv=leg.line, cid=home_id:
                             s.get("competitorId") == cid
                             and _line_eq(s.get("line"), lv))
                away_pred = (lambda s, lv=-leg.line, cid=away_id:
                             s.get("competitorId") == cid
                             and _line_eq(s.get("line"), lv))
                if leg.side == "home":
                    chosen = _od_pick_selection(market, home_pred)
                    opposite = _od_pick_selection(market, away_pred)
                else:
                    chosen = _od_pick_selection(market, away_pred)
                    opposite = _od_pick_selection(market, home_pred)
            elif leg.market_type == "total":
                if leg.side not in ("over", "under"):
                    return None
                other = "under" if leg.side == "over" else "over"
                # Same name-prefix + line predicates as price_sgps'
                # over_sel / under_sel; opposite is the other side at the
                # SAME line.
                def _tot_pred(prefix, lv=leg.line):
                    return lambda s: (
                        (s.get("name") or "").lower().startswith(prefix)
                        and _line_eq(s.get("line"), lv))
                chosen = _od_pick_selection(market, _tot_pred(leg.side))
                opposite = _od_pick_selection(market, _tot_pred(other))
            else:  # "ml"
                if leg.side not in ("home", "away"):
                    return None
                cid = home_id if leg.side == "home" else away_id
                opp_cid = away_id if leg.side == "home" else home_id
                chosen = _od_best_ml_entry(market, cid)
                opposite = _od_best_ml_entry(market, opp_cid)

            if chosen is None:
                return None
            mkt_id = market.get("id")
            out.append(ResolvedLeg(
                ref=_selection_to_leg(event_id, mkt_id, chosen),
                opposite_ref=(_selection_to_leg(event_id, mkt_id, opposite)
                              if opposite is not None else None),
                single_decimal=_od_leg_decimal(chosen),
                opposite_decimal=(_od_leg_decimal(opposite)
                                  if opposite is not None else None),
            ))
        return out
    except Exception as e:
        # Fail-safe contract: malformed structure/legs -> None, never raise.
        # The counter is what tells "book doesn't offer this leg" from "our
        # parser broke".
        if counters is not None:
            counters.record_parse_failure("resolve_legs", logger, e)
        return None


def price_selection_set(client, refs, *,
                        counters: FetchCounters | None = None,
                        ) -> float | None:
    """Price one arbitrary set of PX selections with a single RFQ.

    ``refs`` is a list of ``SelectionLeg`` (from ``resolve_legs`` refs /
    opposite_refs); a 1-element list prices a single. One wire call:
    ``client.submit_parlay_rfq(refs)`` -> ``(chosen_offer, used_fallback)``.
    ``offer["odds"]`` is the AMERICAN parlay odds — converted to decimal
    exactly like the existing orchestrator (``_price_combos_parallel``).
    ``used_fallback=True`` (thin-market teaser tier) still counts as priced.
    Returns None for empty refs, a None offer, or unparseable odds. Never
    raises.
    """
    if not refs:
        return None
    try:
        offer, _used_fallback = client.submit_parlay_rfq(refs)
        if offer is None:
            return None
        am_parlay = offer.get("odds")
        if am_parlay is None:
            return None
        return american_to_decimal(int(am_parlay))
    except BookTransportError:
        # Counted, NOT re-raised — see caesars.price_selection_set.
        if counters is not None:
            counters.bump("transport_errors")
        return None
    except Exception as e:
        if counters is not None:
            counters.record_parse_failure("price_selection_set", logger, e)
        return None
