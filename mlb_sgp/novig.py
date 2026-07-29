"""Novig SGP orchestration library.

Pure function: ``price_sgps(target_lines)`` -> ``list[PricedRow]``.

Consumes ``NovigClient`` (HTTP transport) and reuses the existing
helper functions in ``scraper_novig_sgp.py`` (``match_events``,
``fetch_event_legs``, ``try_integer_fallback_nv``, ``submit_parlay``)
to compose the four-combo SGP pricing per (game, period, spread_line,
total_line).

No DB I/O. Caller (scraper shim or bot sgp_runner) handles persistence.

Source-label tracking
---------------------
Rows produced by the main RFQ path are tagged ``novig_direct``.
Rows produced by the integer-line fallback path (off-main integer
totals, derived from adjacent half-point alts) are tagged
``novig_interpolated``. The dashboard parlay-tab reader filters on
``source IN (..._direct)``, so preserving these labels keeps the
existing dashboard behavior byte-identical post-refactor.

Note: Novig sources its legs from DraftKings (every observed leg
returns ``vendor=DRAFTKINGS``), so its line set is a strict subset
of DK's. Strict line matching at the client level means a missing
leg here typically reflects an off-main target — that's exactly when
the integer-line fallback path is useful.

Sanity filter
-------------
Mirrors the ProphetX defense: drop any combo whose parlay decimal
exceeds ``SANITY_MULT_RATIO`` (1.5x) times the naive independent
leg-product. Legitimate anti-correlated combos top out around ~1.15x;
1.5x leaves headroom while catching systematic mispricings.

Helper-signature deviations from the original plan
--------------------------------------------------
The plan spec referenced ``_select_outcome_ids_for_combo`` as a
function to lift from the legacy scraper. In practice the legacy
scraper doesn't have such a per-combo helper — its leg selection is
already done at the ``fetch_event_legs`` level (returning a per-
period dict ``{home_spread, away_spread, over, under}``). We reuse
``fetch_event_legs`` directly, which is simpler and keeps the
orchestrator behavior byte-identical to the scraper. The integer
fallback path is handled by lifting ``try_integer_fallback_nv``
verbatim.
"""
from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone

from mlb_sgp._shared import (BookTransportError, FetchCounters, PricedRow,
                             TargetLine, decimal_to_american, fetch_counters,
                             price_tally_for)
from mlb_sgp.novig_client import NovigClient

logger = logging.getLogger(__name__)

BOOK_NAME = "novig"
SOURCE_LABEL = "novig_direct"
SOURCE_LABEL_FALLBACK = "novig_interpolated"

# Sanity filter: parlay decimal must be <= naive_product * SANITY_MULT_RATIO.
# Legitimate anti-correlated parlays peak around 1.15x; 1.5x leaves headroom
# while still catching the systematic F5-Over-style mispricings observed on
# ProphetX (and defensively guarded against here in case Novig regresses).
SANITY_MULT_RATIO = 1.5

# Cross-target parallelism. NV's anonymous GraphQL endpoint tolerates
# more concurrency than PX's RFQ, so we run 4 concurrent targets ×
# 4 concurrent combos = 16 in-flight requests per book.
NV_TARGET_PARALLELISM = 4

# Combo names — byte-identical to scraper_novig_sgp.py so the dashboard /
# kalshi_mlb_rfq leg lookups keep matching. The "F5 " prefix is applied to
# F5-period rows in the orchestrator.
_COMBO_DISPLAY = {
    "home_over":  "Home Spread + Over",
    "home_under": "Home Spread + Under",
    "away_over":  "Away Spread + Over",
    "away_under": "Away Spread + Under",
}


def _extract_offered_lines_nv(
    markets: list[dict],
    period_key: str,
    spread_type: str | None = None,
    total_type: str | None = None,
) -> dict[str, set]:
    """Return ``{"spreads": set, "totals": set}`` for one period from NV
    raw GraphQL ``markets`` list.

    Each spread/total market on NV carries one ``strike`` value (the
    line). The home-perspective signed spread set is the union of
    strikes across all SPREAD-type markets for the period; the totals
    set is the union across TOTAL-type markets.

    ``spread_type`` and ``total_type`` are the NV market-type constants
    (``"SPREAD"``/``"TOTAL"`` for FG, ``"SPREAD_1H"``/``"TOTAL_1H"`` for
    F5). They default to ``None`` and are resolved from
    ``scraper_novig_sgp.SPREAD_TYPE`` / ``TOTAL_TYPE`` if not provided —
    that lazy lookup keeps the production call sites unchanged while
    letting tests pass the strings directly (so the helper can be
    exercised without putting ``mlb_sgp/`` on ``sys.path`` and
    accidentally shadowing the project-root ``db`` module).

    Used by Filter A in ``price_sgps`` to drop targets the book doesn't
    offer before the pricing loop runs.
    """
    if spread_type is None or total_type is None:
        from scraper_novig_sgp import SPREAD_TYPE, TOTAL_TYPE
        if spread_type is None:
            spread_type = SPREAD_TYPE[period_key]
        if total_type is None:
            total_type = TOTAL_TYPE[period_key]

    spreads = {
        float(m["strike"]) for m in markets
        if m.get("type") == spread_type and m.get("strike") is not None
    }
    totals = {
        float(m["strike"]) for m in markets
        if m.get("type") == total_type and m.get("strike") is not None
    }
    return {"spreads": spreads, "totals": totals}


def price_sgps(
    target_lines: list[TargetLine],
    periods: tuple[str, ...] = ("FG",),
    client: NovigClient | None = None,
    verbose: bool = False,
) -> list[PricedRow]:
    """Price every target line against the Novig parlay RFQ endpoint.

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
        Optional pre-built ``NovigClient``. If ``None``, a new one is
        created (which opens a real curl_cffi Chrome-TLS session).
        Tests should always pass a mock.
    verbose
        Forwarded to leg-level helpers for debug printing.

    Returns
    -------
    list[PricedRow]
        Up to four PricedRow per (game, period) — one per combo.
        Combos that fail the SANITY_MULT_RATIO check are silently
        dropped. Rows produced via the integer-fallback path are tagged
        with ``source = novig_interpolated`` instead of
        ``novig_direct``.

    This is a thin wrapper over ``_price_sgps`` so issue #35's per-fetch
    counter summary and parse tripwire are emitted on every exit path,
    early returns included.
    """
    with fetch_counters(BOOK_NAME, "sweep", logger) as counters:
        return _price_sgps(target_lines, periods, client, verbose, counters)


def _price_sgps(
    target_lines: list[TargetLine],
    periods: tuple[str, ...],
    client: NovigClient | None,
    verbose: bool,
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
        client = NovigClient(verbose=verbose)

    # Import legacy helpers lazily — keeps the period-filter early-exit
    # test free of any HTTP imports.
    from scraper_novig_sgp import (
        match_events,
        fetch_event_legs,
        try_integer_fallback_nv,
        submit_parlay,
    )

    # Issue #33: the parlay endpoint can be blocked while the GraphQL reads
    # still work. Tally the price calls and raise a "price"-stage transport
    # error if the whole cycle came back empty (see PriceCallTally).
    price_calls = price_tally_for(client, BOOK_NAME)
    # Novig's legacy submit_parlay returns (priced, auth_failed) — a 2-tuple is
    # always truthy, so score on the priced element.
    submit_parlay = price_calls.wrap(submit_parlay, ok=lambda r: bool(r and r[0]))
    tally_start = price_calls.snapshot()

    # ----- Group target lines into the legacy parlay-lines dict shape ----- #
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

    # ----- Phase 1: list Novig events and translate to legacy match shape ----- #
    # match_events consumes the legacy nv_events shape — translate the
    # client's Event dataclasses. Legacy keys: nv_event_id, nv_home,
    # nv_away, nv_home_sym, nv_away_sym, scheduled.
    events = client.list_events()
    nv_events = [
        {
            "nv_event_id": e.event_id,
            "nv_home": e.home_team,
            "nv_away": e.away_team,
            "nv_home_sym": e.home_sym,
            "nv_away_sym": e.away_sym,
            "scheduled": e.start_time,
        }
        for e in events
    ]
    counters.bump("events_seen", len(nv_events))
    matched = match_events(nv_events, parlay_lines)
    counters.bump("events_matched", len(matched or []))
    if not matched:
        return []

    matched_by_gid = {m["game_id"]: m for m in matched}

    out: list[PricedRow] = []
    fetch_now = datetime.now(timezone.utc)

    # Cache fetched legs + raw markets per game_id so two TargetLines on
    # the same game (FG + F5) only hit the GraphQL event-markets endpoint
    # once. fetch_event_legs returns (legs, markets) where `legs` is the
    # per-period dict and `markets` is the raw market list (needed by
    # try_integer_fallback_nv).
    legs_cache: dict[str, tuple[dict, list]] = {}

    # ----- Phase 1.6: group + cache + filter (Filter A) ----- #
    # NV markets are very sparse for alt lines (often only main +/- one
    # half-point neighbor). Drop targets where the spread isn't offered
    # at all OR the total isn't offered AND no half-point neighbors are
    # offered (the latter would let try_integer_fallback_nv rescue an
    # integer total). NV's spread match in fetch_event_legs is strict —
    # if NV's market strike differs from the target, no leg comes back.
    # (SPREAD_TYPE/TOTAL_TYPE are now imported inside the helper.)

    targets_by_game: dict[str, list[TargetLine]] = {}
    for t in targets:
        targets_by_game.setdefault(t.game_id, []).append(t)

    filtered_targets: list[TargetLine] = []
    for game_id, game_targets in targets_by_game.items():
        game = matched_by_gid.get(game_id)
        if game is None:
            continue
        if game_id not in legs_cache:
            legs_cache[game_id] = fetch_event_legs(client.session, game, verbose)
        _legs_dict, markets = legs_cache[game_id]

        # Build offered (spread, total) sets per period via the testable
        # _extract_offered_lines_nv helper.
        offered_per_period: dict[str, dict[str, set]] = {
            period_key: _extract_offered_lines_nv(markets, period_key)
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
            if float(t.total) in offered["totals"]:
                filtered_targets.append(t)
                post += 1
                continue
            # Fallback path requires integer target total + adjacent
            # half-point neighbors both offered.
            if float(t.total).is_integer() and (
                (float(t.total) - 0.5) in offered["totals"]
                and (float(t.total) + 0.5) in offered["totals"]
            ):
                filtered_targets.append(t)
                post += 1
        counters.bump("legs_attempted", pre)
        counters.bump("legs_resolved", post)
        logger.debug("%s: game %s: %d kalshi -> %d offered",
                     BOOK_NAME, game_id, pre, post)
    counters.bump("targets_attempted", len(filtered_targets))

    # ----- Phase 2: per target row, price 4 combos (or fall back) ----- #
    # Layer target-level concurrency on top of the per-combo parallelism
    # inside _price_combos_parallel. 4 concurrent targets × 4 combos =
    # 16 in-flight requests per book. NV's anonymous GraphQL endpoint
    # tolerates more concurrency than PX's RFQ.
    def _price_one_target(t: TargetLine) -> list[PricedRow]:
        """Compute up to 4 PricedRows for a single (game, period, spread, total).

        Pure function over the surrounding closure (legs_cache,
        matched_by_gid, client.session, fetch_now, verbose). Returns an
        empty list when the target can't be priced (any missing leg + no
        fallback hit).
        """
        game = matched_by_gid.get(t.game_id)
        if game is None:
            return []
        period_key = "fg" if t.period == "FG" else "f5"
        legs, markets = legs_cache[t.game_id]

        p = legs.get(period_key, {}) or {}
        home = p.get("home_spread")
        away = p.get("away_spread")
        over = p.get("over")
        under = p.get("under")

        prefix = "" if t.period == "FG" else "F5 "
        target_rows: list[PricedRow] = []

        # ----- Fallback path: any leg missing -> integer-line derivation ----- #
        if not (home and away and over and under):
            fallback = try_integer_fallback_nv(
                client.session, markets,
                game["nv_home_sym"], game["nv_away_sym"],
                t.spread, t.total, period_key,
                verbose=verbose,
            )
            if fallback is None:
                return []
            fair_probs = fallback["fair_probs"]
            for key, base_combo in _COMBO_DISPLAY.items():
                fair_p = fair_probs.get(key)
                if not fair_p or fair_p <= 0:
                    continue
                dec = 1.0 / fair_p
                target_rows.append(PricedRow(
                    game_id=t.game_id,
                    combo=prefix + base_combo,
                    period=t.period,
                    spread_line=t.spread,
                    total_line=t.total,
                    bookmaker=BOOK_NAME,
                    source=SOURCE_LABEL_FALLBACK,
                    sgp_decimal=round(dec, 4),
                    sgp_american=decimal_to_american(dec),
                    fetch_time=fetch_now,
                ))
            return target_rows

        # ----- Main path: price the 4 canonical combos via /parlay/request ----- #
        combos = (
            ("Home Spread + Over",  home, over),
            ("Home Spread + Under", home, under),
            ("Away Spread + Over",  away, over),
            ("Away Spread + Under", away, under),
        )

        priced_by_combo = _price_combos_parallel(
            client.session, combos, submit_parlay, verbose, counters=counters,
        )

        for combo_name, _sp, _to in combos:
            priced = priced_by_combo.get(combo_name)
            if priced is None:
                continue
            dec = priced["decimal"]
            target_rows.append(PricedRow(
                game_id=t.game_id,
                combo=prefix + combo_name,
                period=t.period,
                spread_line=t.spread,
                total_line=t.total,
                bookmaker=BOOK_NAME,
                source=SOURCE_LABEL,
                sgp_decimal=round(dec, 4),
                sgp_american=priced["american"],
                fetch_time=fetch_now,
            ))
        return target_rows

    with ThreadPoolExecutor(max_workers=NV_TARGET_PARALLELISM) as pool:
        futures = [pool.submit(_price_one_target, t) for t in filtered_targets]
        for f in as_completed(futures):
            try:
                out.extend(f.result())
            except Exception as e:
                counters.record_exception("price_target", logger, e)

    # ----- Phase 3: moneyline × total combos (FG only) ----- #
    # Priced ONCE per game (not per spread target) so a game with several
    # alt-spread targets doesn't re-issue identical ML parlay RFQs. Novig
    # resolves over/under at the game's single fg_total_line, so the ML×total
    # is labelled with THAT total (matched_by_gid[gid]["fg_total_line"]) — the
    # one the cached over/under legs actually belong to. spread_line=None.
    fg_games = {t.game_id for t in filtered_targets if t.period == "FG"}

    def _price_ml_game(gid: str) -> list[PricedRow]:
        cached = legs_cache.get(gid)
        if not cached:
            return []
        fg = (cached[0] or {}).get("fg", {})
        home_ml, away_ml = fg.get("home_ml"), fg.get("away_ml")
        over, under = fg.get("over"), fg.get("under")
        total = (matched_by_gid.get(gid) or {}).get("fg_total_line")
        if not (home_ml and away_ml and over and under) or total is None:
            return []
        ml_combos = (
            ("Home ML + Over",  home_ml, over),
            ("Home ML + Under", home_ml, under),
            ("Away ML + Over",  away_ml, over),
            ("Away ML + Under", away_ml, under),
        )
        ml_priced = _price_combos_parallel(
            client.session, ml_combos, submit_parlay, verbose,
            counters=counters)
        rows: list[PricedRow] = []
        for combo_name, _ml, _to in ml_combos:
            priced = ml_priced.get(combo_name)
            if priced is None:
                continue
            rows.append(PricedRow(
                game_id=gid, combo=combo_name, period="FG",
                spread_line=None, total_line=float(total),
                bookmaker=BOOK_NAME, source=SOURCE_LABEL,
                sgp_decimal=round(priced["decimal"], 4),
                sgp_american=priced["american"], fetch_time=fetch_now))
        return rows

    with ThreadPoolExecutor(max_workers=NV_TARGET_PARALLELISM) as pool:
        futs = [pool.submit(_price_ml_game, gid) for gid in fg_games]
        for f in as_completed(futs):
            try:
                out.extend(f.result())
            except Exception as e:
                counters.record_exception("price_ml_target", logger, e)

    counters.bump("prices_returned", len(out))
    price_calls.verdict(tally_start)   # raises if EVERY price call failed
    return out


def _price_combos_parallel(
    session,
    combos: tuple,
    submit_parlay_fn,
    verbose: bool,
    *,
    counters: FetchCounters | None = None,
) -> dict:
    """Price the 4 combo flavors in parallel and return {combo_name: priced}.

    Each combo submits a parlay RFQ via ``submit_parlay`` on the shared
    curl_cffi session. The SANITY_MULT_RATIO filter runs *inside* the
    worker so a combo that fails sanity is omitted from the result dict
    — same behavior as the sequential path, just executed concurrently.

    ``counters`` is keyword-only and optional (direct-call tests).
    """
    counters = counters or FetchCounters(BOOK_NAME, "sweep")
    results: dict = {}

    def _price_one(combo_name, sp, to):
        try:
            priced, _auth_failed = submit_parlay_fn(
                session, [sp["id"], to["id"]], verbose=verbose,
            )
            if priced is None:
                return combo_name, None

            # Sanity filter: drop combos where the parlay decimal exceeds
            # SANITY_MULT_RATIO * naive_leg_product. `available` is implied
            # probability so leg decimal = 1 / available.
            sp_av = sp.get("available")
            to_av = to.get("available")
            if not _passes_sanity_mult_ratio(priced["decimal"], sp_av, to_av):
                counters.bump("sanity_drops")
                # DEBUG, not WARNING: one line per dropped combo across a
                # thread pool. The per-fetch count rides in the summary.
                logger.debug("%s: SANITY-DROP %s parlay=%.2f",
                             BOOK_NAME, combo_name, priced["decimal"])
                return combo_name, None

            return combo_name, priced
        except Exception as e:
            counters.record_exception("price_combo", logger, e)
            return combo_name, None

    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = [
            pool.submit(_price_one, combo_name, sp, to)
            for combo_name, sp, to in combos
        ]
        for fut in as_completed(futures):
            try:
                combo_name, priced = fut.result()
            except Exception as e:
                counters.record_exception("price_combo_future", logger, e)
                continue
            if priced is not None:
                results[combo_name] = priced
    return results


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _passes_sanity_mult_ratio(
    parlay_decimal: float,
    sp_available: float | None,
    to_available: float | None,
) -> bool:
    """Drop combos where parlay decimal blows past naive independent-multiply.

    ``available`` is implied probability from the API. Leg decimal is
    therefore ``1 / available`` and the naive independent-multiply
    decimal is ``1 / (sp_available * to_available)``.

    Returns True (keep) when:
      - both ``available`` values are present and > 0, AND
      - parlay_decimal / naive_decimal <= SANITY_MULT_RATIO.

    Returns True (keep) when either ``available`` is missing — matches
    the legacy scraper's behavior of only running the check when both
    leg implied probs are available (otherwise the combo passes).
    """
    if not (sp_available and to_available):
        return True
    if sp_available <= 0 or to_available <= 0:
        return True
    naive = (1.0 / sp_available) * (1.0 / to_available)
    return parlay_decimal <= naive * SANITY_MULT_RATIO


# ---------------------------------------------------------------------------
# Phase 2 on-demand pricing (kalshi_mlb_mm on-demand engine).
# Append-only: nothing above this line (esp. price_sgps) is modified.
# ---------------------------------------------------------------------------

def _available_to_decimal(leg: dict | None) -> float | None:
    """Novig legs carry `available` = implied probability. Leg decimal is
    1/available, valid only for 0 < available < 1 (available==1 would be a
    decimal of exactly 1.0, i.e. no payout — treat as missing, matching the
    _passes_sanity_mult_ratio convention of only trusting sane probs)."""
    if not isinstance(leg, dict):
        return None
    av = leg.get("available")
    try:
        av = float(av)
    except (TypeError, ValueError):
        return None
    if not (0.0 < av < 1.0):
        return None
    return 1.0 / av


def _lookup_leg(structure: dict, key: str, line: float | None) -> dict | None:
    """Fetch one leg dict from the line-keyed fg bucket. EXACT line match
    only — the integer-line fallback (try_integer_fallback_nv) is
    deliberately NOT applied on the on-demand path (fail-safe: a line the
    book doesn't quote directly resolves to nothing).

    Spread/total keys map ``{line: leg}``; ML keys map straight to a leg
    dict (line-independent), signalled by ``line is None``.
    """
    bucket = structure.get(key)
    if line is None:                       # moneyline: value IS the leg
        return bucket if isinstance(bucket, dict) else None
    if not isinstance(bucket, dict):
        return None
    leg = bucket.get(float(line))
    return leg if isinstance(leg, dict) else None


def resolve_legs(structure, legs, home_team, away_team, *,
                 counters: FetchCounters | None = None):
    """Resolve canonical legs against Novig's FG per-event leg bucket.

    Parameters
    ----------
    structure
        Line-keyed FG bucket built from Novig's EventMarkets tree — the
        per-line generalization of ``fetch_event_legs``'s fg dict
        (scraper_novig_sgp.py:440-459, keys ``home_spread / away_spread /
        over / under / home_ml / away_ml``)::

            {
              "home_spread": {home_line: leg},  # keyed by HOME team's signed line
              "away_spread": {away_line: leg},  # keyed by AWAY team's line (= -home)
              "over":        {total_line: leg},
              "under":       {total_line: leg},
              "home_ml":     leg | None,        # line-independent
              "away_ml":     leg | None,
            }

        where each ``leg`` is the ``_find_outcome_in_spread`` /
        ``_find_outcome_in_total`` shape: ``{"id": <outcome uuid>,
        "available": <implied prob or None>}``. One NV SPREAD market at
        home-perspective strike -1.5 therefore lands in the bucket as
        ``home_spread[-1.5]`` + ``away_spread[+1.5]``.
    legs
        list[CanonicalLeg] — lines are SIGNED home-perspective (negative =
        home favored), so an away-side spread lookup mirrors the sign.
    home_team / away_team
        Unused for Novig (the bucket is already side-labelled); kept for
        the cross-book ``resolve_legs`` signature.

    Returns
    -------
    list[ResolvedLeg] | None
        None when ANY leg's chosen side is unresolvable (or on any
        unexpected input — never raises). A missing OPPOSITE side yields
        ``opposite_ref=None`` on that leg (one-sided → Route B).
        ``single_decimal`` / ``opposite_decimal`` = 1/available when
        0 < available < 1, else None.
    """
    # ResolvedLeg lives in _shared; import matches the module's lazy style.
    from mlb_sgp._shared import ResolvedLeg

    try:
        if not isinstance(structure, dict):
            return None
        out: list = []
        for leg in legs:
            mt, side, line = leg.market_type, leg.side, leg.line
            if mt == "spread":
                if side == "home":
                    chosen_key, opp_key = "home_spread", "away_spread"
                    chosen_line, opp_line = float(line), -float(line)
                elif side == "away":
                    chosen_key, opp_key = "away_spread", "home_spread"
                    chosen_line, opp_line = -float(line), float(line)
                else:
                    return None
            elif mt == "total":
                if side not in ("over", "under"):
                    return None
                chosen_key = side
                opp_key = "under" if side == "over" else "over"
                chosen_line = opp_line = float(line)
            elif mt == "ml":
                if side not in ("home", "away"):
                    return None
                chosen_key = f"{side}_ml"
                opp_key = "away_ml" if side == "home" else "home_ml"
                chosen_line = opp_line = None
            else:
                return None

            chosen = _lookup_leg(structure, chosen_key, chosen_line)
            if chosen is None or not chosen.get("id"):
                return None                     # chosen side miss → fail all
            opposite = _lookup_leg(structure, opp_key, opp_line)
            opp_ref = opposite.get("id") if isinstance(opposite, dict) else None
            out.append(ResolvedLeg(
                ref=chosen["id"],
                opposite_ref=opp_ref or None,
                single_decimal=_available_to_decimal(chosen),
                opposite_decimal=(
                    _available_to_decimal(opposite) if opp_ref else None
                ),
            ))
        return out
    except Exception as e:
        # Fail-safe stays (a raise here would kill the quote); the counter is
        # what tells "book doesn't offer this leg" from "our parser broke".
        if counters is not None:
            counters.record_parse_failure("resolve_legs", logger, e)
        return None


def price_selection_set(client, refs, *,
                        counters: FetchCounters | None = None,
                        ) -> float | None:
    """Price one arbitrary selection set via Novig's anonymous parlay RFQ.

    ``refs`` is a list of outcome UUID strings (a 1-element list prices a
    single). One wire call: ``client.submit_parlay(refs)``
    (novig_client.py:157) returns ``{"decimal": ..., ...}`` or ``{}`` on
    any failure. Returns the float decimal when > 1.0, else None.
    Never raises.
    """
    try:
        if not refs:
            return None
        priced = client.submit_parlay(list(refs)) or {}
        dec = float(priced.get("decimal"))
        return dec if dec > 1.0 else None
    except BookTransportError:
        # Counted, NOT re-raised — see caesars.price_selection_set.
        if counters is not None:
            counters.bump("transport_errors")
        return None
    except Exception as e:
        if counters is not None:
            counters.record_parse_failure("price_selection_set", logger, e)
        return None


def build_line_structure(markets: list, home_sym: str, away_sym: str) -> dict:
    """Build the per-line FG bucket ``resolve_legs`` consumes from Novig's
    RAW market tree (the ``markets`` list ``fetch_event_legs`` returns as its
    second element).

    Generalizes ``scraper_novig_sgp.fetch_event_legs``'s single-target-line
    parse (scraper_novig_sgp.py:391-443) to ALL offered FG lines, reusing its
    exact outcome extractors:

    * SPREAD markets: ``strike`` is the HOME-perspective line (the scraper
      matches ``m["strike"]`` against the home-perspective target at
      scraper_novig_sgp.py:402-404); outcomes are identified by
      ``competitor.symbol`` via ``_find_outcome_in_spread`` →
      ``home_spread[strike]`` / ``away_spread[-strike]``.
    * TOTAL markets: outcomes split on "Over "/"Under " description prefixes
      via ``_find_outcome_in_total`` → ``over[strike]`` / ``under[strike]``.
    * MONEY market: line-independent, symbol-matched like a spread →
      ``home_ml`` / ``away_ml`` (scraper_novig_sgp.py:438-443).

    Mirrors the scraper's ``is_consensus`` preference (scraper_novig_sgp.py:
    398-406): when several markets share a (type, strike), the consensus
    market wins regardless of list order. DEVIATION from ``fetch_event_legs``
    (documented): a ONE-SIDED line (only one outcome present) is still
    recorded — on-demand routes a missing opposite to Route B via
    ``opposite_ref=None``, whereas the grid sweep required both sides.
    FG only (SPREAD/TOTAL/MONEY; ``*_1H`` and prop types ignored).
    Fail-safe: never raises; malformed markets are skipped.
    """
    from scraper_novig_sgp import (
        _find_outcome_in_spread,
        _find_outcome_in_total,
    )

    out: dict = {"home_spread": {}, "away_spread": {}, "over": {},
                 "under": {}, "home_ml": None, "away_ml": None}

    # Group by (type, strike) so the is_consensus preference is applied
    # per line, exactly like fetch_event_legs' spread_matches/total_matches.
    groups: dict[tuple, list] = {}
    for m in markets or []:
        try:
            typ = m.get("type")
            if typ not in ("SPREAD", "TOTAL", "MONEY"):
                continue
            if typ == "MONEY":
                strike = None
            else:
                strike = float(m.get("strike"))
            groups.setdefault((typ, strike), []).append(m)
        except (AttributeError, TypeError, ValueError):
            continue

    for (typ, strike), cands in groups.items():
        m = next((c for c in cands if c.get("is_consensus") is True), cands[0])
        try:
            if typ == "SPREAD":
                home_leg, away_leg = _find_outcome_in_spread(m, home_sym, away_sym)
                if home_leg:
                    out["home_spread"][strike] = home_leg
                if away_leg:
                    out["away_spread"][-strike] = away_leg
            elif typ == "TOTAL":
                over_leg, under_leg = _find_outcome_in_total(m)
                if over_leg:
                    out["over"][strike] = over_leg
                if under_leg:
                    out["under"][strike] = under_leg
            else:  # MONEY
                home_ml, away_ml = _find_outcome_in_spread(m, home_sym, away_sym)
                if home_ml:
                    out["home_ml"] = home_ml
                if away_ml:
                    out["away_ml"] = away_ml
        except (AttributeError, TypeError, ValueError):
            continue
    return out
