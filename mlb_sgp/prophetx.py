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

import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone

from mlb_sgp._shared import PricedRow, TargetLine, american_to_decimal, decimal_to_american
from mlb_sgp.prophetx_client import ProphetXClient, SelectionLeg


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
    """
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
    matched = match_events(px_events, parlay_lines)
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
        if verbose:
            print(f"  game {game_id}: {pre} kalshi → {post} offered", flush=True)

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
            combos, verbose,
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
                if verbose:
                    print(f"  px target error: {e}", flush=True)

    # ----- Phase 3: moneyline × total combos (FG) ----- #
    out.extend(_price_ml_total_for_games(
        client, markets_cache, matched_by_gid, filtered_targets,
        fetch_now, parallelism, verbose,
    ))

    return out


def _price_ml_total_for_games(
    client, markets_cache, matched_by_gid, filtered_targets,
    fetch_now, parallelism, verbose,
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
                client, event_id, ml_mid, total_mid, ml_combos, verbose)
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
                if verbose:
                    print(f"  px ML target error: {e}", flush=True)
    return out


def _price_combos_parallel(
    client: ProphetXClient,
    px_event_id: str,
    spread_mid,
    total_mid,
    combos: tuple,
    verbose: bool,
) -> dict:
    """Price the 4 combo flavors in parallel and return {combo_name: sgp_decimal}.

    Each combo submits a parlay RFQ via ``client.submit_parlay_rfq`` on
    the shared curl_cffi session. The SANITY_MULT_RATIO filter (F5-Over
    defense) runs *inside* the worker so a combo that fails sanity is
    omitted from the result dict — same behavior as the sequential path,
    just executed concurrently.
    """
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
                if verbose:
                    print(
                        f"      SANITY-DROP {combo_name} "
                        f"parlay={sgp_decimal:.2f} naive={leg1_dec * leg2_dec:.2f}"
                    )
                return combo_name, None

            return combo_name, sgp_decimal
        except Exception:
            return combo_name, None

    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = [
            pool.submit(_price_one, combo_name, sp_sel, tot_sel)
            for combo_name, sp_sel, tot_sel in combos
        ]
        for fut in as_completed(futures):
            try:
                combo_name, sgp_decimal = fut.result()
            except Exception:
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
