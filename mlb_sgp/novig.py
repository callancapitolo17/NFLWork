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

from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone

from mlb_sgp._shared import PricedRow, TargetLine, decimal_to_american
from mlb_sgp.novig_client import NovigClient


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
        client = NovigClient(verbose=verbose)

    # Import legacy helpers lazily — keeps the period-filter early-exit
    # test free of any HTTP imports.
    from scraper_novig_sgp import (
        match_events,
        fetch_event_legs,
        try_integer_fallback_nv,
        submit_parlay,
    )

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
    matched = match_events(nv_events, parlay_lines)
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
    from scraper_novig_sgp import SPREAD_TYPE, TOTAL_TYPE

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

        # Build offered (spread, total) sets per period from the raw
        # market list. Each spread/total market carries one `strike`
        # value; the home-perspective spread set is the union of strikes
        # across all SPREAD-type markets for the period.
        offered_per_period: dict[str, dict[str, set]] = {}
        for period_key in ("fg", "f5"):
            sp_type = SPREAD_TYPE[period_key]
            to_type = TOTAL_TYPE[period_key]
            spreads = {
                float(m["strike"]) for m in markets
                if m.get("type") == sp_type and m.get("strike") is not None
            }
            totals = {
                float(m["strike"]) for m in markets
                if m.get("type") == to_type and m.get("strike") is not None
            }
            offered_per_period[period_key] = {
                "spreads": spreads,
                "totals": totals,
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
        if verbose:
            print(f"  game {game_id}: {pre} kalshi → {post} offered", flush=True)

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
            client.session, combos, submit_parlay, verbose,
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
                if verbose:
                    print(f"  nv target error: {e}", flush=True)

    return out


def _price_combos_parallel(
    session,
    combos: tuple,
    submit_parlay_fn,
    verbose: bool,
) -> dict:
    """Price the 4 combo flavors in parallel and return {combo_name: priced}.

    Each combo submits a parlay RFQ via ``submit_parlay`` on the shared
    curl_cffi session. The SANITY_MULT_RATIO filter runs *inside* the
    worker so a combo that fails sanity is omitted from the result dict
    — same behavior as the sequential path, just executed concurrently.
    """
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
                if verbose:
                    print(
                        f"      SANITY-DROP {combo_name} "
                        f"parlay={priced['decimal']:.2f}"
                    )
                return combo_name, None

            return combo_name, priced
        except Exception:
            return combo_name, None

    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = [
            pool.submit(_price_one, combo_name, sp, to)
            for combo_name, sp, to in combos
        ]
        for fut in as_completed(futures):
            try:
                combo_name, priced = fut.result()
            except Exception:
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
