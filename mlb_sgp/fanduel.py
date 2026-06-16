"""FanDuel SGP orchestration library.

Pure function: `price_sgps(target_lines)` -> `list[PricedRow]`.

Consumes `FanDuelClient` (HTTP transport) and calls the existing helper
functions in `scraper_fanduel_sgp.py` (fetch_fd_events, match_events,
fetch_event_runners, price_combo, try_integer_fallback_fd) to compose
the four-combo SGP pricing per (game, period, spread_line, total_line).

No DB I/O. Caller (scraper shim or bot sgp_runner) handles persistence.

Source-label tracking
---------------------
Rows produced by the main API path are tagged ``fanduel_direct``.
Rows produced by the integer-fallback path (off-main integer totals,
derived from adjacent half-point alts via implyBets) are tagged
``fanduel_interpolated``. The dashboard parlay-tab reader filters on
``source IN (..._direct)``, so preserving these labels keeps existing
dashboard behavior byte-identical post-refactor.

Helper-signature deviations from the original plan
--------------------------------------------------
The plan-spec's example structure didn't match FD's actual helper
signatures. This orchestrator adapts:

* ``fetch_event_runners(session, fd_event_id, fd_home, fd_away)`` —
  the legacy SGP helper takes the team-name strings (used to label
  spread runners as home vs away) and returns the nested SGP-filtered
  shape ``{"fg": {"spreads": {(side, line): (mid, sid)}, "totals":
  {("O"|"U", line): (mid, sid)}}, "f5": {...}}``. This is the
  scraper helper, NOT ``FanDuelClient.fetch_event_runners`` (which
  returns a flat Runner list for the singles scraper). The two share
  a name but emit different shapes.

* ``price_combo(session, spread_market, spread_sel, total_market,
  total_sel, verbose=False)`` returns ``{"decimal": float, "american":
  int}`` or ``None``. Unlike DK's ``calculate_sgp``, FD's helper
  already returns finished decimal+american odds (no extra conversion
  needed beyond rounding).

* ``try_integer_fallback_fd(session, sel, home_line, away_line,
  total_line, verbose=False)`` takes the per-period ``sel`` dict
  (``{"spreads": ..., "totals": ...}``) and returns ``{"fair_probs":
  {"home_over", "home_under", "away_over", "away_under"}}`` or
  ``None``. Decimals are recovered as ``1 / fair_prob``.

* No canonical-market filter (unlike DK). FD doesn't reject
  cross-market combos at the implyBets endpoint, so the orchestrator
  only needs the SGP-tab selections and can price any
  (spread_sel, total_sel) pair directly.

* No ``fetch_main_market_nums`` step. FD's ``fetch_event_runners``
  already labels markets as main vs alt via ``_MARKET_MAP`` internally
  and prefers main when both exist for the same (side, line) key.
"""
from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone

from mlb_sgp._shared import PricedRow, TargetLine, decimal_to_american
from mlb_sgp.fd_client import FanDuelClient


BOOK_NAME = "fanduel"
SOURCE_LABEL = "fanduel_direct"
SOURCE_LABEL_FALLBACK = "fanduel_interpolated"

# Combo names — byte-identical to scraper_fanduel_sgp.py so the
# dashboard / kalshi_mlb_rfq leg lookups keep matching. The "F5 "
# prefix is applied to F5-period rows in the orchestrator.
_FALLBACK_KEY_TO_COMBO = {
    "home_over":  "Home Spread + Over",
    "home_under": "Home Spread + Under",
    "away_over":  "Away Spread + Over",
    "away_under": "Away Spread + Under",
}


# Target-level parallelism. FD is not the cycle's long pole (strict
# line matching filters most targets), so the default is conservative;
# env-overridable for ops tuning.
FD_TARGET_PARALLELISM_DEFAULT = 4


def _resolve_parallelism(parallelism: int | None) -> int:
    if parallelism is not None:
        return parallelism
    return int(os.environ.get("MLB_SGP_FD_PARALLELISM",
                              str(FD_TARGET_PARALLELISM_DEFAULT)))


def _extract_offered_lines_fd(
    sel_ids_per_period: dict, period_key: str,
) -> dict[str, set]:
    """Return ``{"spreads": set, "totals": set}`` for one period from FD
    ``fetch_event_runners`` output.

    FD keys spreads as ``(side, signed_line)`` where ``side`` is
    ``"home"`` or ``"away"``. Each magnitude appears twice (once per
    side, with opposite signs), so we collapse to one signed line per
    market by restricting to ``side == "home"`` — that's already the
    home-perspective signed line.

    Totals are keyed ``(over_under, line)`` — we collect ``line``.

    Used by Filter A in ``price_sgps`` to drop targets the book doesn't
    offer before the pricing loop runs.
    """
    sel = sel_ids_per_period.get(period_key, {"spreads": {}, "totals": {}})
    spreads = {
        line for (side, line) in sel.get("spreads", {}).keys()
        if side == "home"
    }
    totals = {line for (ou, line) in sel.get("totals", {}).keys()}
    return {"spreads": spreads, "totals": totals}


def price_sgps(
    target_lines: list[TargetLine],
    periods: tuple[str, ...] = ("FG",),
    client: FanDuelClient | None = None,
    verbose: bool = False,
    parallelism: int | None = None,
    fetchers: dict | None = None,
) -> list[PricedRow]:
    """Price every target line against the FanDuel SGP API.

    Parameters
    ----------
    target_lines
        List of ``TargetLine`` records — one per (game, period, spread,
        total) tuple the caller wants priced.
    periods
        Which periods to actually price. ``("FG",)`` (bot default) or
        ``("FG", "F5")`` (dashboard). Targets in other periods are
        silently dropped.
    client
        Optional pre-built ``FanDuelClient``. If ``None``, a new one is
        created (which opens a real curl_cffi Chrome-TLS session —
        useful for production, but tests should always pass a mock).
    verbose
        Forwarded to leg-level helpers for debug printing.
    parallelism
        Number of targets priced concurrently. ``None`` resolves to the
        ``MLB_SGP_FD_PARALLELISM`` env var or ``FD_TARGET_PARALLELISM_DEFAULT``;
        total in-flight FD requests is ``parallelism × 4`` (4 combos/target).
    fetchers
        Structure-fetch override hooks; default ``None`` uses the legacy
        direct fetches (dashboard path); SGPService injects TTL-cached
        versions; prices are never affected.

    Returns
    -------
    list[PricedRow]
        Up to four PricedRow per (game, period) — one per combo. Rows
        produced via the integer-fallback path are tagged with
        ``source = fanduel_interpolated`` instead of
        ``fanduel_direct``.
    """
    if not target_lines:
        return []

    # Filter to requested periods early; if nothing remains, bail before
    # touching the network / building a client.
    targets = [t for t in target_lines if t.period in periods]
    if not targets:
        return []

    # Build the client lazily so tests that filter everything out via
    # the `periods` argument don't have to mock HTTP at all.
    if client is None:
        client = FanDuelClient(verbose=verbose)

    # Import legacy helpers lazily — keeps the period-filter early-exit
    # test free of any HTTP imports.
    from scraper_fanduel_sgp import (
        fetch_fd_events,
        match_events,
        fetch_event_runners,
        price_combo,
        try_integer_fallback_fd,
    )

    _f = fetchers or {}
    fetch_events_fn = _f.get("fetch_fd_events", fetch_fd_events)
    fetch_runners_fn = _f.get("fetch_event_runners", fetch_event_runners)

    # ----- Group target lines by game ----- #
    # match_events expects the legacy parlay-lines dict shape: one
    # entry per game_id with fg_/f5_ columns. Collapse per-period
    # TargetLine rows back into that shape, leaving the other period's
    # columns as None.
    target_dict: dict[str, dict] = {}
    for t in targets:
        ent = target_dict.setdefault(t.game_id, {
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

    # ----- Phase 1: list FD events and match to our game_ids ----- #
    fd_events = fetch_events_fn(client.session)
    matched = match_events(fd_events, target_dict)
    if not matched:
        return []

    # Index matched games by game_id for fast lookup during the pricing
    # loop. match_events guarantees one entry per game_id (the legacy
    # scraper post-dedupes; we trust match_events to behave identically).
    matched_by_gid = {m["game_id"]: m for m in matched}

    out: list[PricedRow] = []
    fetch_now = datetime.now(timezone.utc)

    # Cache fetched runners per game_id so two TargetLines on the same
    # game (e.g. FG + F5) only hit FD's event-page endpoint once.
    runners_cache: dict[str, dict] = {}

    # ----- Phase 1.6: group targets by game (Filter A) ----- #
    # We need to know what lines FD offers per game before we can drop
    # targets the book doesn't carry. Fetch runners up-front per matched
    # game, then build a per-period set of offered (spread, total) lines.
    targets_by_game: dict[str, list[TargetLine]] = {}
    for t in targets:
        targets_by_game.setdefault(t.game_id, []).append(t)

    filtered_targets: list[TargetLine] = []
    for game_id, game_targets in targets_by_game.items():
        game = matched_by_gid.get(game_id)
        if game is None:
            continue
        if game_id not in runners_cache:
            runners_cache[game_id] = fetch_runners_fn(
                client.session,
                game["fd_event_id"],
                game["fd_home"],
                game["fd_away"],
            )
        sel_ids_per_period = runners_cache[game_id]
        offered_per_period: dict[str, dict[str, set]] = {
            period_key: _extract_offered_lines_fd(sel_ids_per_period, period_key)
            for period_key in ("fg", "f5")
        }
        # Filter this game's targets.
        pre = 0
        post = 0
        for t in game_targets:
            pre += 1
            offered = offered_per_period.get(t.period.lower())
            if offered is None:
                continue
            if t.spread not in offered["spreads"]:
                continue
            if t.total in offered["totals"]:
                filtered_targets.append(t)
                post += 1
                continue
            # Fallback: integer total + adjacent half-point alts both offered.
            if float(t.total).is_integer() and (
                (t.total - 0.5) in offered["totals"]
                and (t.total + 0.5) in offered["totals"]
            ):
                filtered_targets.append(t)
                post += 1
        if verbose:
            print(f"  game {game_id}: {pre} kalshi → {post} offered", flush=True)

    # ----- Phase 2: per target row, price 4 combos ----- #
    # Targets fan out on a thread pool (mirrors draftkings.py /
    # novig.py). Runners are pre-fetched in the filter loop above, so
    # workers only do implyBets calls — no shared-cache mutation races.
    def _price_one_target(t: TargetLine) -> list[PricedRow]:
        target_rows: list[PricedRow] = []
        game = matched_by_gid.get(t.game_id)
        if game is None:
            return target_rows

        period_key = t.period.lower()  # "FG" -> "fg"

        # Runners cache is guaranteed-populated by the filter loop
        # above (we only added a target if its game's runners were
        # fetched successfully).
        sel_ids_per_period = runners_cache[t.game_id]
        sel = sel_ids_per_period.get(period_key, {"spreads": {}, "totals": {}})

        # Sign convention: TargetLine.spread is the home-perspective
        # signed line, so the away leg is keyed by -spread.
        if not sel.get("spreads"):
            return target_rows
        home_line = t.spread
        away_line = -t.spread

        home_spread = sel["spreads"].get(("home", home_line))
        away_spread = sel["spreads"].get(("away", away_line))
        over = sel["totals"].get(("O", t.total))
        under = sel["totals"].get(("U", t.total))

        # If any of the four legs is missing, the requested (spread,
        # total) is off-main / off-alt. Try the integer-line fallback.
        if not (home_spread and away_spread and over and under):
            fallback = try_integer_fallback_fd(
                client.session, sel, home_line, away_line, t.total,
                verbose=verbose,
            )
            if fallback is None:
                return target_rows
            prefix = "" if t.period == "FG" else "F5 "
            fair_probs = fallback["fair_probs"]
            for key, base_combo in _FALLBACK_KEY_TO_COMBO.items():
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

        # ----- Main path: price the 4 canonical combos in parallel ----- #
        prefix = "" if t.period == "FG" else "F5 "
        combos = (
            ("Home Spread + Over",  home_spread, over),
            ("Home Spread + Under", home_spread, under),
            ("Away Spread + Over",  away_spread, over),
            ("Away Spread + Under", away_spread, under),
        )

        priced_by_combo = _price_combos_parallel(
            client.session, combos, price_combo, verbose,
        )

        for combo_name, _sp_pair, _tot_pair in combos:
            result = priced_by_combo.get(combo_name)
            if not result:
                continue
            dec = float(result["decimal"])
            am = int(result["american"])
            target_rows.append(PricedRow(
                game_id=t.game_id,
                combo=prefix + combo_name,
                period=t.period,
                spread_line=t.spread,
                total_line=t.total,
                bookmaker=BOOK_NAME,
                source=SOURCE_LABEL,
                sgp_decimal=round(dec, 4),
                sgp_american=am,
                fetch_time=fetch_now,
            ))
        return target_rows

    n_workers = max(1, _resolve_parallelism(parallelism))
    with ThreadPoolExecutor(max_workers=n_workers) as pool:
        futures = [pool.submit(_price_one_target, t) for t in filtered_targets]
        for f in as_completed(futures):
            try:
                out.extend(f.result())
            except Exception as e:
                if verbose:
                    print(f"  fd target error: {e}", flush=True)

    return out


def _price_combos_parallel(
    session,
    combos: tuple,
    price_combo_fn,
    verbose: bool,
) -> dict:
    """Price the 4 combo flavors in parallel and return {combo_name: result}.

    Each combo runs a single implyBets call (FD's pricer is leg-pair
    direct, no canonical filter needed). 4-way parallelism on the same
    curl_cffi session — safe because cffi sessions are thread-safe for
    separate requests and we never share request state between threads.

    A combo that fails to price (None / falsy result) is omitted from
    the result dict, mirroring the sequential path's ``if not result:
    continue`` behavior.
    """
    results: dict = {}

    def _price_one(combo_name, sp_pair, tot_pair):
        try:
            return combo_name, price_combo_fn(
                session,
                sp_pair[0], sp_pair[1],
                tot_pair[0], tot_pair[1],
                verbose=verbose,
            )
        except Exception:
            return combo_name, None

    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = [
            pool.submit(_price_one, combo_name, sp_pair, tot_pair)
            for combo_name, sp_pair, tot_pair in combos
        ]
        for fut in as_completed(futures):
            try:
                combo_name, result = fut.result()
            except Exception:
                continue
            if result:
                results[combo_name] = result
    return results
