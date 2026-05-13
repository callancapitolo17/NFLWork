"""DraftKings SGP orchestration library.

Pure function: `price_sgps(target_lines)` → `list[PricedRow]`.

Consumes `DraftKingsClient` (HTTP transport) and calls the existing
helper functions in `scraper_draftkings_sgp.py` (match_events,
fetch_main_market_nums, fetch_selection_ids, calculate_sgp,
try_integer_fallback_dk) to compose the four-combo SGP pricing per
(game, period, spread_line, total_line).

No DB I/O. Caller (scraper shim or bot sgp_runner) handles persistence.

Source-label tracking
---------------------
Rows produced by the main API path are tagged ``draftkings_direct``.
Rows produced by the integer-fallback path (off-main total lines,
derived from adjacent half-point alts) are tagged
``draftkings_interpolated``. The dashboard parlay-tab reader filters
on ``source IN (..._direct)``, so preserving these labels keeps the
existing dashboard behavior byte-identical post-refactor.

Helper-signature deviations from the original plan
--------------------------------------------------
The plan-spec invoked the legacy helpers with different signatures
than what scraper_draftkings_sgp.py actually exposes. This orchestrator
adapts to the real signatures:

* ``fetch_selection_ids(session, dk_event_id, main_market_nums=None,
  verbose=False)`` returns ALL selection IDs for the event, grouped by
  period (``fg`` / ``f5``). It does NOT take a spread/total/period
  kwarg — the orchestrator extracts the relevant (sign, line) entries
  itself.

* ``calculate_sgp(session, spread_sel, total_sel, verbose=False)``
  is leg-level: it takes two selection-id strings and returns
  ``{"trueOdds": float, "displayOdds": str}`` or ``None``. The
  orchestrator iterates spread × total candidate selections (filtered
  to canonical markets only) until calculateBets returns a price.

* ``try_integer_fallback_dk(session, sel_ids, spread_line, total_line,
  canonical, verbose=False)`` takes the per-period ``sel_ids`` dict and
  the canonical-market set, and returns
  ``{"fair_probs": {"home_over", "home_under", "away_over", "away_under"}}``
  or ``None``. Decimals are recovered as ``1 / fair_prob``.

* ``match_events(dk_events, parlay_lines)`` expects the legacy
  parlay-lines dict shape (``{game_id: {home_team, away_team,
  commence_time, fg_spread_line, fg_total_line, f5_spread_line,
  f5_total_line}}``). The orchestrator builds that dict from the input
  ``TargetLine`` list, collapsing per-period rows back into one game-
  level entry per ``game_id``.
"""
from __future__ import annotations

from datetime import datetime, timezone

from mlb_sgp._shared import PricedRow, TargetLine, decimal_to_american
from mlb_sgp.dk_client import DraftKingsClient


BOOK_NAME = "draftkings"
SOURCE_LABEL = "draftkings_direct"
SOURCE_LABEL_FALLBACK = "draftkings_interpolated"

# Combo names (must stay byte-identical to scraper_draftkings_sgp.py so
# the dashboard / kalshi_mlb_rfq leg lookups keep matching). The "F5 "
# prefix is applied for F5-period rows.
_COMBO_NAMES = (
    "Home Spread + Over",
    "Home Spread + Under",
    "Away Spread + Over",
    "Away Spread + Under",
)

# Map from fair_probs keys (try_integer_fallback_dk output) to combo names
_FALLBACK_KEY_TO_COMBO = {
    "home_over":  "Home Spread + Over",
    "home_under": "Home Spread + Under",
    "away_over":  "Away Spread + Over",
    "away_under": "Away Spread + Under",
}


def price_sgps(
    target_lines: list[TargetLine],
    periods: tuple[str, ...] = ("FG",),
    client: DraftKingsClient | None = None,
    verbose: bool = False,
) -> list[PricedRow]:
    """Price every target line against the DraftKings SGP API.

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
        Optional pre-built ``DraftKingsClient``. If ``None``, a new one
        is created (which opens a real curl_cffi session — useful for
        production, but tests should always pass a mock).
    verbose
        Forwarded to the leg-level helpers for debug printing.

    Returns
    -------
    list[PricedRow]
        Up to four PricedRow per (game, period) — one per combo.
        Rows produced via the integer-fallback path are tagged with
        ``source = draftkings_interpolated`` instead of
        ``draftkings_direct``.
    """
    if not target_lines:
        return []

    # Filter to requested periods early; if nothing remains, bail before
    # touching the network / building a client.
    targets = [t for t in target_lines if t.period in periods]
    if not targets:
        return []

    # Build the client lazily so tests that filter everything out via the
    # `periods` argument don't have to mock HTTP at all.
    if client is None:
        client = DraftKingsClient(verbose=verbose)

    # Import legacy helpers lazily. Keeping the import inside the function
    # avoids a hard dependency at module-import time (useful for the
    # period-filter early-exit test which never touches HTTP).
    from scraper_draftkings_sgp import (
        _market_num,
        calculate_sgp,
        fetch_dk_events,
        fetch_main_market_nums,
        fetch_selection_ids,
        match_events,
        try_integer_fallback_dk,
    )

    # ----- Group target lines by game ----- #
    # match_events expects a dict shaped like the legacy mlb_parlay_lines
    # table: one entry per game_id with fg_/f5_ columns. We collapse the
    # per-period TargetLine rows back into that shape, leaving the other
    # period's columns as None.
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

    # ----- Phase 1: list DK events and match to our game_ids ----- #
    dk_events = fetch_dk_events(client.session)
    matched = match_events(dk_events, target_dict)
    if not matched:
        return []

    # Index matched games by game_id for fast lookup during the pricing
    # loop. match_events guarantees one entry per game_id.
    matched_by_gid = {m["game_id"]: m for m in matched}

    out: list[PricedRow] = []
    fetch_now = datetime.now(timezone.utc)

    # ----- Phase 2: per game, per period, build and price 4 combos ----- #
    # Iterate target rows directly (not by_game) so a single TargetLine
    # produces exactly one (period, spread, total) priced call. This
    # mirrors the bot use-case where each Kalshi MVE enumeration becomes
    # one target row.
    for t in targets:
        game = matched_by_gid.get(t.game_id)
        if game is None:
            continue

        period_key = t.period.lower()  # "FG" -> "fg" for legacy dict keys

        # Fetch market nums and per-period selection IDs for this game.
        # Note: in the legacy scraper these are cached per-game across
        # both periods. We re-fetch here for simplicity; the bot calls
        # this once per cycle (~minute cadence) so the extra hits are
        # negligible. If profiling shows this matters, we can memoize
        # by (game_id, dk_event_id).
        main_nums = fetch_main_market_nums(client.session, game["dk_event_id"])
        sel_ids_all = fetch_selection_ids(
            client.session, game["dk_event_id"], main_nums, verbose,
        )
        sel_ids = sel_ids_all.get(period_key, {"spreads": {}, "totals": {}})
        canonical = sel_ids.get("canonical", set())

        # Sign convention: home favored (spread < 0) means home is "N"
        # (negative-line side of the spread market), away is "P". Flip
        # when home is the underdog. ``spread`` itself is stored as
        # ``abs(t.spread)`` because DK selection IDs encode unsigned
        # magnitudes (the sign letter carries direction).
        if t.spread < 0:
            home_sign, away_sign = "N", "P"
        else:
            home_sign, away_sign = "P", "N"
        spread_abs = abs(t.spread)

        # _1 = home participant, _3 = away participant in DK sel_id suffix.
        home_spread_sels = sel_ids["spreads"].get((home_sign, spread_abs, "1")) or []
        away_spread_sels = sel_ids["spreads"].get((away_sign, spread_abs, "3")) or []
        over_sels = sel_ids["totals"].get(("O", t.total)) or []
        under_sels = sel_ids["totals"].get(("U", t.total)) or []

        # If any of the four legs is missing, the requested (spread,
        # total) is off-main. Fall back to integer-line derivation.
        if not (home_spread_sels and away_spread_sels and over_sels and under_sels):
            fallback = try_integer_fallback_dk(
                client.session, sel_ids, t.spread, t.total, canonical,
                verbose=verbose,
            )
            if fallback is None:
                continue
            # 4 derived combos at once, all tagged _interpolated
            prefix = "" if t.period == "FG" else "F5 "
            fair_probs = fallback["fair_probs"]
            for key, base_combo in _FALLBACK_KEY_TO_COMBO.items():
                fair_p = fair_probs.get(key)
                if not fair_p or fair_p <= 0:
                    continue
                dec = 1.0 / fair_p
                out.append(PricedRow(
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
            continue

        # ----- Main path: price the 4 canonical combos ----- #
        prefix = "" if t.period == "FG" else "F5 "
        combos = (
            ("Home Spread + Over",  home_spread_sels, over_sels),
            ("Home Spread + Under", home_spread_sels, under_sels),
            ("Away Spread + Over",  away_spread_sels, over_sels),
            ("Away Spread + Under", away_spread_sels, under_sels),
        )

        for combo_name, sp_sels, tot_sels in combos:
            sgp = _price_combo(
                client.session, sp_sels, tot_sels, canonical,
                calculate_sgp, _market_num, verbose,
            )
            if not sgp:
                continue
            dec = sgp["trueOdds"]
            out.append(PricedRow(
                game_id=t.game_id,
                combo=prefix + combo_name,
                period=t.period,
                spread_line=t.spread,
                total_line=t.total,
                bookmaker=BOOK_NAME,
                source=SOURCE_LABEL,
                sgp_decimal=round(dec, 4),
                sgp_american=decimal_to_american(dec),
                fetch_time=fetch_now,
            ))

    return out


def _price_combo(
    session,
    sp_sels: list[str],
    tot_sels: list[str],
    canonical: set,
    calculate_sgp_fn,
    market_num_fn,
    verbose: bool,
):
    """Try every canonical (spread_sel, total_sel) pair until one prices.

    DK rejects "non-canonical" market pairings (e.g. main spread paired
    with a secondary totals-only market) with either a 422 or a
    combinabilityRestrictions response. The legacy scraper filters both
    legs to canonical markets only before calling calculateBets — we
    mirror that here.
    """
    for sp in sp_sels:
        if market_num_fn(sp) not in canonical:
            continue
        for to in tot_sels:
            if market_num_fn(to) not in canonical:
                continue
            sgp = calculate_sgp_fn(session, sp, to, verbose=verbose)
            if sgp:
                return sgp
    return None
