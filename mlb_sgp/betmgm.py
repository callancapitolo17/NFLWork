"""BetMGM SGP orchestration library.

Pure function: ``price_sgps(target_lines)`` -> ``list[PricedRow]``.

Consumes ``BetMGMClient`` (HTTP transport) and prices each (game, period,
spread, total) target's four spread×total corners against BetMGM's Angstrom
bet-builder endpoint (``/cds-api/bettingoffer/picks``), which returns TRUE
correlated odds (not the naive leg product).

No DB I/O. Caller (scraper shim or bot sgp_runner) handles persistence.

Pattern mirrors ``novig.py``:
  * one ``Event`` list call,
  * canonical event matching (teams + UTC-hour bucket for doubleheaders),
  * Filter A — fetch each matched game's full market tree once, build the
    offered (spread, total) sets, and drop targets BetMGM doesn't offer
    BEFORE the pricing loop (BetMGM exposes a full alt grid, so without
    this filter we'd fire thousands of futile combo POSTs per cycle),
  * target-level ThreadPoolExecutor over the 4-combo pricing.

Source label: ``betmgm_direct``. BetMGM offers only half-point lines and the
bot's target lines are all half-points, so there is no integer-line fallback
path (unlike DK/FD/NV) — an off-grid target simply yields no row (honest
"no price" rather than an approximated one).

Sanity filter: mirrors Novig/ProphetX — drop any combo whose correlated
decimal exceeds ``SANITY_MULT_RATIO`` (1.5x) the naive independent leg
product. Legitimate anti-correlated corners top out ~1.15x.
"""
from __future__ import annotations

import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone

from mlb_sgp._shared import PricedRow, TargetLine, decimal_to_american
from mlb_sgp.betmgm_client import BetMGMClient, Event

BOOK_NAME = "betmgm"
SOURCE_LABEL = "betmgm_direct"

SANITY_MULT_RATIO = 1.5

# BetMGM tolerates moderate concurrency on the CDS endpoints. 4 targets ×
# 4 combos = 16 in-flight POSTs per book — same envelope as Novig.
MGM_TARGET_PARALLELISM = 4

# Market-name → (period, kind). Matched case-insensitively by substring so
# minor casing variants ("1st 5 innings - ...") still resolve.
_FG_SPREAD = "run line spread"
_FG_TOTAL = "totals"
_F5_SPREAD = "1st 5 innings - run line spread"
_F5_TOTAL = "1st 5 innings - totals"
# Moneyline (for ML×total combos). BetMGM names it "Money Line" (FG) and
# "1st 5 innings - Money Line" (F5); both are isBetBuilder-eligible.
_FG_ML = "money line"
_F5_ML = "1st 5 innings - money line"


def _to_float(s) -> float | None:
    """Parse a BetMGM number string, which uses a comma decimal ('+1,5')."""
    if s is None:
        return None
    try:
        return float(str(s).replace(",", ".").replace("+", "").strip())
    except (TypeError, ValueError):
        return None


def _classify_market(name: str) -> tuple[str, str] | None:
    """Return (period, kind) for a spread/total market name, else None.

    kind is 'spread' or 'total'; period is 'FG' or 'F5'. F5 is checked
    first because its names contain the FG substrings.
    """
    low = name.lower()
    if _F5_SPREAD in low:
        return ("F5", "spread")
    if _F5_TOTAL in low:
        return ("F5", "total")
    if _F5_ML in low:
        return ("F5", "moneyline")
    # Guard against other-period markets (1st 3 / 1st 7 innings, 1X2, etc.)
    if "innings" in low and "1st 5" not in low:
        return None
    if low == _FG_SPREAD:
        return ("FG", "spread")
    if low == _FG_TOTAL:
        return ("FG", "total")
    if low == _FG_ML:
        return ("FG", "moneyline")
    return None


def parse_markets(markets: list[dict], home_team: str, away_team: str) -> dict:
    """Index BetMGM's raw optionMarkets into a pricing-ready structure.

    Returns::

        {
          "FG": {
            "spreads": {home_signed_line: {"home": leg, "away": leg}},
            "totals":  {line: {"over": leg, "under": leg}},
          },
          "F5": { ... },
        }

    where ``leg`` is ``(market_id, option_id, decimal_odds)``. Spread lines
    are keyed by the HOME-perspective signed line (matching TargetLine.spread).
    Only markets flagged SGP-eligible (``isBetBuilder`` not False) contribute.
    """
    out = {
        "FG": {"spreads": {}, "totals": {}, "moneyline": None},
        "F5": {"spreads": {}, "totals": {}, "moneyline": None},
    }
    home_low = home_team.lower()
    away_low = away_team.lower()

    for om in markets:
        nm = (om.get("name") or {}).get("value", "") if isinstance(om.get("name"), dict) else ""
        cls = _classify_market(nm)
        if cls is None:
            continue
        # SGP eligibility: skip markets explicitly not bet-builder eligible.
        if om.get("isBetBuilder") is False:
            continue
        period, kind = cls
        mid = om.get("id")
        options = om.get("options") or []

        if kind == "spread":
            home_leg = away_leg = None
            home_line = None
            for o in options:
                onm = ((o.get("name") or {}).get("value", "")
                       if isinstance(o.get("name"), dict) else "")
                line = _to_float(o.get("attr"))
                odds = ((o.get("price") or {}).get("odds"))
                oid = o.get("id")
                if line is None or odds is None or oid is None:
                    continue
                onm_low = onm.lower()
                if home_low and home_low in onm_low:
                    home_leg = (mid, oid, float(odds))
                    home_line = line
                elif away_low and away_low in onm_low:
                    away_leg = (mid, oid, float(odds))
            if home_leg and away_leg and home_line is not None:
                out[period]["spreads"][home_line] = {"home": home_leg, "away": away_leg}

        elif kind == "total":
            over_leg = under_leg = None
            line_val = None
            for o in options:
                onm = ((o.get("name") or {}).get("value", "")
                       if isinstance(o.get("name"), dict) else "")
                odds = ((o.get("price") or {}).get("odds"))
                oid = o.get("id")
                if odds is None or oid is None:
                    continue
                m = re.search(r"(\d+[.,]\d+)", onm)
                ln = _to_float(m.group(1)) if m else None
                if ln is None:
                    continue
                line_val = ln
                if onm.lower().startswith("over"):
                    over_leg = (mid, oid, float(odds))
                elif onm.lower().startswith("under"):
                    under_leg = (mid, oid, float(odds))
            if over_leg and under_leg and line_val is not None:
                out[period]["totals"][line_val] = {"over": over_leg, "under": under_leg}

        elif kind == "moneyline":
            # ML option names use SHORT team names ("Giants") while home_team is
            # the full canonical name ("San Francisco Giants"), so match
            # bidirectionally (either string contained in the other). Both
            # options carry decimal `odds`, so the naive-multiply sanity check
            # downstream works unchanged.
            home_leg = away_leg = None
            for o in options:
                onm = ((o.get("name") or {}).get("value", "")
                       if isinstance(o.get("name"), dict) else "")
                odds = ((o.get("price") or {}).get("odds"))
                oid = o.get("id")
                if odds is None or oid is None:
                    continue
                onm_low = onm.lower()
                if home_low and onm_low and (home_low in onm_low or onm_low in home_low):
                    home_leg = (mid, oid, float(odds))
                elif away_low and onm_low and (away_low in onm_low or onm_low in away_low):
                    away_leg = (mid, oid, float(odds))
            if home_leg and away_leg:
                out[period]["moneyline"] = {"home": home_leg, "away": away_leg}

    return out


def _match_events(events: list[Event], targets: list[TargetLine]) -> dict[str, Event]:
    """Map canonical game_id -> BetMGM Event via team-name resolution + UTC
    hour bucket (doubleheader-safe). Imported canonical helpers lazily so the
    early period-filter exit needs no Answer Keys import.

    BetMGM lists both today's and tomorrow's copy of a matchup, so matching
    is bucket-strict: when both the target and event carry a UTC-hour bucket
    they must agree, and a different-day duplicate can never overwrite an
    exact-bucket match. The any-bucket fallback only fires when a side lacks
    a usable timestamp, and never clobbers an exact match.
    """
    from canonical_match import load_team_dict, load_canonical_games, resolve_team_names
    from mlb_sgp._shared import _utc_bucket

    # Target games: game_id -> (canonical home, away, bucket).
    tgames: dict[str, tuple[str, str, str]] = {
        t.game_id: (t.home_team, t.away_team, _utc_bucket(t.commence_time))
        for t in targets
    }

    team_dict = load_team_dict("mlb")
    canonical_games = load_canonical_games("mlb")

    matched: dict[str, Event] = {}
    exact: set[str] = set()
    for ev in events:
        resolved = resolve_team_names(ev.away_team, ev.home_team, team_dict, canonical_games)
        if not resolved or not resolved[0] or not resolved[1]:
            continue
        canon_away, canon_home = resolved
        ev_bucket = _utc_bucket(ev.start_time)
        for gid, (h, a, tb) in tgames.items():
            if h != canon_home or a != canon_away:
                continue
            if tb and ev_bucket:
                if tb == ev_bucket:
                    matched[gid] = ev
                    exact.add(gid)
                # both buckets present but differ → not this event
            elif gid not in exact:
                # one side lacks a timestamp — accept, but never clobber exact
                matched[gid] = ev
    return matched


def price_sgps(
    target_lines: list[TargetLine],
    periods: tuple[str, ...] = ("FG",),
    client: BetMGMClient | None = None,
    verbose: bool = False,
) -> list[PricedRow]:
    """Price every target line against BetMGM's Angstrom bet-builder endpoint.

    Up to four PricedRow per (game, period) — one per spread×total corner.
    Combos BetMGM won't price, or that fail the sanity filter, are dropped.
    """
    if not target_lines:
        return []
    targets = [t for t in target_lines if t.period in periods]
    if not targets:
        return []

    if client is None:
        client = BetMGMClient(verbose=verbose)
    if not client.accessid():
        if verbose:
            print("  [mgm] no accessid — aborting (wrong state / harvest failed)")
        return []

    events = client.list_events()
    if not events:
        return []
    matched_by_gid = _match_events(events, targets)
    if not matched_by_gid:
        return []

    # ----- Filter A: fetch each matched game's market tree once, keep only
    # targets BetMGM actually offers (both spread AND total present). ----- #
    targets_by_game: dict[str, list[TargetLine]] = {}
    for t in targets:
        if t.game_id in matched_by_gid:
            targets_by_game.setdefault(t.game_id, []).append(t)

    parsed_cache: dict[str, dict] = {}
    filtered: list[TargetLine] = []
    for gid, game_targets in targets_by_game.items():
        ev = matched_by_gid[gid]
        markets = client.fetch_markets(ev.event_id)
        if not markets:
            continue
        parsed = parse_markets(markets, ev.home_team, ev.away_team)
        parsed_cache[gid] = parsed
        kept = 0
        for t in game_targets:
            per = parsed.get(t.period)
            if not per:
                continue
            if float(t.spread) in per["spreads"] and float(t.total) in per["totals"]:
                filtered.append(t)
                kept += 1
        if verbose:
            print(f"  [mgm] game {gid[:8]}: {len(game_targets)} targets → {kept} offered")

    if not filtered:
        return []

    fetch_now = datetime.now(timezone.utc)

    def _price_one_target(t: TargetLine) -> list[PricedRow]:
        ev = matched_by_gid.get(t.game_id)
        parsed = parsed_cache.get(t.game_id)
        if ev is None or parsed is None:
            return []
        per = parsed[t.period]
        spread_legs = per["spreads"].get(float(t.spread))
        total_legs = per["totals"].get(float(t.total))
        if not spread_legs or not total_legs:
            return []

        home, away = spread_legs["home"], spread_legs["away"]
        over, under = total_legs["over"], total_legs["under"]
        prefix = "" if t.period == "FG" else "F5 "
        combos = (
            ("Home Spread + Over", home, over),
            ("Home Spread + Under", home, under),
            ("Away Spread + Over", away, over),
            ("Away Spread + Under", away, under),
        )

        def _price_combo(combo_name, sp_leg, to_leg):
            priced = client.price_picks(
                ev.event_id, [(sp_leg[0], sp_leg[1]), (to_leg[0], to_leg[1])]
            )
            if priced is None:
                return combo_name, None
            dec = priced["decimal"]
            naive = sp_leg[2] * to_leg[2]
            if naive > 0 and dec > naive * SANITY_MULT_RATIO:
                if verbose:
                    print(f"      [mgm] SANITY-DROP {combo_name} dec={dec:.2f} naive={naive:.2f}")
                return combo_name, None
            return combo_name, dec

        rows: list[PricedRow] = []
        with ThreadPoolExecutor(max_workers=4) as pool:
            futs = [pool.submit(_price_combo, c, sp, to) for c, sp, to in combos]
            for f in as_completed(futs):
                try:
                    combo_name, dec = f.result()
                except Exception:
                    continue
                if dec is None:
                    continue
                rows.append(PricedRow(
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
        return rows

    out: list[PricedRow] = []
    with ThreadPoolExecutor(max_workers=MGM_TARGET_PARALLELISM) as pool:
        futures = [pool.submit(_price_one_target, t) for t in filtered]
        for f in as_completed(futures):
            try:
                out.extend(f.result())
            except Exception as e:
                if verbose:
                    print(f"  [mgm] target error: {e}")

    # ----- Phase 3: moneyline × total combos (FG only) ----- #
    out.extend(_price_ml_total_for_games(
        client, matched_by_gid, parsed_cache, filtered, fetch_now, verbose))
    return out


def _price_ml_total_for_games(client, matched_by_gid, parsed_cache, filtered,
                              fetch_now, verbose) -> list[PricedRow]:
    """Price the 4 FG ML×total combos per (game, distinct FG total) — ONCE.

    ML_TOTAL_FAMILY is FG-only. Priced per distinct (game, total) rather than
    per spread target so a game with several alt-spread targets at the same
    total doesn't re-issue identical ML price_picks calls (footprint + rate
    limits). spread_line=None marks "no spread leg".
    """
    # Distinct FG totals to price per game (from the targets BetMGM offered).
    totals_by_game: dict[str, set] = {}
    for t in filtered:
        if t.period == "FG":
            totals_by_game.setdefault(t.game_id, set()).add(float(t.total))

    def _price_one_game(gid: str, totals: set) -> list[PricedRow]:
        ev = matched_by_gid.get(gid)
        parsed = parsed_cache.get(gid)
        if ev is None or parsed is None:
            return []
        ml_legs = parsed["FG"].get("moneyline")
        if not ml_legs:
            return []
        ml_home, ml_away = ml_legs["home"], ml_legs["away"]
        rows: list[PricedRow] = []
        jobs = []  # (combo_name, total, ml_leg, to_leg)
        for total in totals:
            tl = parsed["FG"]["totals"].get(total)
            if not tl:
                continue
            over, under = tl["over"], tl["under"]
            jobs += [
                ("Home ML + Over", total, ml_home, over),
                ("Home ML + Under", total, ml_home, under),
                ("Away ML + Over", total, ml_away, over),
                ("Away ML + Under", total, ml_away, under),
            ]
        if not jobs:
            return rows

        def _price(combo_name, total, ml_leg, to_leg):
            priced = client.price_picks(
                ev.event_id, [(ml_leg[0], ml_leg[1]), (to_leg[0], to_leg[1])])
            if priced is None:
                return None
            dec = priced["decimal"]
            naive = ml_leg[2] * to_leg[2]
            if naive > 0 and dec > naive * SANITY_MULT_RATIO:
                if verbose:
                    print(f"      [mgm] SANITY-DROP {combo_name} dec={dec:.2f} naive={naive:.2f}")
                return None
            return (combo_name, total, dec)

        with ThreadPoolExecutor(max_workers=4) as pool:
            futs = [pool.submit(_price, c, tot, ml, to) for c, tot, ml, to in jobs]
            for f in as_completed(futs):
                try:
                    res = f.result()
                except Exception:
                    continue
                if res is None:
                    continue
                combo_name, total, dec = res
                rows.append(PricedRow(
                    game_id=gid, combo=combo_name, period="FG",
                    spread_line=None, total_line=total,
                    bookmaker=BOOK_NAME, source=SOURCE_LABEL,
                    sgp_decimal=round(dec, 4),
                    sgp_american=decimal_to_american(dec),
                    fetch_time=fetch_now))
        return rows

    out: list[PricedRow] = []
    if not totals_by_game:
        return out
    with ThreadPoolExecutor(max_workers=MGM_TARGET_PARALLELISM) as pool:
        futs = [pool.submit(_price_one_game, gid, tots)
                for gid, tots in totals_by_game.items()]
        for f in as_completed(futs):
            try:
                out.extend(f.result())
            except Exception as e:
                if verbose:
                    print(f"  [mgm] ML target error: {e}")
    return out
