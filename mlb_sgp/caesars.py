"""Caesars SGP orchestration library.

Pure function: ``price_sgps(target_lines)`` -> ``list[PricedRow]``.

Consumes ``CaesarsClient`` (token-broker + REST). For each (game, period,
spread, total) target it prices the four spread×total corners against the
Caesars ZeroFlucs engine via ``POST /sb/v2/bets/details`` (correlated price,
not naive product). No DB I/O — caller persists.

Pattern mirrors ``betmgm.py`` / ``novig.py``:
  * one Event list call (tabs feed),
  * canonical event matching (teams + UTC-hour bucket, doubleheader-safe),
  * Filter A — fetch each matched game's full event tree once (alts live in
    ``event.keyMarketGroups[].markets[]``), keep only targets Caesars offers,
  * target-level ThreadPoolExecutor over the 4-combo pricing.

Caesars conveniences vs BetMGM: each selection carries ``type``
(home/away/over/under) directly, and ``market.line`` is already the
HOME-perspective line for run lines (home −1.5 ⇒ ``line:-1.5``), so no
team-name side-matching or sign handling is needed.

Source label ``caesars_direct``. Caesars offers only half-point lines (bot
targets are all half-points) so there is no integer-line fallback — an
off-grid target yields no row (honest "no price"). Sanity filter mirrors the
other books: drop a combo whose decimal exceeds ``SANITY_MULT_RATIO`` (1.5x)
the naive leg product.
"""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone

from mlb_sgp._shared import PricedRow, TargetLine, decimal_to_american
from mlb_sgp.caesars_client import CaesarsClient, Event

BOOK_NAME = "caesars"
SOURCE_LABEL = "caesars_direct"
SANITY_MULT_RATIO = 1.5
CZR_TARGET_PARALLELISM = 3   # gentle — the WAF rate-limits aggressive hits


def _leg(sel: dict, market: dict, event_id: str, competition_id: str,
         line: float | None = None) -> dict:
    """Build a /bets/details leg from a selection + its market (blueprint shape).

    `line` is the selection's own requested line. Caesars run-line markets carry
    a single home-perspective `market.line`; the AWAY selection's line is its
    negation. Totals share the same line for over/under. Passing the wrong sign
    for the away side makes Caesars decline the combo (drops the away corners),
    so the caller supplies the correct per-selection line.
    """
    price = sel.get("price") or {}
    return {
        "selectionId": sel.get("id"),
        "eventId": event_id,
        "marketId": market.get("id"),
        "competitionId": competition_id,
        "priceType": "fp",
        "related": True,
        "eachWay": False,
        "stakePerLine": 0,
        "line": market.get("line") if line is None else line,
        "selectionType": sel.get("type"),
        "name": sel.get("name"),
        "price": price,
        "origPrice": price,
    }


def _classify(name: str) -> tuple[str, str] | None:
    """(period, kind) for the game run-line / total-runs market, else None.

    EXACT-name match (not substring) — Caesars' event-detail tree is full of
    markets whose names contain "total"/"run line" but are NOT the game line:
    player props ("Jacob Wilson - Total Bases"), team totals ("Home Total
    Runs", "Colorado Rockies Total Runs"), other periods ("1st 7 Innings Total
    Runs"), and in-play variants ("Total Runs Live"). Only the canonical
    full-game and 1st-5-innings run-line/total-runs markets (main + alternate)
    qualify. Names are matched after stripping pipes/whitespace, lowercased.
    """
    low = (name or "").replace("|", " ").strip().lower()
    low = " ".join(low.split())  # collapse internal whitespace
    fg_spread = {"run line", "alternate run line"}
    fg_total = {"total runs", "alternate total runs"}
    f5_spread = {"1st 5 innings run line", "alternate 1st 5 innings run line",
                 "1st 5 innings - run line"}
    f5_total = {"1st 5 innings total runs", "alternate 1st 5 innings total runs",
                "1st 5 innings - total runs"}
    if low in fg_spread:
        return ("FG", "spread")
    if low in fg_total:
        return ("FG", "total")
    if low in f5_spread:
        return ("F5", "spread")
    if low in f5_total:
        return ("F5", "total")
    return None


def _price_d(leg: dict) -> float | None:
    p = leg.get("price") or {}
    d = p.get("d") if p.get("d") is not None else p.get("decimal")
    try:
        return float(d) if d else None
    except (TypeError, ValueError):
        return None


def parse_markets(event: dict) -> dict:
    """Index event.keyMarketGroups[].markets[] into pricing-ready legs.

    Returns {"FG": {"spreads": {home_line: {"home": leg, "away": leg}},
                    "totals":  {line: {"over": leg, "under": leg}}},
             "F5": {...}}  where each leg is a /bets/details leg dict.
    """
    out = {"FG": {"spreads": {}, "totals": {}}, "F5": {"spreads": {}, "totals": {}}}
    eid = event.get("id")
    cid = event.get("competitionId")
    for grp in event.get("keyMarketGroups", []) or []:
        for m in grp.get("markets", []) or []:
            cls = _classify(m.get("name", ""))
            if cls is None:
                continue
            period, kind = cls
            sels = m.get("selections") or []
            if not sels:
                continue
            line = m.get("line")
            if line is None:
                continue
            if kind == "spread":
                home = next((s for s in sels if s.get("type") == "home"), None)
                away = next((s for s in sels if s.get("type") == "away"), None)
                if home and away:
                    ml = float(line)
                    out[period]["spreads"][ml] = {
                        # home line is the market line (home perspective);
                        # away selection's line is its negation.
                        "home": _leg(home, m, eid, cid, line=ml),
                        "away": _leg(away, m, eid, cid, line=-ml),
                    }
            else:  # total — line is the total value
                over = next((s for s in sels if s.get("type") == "over"), None)
                under = next((s for s in sels if s.get("type") == "under"), None)
                if over and under:
                    out[period]["totals"][float(line)] = {
                        "over": _leg(over, m, eid, cid),
                        "under": _leg(under, m, eid, cid),
                    }
    return out


def _match_events(events: list[Event], targets: list[TargetLine]) -> dict[str, Event]:
    """game_id -> Event via canonical team resolution + UTC-hour bucket (strict)."""
    from canonical_match import load_team_dict, load_canonical_games, resolve_team_names
    from mlb_sgp._shared import _utc_bucket

    tgames = {t.game_id: (t.home_team, t.away_team, _utc_bucket(t.commence_time))
              for t in targets}
    team_dict = load_team_dict("mlb")
    canonical_games = load_canonical_games("mlb")
    matched: dict[str, Event] = {}
    exact: set[str] = set()
    for ev in events:
        resolved = resolve_team_names(ev.away_team, ev.home_team, team_dict, canonical_games)
        if not resolved or not resolved[0] or not resolved[1]:
            continue
        canon_away, canon_home = resolved
        eb = _utc_bucket(ev.start_time)
        for gid, (h, a, tb) in tgames.items():
            if h != canon_home or a != canon_away:
                continue
            if tb and eb:
                if tb == eb:
                    matched[gid] = ev
                    exact.add(gid)
            elif gid not in exact:
                matched[gid] = ev
    return matched


def price_sgps(
    target_lines: list[TargetLine],
    periods: tuple[str, ...] = ("FG",),
    client: CaesarsClient | None = None,
    verbose: bool = False,
) -> list[PricedRow]:
    if not target_lines:
        return []
    targets = [t for t in target_lines if t.period in periods]
    if not targets:
        return []

    if client is None:
        client = CaesarsClient(verbose=verbose)
    events = client.list_events()
    if not events:
        if verbose:
            print("  [czr] no events (token unverified / off-day) — no rows")
        return []
    matched_by_gid = _match_events(events, targets)
    if not matched_by_gid:
        return []

    targets_by_game: dict[str, list[TargetLine]] = {}
    for t in targets:
        if t.game_id in matched_by_gid:
            targets_by_game.setdefault(t.game_id, []).append(t)

    parsed_cache: dict[str, dict] = {}
    filtered: list[TargetLine] = []
    for gid, game_targets in targets_by_game.items():
        ev = matched_by_gid[gid]
        event = client.fetch_event(ev.event_id)
        if not event:
            continue
        parsed = parse_markets(event)
        parsed_cache[gid] = parsed
        kept = 0
        for t in game_targets:
            per = parsed.get(t.period)
            if per and float(t.spread) in per["spreads"] and float(t.total) in per["totals"]:
                filtered.append(t)
                kept += 1
        if verbose:
            print(f"  [czr] game {gid[:8]}: {len(game_targets)} targets → {kept} offered")

    if not filtered:
        return []

    fetch_now = datetime.now(timezone.utc)

    def _price_one_target(t: TargetLine) -> list[PricedRow]:
        parsed = parsed_cache.get(t.game_id)
        if not parsed:
            return []
        per = parsed[t.period]
        sp = per["spreads"].get(float(t.spread))
        to = per["totals"].get(float(t.total))
        if not sp or not to:
            return []
        prefix = "" if t.period == "FG" else "F5 "
        combos = (
            ("Home Spread + Over", sp["home"], to["over"]),
            ("Home Spread + Under", sp["home"], to["under"]),
            ("Away Spread + Over", sp["away"], to["over"]),
            ("Away Spread + Under", sp["away"], to["under"]),
        )

        def _price_combo(combo_name, sp_leg, to_leg):
            priced = client.price_combo([sp_leg, to_leg])
            if priced is None:
                return combo_name, None
            dec = priced["decimal"]
            sd, td = _price_d(sp_leg), _price_d(to_leg)
            if sd and td:
                naive = sd * td
                if naive > 0 and dec > naive * SANITY_MULT_RATIO:
                    if verbose:
                        print(f"      [czr] SANITY-DROP {combo_name} dec={dec:.2f} naive={naive:.2f}")
                    return combo_name, None
            return combo_name, (dec, priced.get("american"))

        rows: list[PricedRow] = []
        with ThreadPoolExecutor(max_workers=4) as pool:
            futs = [pool.submit(_price_combo, c, s, o) for c, s, o in combos]
            for f in as_completed(futs):
                try:
                    combo_name, res = f.result()
                except Exception:
                    continue
                if res is None:
                    continue
                dec, am = res
                rows.append(PricedRow(
                    game_id=t.game_id,
                    combo=prefix + combo_name,
                    period=t.period,
                    spread_line=t.spread,
                    total_line=t.total,
                    bookmaker=BOOK_NAME,
                    source=SOURCE_LABEL,
                    sgp_decimal=round(dec, 4),
                    sgp_american=am if am is not None else decimal_to_american(dec),
                    fetch_time=fetch_now,
                ))
        return rows

    out: list[PricedRow] = []
    with ThreadPoolExecutor(max_workers=CZR_TARGET_PARALLELISM) as pool:
        futures = [pool.submit(_price_one_target, t) for t in filtered]
        for f in as_completed(futures):
            try:
                out.extend(f.result())
            except Exception as e:
                if verbose:
                    print(f"  [czr] target error: {e}")
    return out
