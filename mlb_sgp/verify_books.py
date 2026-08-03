"""Six-book live verification harness — the bot-restart gate (issue #42).

Answers one question per book, on the path a quote will actually use
(``SGPService.price_on_demand``): **can this book price a real combo right
now, and how fast?**

Why this exists rather than a checklist of curl commands: a book returning
HTTP 200 proves nothing (issue #32's meta-bug), and a book returning no price
is ambiguous — it may be down, or it may simply not offer the rung you asked
for. Every one of the three Phase-3 repairs (#39 DK, #40 Novig, #41 Caesars)
had a root cause that differed from the filed diagnosis, so this harness
reports those two cases as DIFFERENT outcomes and never collapses them.

Outcomes per book:

  ``priced``       resolved legs and returned a fair
  ``not_offered``  structure is healthy but carries none of the candidate
                   lines — a coverage fact, NOT a failure
  ``no_price``     legs resolved, book declined to price them
  ``down``         event match or structure fetch failed — a real failure

Coverage is probed through each book's own ``resolve`` hook, which is a pure
lookup against an already-fetched structure: it costs ZERO extra wire calls
and exercises exactly the function the quote path calls. That matters because
Novig rate-limits (26 rapid price calls trip a 403, issue #40).

Usage
-----
    python3 mlb_sgp/verify_books.py                 # all books, human report
    python3 mlb_sgp/verify_books.py --json out.json # machine-readable too
    python3 mlb_sgp/verify_books.py --books novig,caesars --repeat 5

Inputs:  live book APIs; the MLB canonical team dictionary.
Outputs: a report on stdout; optionally a JSON blob.
Side effects: **none** — no DB is opened for write (``health_db_path=None``),
no bot is started or stopped, no file in the repo is modified.
"""
from __future__ import annotations

import argparse
import json
import statistics
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent
for _p in (str(_REPO_ROOT), str(_REPO_ROOT / "mlb_sgp"),
           str(_REPO_ROOT / "Answer Keys")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from kalshi_common.legset import CanonicalLeg          # noqa: E402
from kalshi_common.sgp_service import SGPService       # noqa: E402
from mlb_sgp._shared import FetchCounters, GameRef     # noqa: E402

ALL_BOOKS = ("draftkings", "fanduel", "prophetx", "novig", "betmgm", "caesars")

# Candidate ladder probed against each book's structure. Deliberately spans
# BOTH half-point and whole-number totals: the books split on this — Novig and
# BetMGM quote .5 ladders, FanDuel and Caesars carry only the game's main
# line, which is often a whole number.
CANDIDATE_SPREADS = (-1.5, 1.5, -2.5, 2.5)
CANDIDATE_TOTALS = (7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5)

# Novig 403s on bursts (26 rapid price calls, issue #40) and every shape here
# fans out into a 2^N partition, so a full run is easily 60+ wire price calls
# per book. Pace between calls; --pacing raises it when a run needs to
# distinguish "this book is slow" from "this book rate-limited us".
DEFAULT_PACING_SEC = 1.0

GAME_ID = "verify"


def _leg(market_type: str, line, side) -> CanonicalLeg:
    return CanonicalLeg(game_id=GAME_ID, market_type=market_type,
                        line=line, side=side)


def shapes_for(spread: float, total: float) -> dict:
    """The shapes issue #42 requires, at one (spread, total).

    ``three_leg`` exercises Route A's leg ceiling
    (ON_DEMAND_MAX_PARTITION_LEGS = 3) and the partition/transfer choice.

    ``spread_x_ml`` is measured separately because it is the shape books
    disagree about: DK, FanDuel and BetMGM decline it, which is also what
    makes ``three_leg`` unpriceable at those books (an MLB same-game 3-leg
    has only spread / total / ml to draw on, so it necessarily contains the
    spread+ml pair). Without this row, "3-leg fails at DK" reads like a
    transport problem instead of a book rule.
    """
    return {
        "spread_x_total": [_leg("spread", spread, "home"),
                           _leg("total", total, "over")],
        "ml_x_total": [_leg("ml", None, "home"),
                       _leg("total", total, "over")],
        "spread_x_ml": [_leg("spread", spread, "home"),
                        _leg("ml", None, "home")],
        "three_leg": [_leg("spread", spread, "home"),
                      _leg("total", total, "over"),
                      _leg("ml", None, "home")],
    }


def pick_game() -> GameRef:
    """A live, not-yet-started MLB game, resolved to canonical team names.

    Sourced from Novig's event list purely because it is the cheapest
    anonymous listing available; the GameRef is book-agnostic from there.
    """
    from canonical_match import (load_canonical_games, load_team_dict,
                                 resolve_team_names)
    from mlb_sgp.novig_client import NovigClient

    team_dict = load_team_dict("mlb")
    canonical = load_canonical_games("mlb")
    now = datetime.now(timezone.utc)
    for event in NovigClient(verbose=False).list_events():
        resolved = resolve_team_names(event.away_team, event.home_team,
                                      team_dict, canonical)
        if not resolved or not resolved[0] or not resolved[1]:
            continue
        start = datetime.fromisoformat(event.start_time)
        if start <= now:
            continue
        return GameRef(game_id=GAME_ID, home_team=resolved[1],
                       away_team=resolved[0], commence_time=start)
    raise SystemExit("no upcoming MLB game could be resolved")


def probe_coverage(service: SGPService, book: str, game: GameRef) -> dict:
    """Which candidate lines this book's structure actually carries.

    Returns ``{"status": ..., "structure_sec": float, "spreads": [...],
    "totals": [...]}``. Only ONE structure fetch happens here; every line
    check after it is a pure in-memory resolve.
    """
    counters = FetchCounters(book, "verify")
    hooks = service._book_on_demand_hooks(book, counters=counters)
    state = service._ensure_client(book)

    started = time.time()
    try:
        event = hooks["match_event"](state.client, game)
    except Exception as exc:                       # transport / auth failure
        return {"status": "down", "error": f"match_event: {exc}",
                "structure_sec": time.time() - started,
                "spreads": [], "totals": []}
    if event is None:
        return {"status": "down", "error": "event did not match",
                "structure_sec": time.time() - started,
                "spreads": [], "totals": []}

    try:
        structure = hooks["build_structure"](state.client, event, game)
    except Exception as exc:
        return {"status": "down", "error": f"build_structure: {exc}",
                "structure_sec": time.time() - started,
                "spreads": [], "totals": []}
    structure_sec = time.time() - started
    if structure is None:
        return {"status": "down", "error": "structure fetch returned nothing",
                "structure_sec": structure_sec, "spreads": [], "totals": []}

    def offers(legs) -> bool:
        try:
            return hooks["resolve"](structure, legs, game.home_team,
                                    game.away_team) is not None
        except Exception:
            return False

    spreads = [s for s in CANDIDATE_SPREADS
               if offers([_leg("spread", s, "home")])]
    totals = [t for t in CANDIDATE_TOTALS
              if offers([_leg("total", t, "over")])]
    status = "ok" if (spreads and totals) else "not_offered"
    return {"status": status, "structure_sec": structure_sec,
            "spreads": spreads, "totals": totals}


def price_shape(service: SGPService, book: str, game: GameRef, legs) -> dict:
    """One timed ``price_on_demand`` call. Never raises."""
    started = time.time()
    result = service.price_on_demand(book, game, legs)
    elapsed = time.time() - started
    if result is None or getattr(result, "fair", None) is None:
        return {"outcome": "no_price", "sec": elapsed}
    return {"outcome": "priced", "fair": result.fair, "route": result.route,
            "cells": result.n_cells_priced, "sec": elapsed}


def verify_book(book: str, game: GameRef, repeat: int,
                pacing: float = DEFAULT_PACING_SEC) -> dict:
    """Cold pass on a fresh service, then ``repeat`` warm passes.

    Cold is measured on a brand-new ``SGPService`` so it includes event
    discovery, structure fetch and (for Caesars) the AWS-WAF token mint —
    the latency a maker actually pays on the first RFQ of a game.
    """
    report = {"book": book}
    cold_service = SGPService(books=(book,), health_db_path=None)
    coverage = probe_coverage(cold_service, book, game)
    report["coverage"] = coverage
    if coverage["status"] == "down":
        report["outcome"] = "down"
        return report
    if coverage["status"] == "not_offered":
        report["outcome"] = "not_offered"
        return report

    spread = coverage["spreads"][0]
    total = coverage["totals"][0]
    report["priced_at"] = {"spread": spread, "total": total}
    shapes = shapes_for(spread, total)

    # Cold: a fresh service so nothing is cached. The coverage probe above
    # already warmed THIS service's structure cache, so cold latency is
    # measured on yet another fresh one.
    report["cold"] = {}
    for name, legs in shapes.items():
        fresh = SGPService(books=(book,), health_db_path=None)
        report["cold"][name] = price_shape(fresh, book, game, legs)
        time.sleep(pacing)

    report["warm"] = {}
    for name, legs in shapes.items():
        samples = []
        for _ in range(repeat):
            samples.append(price_shape(cold_service, book, game, legs))
            time.sleep(pacing)
        report["warm"][name] = samples

    outcomes = {s["outcome"] for s in report["cold"].values()}
    report["outcome"] = "priced" if "priced" in outcomes else "no_price"
    if hasattr(cold_service._state.get(book), "client"):
        mint = getattr(cold_service._state[book].client, "last_mint_sec", None)
        if mint:
            report["token_mint_sec"] = mint
    return report


def _percentile(values, pct):
    if not values:
        return None
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    return statistics.quantiles(ordered, n=100)[min(pct, 99) - 1]


def _secs(value) -> str:
    """Format a latency for the table; '-' when there is no sample."""
    return "-" if value is None else f"{value:.2f}"


def summarize(reports: list[dict]) -> None:
    print("\n" + "=" * 78)
    print("PER-BOOK OUTCOME")
    print("=" * 78)
    print(f"{'book':11s} {'outcome':12s} {'struct s':>9s} "
          f"{'spreads':>8s} {'totals':>7s}  note")
    for r in reports:
        cov = r.get("coverage", {})
        note = cov.get("error", "")
        if r["outcome"] == "priced":
            at = r["priced_at"]
            note = f"priced at spread={at['spread']} total={at['total']}"
        print(f"{r['book']:11s} {r['outcome']:12s} "
              f"{cov.get('structure_sec', 0):9.2f} "
              f"{len(cov.get('spreads', [])):8d} {len(cov.get('totals', [])):7d}"
              f"  {note}")

    print("\n" + "=" * 78)
    print("SHAPE COVERAGE (cold call)")
    print("=" * 78)
    shape_names = ("spread_x_total", "ml_x_total", "spread_x_ml",
                   "three_leg")
    print(f"{'book':11s} " + " ".join(f"{n:>16s}" for n in shape_names))
    for r in reports:
        if "cold" not in r:
            print(f"{r['book']:11s} " + " ".join(f"{'-':>16s}"
                                                 for _ in shape_names))
            continue
        cells = []
        for name in shape_names:
            s = r["cold"].get(name, {})
            if s.get("outcome") == "priced":
                cells.append(f"{s['fair']:.4f}/{s['route'][:4]}")
            else:
                cells.append("no_price")
        print(f"{r['book']:11s} " + " ".join(f"{c:>16s}" for c in cells))

    print("\n" + "=" * 78)
    print("LATENCY — price_on_demand, seconds")
    print("=" * 78)
    print(f"{'book':11s} {'shape':16s} {'cold':>7s} {'warm p50':>9s} "
          f"{'warm p95':>9s} {'n':>3s}")
    for r in reports:
        if "warm" not in r:
            continue
        for name in shape_names:
            cold = r["cold"].get(name, {}).get("sec")
            warm = [s["sec"] for s in r["warm"].get(name, [])
                    if s["outcome"] == "priced"]
            p50 = statistics.median(warm) if warm else None
            p95 = _percentile(warm, 95)
            print(f"{r['book']:11s} {name:16s} {_secs(cold):>7s} "
                  f"{_secs(p50):>9s} {_secs(p95):>9s} {len(warm):>3d}")
        if r.get("token_mint_sec"):
            print(f"{r['book']:11s} {'(token mint)':16s} "
                  f"{r['token_mint_sec']:>7.2f}")

    cross_book_agreement(reports)


def cross_book_agreement(reports: list[dict]) -> None:
    """Books that priced the SAME rung must agree — the strongest correctness
    signal available without a ground truth."""
    by_line: dict[tuple, list] = {}
    for r in reports:
        if r.get("outcome") != "priced":
            continue
        cell = r["cold"].get("spread_x_total", {})
        if cell.get("outcome") != "priced":
            continue
        at = r["priced_at"]
        by_line.setdefault((at["spread"], at["total"]), []).append(
            (r["book"], cell["fair"]))

    print("\n" + "=" * 78)
    print("CROSS-BOOK AGREEMENT — spread x total, by rung")
    print("=" * 78)
    for line, entries in sorted(by_line.items()):
        fairs = [f for _b, f in entries]
        spread_pct = ((max(fairs) - min(fairs)) / min(fairs) * 100
                      if len(fairs) > 1 else 0.0)
        print(f"  spread={line[0]} total={line[1]}: " +
              ", ".join(f"{b} {f:.4f}" for b, f in sorted(entries)) +
              (f"   range {max(fairs) - min(fairs):.4f} ({spread_pct:.1f}%)"
               if len(fairs) > 1 else "   (single book)"))
    if not by_line:
        print("  no book priced a spread x total combo")


def price_common_rung(reports: list[dict], game: GameRef,
                      pacing: float = DEFAULT_PACING_SEC) -> dict:
    """Re-price ONE rung — the one the most books carry — at every book that
    carries it.

    The per-book pass lets each book pick its own first offered line, which
    makes those fairs incomparable. Cross-book agreement is the strongest
    correctness signal available without a ground truth (#39 closed with it
    unverified at n=1), so it needs every book answering the SAME question.
    """
    offered_by: dict[tuple, list] = {}
    for r in reports:
        cov = r.get("coverage", {})
        for spread in cov.get("spreads", ()):
            for total in cov.get("totals", ()):
                offered_by.setdefault((spread, total), []).append(r["book"])
    if not offered_by:
        return {}
    rung, books = max(offered_by.items(), key=lambda kv: len(kv[1]))

    print(f"\n--- common rung spread={rung[0]} total={rung[1]} "
          f"({len(books)} books) ...", flush=True)
    legs = shapes_for(rung[0], rung[1])["spread_x_total"]
    fairs = {}
    for book in books:
        service = SGPService(books=(book,), health_db_path=None)
        cell = price_shape(service, book, game, legs)
        if cell["outcome"] == "priced":
            fairs[book] = cell["fair"]
        time.sleep(pacing)
    return {"rung": rung, "books_offering": books, "fairs": fairs}


def report_common_rung(common: dict) -> None:
    print("\n" + "=" * 78)
    print("CROSS-BOOK AGREEMENT — same rung, every book that carries it")
    print("=" * 78)
    if not common or not common.get("fairs"):
        print("  no rung was carried by any book")
        return
    rung = common["rung"]
    fairs = common["fairs"]
    print(f"  spread={rung[0]} total={rung[1]}  "
          f"({len(common['books_offering'])} books carry it, "
          f"{len(fairs)} priced)")
    for book, fair in sorted(fairs.items(), key=lambda kv: kv[1]):
        print(f"      {book:11s} {fair:.4f}")
    if len(fairs) > 1:
        lo, hi = min(fairs.values()), max(fairs.values())
        median = statistics.median(fairs.values())
        print(f"  range {lo:.4f}..{hi:.4f}  spread {hi - lo:.4f} "
              f"({(hi - lo) / lo * 100:.1f}% of low)  median {median:.4f}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--books", default=",".join(ALL_BOOKS))
    parser.add_argument("--repeat", type=int, default=3,
                        help="warm samples per shape (default 3)")
    parser.add_argument("--json", dest="json_path", default=None)
    parser.add_argument("--pacing", type=float,
                        default=DEFAULT_PACING_SEC,
                        help="seconds between price calls "
                             "(raise to clear Novig's rate limit)")
    args = parser.parse_args()

    books = tuple(b.strip() for b in args.books.split(",") if b.strip())
    game = pick_game()
    print(f"game: {game.away_team} @ {game.home_team}  "
          f"start={game.commence_time}")
    print(f"books: {', '.join(books)}   warm samples: {args.repeat}")

    reports = []
    for book in books:
        print(f"\n--- {book} ...", flush=True)
        try:
            reports.append(verify_book(book, game, args.repeat,
                                       args.pacing))
        except Exception as exc:
            # A harness bug must not be reported as a dead book.
            reports.append({"book": book, "outcome": "harness_error",
                            "coverage": {"error": repr(exc)}})
    common = price_common_rung(reports, game, args.pacing)
    summarize(reports)
    report_common_rung(common)

    if args.json_path:
        blob = {"game": {"home": game.home_team, "away": game.away_team,
                         "start": str(game.commence_time)},
                "captured_at": datetime.now(timezone.utc).isoformat(),
                "common_rung": common,
                "reports": reports}
        Path(args.json_path).write_text(json.dumps(blob, indent=2, default=str))
        print(f"\nwrote {args.json_path}")

    down = [r["book"] for r in reports if r["outcome"] in
            ("down", "harness_error")]
    if down:
        print(f"\nFAIL — books down: {', '.join(down)}")
        return 1
    print("\nAll books reachable and structurally healthy.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
