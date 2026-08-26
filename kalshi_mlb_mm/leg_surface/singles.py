"""Singles-route ingest: one whole-slate scrape per pass (issue #96).

Inputs: the DK / FD singles scrapers (in-process) and the current slate.
Outputs: SurfaceRows + counts. Side effects: live HTTP to one book.

This route exists because two books cannot serve the structure route:
  * DraftKings publishes no structure odds at all (#95: 21/21
    ``no_structure_odds``), and its ``calculateBets`` price host returns 403
    for every set size as of 2026-08-25 — so DK's read endpoints, which the
    singles scraper uses, are its only working path.
  * FanDuel's SGP structure carries exactly ONE total line per period (its own
    main) while its singles scraper has the full ladder, so FD's FG/F5 TOTALS
    come from here and everything else from the structure.

It calls ``collect_singles_rows``, never ``scrape_singles``: the latter does a
CREATE OR REPLACE of the production dk_odds / fd_odds snapshot that MLB.R and
the dashboard read, and the maker's ingest cadence has no business becoming
the dashboard's refresh rate. Reading those DuckDBs instead was also rejected
— neither scraper has a cron and both files were last written 2026-08-17.

Neither scraper emits single-inning markets (``_SINGLE_INNING_RE`` in DK's
``classify_market``, mirrored in FD's), so I1 legs never come from this route.
"""
from __future__ import annotations

import logging
import sys
from datetime import datetime, timezone
from pathlib import Path

from mlb_sgp._shared import american_to_decimal

from kalshi_mlb_mm.leg_surface.devig import ExclusionCounts, devig_rung
from kalshi_mlb_mm.leg_surface.slate import rungs

log = logging.getLogger(__name__)

ROUTE = "singles"

# The scrapers' wide rows carry FG/F3/F5/F7; the surface prices FG and F5.
_SUPPORTED_PERIODS = ("FG", "F5")


_REPO_ROOT = Path(__file__).resolve().parents[2]


def _ensure_scraper_path() -> None:
    """Make the scrapers' top-level imports resolvable.

    They import their clients by bare name (``from dk_client import ...``),
    which only works with ``mlb_sgp/`` itself on sys.path — true for CLI runs
    (cwd=mlb_sgp/) but not for a bot at the repo root. SGPService does the
    same insert for the structure route; this route never builds one, so it
    cannot rely on that having happened.
    """
    for path in (_REPO_ROOT / "mlb_sgp", _REPO_ROOT / "Answer Keys"):
        if str(path) not in sys.path:
            sys.path.insert(0, str(path))


def _scrape(book: str) -> list[dict]:
    """One book's whole-slate singles rows, no production DB write."""
    _ensure_scraper_path()
    if book == "draftkings":
        from scraper_draftkings_singles import collect_singles_rows
    elif book == "fanduel":
        from scraper_fanduel_singles import collect_singles_rows
    else:
        raise ValueError(f"no singles scraper for book {book!r}")
    return collect_singles_rows()


def _as_utc(value) -> datetime | None:
    """Scraper start times arrive as ISO strings or datetimes; normalize to
    aware UTC so the tolerance comparison is never tz-mixed."""
    if value is None:
        return None
    if isinstance(value, str):
        text = value.replace("Z", "+00:00")
        try:
            value = datetime.fromisoformat(text)
        except ValueError:
            return None
    if not isinstance(value, datetime):
        return None
    return (value.replace(tzinfo=timezone.utc) if value.tzinfo is None
            else value.astimezone(timezone.utc))


def match_book_game(game, rows: list[dict], *, tolerance_min: float):
    """The book's rows for ONE Kalshi game -> (rows, reason).

    Matched on canonical team names PLUS start time within tolerance, and
    AMBIGUITY FAILS CLOSED. Team names alone silently return the wrong game of
    a doubleheader: #95's FD slate carried two "Philadelphia Phillies @ Seattle
    Mariners" rows a day apart, and matching on teams took tomorrow's game,
    producing fairs off by 0.05-0.11 in probability. That is a wrong number,
    not a decline — the one failure a maker must never have.
    """
    kalshi_start = _as_utc(game.start_utc)
    if kalshi_start is None:
        return [], "game_unmatched"
    tolerance_sec = tolerance_min * 60.0
    by_book_game: dict = {}
    for row in rows:
        if row.get("home_team") != game.home_team:
            continue
        if row.get("away_team") != game.away_team:
            continue
        start = _as_utc(row.get("game_start_time"))
        if start is None:
            continue
        if abs((start - kalshi_start).total_seconds()) > tolerance_sec:
            continue
        by_book_game.setdefault(row.get("game_id"), []).append(row)
    if not by_book_game:
        return [], "game_unmatched"
    if len(by_book_game) > 1:
        log.warning("surface singles: %s matched %d book games within %.0fmin "
                    "— declining", game.game_id, len(by_book_game),
                    tolerance_min)
        return [], "game_ambiguous"
    return next(iter(by_book_game.values())), None


def book_decimals(book_rows: list[dict]) -> dict:
    """The book's wide rows -> {(period, market_type, line): {side: decimal}}.

    The scrapers emit one wide row per (period, market kind, line) with both
    sides' American prices in it, which is already a rung. Spread lines are
    stored home-perspective to match ``CanonicalLeg.line``.
    """
    out: dict = {}

    def put(key, side, american):
        if american is None:
            return
        try:
            decimal = american_to_decimal(int(american))
        except (TypeError, ValueError):
            return
        out.setdefault(key, {})[side] = decimal

    for row in book_rows:
        period = row.get("period")
        if period not in _SUPPORTED_PERIODS:
            continue
        if row.get("home_ml") is not None or row.get("away_ml") is not None:
            put((period, "ml", None), "home", row.get("home_ml"))
            put((period, "ml", None), "away", row.get("away_ml"))
        home_spread = row.get("home_spread")
        if home_spread is not None:
            key = (period, "spread", round(float(home_spread), 2))
            put(key, "home", row.get("home_spread_price"))
            put(key, "away", row.get("away_spread_price"))
        total = row.get("total")
        if total is not None:
            key = (period, "total", round(float(total), 2))
            put(key, "over", row.get("over_price"))
            put(key, "under", row.get("under_price"))
    return out


def price_game(book: str, game, book_rows: list[dict], *, owned_markets,
               built_at: datetime, band_min: float,
               band_max: float) -> tuple[list, ExclusionCounts]:
    """One (book, game) from already-scraped rows. No network here."""
    counts = ExclusionCounts()
    decimals_by_rung = book_decimals(book_rows)
    rows = []
    for (period, market_type, line), rung_legs in rungs(game.legs).items():
        if (market_type, period) not in owned_markets:
            continue          # this route does not own the rung; not a miss
        key = (period, market_type,
               None if line is None else round(float(line), 2))
        outcome = devig_rung(book=book, route=ROUTE, game=game,
                             legs=rung_legs,
                             decimals=decimals_by_rung.get(key, {}),
                             built_at=built_at, band_min=band_min,
                             band_max=band_max)
        if outcome.reason is not None:
            counts.bump(outcome.reason)
            continue
        rows.extend(outcome.rows)
    return rows, counts


def run_pass(book: str, games, *, owned_markets, band_min: float,
             band_max: float, tolerance_min: float
             ) -> tuple[list, ExclusionCounts, int]:
    """One book's pass: scrape the whole slate once, then match per game.

    ``built_at`` is the scrape's own fetch_time, so every row carries the age
    of the payload it came from rather than the age of the local match.
    """
    scraped = _scrape(book)
    built_at = _as_utc(scraped[0].get("fetch_time")) if scraped else None
    built_at = built_at or datetime.now(timezone.utc)
    rows: list = []
    counts = ExclusionCounts()
    games_priced = 0
    for game in games:
        book_rows, reason = match_book_game(game, scraped,
                                            tolerance_min=tolerance_min)
        if reason is not None:
            counts.bump(reason)
            continue
        game_rows, game_counts = price_game(
            book, game, book_rows, owned_markets=owned_markets,
            built_at=built_at, band_min=band_min, band_max=band_max)
        counts.add(game_counts)
        if game_rows:
            games_priced += 1
            rows.extend(game_rows)
    return rows, counts, games_priced
