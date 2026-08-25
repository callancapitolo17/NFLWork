"""Book-consensus fair value for P(YRFI) — one game, 4 books, live.

Reuses the issue-#87 pipeline: a KXMLBRFI market is exactly a 1st-inning
total at 0.5 (period "I1"), and SGPService.price_on_demand on a single such
leg runs the 2-cell partition — an exact two-way probit devig of the book's
YRFI/NRFI prices. No new book code.

Inputs: RfiGame. Outputs: FairResult or None. Side effects: live HTTP to the
books via the shared per-book clients (no DB writes).
"""
import logging
import statistics
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from datetime import datetime, timezone

from scipy.stats import norm

from kalshi_common.legset import CanonicalLeg
from kalshi_common.sgp_service import SGPService
from mlb_sgp._shared import GameRef

from kalshi_rfi import config
from kalshi_rfi.discovery import RfiGame

log = logging.getLogger(__name__)


@dataclass
class FairResult:
    consensus_yes: float            # median devigged P(YRFI) of gated books
    book_fairs: dict = field(default_factory=dict)  # book -> devigged fair
    sigma_z: float | None = None    # sample stddev of probit z's (None if n<2)
    fetched_at: datetime = field(
        default_factory=lambda: datetime.now(timezone.utc))


def consensus(book_fairs: dict) -> tuple[float, float] | None:
    """(median fair, sigma_z) if the dispersion gate passes, else None.

    Same semantics as kalshi_mlb_mm issue #20: >= MIN_BOOKS books, sample
    stddev (ddof=1) of probit-transformed fairs <= SIGMA_Z_MAX, no outlier
    removal — a loud dissenter may be the informed book, so decline.
    """
    fairs = {b: f for b, f in book_fairs.items() if f is not None}
    if len(fairs) < max(config.MIN_BOOKS, 2):
        return None
    zs = [float(norm.ppf(min(max(f, 1e-6), 1.0 - 1e-6)))
          for f in fairs.values()]
    sigma_z = statistics.stdev(zs)
    if sigma_z > config.SIGMA_Z_MAX:
        return None
    return statistics.median(fairs.values()), sigma_z


class FairService:
    """Holds one SGPService (persistent per-book clients) for the daemon."""

    def __init__(self, service: SGPService | None = None):
        # structure_ttl_sec=0.0: our fair IS the structure's odds, so every
        # fetch must hit the book's wire — the TTL cache exists for the MM's
        # cadence, and serving cached odds here would be stale fair sold as
        # live (adversarial review 2026-08-25, finding 1).
        # single_leg_structure_fair=True: opt in to the n==1 fast path;
        # default-off keeps the MM/taker SGP pipeline untouched.
        self._service = service or SGPService(
            books=config.BOOKS, health_db_path=None,
            structure_ttl_sec=0.0, single_leg_structure_fair=True)
        self._pool = ThreadPoolExecutor(max_workers=len(config.BOOKS))

    def fetch(self, game: RfiGame) -> FairResult | None:
        """Live 4-book fetch for one game. None if the gate declines."""
        leg = CanonicalLeg(game.suffix, "total", 0.5, "over", "I1")
        ref = GameRef(game_id=game.suffix, home_team=game.home_team,
                      away_team=game.away_team,
                      commence_time=game.commence_utc)

        def one(book: str):
            res = self._service.price_on_demand(book, ref, [leg])
            return book, (res.fair if res is not None else None)

        book_fairs = dict(self._pool.map(one, config.BOOKS))
        gated = consensus(book_fairs)
        priced = {b: f for b, f in book_fairs.items() if f is not None}
        if gated is None:
            log.info("fair: %s declined (books=%d fairs=%s)",
                     game.ticker, len(priced),
                     {b: round(f, 3) for b, f in priced.items()})
            return None
        fair, sigma_z = gated
        return FairResult(consensus_yes=fair, book_fairs=priced,
                          sigma_z=sigma_z)

    def close(self):
        self._pool.shutdown(wait=False)
        self._service.close()
