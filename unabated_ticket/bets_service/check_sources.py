"""One-shot login check: can this machine reach every venue the bets service polls?

    python -m unabated_ticket.bets_service.check_sources     (from the repo root)

Built for moving the service off the Mac (phone page plan, step 0): run it once
on the new server and it says, per venue, whether the login and one full read
worked from there. Offshore books may refuse or flag logins from a data-center
address, and this answers that before anything is built on top.

Inputs:  the same credentials the service reads (config.py: env, bets_service/.env,
         kalshi_draft/.env, bet_logger/.env, the Novig token file, the BetOnline
         cookie file), and the same source list (service.OPTIONAL_SOURCE_FACTORIES).
         Also GETs one public Unabated league snapshot, the feed the Edges scan reads.
Outputs: one line per venue on stdout: ok with a record count and seconds, not
         configured, or FAILED with the error. Record contents are never printed.
         Exit code 0 when nothing FAILED (not-configured venues don't fail it), 1 otherwise.
Side effects: one fetch() per configured source, i.e. the same read-only venue
         requests one service poll makes (logins included). BetOnline's refresh
         can rewrite its cookie file, exactly as a poll does. No DuckDB read or write.
"""
import logging
import sys
import time
import urllib.request
from collections.abc import Callable
from dataclasses import dataclass

from unabated_ticket.bets_service import service
from unabated_ticket.bets_service.sources import Source
from unabated_ticket.bets_service.sources.kalshi import KalshiSource

# NFL's league snapshot: public, no login. Proves the Edges scan can run here.
UNABATED_SNAPSHOT_URL = "https://content.unabated.com/markets/v2/league/1/odds.json"
UNABATED_TIMEOUT_SEC = 30
# CloudFront answers some clients with a block page; a real snapshot is far larger.
MIN_SNAPSHOT_BYTES = 1000
MAX_ERROR_CHARS = 300


@dataclass(frozen=True)
class CheckResult:
    venue: str
    status: str  # "ok" | "not configured" | "FAILED"
    detail: str


def check_source(venue: str, build: Callable[[], Source | None],
                 clock: Callable[[], float] = time.monotonic) -> CheckResult:
    """Build the venue's source and run one fetch(). A factory returning None
    means its credentials are absent; a factory raising means they are present
    but unusable, which is a failure worth seeing."""
    try:
        source = build()
    except Exception as exc:  # noqa: BLE001 — reported per venue, never swallowed
        return CheckResult(venue, "FAILED", f"setup: {_describe(exc)}")
    if source is None:
        return CheckResult(venue, "not configured", "no credentials on this machine")
    started = clock()
    try:
        records = source.fetch()
    except Exception as exc:  # noqa: BLE001
        return CheckResult(venue, "FAILED", f"{_describe(exc)} after {clock() - started:.1f}s")
    return CheckResult(venue, "ok", f"{len(records)} records in {clock() - started:.1f}s")


def check_unabated_feed(open_url: Callable[..., object] = urllib.request.urlopen) -> CheckResult:
    """GET the public NFL snapshot the Edges scan starts from."""
    try:
        with open_url(UNABATED_SNAPSHOT_URL, timeout=UNABATED_TIMEOUT_SEC) as response:
            size = len(response.read())
    except Exception as exc:  # noqa: BLE001
        return CheckResult("unabated feed", "FAILED", _describe(exc))
    if size < MIN_SNAPSHOT_BYTES:
        return CheckResult("unabated feed", "FAILED",
                           f"expected a snapshot of >= {MIN_SNAPSHOT_BYTES} bytes, got {size}")
    return CheckResult("unabated feed", "ok", f"NFL snapshot {size:,} bytes")


def run_checks() -> list[CheckResult]:
    """Every venue the service polls, in the order it registers them, then the feed."""
    results = [check_source("kalshi", KalshiSource)]
    for venue, factory in service.OPTIONAL_SOURCE_FACTORIES:
        results.append(check_source(venue, factory))
    results.append(check_unabated_feed())
    return results


def format_results(results: list[CheckResult]) -> str:
    width = max(len(result.venue) for result in results)
    return "\n".join(f"{result.venue:<{width}}  {result.status:<14}  {result.detail}" for result in results)


def _describe(exc: Exception) -> str:
    text = f"{type(exc).__name__}: {exc}"
    return text if len(text) <= MAX_ERROR_CHARS else text[:MAX_ERROR_CHARS] + "..."


def main() -> int:
    # Sources log their own warnings (e.g. a missing login and its fix); keep them
    # on stderr so stdout is just the table.
    logging.basicConfig(level=logging.WARNING, stream=sys.stderr, format="%(levelname)s %(name)s: %(message)s")
    results = run_checks()
    print(format_results(results))
    return 1 if any(result.status == "FAILED" for result in results) else 0


if __name__ == "__main__":
    sys.exit(main())
