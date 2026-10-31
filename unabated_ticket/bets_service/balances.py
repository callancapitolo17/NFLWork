"""Live account balance per venue, for the tracker's "Balance by book" card.

A source that can read its account balance defines `balance() -> dict` with
    cash       dollars on account, money in open bets left out
    inBets     dollars in open bets (Kalshi: what the open positions are worth now)
    available  dollars the book lets you bet now (credit included), or None
    creditLimit dollars of credit, or None
and raises on any failure. The poll loop calls it for every source at startup,
then after each fetch(). A source holding several wallets (Kalshi's primary
and the bot's subaccount) instead defines `balance_reads() -> {venue: read}`,
one balance() per wallet, each recorded (and failing) on its own.

Inputs:  the sources' own sessions (no new logins).
Outputs: GET /balances.json -> {generatedAt, balances: {venue: reading}}.
Side effects: none on disk. Readings live in memory only (a restart re-reads
         every venue at once); a failed read keeps the last good
         numbers and says why, and a venue never read reads "unavailable", never $0.
"""
import logging
import threading
from collections.abc import Callable
from datetime import datetime, timezone

log = logging.getLogger(__name__)

BALANCE_FIELDS = ("cash", "inBets", "available", "creditLimit")


def _iso_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def validate_reading(raw: object) -> dict:
    """The four fields as floats (or None); raises when `cash` is missing,
    since a balance with no cash figure is no balance at all."""
    if not isinstance(raw, dict):
        raise RuntimeError(f"balance() must return a dict, got {type(raw).__name__}")
    reading = {}
    for field in BALANCE_FIELDS:
        value = raw.get(field)
        if value is not None and (isinstance(value, bool) or not isinstance(value, (int, float))):
            raise RuntimeError(f"balance field {field} must be a number or None, got {value!r}")
        reading[field] = float(value) if value is not None else None
    if reading["cash"] is None:
        raise RuntimeError("balance() returned no cash figure")
    return reading


class BalanceBook:
    """The latest balance per venue, shared by the poll thread and the HTTP handlers."""

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._by_venue: dict[str, dict] = {}

    def record_success(self, venue: str, reading: dict) -> None:
        with self._lock:
            self._by_venue[venue] = {**reading, "fetchedAt": _iso_now(), "error": None}

    def record_failure(self, venue: str, error: str) -> None:
        """Keep the last good numbers (and their fetchedAt) so the card can show
        how stale they are; a venue never read keeps cash None."""
        with self._lock:
            held = self._by_venue.get(venue) or {field: None for field in BALANCE_FIELDS} | {"fetchedAt": None}
            self._by_venue[venue] = {**held, "error": error, "failedAt": _iso_now()}

    def snapshot(self) -> dict[str, dict]:
        with self._lock:
            return {venue: dict(reading) for venue, reading in self._by_venue.items()}


def balance_reads_of(source: object) -> dict[str, Callable[[], dict]]:
    """{venue: read} for a source: its balance_reads() when it has one, else
    {source.name: balance} when it has balance(), else nothing."""
    reads = getattr(source, "balance_reads", None)
    if reads is not None:
        return reads()
    read = getattr(source, "balance", None)
    return {source.name: read} if read is not None else {}


def read_balance_once(source: object, book: BalanceBook) -> None:
    """One read of each of the source's balances; never raises."""
    for venue, read in balance_reads_of(source).items():
        try:
            book.record_success(venue, validate_reading(read()))
        except Exception as error:  # noqa: BLE001 — a failed balance is shown, never fatal
            log.warning("balance %s failed: %s: %s", venue, type(error).__name__, error)
            book.record_failure(venue, f"{type(error).__name__}: {error}")


def balances_payload(book: BalanceBook) -> dict:
    return {"generatedAt": _iso_now(), "balances": book.snapshot()}
