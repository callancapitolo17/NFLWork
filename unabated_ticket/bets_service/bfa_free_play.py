"""BFA free play for the Place buttons (Cal 2026-11-09): a straight or a teaser bet with the
account's free play instead of its cash. A free play wager is the cash wager with UseFreePlay
true on every pick; what BFA's site code and account reads showed is in the README (BFA free
play on Place).

Every check here is a refusal before anything is sent.
Inputs: the account's metadata and the BFA source's bearer headers. Outputs: the balance, the
rules for one bet type, and why a free play would be refused. Side effects: one GET (the
balance); nothing is written.
"""
import threading
import time
from collections.abc import Callable


from unabated_ticket.bets_service.sources.bfa import BALANCE_URL, HTTP_TIMEOUT_SEC, is_number

# How long GET /bfa_free_play.json serves one read (the Ticket opens often); a Place reads afresh.
BALANCE_CACHE_SEC = 30.0


def read_free_play_balance(session: object, headers: dict, player_id: str) -> float:
    """The account's free play balance in dollars: freePlayAmount off BALANCE_URL. Raises on
    anything but a 200 carrying a number there."""
    response = session.get(BALANCE_URL, params={"playerId": player_id}, headers=headers, timeout=HTTP_TIMEOUT_SEC)
    if response.status_code != 200:
        raise RuntimeError(f"BFA balance: HTTP {response.status_code}")
    body = response.json()
    amount = body.get("freePlayAmount") if isinstance(body, dict) else None
    if not is_number(amount) or amount < 0:
        raise RuntimeError(f"BFA balance: expected a freePlayAmount of 0 or more, got {amount!r}")
    return float(amount)


def free_play_rules_of(metadata: dict, wager_type_id: int, sport_codes: set[str]) -> dict | str:
    """{oddsLimit, maxPayout} for a free play on the account's bet type `wager_type_id` in
    every sport of `sport_codes`, or why free play is not offered there. Every metadata
    betType row of that type and those sports must offer free play; the tightest limits win."""
    rows = [row for row in metadata.get("betType") or [] if isinstance(row, dict)
            and row.get("wagerTypeId") == wager_type_id and row.get("sportId") in sport_codes]
    if {row.get("sportId") for row in rows} != sport_codes:
        return f"BFA's account metadata lists bet type {wager_type_id} for {sorted({row.get('sportId') for row in rows})}, not {sorted(sport_codes)}"
    if not all(row.get("freePlay") is True for row in rows):
        return "BFA does not offer free play on this bet type"
    odds_limits = [row.get("fpOddsLimit") for row in rows]
    max_payouts = [row.get("fpMaxPayout") for row in rows]
    # An American limit may be negative (favorites only); a payout cap must be above 0.
    if not (all(is_number(value) and value != 0 for value in odds_limits)
            and all(is_number(value) and value > 0 for value in max_payouts)):
        return f"BFA's account metadata carries no free play limits for this bet type ({odds_limits}, {max_payouts})"
    return {"oddsLimit": min(odds_limits), "maxPayout": float(min(max_payouts))}


def free_play_problem(rules: dict, balance: float, stake: int, to_win: float,
                      straight_price: float | None) -> str | None:
    """Why this free play cannot go out, or None. `straight_price` is the straight's American
    price (None on a teaser: the site checks no odds limit on teaser legs)."""
    if stake > balance + 1e-9:
        return f"${stake} is more than your ${balance:,.2f} of BFA free play (BFA never covers the rest with cash)"
    if straight_price is not None and straight_price > rules["oddsLimit"]:
        return f"BFA takes free play at {rules['oddsLimit']:+g} or shorter, this is {straight_price:+g}"
    if to_win > rules["maxPayout"] + 0.005:
        return f"a free play may win at most ${rules['maxPayout']:,.0f} at BFA, ${stake} wins ${to_win:,.2f}"
    return None


def free_play_refusal(session: object, headers: dict, player_id: str, metadata: dict, wager_type_id: int,
                      sport_codes: set[str], stake: int, to_win: float, straight_price: float | None) -> str | None:
    """Why this bet cannot go out on free play, or None: the type's rules, then the balance read
    NOW (free play spent since the Ticket last read it must not go out twice). Raises when the
    balance cannot be read: the caller refuses, nothing is sent."""
    rules = free_play_rules_of(metadata, wager_type_id, sport_codes)
    if isinstance(rules, str):
        return rules
    balance = read_free_play_balance(session, headers, player_id)
    return free_play_problem(rules, balance, stake, to_win, straight_price)


class FreePlayReader:
    """GET /bfa_free_play.json: the free play balance, held BALANCE_CACHE_SEC, and the tightest
    limits across the bet types that offer free play (the Ticket's default stake; each Place
    checks its own type afresh).

    `session_factory()` builds an HTTP session (tests inject a fake); `read_metadata(headers,
    player_id)` is the placer's cached metadata read; `clock()` paces the cache."""

    def __init__(self, bfa_source: object, session_factory: Callable[[], object],
                 read_metadata: Callable[[dict, str], dict], clock: Callable[[], float] = time.monotonic):
        self._bfa = bfa_source
        self._session_factory = session_factory
        self._read_metadata = read_metadata
        self._clock = clock
        self._lock = threading.Lock()
        self._held: dict | None = None
        self._held_at = 0.0

    def read(self) -> dict:
        """{ok: true, freePlay, maxPayout, oddsLimit} or {ok: false, error}. Never raises."""
        with self._lock:
            if self._held is not None and self._clock() - self._held_at < BALANCE_CACHE_SEC:
                return self._held
            try:
                headers, player_id = self._bfa.authorized()
                balance = read_free_play_balance(self._session_factory(), headers, player_id)
                metadata = self._read_metadata(headers, player_id)
            except Exception as error:  # noqa: BLE001 — the Ticket shows the read failing; Place re-reads
                return {"ok": False, "error": f"{type(error).__name__}: {error}"}
            rows = [row for row in metadata.get("betType") or [] if isinstance(row, dict) and row.get("freePlay") is True]
            payouts = [row["fpMaxPayout"] for row in rows if is_number(row.get("fpMaxPayout")) and row["fpMaxPayout"] > 0]
            limits = [row["fpOddsLimit"] for row in rows if is_number(row.get("fpOddsLimit")) and row["fpOddsLimit"] > 0]
            self._held = {"ok": True, "freePlay": balance, "maxPayout": float(min(payouts)) if payouts else None,
                          "oddsLimit": min(limits) if limits else None}
            self._held_at = self._clock()
            return self._held

    def forget(self) -> None:
        """Drop the held read: a free play just went out, so the next read asks BFA."""
        with self._lock:
            self._held = None
