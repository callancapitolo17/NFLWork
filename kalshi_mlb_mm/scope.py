"""Decide whether an incoming RFQ's combo is one we can price.

Scope is now a cheap structural gate: every leg must be a supported MLB market
type (spread / total / moneyline). It does NOT require a specific shape or count
— cross-game and N-leg combos are in scope. Whether a particular combo is
actually priceable (enough agreeing books, a known same-game shape) is decided
later by `combo_pricer.combo_fair`, which returns no_fair when it can't price.

Player props (KXMLBHR / KXMLBHIT / KXMLBTB / KXMLBKS / KXMLBHRR) and all
non-MLB legs are out of scope here — a single unsupported leg drops the combo,
since the maker can only quote a combo whose EVERY leg it can price.
"""

# Market-ticker prefixes the maker can price (Phase 0-3). Props are added in a
# later phase; until then a prop leg makes the whole combo out-of-scope.
SUPPORTED_PREFIXES = ("KXMLBSPREAD-", "KXMLBTOTAL-", "KXMLBGAME-")


def decode_legs(market: dict) -> list[dict] | None:
    """Pull the constituent legs from a GET /markets/{ticker} 'market' dict.
    The API returns them under mve_selected_legs in the exact
    {event_ticker, market_ticker, side} shape the leg-typer consumes."""
    legs = market.get("mve_selected_legs")
    if not legs or not isinstance(legs, list):
        return None
    return legs


def is_in_scope(legs: list[dict]) -> bool:
    """True iff every leg is a supported MLB market type (spread/total/moneyline).
    Empty leg lists are out of scope."""
    if not legs:
        return False
    return all(
        str(l.get("market_ticker", "")).startswith(SUPPORTED_PREFIXES)
        for l in legs
    )


def is_spread_total_2leg(legs: list[dict]) -> bool:
    """True iff exactly one KXMLBSPREAD- leg and one KXMLBTOTAL- leg.
    Retained for the existing same-game-shape checks/tests; `is_in_scope` is the
    general gate the discovery loop now uses."""
    if len(legs) != 2:
        return False
    tickers = [str(l.get("market_ticker", "")) for l in legs]
    has_spread = sum(t.startswith("KXMLBSPREAD-") for t in tickers) == 1
    has_total = sum(t.startswith("KXMLBTOTAL-") for t in tickers) == 1
    return has_spread and has_total
