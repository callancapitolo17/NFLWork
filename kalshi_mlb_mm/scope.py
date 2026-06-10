"""Decide whether an incoming RFQ's combo is one we model (2-leg spread×total)."""


def decode_legs(market: dict) -> list[dict] | None:
    """Pull the constituent legs from a GET /markets/{ticker} 'market' dict.
    The API returns them under mve_selected_legs in the exact
    {event_ticker, market_ticker, side} shape the leg-typer consumes."""
    legs = market.get("mve_selected_legs")
    if not legs or not isinstance(legs, list):
        return None
    return legs


def is_spread_total_2leg(legs: list[dict]) -> bool:
    """True iff exactly one KXMLBSPREAD- leg and one KXMLBTOTAL- leg."""
    if len(legs) != 2:
        return False
    tickers = [str(l.get("market_ticker", "")) for l in legs]
    has_spread = sum(t.startswith("KXMLBSPREAD-") for t in tickers) == 1
    has_total = sum(t.startswith("KXMLBTOTAL-") for t in tickers) == 1
    return has_spread and has_total
