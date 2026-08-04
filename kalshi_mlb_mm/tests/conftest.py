"""Shared fixture helpers for the maker's tests.

`kalshi_mlb_mm/tests` has no other conftest; this exists so the Kalshi leg-dict
shape is built in ONE place. Issue #71 was hidden for months by fixtures that
gave a game's spread, total and moneyline legs a single shared `event_ticker` —
input Kalshi never emits — so the same-game grid tests passed while the grids
were unreachable in production. Building leg dicts through `leg()` keeps every
fixture on the real shape by construction.
"""


def event_ticker_for(market_ticker: str) -> str:
    """The family-specific event a Kalshi market ticker belongs to.

    Kalshi nests one market under one event, and each market FAMILY gets its own
    event, so the event ticker is the market ticker minus its final segment:

        KXMLBSPREAD-25JUN271905NYYBOS-BOS2  ->  KXMLBSPREAD-25JUN271905NYYBOS
        KXMLBTOTAL-25JUN271905NYYBOS-9      ->  KXMLBTOTAL-25JUN271905NYYBOS

    Both reduce to the game code `25JUN271905NYYBOS` via `legset.game_id_of` —
    which is what `partition_by_game` keys on.
    """
    return market_ticker.rsplit("-", 1)[0]


def leg(market_ticker: str, side: str = "yes") -> dict:
    """One Kalshi leg dict, with the event ticker derived from the market."""
    return {"market_ticker": market_ticker,
            "event_ticker": event_ticker_for(market_ticker),
            "side": side}
