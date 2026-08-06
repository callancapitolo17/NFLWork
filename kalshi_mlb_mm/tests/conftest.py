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


class FakeLiveEngine:
    """Duck-typed OnDemandEngine serving canned live fairs for EVERY hash.

    Post-#54 the quote path is live-only: every discovery/confirm test that
    reaches pricing needs an engine, even when the test is about something
    else entirely (caps, cooldowns, hysteresis, corr-sanity). Install this
    via `monkeypatch.setattr(main, "_ENGINE", FakeLiveEngine(...))` so those
    tests exercise their own subject instead of stalling on_demand_pending.
    Tests about the engine/feed itself build their own stubs instead.
    """

    def __init__(self, fairs=None):
        self.fairs = dict(fairs) if fairs else {"draftkings": 0.55,
                                                "fanduel": 0.55}
        self.ensure_calls = []
        self.refetch_calls = []

    def lookup(self, h):
        return dict(self.fairs)

    def lookup_results(self, h):
        from mlb_sgp._shared import OnDemandBookResult
        return {b: OnDemandBookResult(book=b, fair=f, route="partition",
                                      n_cells_priced=4, latency_sec=1.0)
                for b, f in self.fairs.items()}

    def landed_at(self, h):
        return 1000.0

    def landed_empty(self, h):
        return False

    def result_age_sec(self, h):
        return 1.0

    def ensure_fetch(self, h, game, legs):
        self.ensure_calls.append((h, game, tuple(legs)))
        return True

    def refetch_now(self, jobs, deadline_sec):
        self.refetch_calls.append((list(jobs), deadline_sec))
        return True
