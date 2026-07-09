from abc import ABC, abstractmethod
from dataclasses import dataclass


@dataclass
class Candidate:
    market_ticker: str
    side: str           # "yes" | "no"
    fair_prob: float
    label: str          # e.g. "over_2.5", "under_2.5"
    meta: dict | None = None   # provenance for the research firehose (book, alt, overround)


class SportAdapter(ABC):
    sport: str
    league_prefix: str      # Unabated league key prefix, e.g. "lg21"

    @property
    def league_id(self) -> int:
        """Numeric Unabated league id (drives the v2 odds URL), e.g. "lg21" -> 21."""
        return int(self.league_prefix[2:])

    @abstractmethod
    def canon_team(self, name: str) -> str: ...

    @abstractmethod
    def kalshi_series(self) -> str: ...

    @abstractmethod
    def event_teams(self, kalshi_event: dict) -> frozenset:
        """Canon team pair parsed from the Kalshi event title, for Unabated↔Kalshi matching."""
        ...

    @abstractmethod
    def price_event(self, state, event_meta, kalshi_event) -> list[Candidate]:
        """Return priced Candidates for the given event pair. Return [] to fail closed."""
        ...
