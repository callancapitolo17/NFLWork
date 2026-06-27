from abc import ABC, abstractmethod

class SportAdapter(ABC):
    sport: str
    league_prefix: str
    outcomes: tuple

    @abstractmethod
    def fair(self, state, event_meta) -> dict | None: ...
    @abstractmethod
    def canon_team(self, name: str) -> str: ...
    @abstractmethod
    def kalshi_series(self) -> str: ...
    @abstractmethod
    def map_outcome_tickers(self, kalshi_event: dict) -> dict: ...
