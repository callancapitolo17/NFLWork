"""Source protocol — the plug-in point for every venue (#115 BetOnline,
#116 Novig, #117 ProphetX add a module here).

A source turns one venue's account data into normalised bet records
(contract in docs/2026-09-11-issue-114-bet-history-plan.md). The service
loop calls `fetch()` every `poll_sec`; a raised exception is a failed run —
the service logs it to `source_runs`, keeps the source's previous records and
moves on. Records must carry a stable `id` ("<venue>:<native id>") because the
store upserts on it.
"""
from typing import Protocol


class Source(Protocol):
    name: str
    poll_sec: float

    def fetch(self) -> list[dict]:
        """Return every record the venue currently knows (open and settled).
        Raise on any failure; never return a partial list as if complete."""
        ...
