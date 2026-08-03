"""Deterministic golden baselines for all six SGP orchestrators (issue #42).

NOT a test module — imported by ``test_golden_baselines.py`` and used by
``capture_goldens.py``.

Why this replaces the old baselines
-----------------------------------
``test_sgp_regression.py`` shipped two "golden" tests for DK and FD that
shelled out to the real scraper against the LIVE book API, wrote into the
shared production DuckDB (``Answer Keys/mlb_mm.duckdb``), and diffed the
result against a CSV keyed by ``game_id``. Every one of those game_ids expires
within hours of capture, so the tests failed with "no keys in common" and
regenerating the CSV would buy exactly one day. They also could not run
offline, and they mutated production data from a test run.

This harness inverts that: capture each book's raw payloads ONCE into a
frozen fixture, then replay them offline forever. Same idea as the per-book
fixture tests, but pinned at the level the old goldens cared about — the
orchestrator's ``PricedRow`` output.

How it works
------------
Each book declares SEAMS: the (module, attribute) pairs its ``price_sgps``
reaches the network through, plus a ``key`` function turning the call's
arguments into a stable string. In capture mode the real function runs and
its return value is recorded under that key; in replay mode the recorded
value is served and the network is never touched.

A replay miss is a LOUD failure, never a silent empty result — that is the
whole lesson of this epic (issue #33): a book returning nothing must never
look like a book with nothing to say.
"""
from __future__ import annotations

import datetime as _dt
import gzip
import json
from dataclasses import asdict, is_dataclass
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

GOLDEN_DIR = Path(__file__).resolve().parent / "golden"

# Rows are compared on these fields. fetch_time is excluded deliberately: it
# is wall-clock and would differ on every run.
ROW_FIELDS = ("game_id", "combo", "period", "spread_line", "total_line",
              "bookmaker", "source", "sgp_decimal", "sgp_american")


# --------------------------------------------------------------------- #
# JSON round-tripping for the payload shapes the seams return            #
# --------------------------------------------------------------------- #

def encode(obj):
    """Make a seam's return value JSON-safe, preserving enough type
    information for ``decode`` to rebuild what the orchestrator expects.

    Tuples matter: FanDuel's runner dict maps to ``(market_id, selection_id)``
    tuples and BetMGM's legs are 3-tuples, both of which JSON would flatten
    into lists and break identity comparisons downstream.
    """
    if is_dataclass(obj) and not isinstance(obj, type):
        cls = type(obj)
        # MODULE-qualified: four of the six book clients define their own
        # dataclass named `Event`, so a bare class name silently rebuilds a
        # Novig Event from ProphetX fields.
        return {"__dataclass__": f"{cls.__module__}.{cls.__qualname__}",
                "fields": {k: encode(v) for k, v in asdict(obj).items()}}
    if isinstance(obj, SimpleNamespace):
        return {"__namespace__": {k: encode(v) for k, v in vars(obj).items()}}
    if isinstance(obj, tuple):
        return {"__tuple__": [encode(v) for v in obj]}
    if isinstance(obj, set):
        return {"__set__": [encode(v) for v in obj]}
    if isinstance(obj, dict):
        # Keys can be tuples (FD runners) or floats (CZR/MGM line maps), so
        # they are encoded alongside the values rather than as JSON keys.
        return {"__dict__": [[encode(k), encode(v)] for k, v in obj.items()]}
    if isinstance(obj, list):
        return [encode(v) for v in obj]
    if isinstance(obj, (_dt.datetime, _dt.date)):
        # Without this a datetime falls through to json.dump(default=str),
        # which stringifies it with a SPACE separator and never becomes a
        # datetime again on replay — a silent type drift in a repo whose
        # timestamp handling is a regression gate (CLAUDE.md housekeeping #6).
        return {"__datetime__": obj.isoformat()}
    return obj


def _resolve_dataclass(path: str):
    """Import a dataclass from its module-qualified name, or None."""
    import importlib
    module_name, _, cls_name = path.rpartition(".")
    if not module_name:
        return None
    try:
        return getattr(importlib.import_module(module_name), cls_name, None)
    except ImportError:
        return None


def decode(obj):
    if isinstance(obj, dict):
        if "__dataclass__" in obj:
            fields = {k: decode(v) for k, v in obj["fields"].items()}
            cls = _resolve_dataclass(obj["__dataclass__"])
            return cls(**fields) if cls else SimpleNamespace(**fields)
        if "__namespace__" in obj:
            return SimpleNamespace(
                **{k: decode(v) for k, v in obj["__namespace__"].items()})
        if "__tuple__" in obj:
            return tuple(decode(v) for v in obj["__tuple__"])
        if "__set__" in obj:
            return {decode(v) for v in obj["__set__"]}
        if "__datetime__" in obj:
            return _dt.datetime.fromisoformat(obj["__datetime__"])
        if "__dict__" in obj:
            return {decode(k): decode(v) for k, v in obj["__dict__"]}
        return {k: decode(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [decode(v) for v in obj]
    return obj


# Argument names that vary run to run without selecting a different payload.
_NOISE_KWARGS = ("verbose", "profile", "on_decline", "counters")


def call_key(args, kwargs, skip: int) -> str:
    """Canonical key for one seam call.

    ``skip`` drops the leading positional arguments that are connections
    rather than selectors — a curl session or a client object. They are
    identity-only, so including them would make every replay a miss.
    """
    payload = [encode(a) for a in args[skip:]]
    extras = {k: encode(v) for k, v in sorted(kwargs.items())
              if k not in _NOISE_KWARGS}
    return json.dumps([payload, extras], sort_keys=True, default=str)


# --------------------------------------------------------------------- #
# Per-book seam declarations                                             #
# --------------------------------------------------------------------- #
# module_patches:  (module, attribute, leading args to ignore)
# client_methods:  {method: leading args to ignore}   (self is already bound)
# fetchers:        {name: leading args to ignore}     passed to price_sgps
#
# Names verified against each orchestrator's own call sites rather than
# assumed — e.g. ProphetX prices through client.submit_parlay_rfq and
# Caesars through client.price_combo, neither of which is a module function.

def _seams(book: str) -> dict:
    if book == "novig":
        return {
            "module_patches": [
                # (session, game, verbose) -> keyed on the game dict
                ("scraper_novig_sgp", "fetch_event_legs", 1),
                # (session, outcome_ids, verbose) -> keyed on the ids
                ("scraper_novig_sgp", "submit_parlay", 1),
            ],
            "client_methods": {"list_events": 0},
        }
    if book == "draftkings":
        return {
            "module_patches": [("scraper_draftkings_sgp", "calculate_sgp", 1)],
            "fetchers": {"fetch_dk_events": 1,
                         "fetch_main_market_nums": 1,
                         "fetch_selection_ids": 1},
        }
    if book == "fanduel":
        return {
            "module_patches": [("scraper_fanduel_sgp", "price_combo", 1)],
            "fetchers": {"fetch_fd_events": 1, "fetch_event_runners": 1},
        }
    if book == "prophetx":
        return {"client_methods": {"list_events": 0,
                                   "fetch_event_markets": 0,
                                   "submit_parlay_rfq": 0}}
    if book == "betmgm":
        # parse_markets runs for real against the frozen fetch_markets
        # payload, so the parser stays under test.
        return {"client_methods": {"accessid": 0, "list_events": 0,
                                   "fetch_markets": 0, "price_picks": 0}}
    if book == "caesars":
        return {"client_methods": {"list_events": 0, "fetch_event": 0,
                                   "price_combo": 0}}
    raise ValueError(f"unknown book {book!r}")


def rows_to_json(rows) -> list:
    return [{f: getattr(r, f) for f in ROW_FIELDS} for r in rows]


def sort_rows(rows: list) -> list:
    """Thread pools return rows in nondeterministic order."""
    return sorted(rows, key=lambda r: (str(r["game_id"]), str(r["combo"]),
                                       str(r["period"]),
                                       str(r["spread_line"]),
                                       str(r["total_line"])))


def golden_path(book: str) -> Path:
    # gzipped: a raw capture is up to 4.5 MB of market tree per book (7.3 MB
    # across six), and every recapture would add another blob. These compress
    # ~35x, so the whole set is ~200 KB.
    return GOLDEN_DIR / f"{book}_golden.json.gz"


def load_golden(book: str) -> dict | None:
    path = golden_path(book)
    if not path.exists():
        return None
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)


def save_golden(book: str, blob: dict) -> Path:
    path = golden_path(book)
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        json.dump(blob, handle, indent=1, sort_keys=True, default=str)
    return path


class ReplayMiss(AssertionError):
    """A seam asked for a payload the fixture does not contain."""


class _Replayer:
    """Serves recorded returns for one seam. A miss raises rather than
    returning an empty result — see the module docstring."""

    def __init__(self, book: str, seam: str, skip: int, recorded: dict):
        self.book, self.seam, self.skip = book, seam, skip
        self.recorded = recorded
        self.calls: list[str] = []

    def __call__(self, *args, **kwargs):
        key = call_key(args, kwargs, self.skip)
        self.calls.append(key)
        if key not in self.recorded:
            raise ReplayMiss(
                f"{self.book}/{self.seam}: no recorded payload for {key}. "
                f"The fixture holds {len(self.recorded)} key(s); recapture "
                f"with mlb_sgp/tests/capture_goldens.py if the orchestrator's "
                f"call pattern changed on purpose.")
        return decode(self.recorded[key])


# Seams whose captured value is auth material rather than market data. It is
# never read on replay (every downstream call is replayed too), so a constant
# keeps the fixture honest without committing a live credential.
_REDACT_SEAMS = {"accessid"}
REDACTED = "REDACTED-BY-CAPTURE"


class _Recorder:
    """Wraps a real seam, remembering every (key -> return) it produced."""

    def __init__(self, real, skip: int, sink: dict, seam: str = ""):
        self.real, self.skip, self.sink = real, skip, sink
        self.seam = seam

    def __call__(self, *args, **kwargs):
        result = self.real(*args, **kwargs)
        stored = REDACTED if self.seam in _REDACT_SEAMS else encode(result)
        try:
            self.sink[call_key(args, kwargs, self.skip)] = stored
        except (TypeError, ValueError):
            # An unencodable payload must not break a live capture. The
            # replay would fail loudly on the missing key, and
            # verify_capture() catches it before the fixture is written.
            pass
        return result


def replay(book: str, monkeypatch, recorded: dict):
    """Patch every seam of ``book`` to serve ``recorded``.

    Returns ``(client, fetchers)`` to hand to that book's ``price_sgps``.
    Nothing here touches the network.
    """
    import importlib

    seams = _seams(book)

    for module_name, attr, skip in seams.get("module_patches", []):
        module = importlib.import_module(module_name)
        monkeypatch.setattr(module, attr,
                            _Replayer(book, attr, skip,
                                      recorded.get(attr, {})))

    client = MagicMock()
    for method, skip in seams.get("client_methods", {}).items():
        setattr(client, method,
                _Replayer(book, method, skip, recorded.get(method, {})))

    fetchers = None
    if "fetchers" in seams:
        fetchers = {name: _Replayer(book, name, skip,
                                    recorded.get(name, {}))
                    for name, skip in seams["fetchers"].items()}
    return client, fetchers


def record(book: str, monkeypatch, client, sink: dict):
    """Wrap every seam of ``book`` so a LIVE run fills ``sink``.

    Returns the fetchers dict (or None). The caller supplies a real client;
    its methods are wrapped in place.
    """
    import importlib

    seams = _seams(book)

    for module_name, attr, skip in seams.get("module_patches", []):
        module = importlib.import_module(module_name)
        monkeypatch.setattr(module, attr,
                            _Recorder(getattr(module, attr), skip,
                                      sink.setdefault(attr, {}), attr))

    for method, skip in seams.get("client_methods", {}).items():
        monkeypatch.setattr(client, method,
                            _Recorder(getattr(client, method), skip,
                                      sink.setdefault(method, {}), method))

    fetchers = None
    if "fetchers" in seams:
        import importlib as _il
        legacy = _il.import_module(_LEGACY_MODULE[book])
        fetchers = {name: _Recorder(getattr(legacy, name), skip,
                                    sink.setdefault(name, {}), name)
                    for name, skip in seams["fetchers"].items()}
    return fetchers


_LEGACY_MODULE = {
    "draftkings": "scraper_draftkings_sgp",
    "fanduel": "scraper_fanduel_sgp",
    "novig": "scraper_novig_sgp",
    "prophetx": "scraper_prophetx_sgp",
}


def price_with(book: str, targets, client, fetchers, periods=("FG",)):
    """Call one book's ``price_sgps`` with the right signature.

    DK and FD take a ``fetchers`` dict and a ``parallelism``; the other four
    do not. Keeping that knowledge here is what lets the golden test body be
    identical for all six books.
    """
    import importlib

    module = importlib.import_module(f"mlb_sgp.{book}")
    if fetchers is not None:
        return module.price_sgps(targets, periods=periods, client=client,
                                 fetchers=fetchers)
    return module.price_sgps(targets, periods=periods, client=client)
