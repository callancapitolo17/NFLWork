"""``price_sgps(..., counters=...)`` — the sweep's counter seam (issue #38).

#35 scoped every sweep with its own ``FetchCounters``, created and discarded
inside ``price_sgps``. #38 needs those numbers on a health row, and the only
caller that can persist them is ``SGPService`` — which sees nothing but
``list[PricedRow] | None``.

So the object becomes an OPTIONAL, KEYWORD-ONLY parameter: the service makes
it, hands it down, and reads it afterwards. Keyword-only is not stylistic —
DK/FD/PX take ``parallelism``/``fetchers`` positionally-capable, and an
argument that silently bound to the wrong slot is the exact failure this
epic already nearly shipped once.

The drives in ``sweep_drives`` run each book's REAL orchestration, so these
tests prove the counters are actually filled, not merely accepted.
"""
from __future__ import annotations

import inspect
import logging

import pytest

from mlb_sgp import betmgm, caesars, draftkings, fanduel, novig, prophetx
from mlb_sgp._shared import FetchCounters, FetchCountersSnapshot
from mlb_sgp.tests.sweep_drives import DRIVES

BOOK_MODULES = {
    "draftkings": draftkings, "fanduel": fanduel, "prophetx": prophetx,
    "novig": novig, "betmgm": betmgm, "caesars": caesars,
}


def test_book_name_matches_the_service_key():
    """``fetch_counters`` REJECTS a counters object whose (book, path)
    disagrees with the fetch it is scoping — a mislabeled health row is
    worse than none. That guard turns a BOOK_NAME rename into an exception
    inside every sweep, swallowed as outcome='error'. Pin the invariant."""
    from kalshi_common.sgp_service import DEFAULT_BOOKS

    assert set(DEFAULT_BOOKS) == set(BOOK_MODULES)
    for key, mod in BOOK_MODULES.items():
        assert mod.BOOK_NAME == key, (
            f"{key}: BOOK_NAME={mod.BOOK_NAME!r} would break the sweep")


def test_mismatched_counters_are_rejected():
    """Proof the guard above is real, not decorative."""
    from mlb_sgp._shared import fetch_counters

    wrong = FetchCounters("novig", "on_demand")
    with pytest.raises(ValueError, match="novig"):
        with fetch_counters("caesars", "sweep", logging.getLogger("t"),
                            counters=wrong):
            pass


@pytest.mark.parametrize("book", sorted(BOOK_MODULES))
def test_price_sgps_counters_is_keyword_only(book):
    param = inspect.signature(BOOK_MODULES[book].price_sgps).parameters["counters"]
    assert param.kind is inspect.Parameter.KEYWORD_ONLY
    assert param.default is None


@pytest.mark.parametrize("book", sorted(BOOK_MODULES))
def test_shim_call_shape_still_binds(book):
    """The dashboard shims call price_sgps(targets, periods=, client=,
    verbose=). Adding the seam must not disturb that."""
    sig = inspect.signature(BOOK_MODULES[book].price_sgps)
    bound = sig.bind(["t"], periods=("FG",), client=object(), verbose=False)
    assert "counters" not in bound.arguments


@pytest.mark.parametrize("book", sorted(DRIVES))
def test_caller_supplied_counters_are_filled_by_the_real_sweep(book,
                                                               monkeypatch):
    counters = FetchCounters(book, "sweep")
    rows = DRIVES[book](monkeypatch, healthy=True, counters=counters)
    assert rows, f"{book}: drive produced no rows — fixture drifted"

    snap = counters.snapshot()
    assert isinstance(snap, FetchCountersSnapshot)
    assert snap.book == book and snap.path == "sweep"
    assert snap.legs_attempted > 0, f"{book}: nothing counted into MY object"
    assert snap.prices_returned > 0


@pytest.mark.parametrize("book", sorted(DRIVES))
def test_broken_parse_trips_the_wire_on_the_caller_object(book, monkeypatch):
    """The tripwire must be readable from the object the CALLER owns —
    that is what ends up in the health row's counters_json."""
    counters = FetchCounters(book, "sweep")
    DRIVES[book](monkeypatch, healthy=False, counters=counters)
    assert counters.tripwires(), f"{book}: no tripwire on the caller's object"


@pytest.mark.parametrize("book", sorted(DRIVES))
def test_omitting_counters_still_works(book, monkeypatch):
    """Every pre-#38 caller (dashboard shims, CLI) passes nothing."""
    assert DRIVES[book](monkeypatch, healthy=True)
