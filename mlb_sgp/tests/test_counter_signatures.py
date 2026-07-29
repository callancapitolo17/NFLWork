"""Signature guards for the `counters` parameter (issue #35).

Every function this ticket touched takes ``counters`` KEYWORD-ONLY. That is
not stylistic: BetMGM's ``price_selection_set(client, refs, fixture_id)``
takes a positional 3rd argument, and an earlier epic ticket nearly bound a
positional arg to the wrong parameter. Keyword-only makes that class of bug
impossible, and ``inspect.signature().bind()`` proves it here rather than by
eye.
"""
from __future__ import annotations

import inspect

import pytest

from mlb_sgp import betmgm, caesars, draftkings, fanduel, novig, prophetx
from mlb_sgp._shared import FetchCounters

BOOK_MODULES = {
    "draftkings": draftkings, "fanduel": fanduel, "prophetx": prophetx,
    "novig": novig, "betmgm": betmgm, "caesars": caesars,
}


@pytest.mark.parametrize("book", sorted(BOOK_MODULES))
@pytest.mark.parametrize("fname", ["resolve_legs", "price_selection_set"])
def test_counters_is_keyword_only(book, fname):
    fn = getattr(BOOK_MODULES[book], fname)
    param = inspect.signature(fn).parameters["counters"]
    assert param.kind is inspect.Parameter.KEYWORD_ONLY, (
        f"{book}.{fname}: counters must be keyword-only")
    assert param.default is None, (
        f"{book}.{fname}: counters must default to None")


@pytest.mark.parametrize("book", sorted(BOOK_MODULES))
def test_resolve_legs_positional_call_still_binds(book):
    """The pre-#35 4-positional-arg call site must still be valid."""
    fn = BOOK_MODULES[book].resolve_legs
    bound = inspect.signature(fn).bind({"structure": 1}, [], "home", "away")
    assert "counters" not in bound.arguments


def test_betmgm_price_selection_set_keeps_fixture_id_positional():
    """The exact hazard: a 3rd positional arg must bind to ``fixture_id``,
    never to ``counters``."""
    sig = inspect.signature(betmgm.price_selection_set)
    bound = sig.bind("client", ["r0"], "fixture-123")
    assert bound.arguments["fixture_id"] == "fixture-123"
    assert "counters" not in bound.arguments

    counters = FetchCounters("betmgm", "on_demand")
    bound2 = sig.bind("client", ["r0"], "fixture-123", counters=counters)
    assert bound2.arguments["counters"] is counters


@pytest.mark.parametrize("book", ["draftkings", "fanduel", "prophetx",
                                  "novig", "caesars"])
def test_non_mgm_price_selection_set_binds_two_positionals(book):
    sig = inspect.signature(BOOK_MODULES[book].price_selection_set)
    bound = sig.bind("client", ["r0"])
    assert "counters" not in bound.arguments


@pytest.mark.parametrize("book", sorted(BOOK_MODULES))
def test_price_sgps_public_signature_is_unchanged(book):
    """The sweep wrapper must keep the shim/dashboard call shape: every
    parameter after ``target_lines`` still has a default, and the names the
    shims pass by keyword still exist."""
    sig = inspect.signature(BOOK_MODULES[book].price_sgps)
    names = list(sig.parameters)
    assert names[0] == "target_lines"
    for name in names[1:]:
        assert sig.parameters[name].default is not inspect.Parameter.empty, (
            f"{book}.price_sgps: {name} lost its default")
    # Shim call shape: price_sgps(targets, periods=..., client=..., verbose=...)
    sig.bind(["t"], periods=("FG",), client=object(), verbose=False)
    assert "counters" not in sig.parameters, (
        f"{book}.price_sgps must not leak the counters seam to callers")


@pytest.mark.parametrize("book", sorted(BOOK_MODULES))
def test_ml_total_helper_takes_counters_keyword_only(book):
    """Novig prices ML×total inline; the other five use a helper."""
    fn = getattr(BOOK_MODULES[book], "_price_ml_total_for_games", None)
    if fn is None:
        pytest.skip(f"{book} has no _price_ml_total_for_games helper")
    param = inspect.signature(fn).parameters["counters"]
    assert param.kind is inspect.Parameter.KEYWORD_ONLY
    assert param.default is None


def test_sgp_service_hook_factory_counters_is_keyword_only():
    from kalshi_common.sgp_service import SGPService

    param = inspect.signature(
        SGPService._book_on_demand_hooks).parameters["counters"]
    assert param.kind is inspect.Parameter.KEYWORD_ONLY
    assert param.default is None
