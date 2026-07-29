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


@pytest.mark.parametrize("book", sorted(BOOK_MODULES))
def test_wrapper_delegates_positionally_in_the_right_order(book):
    """``price_sgps`` forwards to ``_price_sgps`` POSITIONALLY. A parameter
    reordered on one side only would silently swap e.g. `verbose` and
    `parallelism`, so check the names line up rather than trusting the eye."""
    import ast

    mod = BOOK_MODULES[book]
    src = inspect.getsource(mod.price_sgps).lstrip()
    call = next(n for n in ast.walk(ast.parse(src))
                if isinstance(n, ast.Call)
                and getattr(n.func, "id", "") == "_price_sgps")
    passed = [getattr(a, "id", "<expr>") for a in call.args]
    assert not call.keywords, f"{book}: delegation should be all-positional"
    assert passed == list(inspect.signature(mod._price_sgps).parameters), (
        f"{book}: wrapper args do not line up with _price_sgps parameters")


def test_concurrent_fetches_do_not_share_counters():
    """The sweep and an on-demand call for the SAME book run concurrently by
    design (price_on_demand runs on the caller's thread, alongside refresh()).
    Each fetch scope must own its counters, or one fetch's numbers would bleed
    into the other's and both tripwire verdicts would be wrong."""
    import logging
    from concurrent.futures import ThreadPoolExecutor

    from mlb_sgp._shared import fetch_counters

    logger = logging.getLogger("test.concurrent_counters")
    logger.addHandler(logging.NullHandler())
    seen = []

    def one_fetch(n):
        with fetch_counters("caesars", "sweep", logger) as c:
            seen.append(c)          # keep a REFERENCE: id() alone would be
            c.bump("legs_attempted", n)   # ambiguous once CPython recycles
            c.bump("legs_resolved", n)    # a freed object's id.
            return c.snapshot()

    with ThreadPoolExecutor(max_workers=8) as pool:
        snaps = list(pool.map(one_fetch, range(1, 9)))

    assert len({id(c) for c in seen}) == 8, "fetch scopes shared a counters object"
    # Each scope sees only its OWN bumps, never the other seven's.
    assert sorted(s.legs_attempted for s in snaps) == list(range(1, 9))
