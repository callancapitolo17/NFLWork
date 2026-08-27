"""Regression tests for issue #102 — DK's on-demand price 403 was invisible.

CONTEXT. On 2026-08-27 ``POST */api/wager/v1/calculateBets`` began answering
403 AkamaiGHost for every set size, on both wager hosts, with or without the
``/en/`` locale prefix that issue #39 introduced. Recon (see
``mlb_sgp/README.md`` § DraftKings price host) showed DK's OWN betslip fails
identically from this egress, so the request form is not the defect.

THE DEFECT PINNED HERE is the diagnostic one, and it is the same shape issue
#39 fixed for the sweep and issue #40 fixed for Novig — this time on the
on-demand path, which since issue #81 is the ONLY path the maker prices on.

``draftkings.price_selection_set`` is the only book's price hook that issues
its own HTTP call inline; every other book delegates to a client that runs the
response through ``check_response``. It classified nothing::

    if resp.status_code != 200:
        return None          # 422 non-combinable AND 403 blockade, identically

So a blocked DK produced a fetch with ``transport_errors == 0`` and
``prices_returned == 0``, which is precisely the ``prices_empty`` tripwire's
firing condition — a tripwire whose documented meaning (``_shared.py``) is
"priced targets, got nothing back, and transport never complained", i.e. a
PARSER regression. A restart of the maker would therefore have pointed the
next fixer at DK's parser while the real story was a dead endpoint: the exact
misdiagnosis epic #32 exists to remove.

No network: every test injects a fake session.
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

# The DK orchestrator lazily imports its legacy scraper module by TOP-LEVEL
# name, which only resolves with mlb_sgp/ itself on sys.path.
_MLB_SGP_DIR = Path(__file__).resolve().parents[1]
if str(_MLB_SGP_DIR) not in sys.path:
    sys.path.insert(0, str(_MLB_SGP_DIR))

from mlb_sgp import draftkings                                   # noqa: E402
from mlb_sgp._shared import FetchCounters                        # noqa: E402

_LOG = logging.getLogger(__name__)


class _Response:
    def __init__(self, status_code, body=None, text=""):
        self.status_code = status_code
        self._body = body
        self.text = text

    def json(self):
        if self._body is None:
            raise ValueError("Expecting value: line 1 column 1 (char 0)")
        return self._body


class _Session:
    def __init__(self, response):
        self._response = response
        self.posts = 0

    def post(self, url, **kwargs):
        self.posts += 1
        return self._response


class _Client:
    def __init__(self, response):
        self.session = _Session(response)


AKAMAI_DENY = _Response(
    403, text="<HTML><HEAD>\n<TITLE>Access Denied</TITLE>\n</HEAD>")
NON_COMBINABLE = _Response(422, body={"statusCode": "NonCombinable"})
PRICED = _Response(200, body={"bets": [{"trueOdds": "3.72",
                                        "selectionsMapped": ["a", "b"]}]})

_REFS = ["0HC1P150_1", "0OU1O950_1"]


def _price(response):
    """Run one on-demand price call and return (result, counters snapshot)."""
    counters = FetchCounters("draftkings", "on_demand")
    counters.bump("targets_attempted")       # the caller counts the target
    result = draftkings.price_selection_set(_Client(response), _REFS,
                                            counters=counters)
    return result, counters


# --------------------------------------------------------------------------- #
# A blocked endpoint must be counted as a transport failure                    #
# --------------------------------------------------------------------------- #

def test_blocked_price_endpoint_counts_a_transport_error():
    """403 is the book refusing US, not the book refusing this combo."""
    result, counters = _price(AKAMAI_DENY)
    assert result is None                    # unchanged: declines never raise
    assert counters.snapshot().transport_errors == 1, (
        "a 403 blockade must be counted as a transport error, or a dead DK "
        "endpoint is indistinguishable from a run of unpriceable combos")


def test_blocked_price_endpoint_does_not_arm_the_prices_empty_tripwire():
    """``prices_empty`` means "the parser went silent". A dead endpoint must
    not fire it — that is what sent issue #102 looking at the wrong layer."""
    _, counters = _price(AKAMAI_DENY)
    assert "prices_empty" not in counters.tripwires()


def test_rate_limited_price_endpoint_counts_a_transport_error():
    """429 survives ``request_with_retry`` only when retries are exhausted, so
    reaching here means DK is throttling us — a transport verdict, not a
    per-combo one."""
    _, counters = _price(_Response(429, text="Too Many Requests"))
    assert counters.snapshot().transport_errors == 1


# --------------------------------------------------------------------------- #
# A per-combo decline must stay a per-combo decline                            #
# --------------------------------------------------------------------------- #

def test_non_combinable_is_not_a_transport_error():
    """422 is normal and high-volume — most combos decline this way. Counting
    it as transport would mask a real blockade behind constant noise."""
    result, counters = _price(NON_COMBINABLE)
    assert result is None
    assert counters.snapshot().transport_errors == 0


def test_non_combinable_still_arms_the_prices_empty_tripwire():
    """The mirror of the 403 case: when transport is healthy and nothing
    priced, the tripwire SHOULD fire. Silencing every non-200 would trade one
    blind spot for another."""
    _, counters = _price(NON_COMBINABLE)
    assert "prices_empty" in counters.tripwires()


def test_successful_price_is_unchanged():
    result, counters = _price(PRICED)
    assert result == 3.72
    assert counters.snapshot().transport_errors == 0


def test_blocked_endpoint_makes_exactly_one_wire_call():
    """A 403 must not be retried: issue #90 found our own retry volume is what
    holds a reputation block open."""
    client = _Client(AKAMAI_DENY)
    draftkings.price_selection_set(client, _REFS,
                                   counters=FetchCounters("draftkings",
                                                          "on_demand"))
    assert client.session.posts == 1


def test_price_selection_set_still_never_raises():
    """The on-demand contract: one unpriceable combo must never abort a
    flight, whatever the status."""
    assert draftkings.price_selection_set(_Client(AKAMAI_DENY), _REFS) is None
