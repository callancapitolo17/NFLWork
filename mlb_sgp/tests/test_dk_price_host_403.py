"""Regression tests for issue #39 — DraftKings' price host returned 403.

Two distinct defects, both pinned here.

1. THE OUTAGE. DK reads from ``sportsbook-nash.draftkings.com`` but prices on
   ``gaming-us-nj.draftkings.com``. On 2026-06-24 an Akamai edge rule started
   denying ``POST /api/wager/v1/calculateBets`` — events and structure stayed
   green (200, full selection lists) while every price call died, so DK's rows
   went from ~18k/day to zero. The rule matches that path EXACTLY; DK's own
   betslip builds the URL as ``{locale}/api/wager/v1/calculateBets`` (English
   maps to an empty prefix), and the origin still routes the ``/en/`` form,
   which the WAF rule does not cover. Verified live 2026-07-30:
       POST /api/wager/v1/calculateBets     -> 403 AkamaiGHost "Access Denied"
       POST /en/api/wager/v1/calculateBets  -> 200 correlated SGP price
   The URL constant is therefore load-bearing, and the ``/en/`` segment must
   never be "tidied away".

2. THE DIAGNOSIS COST. ``calculate_sgp`` swallows a non-200 into ``None``, so
   the status never reached ``PriceCallTally``'s all-failed verdict. The
   verdict raised ``BookTransportError:price`` with ``status_code=None``, and
   ``sgp_fetch_health.error_class`` lost the ``:403`` — the single field that
   would have said "the price HOST is blocked" instead of costing a session of
   host-by-host bisection. The tally now carries the last decline status.

No network: every test injects a fake session.
"""
from __future__ import annotations

import inspect
import sys
from pathlib import Path

import pytest

# The DK client lazily imports its legacy scraper module by TOP-LEVEL name,
# which only resolves with mlb_sgp/ itself on sys.path (see
# test_book_transport_errors.py for the same insert).
_MLB_SGP_DIR = Path(__file__).resolve().parents[1]
if str(_MLB_SGP_DIR) not in sys.path:
    sys.path.insert(0, str(_MLB_SGP_DIR))

from mlb_sgp._shared import BookTransportError, PriceCallTally  # noqa: E402


# --------------------------------------------------------------------------- #
# Fake transport                                                               #
# --------------------------------------------------------------------------- #

class _Response:
    def __init__(self, status_code, body=None, text=""):
        self.status_code = status_code
        self._body = body
        self.text = text

    def json(self):
        if self._body is None:
            raise ValueError("Expecting value: line 1 column 1 (char 0)")
        return self._body


class _RecordingSession:
    """Captures the URL of every POST and replies with a canned response."""

    def __init__(self, response):
        self._response = response
        self.posted_urls: list[str] = []

    def post(self, url, **kwargs):
        self.posted_urls.append(url)
        return self._response

    def get(self, url, **kwargs):
        return self._response


AKAMAI_DENY = _Response(403, body=None, text="<HTML><HEAD>\n<TITLE>Access Denied</TITLE>")


# --------------------------------------------------------------------------- #
# 1. The URL must carry DK's locale prefix                                     #
# --------------------------------------------------------------------------- #

def test_dk_calculate_bets_url_carries_the_en_locale_prefix():
    """The bare ``/api/wager/v1/calculateBets`` path is WAF-denied; only the
    locale-prefixed form reaches the origin. Removing ``/en/`` re-breaks DK."""
    from scraper_draftkings_sgp import DK_CALCULATE_BETS_URL
    assert DK_CALCULATE_BETS_URL.endswith("/en/api/wager/v1/calculateBets"), (
        "DK's price URL lost its /en/ locale prefix — the bare wager path is "
        "denied at Akamai's edge (issue #39)")


def test_dk_calculate_bets_url_still_targets_the_wager_host():
    """Guard the other half: the prefix fix must not have moved the host."""
    from scraper_draftkings_sgp import DK_CALCULATE_BETS_URL
    assert DK_CALCULATE_BETS_URL.startswith("https://gaming-us-nj.draftkings.com/")


def test_dk_price_calls_hit_the_locale_prefixed_url():
    """The constant is what both price paths actually POST to."""
    from scraper_draftkings_sgp import DK_CALCULATE_BETS_URL, calculate_sgp
    session = _RecordingSession(AKAMAI_DENY)
    calculate_sgp(session, "0HC123P150_1", "0OU123O950_1")
    assert session.posted_urls == [DK_CALCULATE_BETS_URL]
    assert "/en/api/wager/v1/calculateBets" in session.posted_urls[0]


def test_no_second_hardcoded_copy_of_the_wager_url():
    """The trifecta scraper used to redefine this URL, so the fix would have
    repaired the SGP scraper and left that one 403ing. Every DK price caller
    must import the one constant."""
    pkg = Path(__file__).resolve().parents[1]
    offenders = []
    for path in pkg.glob("*.py"):
        # scraper_draftkings_sgp.py owns the one true definition.
        if path.name == "scraper_draftkings_sgp.py":
            continue
        for lineno, line in enumerate(path.read_text().splitlines(), 1):
            # A URL LITERAL, not prose: docstrings name the host when
            # explaining why the price stage is separate.
            if "https://gaming-us-nj" not in line:
                continue
            if line.lstrip().startswith("#"):
                continue
            offenders.append(f"{path.name}:{lineno}")
    assert not offenders, (
        f"hardcoded DK wager host outside scraper_draftkings_sgp.py: {offenders}")


def test_trifecta_imports_the_shared_price_url():
    import scraper_draftkings_trifecta as trifecta
    from scraper_draftkings_sgp import DK_CALCULATE_BETS_URL
    assert trifecta.DK_CALCULATE_BETS_URL == DK_CALCULATE_BETS_URL


def test_on_demand_price_path_shares_the_same_url_constant():
    """``price_selection_set`` is the on-demand (live quote) price call. It
    must be fixed by the same constant, or on-demand stays dark while the
    sweep recovers."""
    from mlb_sgp import draftkings
    src = inspect.getsource(draftkings.price_selection_set)
    assert "DK_CALCULATE_BETS_URL" in src


# --------------------------------------------------------------------------- #
# 2. The all-failed verdict must carry the decline status                      #
# --------------------------------------------------------------------------- #

def test_price_tally_verdict_reports_the_last_decline_status():
    tally = PriceCallTally("draftkings")
    start = tally.snapshot()
    for _ in range(PriceCallTally.MIN_ATTEMPTS_FOR_VERDICT):
        tally.record(False, status=403)
    with pytest.raises(BookTransportError) as exc:
        tally.verdict(start)
    assert exc.value.stage == "price"
    assert exc.value.status_code == 403
    assert "403" in str(exc.value)


def test_price_tally_verdict_status_is_none_when_no_status_was_seen():
    """A decline with no HTTP status (parse failure, None return) must not
    invent one — status stays absent, as before issue #39."""
    tally = PriceCallTally("draftkings")
    start = tally.snapshot()
    for _ in range(PriceCallTally.MIN_ATTEMPTS_FOR_VERDICT):
        tally.record(False)
    with pytest.raises(BookTransportError) as exc:
        tally.verdict(start)
    assert exc.value.status_code is None


def test_price_tally_status_is_keyword_only():
    """Positional drift here would silently rebind ``ok``."""
    sig = inspect.signature(PriceCallTally.record)
    with pytest.raises(TypeError):
        sig.bind(PriceCallTally("dk"), False, 403)
    sig.bind(PriceCallTally("dk"), False, status=403)


def test_price_tally_success_clears_the_recorded_status():
    """A cycle that priced anything is healthy — no stale status should leak
    into a later cycle's verdict."""
    tally = PriceCallTally("draftkings")
    for _ in range(PriceCallTally.MIN_ATTEMPTS_FOR_VERDICT):
        tally.record(False, status=403)
    start = tally.snapshot()                 # a NEW cycle begins
    for _ in range(PriceCallTally.MIN_ATTEMPTS_FOR_VERDICT):
        tally.record(False)                  # declines with no status
    with pytest.raises(BookTransportError) as exc:
        tally.verdict(start)
    assert exc.value.status_code is None


# --------------------------------------------------------------------------- #
# 3. calculate_sgp must report the decline status                              #
# --------------------------------------------------------------------------- #

def test_calculate_sgp_reports_decline_status_to_callback():
    from scraper_draftkings_sgp import calculate_sgp
    seen: list[int] = []
    session = _RecordingSession(AKAMAI_DENY)
    result = calculate_sgp(session, "0HC1P150_1", "0OU1O950_1",
                           on_decline=seen.append)
    assert result is None            # unchanged contract: declines are not fatal
    assert seen == [403]


def test_calculate_sgp_reports_422_as_a_decline_status():
    """422 (non-combinable) is a normal, per-combo outcome. It must still be
    reported so a cycle of pure 422s is not mistaken for a 403 blockade."""
    from scraper_draftkings_sgp import calculate_sgp
    seen: list[int] = []
    session = _RecordingSession(_Response(422, body={"statusCode": "NonCombinable"}))
    assert calculate_sgp(session, "0HC1P150_1", "0OU1O950_1",
                         on_decline=seen.append) is None
    assert seen == [422]


def test_calculate_sgp_does_not_report_on_success():
    from scraper_draftkings_sgp import calculate_sgp
    seen: list[int] = []
    ok = _Response(200, body={"bets": [{"trueOdds": 3.72, "displayOdds": "+272",
                                        "selectionsMapped": ["a", "b"]}]})
    result = calculate_sgp(_RecordingSession(ok), "0HC1P150_1", "0OU1O950_1",
                           on_decline=seen.append)
    assert result == {"trueOdds": 3.72, "displayOdds": "+272"}
    assert seen == []


def test_calculate_sgp_on_decline_is_keyword_only():
    sig = inspect.signature(
        __import__("scraper_draftkings_sgp").calculate_sgp)
    with pytest.raises(TypeError):
        sig.bind(object(), "a", "b", False, print)
    sig.bind(object(), "a", "b", on_decline=print)


def test_calculate_sgp_without_callback_is_unchanged():
    """Every existing call site passes no callback — that path must not raise."""
    from scraper_draftkings_sgp import calculate_sgp
    assert calculate_sgp(_RecordingSession(AKAMAI_DENY), "0HC1P150_1",
                         "0OU1O950_1") is None


# --------------------------------------------------------------------------- #
# 4. End-to-end: a blocked price host yields error_class price:403             #
# --------------------------------------------------------------------------- #

def test_blocked_price_host_produces_error_class_with_status():
    """The whole point of issue #39's diagnostic half: one string in
    ``sgp_fetch_health`` that names the stage AND the status."""
    from kalshi_common.sgp_health import transport_error_class
    from scraper_draftkings_sgp import calculate_sgp

    tally = PriceCallTally("draftkings")
    priced = tally.wrap(calculate_sgp)
    start = tally.snapshot()
    session = _RecordingSession(AKAMAI_DENY)
    for _ in range(PriceCallTally.MIN_ATTEMPTS_FOR_VERDICT):
        priced(session, "0HC1P150_1", "0OU1O950_1", on_decline=tally.note_status)

    with pytest.raises(BookTransportError) as exc:
        tally.verdict(start)
    assert transport_error_class(exc.value) == "BookTransportError:price:403"


# --------------------------------------------------------------------------- #
# 5. The orchestrator must actually WIRE the callback                          #
# --------------------------------------------------------------------------- #

def test_orchestrator_binds_the_decline_callback_to_the_real_price_fn():
    """``_accepts_on_decline`` degrades silently by design, so pin the happy
    path: the real ``calculate_sgp`` must be recognised as accepting the
    callback. A rename that breaks the binding fails here rather than quietly
    dropping the status from every future verdict."""
    from mlb_sgp.draftkings import _accepts_on_decline
    from scraper_draftkings_sgp import calculate_sgp
    assert _accepts_on_decline(calculate_sgp) is True


def test_orchestrator_tolerates_a_substitute_without_the_callback():
    """Tests and the dashboard monkeypatch this seam with pre-#39 stand-ins;
    those must keep pricing rather than raising TypeError on every combo."""
    from mlb_sgp.draftkings import _accepts_on_decline

    def legacy_stub(session, spread_sel, total_sel, verbose=False):
        return {"trueOdds": 2.0, "displayOdds": "+100"}

    assert _accepts_on_decline(legacy_stub) is False


def test_accepts_on_decline_survives_an_uninspectable_callable():
    from unittest.mock import MagicMock

    from mlb_sgp.draftkings import _accepts_on_decline
    assert _accepts_on_decline(MagicMock()) in (True, False)   # must not raise
    assert _accepts_on_decline(len) is False
