"""Regression tests for issue #33 — "book down" must not look like "no markets".

The meta-bug: every client swallowed a 403 / DNS failure / malformed JSON into
an empty list. The orchestrator then returned ``[]`` (not ``None``), which
``SGPService._book_done`` scores as SUCCESS — resetting the failure counter so
the 3-strike client reinit never fires — and which ``sgp_cycle`` treats as "the
book has no markets today", running ``clear_source`` and DELETING the book's
rows. A dead book was indistinguishable from an empty slate.

The contract these tests pin, per book:

  transport failure (non-200 other than 404, connection error, timeout,
  malformed JSON, auth gate)          -> raise BookTransportError
  valid-but-empty payload             -> return [] / empty structure
  404 on ONE event's structure fetch  -> skip that event, keep the cycle

No network: every test injects a fake session.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

# The clients lazily import their legacy scraper module by TOP-LEVEL name
# (``from scraper_novig_sgp import ...``), which only resolves with mlb_sgp/
# itself on sys.path. SGPService does this insert at construction time; here we
# do it explicitly so these tests never depend on import order.
_MLB_SGP_DIR = Path(__file__).resolve().parents[1]
if str(_MLB_SGP_DIR) not in sys.path:
    sys.path.insert(0, str(_MLB_SGP_DIR))

from mlb_sgp._shared import BookTransportError  # noqa: E402

FIX = Path(__file__).parent / "fixtures"


# --------------------------------------------------------------------------- #
# Fake transport helpers                                                       #
# --------------------------------------------------------------------------- #

class FakeResponse:
    """Minimal stand-in for a curl_cffi response."""

    def __init__(self, status_code=200, body=None, text=None):
        self.status_code = status_code
        self._body = body
        self.text = text if text is not None else json.dumps(body or {})

    def json(self):
        if self._body is None:
            raise ValueError("Expecting value: line 1 column 1 (char 0)")
        return self._body

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(f"HTTP {self.status_code}")


class FakeSession:
    """Returns the same response for every request (get or post)."""

    def __init__(self, response):
        self._response = response

    def get(self, url, **kwargs):
        return self._response

    def post(self, url, **kwargs):
        return self._response


class BoomSession:
    """Every request raises — models a DNS failure or connection reset."""

    def __init__(self, exc=None):
        self._exc = exc or ConnectionError("[Errno 8] nodename nor servname provided")

    def get(self, url, **kwargs):
        raise self._exc

    def post(self, url, **kwargs):
        raise self._exc


def _client(cls, session):
    """Build a client with __init__ bypassed (no real HTTP session)."""
    c = cls.__new__(cls)
    c.session = session
    c.verbose = False
    return c


FORBIDDEN = FakeResponse(status_code=403, body={"error": "forbidden"},
                         text="Access Denied")
SERVER_ERROR = FakeResponse(status_code=500, body={}, text="oops")
MALFORMED = FakeResponse(status_code=200, body=None, text="<html>WAF challenge</html>")
NOT_FOUND = FakeResponse(status_code=404, body={}, text="not found")


# --------------------------------------------------------------------------- #
# ProphetX                                                                     #
# --------------------------------------------------------------------------- #

def test_px_list_events_raises_on_403():
    from mlb_sgp.prophetx_client import ProphetXClient
    client = _client(ProphetXClient, FakeSession(FORBIDDEN))
    with pytest.raises(BookTransportError) as exc:
        client.list_events()
    assert exc.value.book == "prophetx"
    assert exc.value.stage == "events"
    assert exc.value.status_code == 403


def test_px_list_events_raises_on_connection_error():
    from mlb_sgp.prophetx_client import ProphetXClient
    client = _client(ProphetXClient, BoomSession())
    with pytest.raises(BookTransportError):
        client.list_events()


def test_px_list_events_raises_on_malformed_json():
    from mlb_sgp.prophetx_client import ProphetXClient
    client = _client(ProphetXClient, FakeSession(MALFORMED))
    with pytest.raises(BookTransportError):
        client.list_events()


def test_px_list_events_valid_empty_returns_empty_list():
    from mlb_sgp.prophetx_client import ProphetXClient
    ok = FakeResponse(200, {"data": {"tournaments": []}})
    client = _client(ProphetXClient, FakeSession(ok))
    assert client.list_events() == []


def test_px_fetch_event_markets_raises_on_500():
    from mlb_sgp.prophetx_client import ProphetXClient
    client = _client(ProphetXClient, FakeSession(SERVER_ERROR))
    with pytest.raises(BookTransportError) as exc:
        client.fetch_event_markets("10077925")
    assert exc.value.stage == "structure"


def test_px_fetch_event_markets_404_skips_event_without_raising():
    """A single removed/postponed game must not abort the whole book."""
    from mlb_sgp.prophetx_client import ProphetXClient
    client = _client(ProphetXClient, FakeSession(NOT_FOUND))
    assert client.fetch_event_markets("gone") == []


# --------------------------------------------------------------------------- #
# Novig                                                                        #
# --------------------------------------------------------------------------- #

def test_nv_list_events_raises_on_dns_failure():
    """api.novig.us stopped resolving in production — it must be loud."""
    from mlb_sgp.novig_client import NovigClient
    client = _client(NovigClient, BoomSession())
    with pytest.raises(BookTransportError) as exc:
        client.list_events()
    assert exc.value.book == "novig"
    assert exc.value.stage == "events"


def test_nv_list_events_raises_on_500():
    from mlb_sgp.novig_client import NovigClient
    client = _client(NovigClient, FakeSession(SERVER_ERROR))
    with pytest.raises(BookTransportError):
        client.list_events()


def test_nv_list_events_raises_on_malformed_json():
    from mlb_sgp.novig_client import NovigClient
    client = _client(NovigClient, FakeSession(MALFORMED))
    with pytest.raises(BookTransportError):
        client.list_events()


def test_nv_list_events_valid_empty_returns_empty_list():
    from mlb_sgp.novig_client import NovigClient
    ok = FakeResponse(200, {"data": {"event": []}})
    client = _client(NovigClient, FakeSession(ok))
    assert client.list_events() == []


def test_nv_fetch_event_legs_raises_on_500():
    from mlb_sgp.novig_client import NovigClient
    client = _client(NovigClient, FakeSession(SERVER_ERROR))
    with pytest.raises(BookTransportError) as exc:
        client.fetch_event_legs("uuid-1")
    assert exc.value.stage == "structure"


def test_nv_fetch_event_legs_404_skips_event_without_raising():
    from mlb_sgp.novig_client import NovigClient
    client = _client(NovigClient, FakeSession(NOT_FOUND))
    legs = client.fetch_event_legs("gone")
    assert legs.spread_legs == [] and legs.total_legs == []


# --------------------------------------------------------------------------- #
# DraftKings (transport lives in the legacy scraper helpers + dk_client)        #
# --------------------------------------------------------------------------- #

def test_dk_fetch_events_raises_on_403():
    """DK has been returning 403 in production since ~2026-06."""
    from scraper_draftkings_sgp import fetch_dk_events
    with pytest.raises(BookTransportError) as exc:
        fetch_dk_events(FakeSession(FORBIDDEN))
    assert exc.value.book == "draftkings"
    assert exc.value.stage == "events"


def test_dk_fetch_events_valid_empty_returns_empty_list():
    from scraper_draftkings_sgp import fetch_dk_events
    ok = FakeResponse(200, {"events": []})
    assert fetch_dk_events(FakeSession(ok)) == []


def test_dk_subcat_markets_raises_on_403():
    from scraper_draftkings_sgp import _fetch_subcat_markets
    with pytest.raises(BookTransportError) as exc:
        _fetch_subcat_markets(FakeSession(FORBIDDEN), "evt1", "4519")
    assert exc.value.stage == "structure"


def test_dk_subcat_markets_404_skips_event_without_raising():
    from scraper_draftkings_sgp import _fetch_subcat_markets
    assert _fetch_subcat_markets(FakeSession(NOT_FOUND), "gone", "4519") == []


def test_dk_selection_ids_raises_on_403():
    from scraper_draftkings_sgp import fetch_selection_ids
    with pytest.raises(BookTransportError) as exc:
        fetch_selection_ids(FakeSession(FORBIDDEN), "evt1")
    assert exc.value.stage == "structure"


def test_dk_selection_ids_404_skips_event_without_raising():
    from scraper_draftkings_sgp import fetch_selection_ids
    out = fetch_selection_ids(FakeSession(NOT_FOUND), "gone")
    assert out["fg"]["spreads"] == {} and out["f5"]["spreads"] == {}


def test_dk_client_fetch_event_selections_raises_on_403():
    from mlb_sgp.dk_client import DraftKingsClient
    client = _client(DraftKingsClient, FakeSession(FORBIDDEN))
    with pytest.raises(BookTransportError):
        client.fetch_event_selections("evt1")


# --------------------------------------------------------------------------- #
# FanDuel                                                                      #
# --------------------------------------------------------------------------- #

def test_fd_fetch_events_raises_on_403():
    from scraper_fanduel_sgp import fetch_fd_events
    with pytest.raises(BookTransportError) as exc:
        fetch_fd_events(FakeSession(FORBIDDEN))
    assert exc.value.book == "fanduel"
    assert exc.value.stage == "events"


def test_fd_fetch_events_raises_on_connection_error():
    from scraper_fanduel_sgp import fetch_fd_events
    with pytest.raises(BookTransportError):
        fetch_fd_events(BoomSession())


def test_fd_event_runners_raises_on_403():
    from scraper_fanduel_sgp import fetch_event_runners
    with pytest.raises(BookTransportError) as exc:
        fetch_event_runners(FakeSession(FORBIDDEN), "evt1", "Home", "Away")
    assert exc.value.stage == "structure"


def test_fd_event_runners_404_skips_event_without_raising():
    from scraper_fanduel_sgp import fetch_event_runners
    out = fetch_event_runners(FakeSession(NOT_FOUND), "gone", "Home", "Away")
    assert out["fg"]["spreads"] == {} and out["f5"]["totals"] == {}


def test_fd_client_event_page_raises_on_403():
    from mlb_sgp.fd_client import FanDuelClient
    client = _client(FanDuelClient, FakeSession(FORBIDDEN))
    with pytest.raises(BookTransportError):
        client.fetch_event_page("evt1")


def test_fd_client_event_page_valid_empty_returns_empty():
    from mlb_sgp.fd_client import FanDuelClient
    client = _client(FanDuelClient, FakeSession(FakeResponse(200, {})))
    markets, runners = client.fetch_event_page("evt1")
    assert markets == [] and runners == []


# --------------------------------------------------------------------------- #
# BetMGM                                                                       #
# --------------------------------------------------------------------------- #

def _mgm_client(session, accessid="abc123"):
    from mlb_sgp.betmgm_client import BetMGMClient
    c = BetMGMClient.__new__(BetMGMClient)
    c.state = "pa"
    c.verbose = False
    c.session = session
    c._accessid = accessid
    return c


def test_mgm_accessid_raises_when_config_carries_no_id():
    """A clientconfig response with no publicAccessId is an AUTH failure,
    not an empty slate — every downstream call is gated on it."""
    from mlb_sgp.betmgm_client import BetMGMClient
    client = _mgm_client(FakeSession(FakeResponse(200, {}, text="{}")), accessid="")
    with pytest.raises(BookTransportError) as exc:
        client.accessid()
    assert exc.value.book == "betmgm"
    assert exc.value.stage == "auth"


def test_mgm_accessid_raises_on_connection_error():
    client = _mgm_client(BoomSession(), accessid="")
    with pytest.raises(BookTransportError):
        client.accessid()


def test_mgm_list_events_raises_on_500():
    client = _mgm_client(FakeSession(SERVER_ERROR))
    with pytest.raises(BookTransportError) as exc:
        client.list_events()
    assert exc.value.stage == "events"


def test_mgm_list_events_raises_on_connection_error():
    client = _mgm_client(BoomSession())
    with pytest.raises(BookTransportError):
        client.list_events()


def test_mgm_list_events_valid_empty_returns_empty_list():
    client = _mgm_client(FakeSession(FakeResponse(200, {"fixtures": []})))
    assert client.list_events() == []


def test_mgm_fetch_markets_raises_on_500():
    client = _mgm_client(FakeSession(SERVER_ERROR))
    with pytest.raises(BookTransportError) as exc:
        client.fetch_markets("fx1")
    assert exc.value.stage == "structure"


def test_mgm_fetch_markets_404_skips_event_without_raising():
    client = _mgm_client(FakeSession(NOT_FOUND))
    assert client.fetch_markets("gone") == []


# --------------------------------------------------------------------------- #
# Caesars                                                                      #
# --------------------------------------------------------------------------- #

def _czr_client(session, token="tok"):
    from mlb_sgp.caesars_client import CaesarsClient
    c = CaesarsClient.__new__(CaesarsClient)
    c.state = "nj"
    c.verbose = False
    c.session = session
    c._token = token
    c._device = "dev"
    c._cookies = {}
    c._minted_at = 1e12          # far future -> ensure_token short-circuits True
    return c


def test_czr_list_events_raises_when_token_cannot_be_minted(monkeypatch):
    """The WAF broker's quiet skip made Caesars look like a permanent off-day."""
    from mlb_sgp.caesars_client import CaesarsClient
    client = _czr_client(FakeSession(FakeResponse(200, {"competitions": []})), token="")
    monkeypatch.setattr(client, "ensure_token", lambda: False)
    with pytest.raises(BookTransportError) as exc:
        client.list_events()
    assert exc.value.book == "caesars"
    assert exc.value.stage == "auth"


def test_czr_list_events_raises_when_waf_serves_a_challenge(monkeypatch):
    import mlb_sgp.caesars_client as czr
    monkeypatch.setattr(czr.time, "sleep", lambda *_: None)
    client = _czr_client(FakeSession(MALFORMED))
    with pytest.raises(BookTransportError) as exc:
        client.list_events()
    assert exc.value.stage == "events"


def test_czr_list_events_raises_on_connection_error(monkeypatch):
    import mlb_sgp.caesars_client as czr
    monkeypatch.setattr(czr.time, "sleep", lambda *_: None)
    client = _czr_client(BoomSession())
    with pytest.raises(BookTransportError):
        client.list_events()


def test_czr_list_events_valid_empty_returns_empty_list():
    client = _czr_client(FakeSession(FakeResponse(200, {"competitions": []})))
    assert client.list_events() == []


def test_czr_fetch_event_raises_on_waf_challenge(monkeypatch):
    import mlb_sgp.caesars_client as czr
    monkeypatch.setattr(czr.time, "sleep", lambda *_: None)
    client = _czr_client(FakeSession(MALFORMED))
    with pytest.raises(BookTransportError) as exc:
        client.fetch_event("evt1")
    assert exc.value.stage == "structure"


def test_czr_fetch_event_404_skips_event_without_raising():
    client = _czr_client(FakeSession(NOT_FOUND))
    assert client.fetch_event("gone") is None


def test_mgm_fetch_markets_retries_once_after_a_stale_accessid_gate():
    """BetMGM's 400 "Access id ..." is the one non-200 worth retrying."""
    from mlb_sgp.betmgm_client import BetMGMClient

    gate = FakeResponse(400, {}, text='{"message":"Access id is invalid"}')
    ok = FakeResponse(200, {"fixtures": [{"optionMarkets": [{"id": "m1"}]}]})

    class TwoStepSession:
        def __init__(self):
            self.calls = 0

        def get(self, url, **kwargs):
            self.calls += 1
            # 1: gated markets fetch, 2: accessid re-harvest, 3: retried fetch
            if self.calls == 1:
                return gate
            if self.calls == 2:
                return FakeResponse(200, {}, text='{"publicAccessId":"newid"}')
            return ok

    session = TwoStepSession()
    client = _mgm_client(session)
    assert client.fetch_markets("fx1") == [{"id": "m1"}]
    assert session.calls == 3


def test_mgm_fetch_markets_raises_when_the_gate_persists_after_refresh():
    from mlb_sgp.betmgm_client import BetMGMClient
    gate = FakeResponse(400, {}, text='{"message":"Access id is invalid"}')

    class AlwaysGated:
        def get(self, url, **kwargs):
            if "clientconfig" in url:
                return FakeResponse(200, {}, text='{"publicAccessId":"newid"}')
            return gate

    client = _mgm_client(AlwaysGated())
    with pytest.raises(BookTransportError) as exc:
        client.fetch_markets("fx1")
    assert exc.value.stage == "structure"
    assert exc.value.status_code == 400


def test_czr_404_on_the_tabs_feed_raises_rather_than_looking_empty(monkeypatch):
    """A 404 on ONE event is a skip; a 404 on the FEED means it moved."""
    import mlb_sgp.caesars_client as czr
    monkeypatch.setattr(czr.time, "sleep", lambda *_: None)
    client = _czr_client(FakeSession(NOT_FOUND))
    with pytest.raises(BookTransportError) as exc:
        client.list_events()
    assert exc.value.stage == "events"
