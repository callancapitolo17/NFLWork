"""Tests for the #42 verification gate and its golden record/replay harness.

The gate exists to tell a dead book from a book that simply does not offer the
line you asked for. If IT is wrong, every conclusion drawn from it is wrong —
so the classification logic and the seam declarations are pinned here.
"""
from __future__ import annotations

import importlib
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

MLB_SGP_DIR = Path(__file__).resolve().parents[1]
if str(MLB_SGP_DIR) not in sys.path:
    sys.path.insert(0, str(MLB_SGP_DIR))

from mlb_sgp._shared import GameRef                    # noqa: E402
from mlb_sgp.tests import golden_replay as gr          # noqa: E402
from mlb_sgp import verify_books as vb                 # noqa: E402

GAME = GameRef(game_id="verify", home_team="Philadelphia Phillies",
               away_team="Washington Nationals", commence_time=None)


# --------------------------------------------------------------------- #
# Shape construction                                                     #
# --------------------------------------------------------------------- #

def test_shapes_cover_the_four_required_forms():
    shapes = vb.shapes_for(-1.5, 8.5)
    assert set(shapes) == {"spread_x_total", "ml_x_total", "spread_x_ml",
                           "three_leg"}
    assert len(shapes["three_leg"]) == 3
    assert all(len(v) == 2 for k, v in shapes.items() if k != "three_leg")


def test_three_leg_contains_the_pair_no_book_prices():
    """The 3-leg's ❌ column is caused by the spread+ml pair, so the two
    shapes must stay coupled — if this ever stops holding, the README's
    explanation of the capability matrix is wrong."""
    shapes = vb.shapes_for(-1.5, 8.5)
    pair = {(leg.market_type, leg.side) for leg in shapes["spread_x_ml"]}
    three = {(leg.market_type, leg.side) for leg in shapes["three_leg"]}
    assert pair <= three


def test_shapes_use_the_lines_they_were_given():
    shapes = vb.shapes_for(-2.5, 10.0)
    spread, total = shapes["spread_x_total"]
    assert spread.line == -2.5 and total.line == 10.0
    ml, _ = shapes["ml_x_total"]
    assert ml.line is None, "a moneyline leg carries no line"


# --------------------------------------------------------------------- #
# Outcome classification — the whole point of the harness                #
# --------------------------------------------------------------------- #

class _FakeService:
    """Stands in for SGPService with scripted hooks."""

    def __init__(self, *, match=None, structure=None, offers=(),
                 match_raises=None, structure_raises=None):
        self._match, self._structure = match, structure
        self._offers, self._match_raises = offers, match_raises
        self._structure_raises = structure_raises

    def _ensure_client(self, book):
        return SimpleNamespace(client=object())

    def _book_on_demand_hooks(self, book, *, counters=None):
        def match_event(client, game):
            if self._match_raises:
                raise self._match_raises
            return self._match

        def build_structure(client, event, game):
            if self._structure_raises:
                raise self._structure_raises
            return self._structure

        def resolve(structure, legs, home, away):
            leg = legs[0]
            return legs if (leg.market_type, leg.line) in self._offers else None

        return {"match_event": match_event, "build_structure": build_structure,
                "resolve": resolve, "price": lambda *a, **k: None}


def test_transport_failure_is_down_not_not_offered():
    """A raising client is a REAL failure and must never be softened into a
    coverage fact — that conflation is issue #32's meta-bug."""
    service = _FakeService(match_raises=RuntimeError("403 blocked"))
    result = vb.probe_coverage(service, "novig", GAME)
    assert result["status"] == "down"
    assert "403 blocked" in result["error"]


def test_unmatched_event_is_down():
    result = vb.probe_coverage(_FakeService(match=None), "novig", GAME)
    assert result["status"] == "down"


def test_empty_structure_is_down():
    service = _FakeService(match={"id": "e1"}, structure=None)
    assert vb.probe_coverage(service, "novig", GAME)["status"] == "down"


def test_structure_exception_is_down():
    service = _FakeService(match={"id": "e1"},
                           structure_raises=ValueError("bad tree"))
    result = vb.probe_coverage(service, "novig", GAME)
    assert result["status"] == "down"
    assert "bad tree" in result["error"]


def test_healthy_structure_with_no_candidate_lines_is_not_offered():
    """Caesars and FanDuel carry only the main line. That must read as
    coverage, not as a dead book."""
    service = _FakeService(match={"id": "e1"}, structure={"ok": True},
                           offers=())
    result = vb.probe_coverage(service, "caesars", GAME)
    assert result["status"] == "not_offered"
    assert result["spreads"] == [] and result["totals"] == []


def test_partial_coverage_reports_only_the_offered_lines():
    offers = {("spread", -1.5), ("total", 9.0)}
    service = _FakeService(match={"id": "e1"}, structure={"ok": True},
                           offers=offers)
    result = vb.probe_coverage(service, "fanduel", GAME)
    assert result["status"] == "ok"
    assert result["spreads"] == [-1.5]
    assert result["totals"] == [9.0]


def test_candidate_ladder_spans_whole_and_half_totals():
    """The books split on this; a ladder of only .5 lines would report
    FanDuel and Caesars as offering nothing."""
    assert any(float(t).is_integer() for t in vb.CANDIDATE_TOTALS)
    assert any(not float(t).is_integer() for t in vb.CANDIDATE_TOTALS)


def test_pick_game_falls_through_an_unreachable_book(monkeypatch):
    """A gate for six books must not fail to start because ONE book is down.

    This is not hypothetical — pick_game read Novig only, and a Novig TLS
    drop made the whole harness unrunnable mid-review.
    """
    from datetime import datetime, timedelta, timezone
    soon = datetime.now(timezone.utc) + timedelta(hours=3)
    good = SimpleNamespace(home_team="Philadelphia Phillies",
                           away_team="Washington Nationals",
                           start_time=soon.isoformat())

    calls = []

    def fake_events(book):
        calls.append(book)
        if book == "novig":
            return []          # unreachable
        return [good]

    monkeypatch.setattr(vb, "_events_from", fake_events)
    game = vb.pick_game(source_books=("novig", "caesars"))
    assert calls == ["novig", "caesars"], "must move on to the next book"
    assert game.home_team == "Philadelphia Phillies"


def test_pick_game_skips_games_already_started(monkeypatch):
    from datetime import datetime, timedelta, timezone
    started = datetime.now(timezone.utc) - timedelta(minutes=30)
    upcoming = datetime.now(timezone.utc) + timedelta(hours=2)
    monkeypatch.setattr(vb, "_events_from", lambda book: [
        SimpleNamespace(home_team="Philadelphia Phillies",
                        away_team="Washington Nationals",
                        start_time=started.isoformat()),
        SimpleNamespace(home_team="New York Yankees",
                        away_team="Boston Red Sox",
                        start_time=upcoming.isoformat()),
    ])
    game = vb.pick_game(source_books=("novig",))
    assert game.home_team == "New York Yankees"


def test_pick_game_raises_when_every_book_is_unreachable(monkeypatch):
    monkeypatch.setattr(vb, "_events_from", lambda book: [])
    with pytest.raises(SystemExit) as exc:
        vb.pick_game(source_books=("novig", "caesars"))
    assert "novig" in str(exc.value) and "caesars" in str(exc.value)


def test_secs_formats_a_missing_sample():
    assert vb._secs(None) == "-"
    assert vb._secs(1.234) == "1.23"


# --------------------------------------------------------------------- #
# Golden harness                                                         #
# --------------------------------------------------------------------- #

@pytest.mark.parametrize("value", [
    {"a": (1, 2)},
    {("home", -1.5): ("M0", "S0")},          # FanDuel runner shape
    [("m1", "o1", 1.9)],                     # BetMGM leg shape
    {"nested": {"deep": [{"x": (1,)}]}},
    {7.5: {"over": {"id": "o"}}},            # float-keyed line map
    SimpleNamespace(a=1, b=(2, 3)),
    {1, 2, 3},
])
def test_encode_decode_round_trips(value):
    import json
    restored = gr.decode(json.loads(json.dumps(gr.encode(value))))
    assert restored == value
    assert type(restored) is type(value)


def test_encode_preserves_tuple_vs_list():
    """JSON would flatten both to a list; FD keys its runner dict on tuples,
    so the distinction has to survive."""
    import json
    restored = gr.decode(json.loads(json.dumps(gr.encode({"t": (1, 2),
                                                          "l": [1, 2]}))))
    assert isinstance(restored["t"], tuple)
    assert isinstance(restored["l"], list)


def test_call_key_skips_the_connection_argument():
    """args[0] is a session or client — identity-only, and including it would
    make every replay a miss."""
    session_a, session_b = object(), object()
    assert (gr.call_key((session_a, "e1"), {}, skip=1)
            == gr.call_key((session_b, "e1"), {}, skip=1))


def test_call_key_distinguishes_real_selectors():
    assert (gr.call_key((None, "e1"), {}, skip=1)
            != gr.call_key((None, "e2"), {}, skip=1))


def test_call_key_ignores_noise_kwargs():
    assert (gr.call_key((None, "e1"), {"verbose": True}, skip=1)
            == gr.call_key((None, "e1"), {"verbose": False}, skip=1))


def test_replay_miss_raises_loudly():
    """A missing payload must never be served as an empty result — that is
    exactly the silent-failure mode this epic exists to remove."""
    replayer = gr._Replayer("novig", "submit_parlay", 1, {})
    with pytest.raises(gr.ReplayMiss) as exc:
        replayer(object(), ["a", "b"])
    assert "novig/submit_parlay" in str(exc.value)


def test_recorder_stores_under_the_key_replay_will_ask_for():
    sink = {}
    recorder = gr._Recorder(lambda _s, ids: {"decimal": 2.0}, 1, sink)
    recorder(object(), ["x", "y"])
    replayer = gr._Replayer("novig", "submit_parlay", 1, sink)
    assert replayer(object(), ["x", "y"]) == {"decimal": 2.0}


def test_datetimes_survive_the_json_round_trip():
    """encode() used to pass a datetime straight through, so json's
    default=str stringified it with a SPACE separator and it never became a
    datetime again — a silent type drift, in a repo where timestamp handling
    is a regression gate."""
    from datetime import datetime, timezone
    value = {"start": datetime(2026, 8, 3, 22, 40, tzinfo=timezone.utc)}
    restored = gr.decode(json.loads(json.dumps(gr.encode(value))))
    assert restored == value
    assert isinstance(restored["start"], datetime)


def test_auth_seams_are_redacted_at_capture():
    """A captured accessid is auth material, is never read on replay, and
    must not be committed."""
    sink = {}
    recorder = gr._Recorder(lambda: "a-real-live-accessid", 0, sink,
                            "accessid")
    assert recorder() == "a-real-live-accessid", "the live call is unchanged"
    assert list(sink.values()) == [gr.REDACTED]


def test_encode_qualifies_dataclasses_by_module():
    """Four of the six clients define a dataclass named `Event`; a bare class
    name silently rebuilds the wrong one."""
    from mlb_sgp.novig_client import Event
    encoded = gr.encode(Event("e1", "Home", "Away", "H", "A", "2026-08-03"))
    assert encoded["__dataclass__"] == "mlb_sgp.novig_client.Event"
    assert gr.decode(encoded) == Event("e1", "Home", "Away", "H", "A",
                                       "2026-08-03")


# --------------------------------------------------------------------- #
# Seam declarations must match reality                                   #
# --------------------------------------------------------------------- #

_CLIENT_MODULES = {
    "draftkings": ("mlb_sgp.dk_client", "DraftKingsClient"),
    "fanduel": ("mlb_sgp.fd_client", "FanDuelClient"),
    "prophetx": ("mlb_sgp.prophetx_client", "ProphetXClient"),
    "novig": ("mlb_sgp.novig_client", "NovigClient"),
    "betmgm": ("mlb_sgp.betmgm_client", "BetMGMClient"),
    "caesars": ("mlb_sgp.caesars_client", "CaesarsClient"),
}


@pytest.mark.parametrize("book", vb.ALL_BOOKS)
def test_declared_seams_exist(book):
    """Every seam the golden harness patches must actually be there.

    Guessed seam names are the failure mode this catches: patching a name a
    book does not have leaves the real function live, so a "golden" would
    silently hit the network.
    """
    seams = gr._seams(book)
    for module_name, attr, _skip in seams.get("module_patches", []):
        module = importlib.import_module(module_name)
        assert hasattr(module, attr), f"{module_name}.{attr} does not exist"

    if seams.get("client_methods"):
        mod_name, cls_name = _CLIENT_MODULES[book]
        client_cls = getattr(importlib.import_module(mod_name), cls_name)
        for method in seams["client_methods"]:
            assert hasattr(client_cls, method), (
                f"{cls_name}.{method} does not exist")

    for name in seams.get("fetchers", {}):
        legacy = importlib.import_module(gr._LEGACY_MODULE[book])
        assert hasattr(legacy, name), f"{book}: no legacy fetcher {name}"


@pytest.mark.parametrize("book", vb.ALL_BOOKS)
def test_every_book_declares_a_price_seam(book):
    """A golden that froze structure but not pricing would still hit the wire
    for every combo."""
    seams = gr._seams(book)
    names = ([a for _m, a, _s in seams.get("module_patches", [])]
             + list(seams.get("client_methods", {})))
    assert any(k in n for n in names
               for k in ("price", "parlay", "calculate", "picks")), (
        f"{book}: no pricing seam declared, goldens would call the book")
