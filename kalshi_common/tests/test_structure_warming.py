"""Issue #50 item 2: background structure warming — no pricing calls, ever.

``SGPService.warm_structures(games)`` runs match_event + build_structure for
every slate game at every book so the events/structure TTL caches (and
Caesars' WAF token) are warm before an RFQ arrives. The contract under test:

  1. the ``price`` hook is NEVER invoked (warming is structure-only);
  2. a warmed game's next ``price_on_demand`` pays zero wire structure calls
     (the cold-start penalty this ticket removes);
  3. each book writes one ``sgp_fetch_health`` row with ``path='warming'``;
  4. a transport failure during warming feeds the same 3-strike client
     accounting the sweep and on-demand paths share;
  5. a Caesars-style client gets its WAF token pre-warmed with enough
     remaining TTL to cover the warming cadence;
  6. ``sgp_runner.warm_cycle`` derives the slate from ``mlb_target_lines``.
"""
import inspect
from datetime import datetime, timezone
from types import SimpleNamespace

import duckdb
import pytest

from kalshi_common.sgp_service import SGPService
from mlb_sgp._shared import BookTransportError, GameRef, TargetLine

GAME = GameRef(game_id="g1", home_team="Boston Red Sox",
               away_team="New York Yankees", commence_time=None)
GAME2 = GameRef(game_id="g2", home_team="Chicago Cubs",
                away_team="St. Louis Cardinals", commence_time=None)


def _seam_hooks(price_fn=None, match=lambda c, g: "EV",
                build=lambda c, ev, g: {"fake": "structure"}):
    if price_fn is None:
        def price_fn(client, refs, event):     # pragma: no cover — guard
            raise AssertionError("warming must never price")
    return {"match_event": match,
            "build_structure": build,
            "resolve": lambda structure, legs, home, away: None,
            "price": price_fn}


def _seam_service(book="betmgm", hooks=None, **kwargs):
    svc = SGPService(books=(book,),
                     on_demand_hooks={book: hooks or _seam_hooks()},
                     **kwargs)
    # the hook seam never builds clients; warming must tolerate that
    svc._state[book].client = SimpleNamespace()
    return svc


# --------------------------------------------------------------------- #
# 1. structure-only                                                      #
# --------------------------------------------------------------------- #

def test_warm_structures_exists_and_never_calls_price():
    priced = []

    def price(client, refs, event):
        priced.append(refs)
        return 2.0

    svc = _seam_service(hooks=_seam_hooks(price_fn=price))
    warmed = svc.warm_structures([GAME, GAME2])
    assert priced == []                        # the whole point of the ticket
    assert warmed == {"betmgm": 2}


def test_warm_structures_skips_unlisted_games_without_failing():
    svc = _seam_service(hooks=_seam_hooks(
        match=lambda c, g: ("EV" if g.game_id == "g1" else None)))
    warmed = svc.warm_structures([GAME, GAME2])
    assert warmed == {"betmgm": 1}


def test_warm_structures_signature_is_keyword_safe():
    sig = inspect.signature(SGPService.warm_structures)
    extras = [p for name, p in sig.parameters.items()
              if name not in ("self", "games")]
    assert extras, "expected at least one tuning knob beyond games"
    assert all(p.kind is inspect.Parameter.KEYWORD_ONLY for p in extras)


# --------------------------------------------------------------------- #
# 2. warming removes the cold-structure penalty                          #
# --------------------------------------------------------------------- #

def test_price_on_demand_after_warming_is_warm(tmp_path, monkeypatch):
    from kalshi_common.tests.test_on_demand_cache_state import (
        _legs, _patched_mgm_service)
    db = str(tmp_path / "health.duckdb")
    svc, st = _patched_mgm_service(monkeypatch, db)

    svc.warm_structures([GAME])
    paid_by_warming = st.client.structure_calls
    assert paid_by_warming > 0                 # warming did the wire work

    res = svc.price_on_demand("betmgm", GAME, _legs())
    assert res is not None
    assert st.client.structure_calls == paid_by_warming   # RFQ paid nothing

    svc.health.flush()
    con = duckdb.connect(db, read_only=True)
    try:
        row = con.execute(
            "SELECT cache_state FROM sgp_fetch_health "
            "WHERE path = 'on_demand'").fetchone()
    finally:
        con.close()
    assert row == ("warm",)


# --------------------------------------------------------------------- #
# 3. health rows                                                         #
# --------------------------------------------------------------------- #

def test_warm_structures_records_warming_health_rows(tmp_path):
    db = str(tmp_path / "health.duckdb")
    svc = _seam_service(health_db_path=db)
    svc.warm_structures([GAME, GAME2])
    assert svc.health.flush() == 1             # one row per book per pass
    con = duckdb.connect(db, read_only=True)
    try:
        rows = con.execute(
            "SELECT book, path, outcome FROM sgp_fetch_health").fetchall()
    finally:
        con.close()
    assert rows == [("betmgm", "warming", "ok")]


# --------------------------------------------------------------------- #
# 4. transport failure -> shared strike accounting                       #
# --------------------------------------------------------------------- #

def test_warming_transport_error_counts_a_strike(tmp_path):
    def match(client, game):
        raise BookTransportError("events feed 403", stage="events",
                                 status_code=403)

    db = str(tmp_path / "health.duckdb")
    svc = _seam_service(hooks=_seam_hooks(match=match), health_db_path=db)
    warmed = svc.warm_structures([GAME])
    assert warmed == {"betmgm": 0}
    assert svc._state["betmgm"].failures == 1  # same path as sweep/on-demand
    svc.health.flush()
    con = duckdb.connect(db, read_only=True)
    try:
        row = con.execute(
            "SELECT outcome, error_class FROM sgp_fetch_health").fetchone()
    finally:
        con.close()
    assert row[0] == "transport_error"
    assert "403" in row[1]


# --------------------------------------------------------------------- #
# 5. Caesars token pre-warm                                              #
# --------------------------------------------------------------------- #

def test_warming_prewarms_a_token_bearing_client():
    asked = []

    class FakeCZRClient:
        def ensure_fresh_token(self, *, min_remaining_sec):
            asked.append(min_remaining_sec)
            return True

    svc = SGPService(books=("caesars",),
                     on_demand_hooks={"caesars": _seam_hooks()})
    svc._state["caesars"].client = FakeCZRClient()
    svc.warm_structures([GAME])
    assert len(asked) == 1
    # remaining TTL must outlive a warming cadence gap, not just be nonzero
    assert asked[0] >= 120.0


def test_warming_tolerates_clients_without_token_prewarm():
    svc = _seam_service()                      # SimpleNamespace client
    assert svc.warm_structures([GAME]) == {"betmgm": 1}


# --------------------------------------------------------------------- #
# 6. warm_cycle: slate from mlb_target_lines                             #
# --------------------------------------------------------------------- #

def test_warm_cycle_derives_games_from_target_lines(tmp_path):
    from kalshi_common import sgp_runner
    db = str(tmp_path / "market.duckdb")
    ct = datetime(2026, 8, 5, 23, 5, tzinfo=timezone.utc)
    targets = [
        TargetLine(game_id="g1", home_team="Boston Red Sox",
                   away_team="New York Yankees", commence_time=ct,
                   period="FG", spread=-1.5, total=8.5),
        TargetLine(game_id="g1", home_team="Boston Red Sox",
                   away_team="New York Yankees", commence_time=ct,
                   period="FG", spread=-2.5, total=9.5),
        TargetLine(game_id="g2", home_team="Chicago Cubs",
                   away_team="St. Louis Cardinals", commence_time=ct,
                   period="F5", spread=-0.5, total=4.5),
    ]
    sgp_runner.write_target_lines(targets, db_path=db)

    seen = []

    class FakeService:
        def warm_structures(self, games):
            seen.extend(games)
            return {"betmgm": len(games)}

    out = sgp_runner.warm_cycle(bot_market_db=db, service=FakeService())
    assert out == {"betmgm": 2}
    # one GameRef per game, not per target line
    assert sorted(g.game_id for g in seen) == ["g1", "g2"]
    g1 = next(g for g in seen if g.game_id == "g1")
    assert g1.home_team == "Boston Red Sox"
    assert g1.away_team == "New York Yankees"
    # commence_time round-trips through the naive-UTC mlb_target_lines
    # convention (what the legacy hour-bucket matchers expect)
    assert g1.commence_time is not None
    assert g1.commence_time.replace(tzinfo=None).hour == 23


def test_warm_cycle_params_are_keyword_only():
    from kalshi_common import sgp_runner
    sig = inspect.signature(sgp_runner.warm_cycle)
    assert all(p.kind is inspect.Parameter.KEYWORD_ONLY
               for p in sig.parameters.values())
    with pytest.raises(TypeError):
        sig.bind("db_path", object())


# --------------------------------------------------------------------- #
# 7. no cold window: warming REFRESHES entries older than one cadence     #
# --------------------------------------------------------------------- #
# get_or_fetch does not extend TTL on a hit, so warming that only fills
# missing entries leaves a cold window: entry fetched at t=0, pass at 120
# hits (age 120 < TTL), entry expires at TTL, next pass refreshes at 240 —
# every RFQ in [TTL, 240) pays cold. The pass must evict-then-refetch
# entries older than the cadence so entry age never reaches the TTL.

def test_ttlcache_evict_older_than():
    from mlb_sgp._shared import TTLCache
    clock = SimpleNamespace(t=1000.0)
    cache = TTLCache(ttl_sec=420, now_fn=lambda: clock.t)
    cache.get_or_fetch("old", lambda: "a")
    clock.t += 130
    cache.get_or_fetch("young", lambda: "b")
    assert cache.evict_older_than(120) == 1          # only "old"
    fetches = []
    cache.get_or_fetch("old", lambda: fetches.append(1) or "a2")
    cache.get_or_fetch("young", lambda: fetches.append(1) or "b2")
    assert fetches == [1]                            # "young" still cached


def test_warming_pass_refetches_structures_older_than_refresh_age(
        tmp_path, monkeypatch):
    from mlb_sgp._shared import TTLCache
    from kalshi_common.tests.test_on_demand_cache_state import (
        _legs, _patched_mgm_service)
    clock = SimpleNamespace(t=1000.0)
    svc, st = _patched_mgm_service(monkeypatch, str(tmp_path / "h.duckdb"))
    # inject a controllable clock into this book's caches (pre-seed wins:
    # _structure_caches uses setdefault)
    st.caches = {"events": TTLCache(900, now_fn=lambda: clock.t),
                 "structure": TTLCache(420, now_fn=lambda: clock.t)}

    svc.warm_structures([GAME])
    after_first = st.client.structure_calls
    assert after_first > 0

    clock.t += 60                                    # young: pass is a no-op
    svc.warm_structures([GAME])
    assert st.client.structure_calls == after_first

    clock.t += 70                                    # age 130 > cadence 120
    svc.warm_structures([GAME])
    assert st.client.structure_calls > after_first   # refreshed, not hit

    # and the RFQ right after still pays nothing
    before = st.client.structure_calls
    res = svc.price_on_demand("betmgm", GAME, _legs())
    assert res is not None
    assert st.client.structure_calls == before


def test_warm_cycle_empty_target_table_is_a_noop(tmp_path):
    from kalshi_common import sgp_runner
    db = str(tmp_path / "market.duckdb")
    sgp_runner.write_target_lines([], db_path=db)

    class NeverCalled:
        def warm_structures(self, games):      # pragma: no cover — guard
            raise AssertionError("no games -> no warming pass")

    assert sgp_runner.warm_cycle(bot_market_db=db, service=NeverCalled()) == {}
