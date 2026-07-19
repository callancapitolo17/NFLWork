import datetime
import pytest
from unabated_edge import config
from unabated_edge.feed import EventMeta
from unabated_edge.maker import engine as mengine, state as mstate
from unabated_edge.sports.soccer import Soccer

_NOW = datetime.datetime(2026, 7, 11, 12, 0, tzinfo=datetime.timezone.utc)
_KICK = _NOW + datetime.timedelta(hours=3)
_EM = EventMeta(event_id=1, league_key="lg21", start_utc=_KICK,
                home_id=10, away_id=20, home="TeamA", away="TeamB")
_KEV = {"title": "A vs B: Regulation Time Total Goals",
        "markets": [{"ticker": "T-O25", "strike_type": "greater", "floor_strike": 2.5}]}
_LADDER = {2.5: {"p_over": 0.42, "book": 7, "alt": False, "overround": 1.05,
                "modified_on": _NOW.isoformat()}}
# wide crowd: 0.30 yes bid / 0.55 yes ask (no bid 0.45) -> our quotes fit inside
_BOOK = {"yes_bids": [(0.30, 500.0)], "no_bids": [(0.45, 500.0)]}
_BOOKS = {"T-O25": _BOOK}


class FakeGateway:
    is_live = False
    def __init__(self):
        self.placed, self.cancelled, self._n = [], [], 0
    def place(self, ticker, side, price_cents, count, client_order_id):
        self._n += 1
        self.placed.append((ticker, side, price_cents, count))
        return f"o-{self._n}"
    def cancel(self, order_id):
        self.cancelled.append(order_id)
        return True


@pytest.fixture
def eng(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "MAKER_DB_PATH", tmp_path / "mk.duckdb")
    from unabated_edge.maker import store
    store.init()
    return mengine.MakerEngine(FakeGateway(), mstate.MakerState())


def test_rests_two_sides_at_fair_minus_margin(eng):
    eng.on_match(Soccer(), _EM, _KEV, _LADDER, _BOOKS, _NOW)
    sides = {(s, p) for _t, s, p, _n in eng.gateway.placed}
    # fee(0.42)*100+1 ≈ 1.44c, roi 3%*42 = 1.26c -> m=ceil(1.44)=2c -> yes 40c
    # no side: fair 58c, fee≈1.43c+1, roi 1.74c -> m=2c -> 56c, capped by ask? no ask=1-0.30=0.70 -> cap 69 -> 56c
    assert ("yes", 40) in sides and ("no", 56) in sides


def test_holds_same_price_no_churn(eng):
    eng.on_match(Soccer(), _EM, _KEV, _LADDER, _BOOKS, _NOW)
    n = len(eng.gateway.placed)
    eng.on_match(Soccer(), _EM, _KEV, _LADDER, _BOOKS, _NOW + datetime.timedelta(seconds=5))
    assert len(eng.gateway.placed) == n and eng.gateway.cancelled == []


def test_requotes_on_fair_move(eng):
    eng.on_match(Soccer(), _EM, _KEV, _LADDER, _BOOKS, _NOW)
    moved = {2.5: {**_LADDER[2.5], "p_over": 0.47}}
    eng.on_match(Soccer(), _EM, _KEV, moved, _BOOKS, _NOW + datetime.timedelta(seconds=5))
    assert len(eng.gateway.cancelled) == 2          # both sides replaced
    assert any(p == 45 for _t, s, p, _n in eng.gateway.placed if s == "yes")  # 47 - 2c margin


def test_never_crosses_the_ask(eng):
    tight = {"T-O25": {"yes_bids": [(0.30, 500.0)], "no_bids": [(0.59, 500.0)]}}  # yes ask 41c
    eng.on_match(Soccer(), _EM, _KEV, _LADDER, tight, _NOW)
    yes = [p for _t, s, p, _n in eng.gateway.placed if s == "yes"]
    assert yes == [40]                              # min(40, 41-1) = 40 still fine
    tighter = {"T-O25": {"yes_bids": [(0.30, 500.0)], "no_bids": [(0.64, 500.0)]}}  # yes ask 36c
    eng.on_match(Soccer(), _EM, _KEV, _LADDER, tighter, _NOW + datetime.timedelta(seconds=5))
    # 36-1=35 < fair 42 - MAX_MARGIN 5 = 37 -> crowd tighter than we'll quote: cancel, don't rest
    assert any(o for o in eng.gateway.cancelled)
    assert not any(p == 35 for _t, s, p, _n in eng.gateway.placed if s == "yes")


def test_alt_rung_gated_on_overround(eng):
    bad_alt = {2.5: {"p_over": 0.42, "book": 7, "alt": True, "overround": 1.30,
                     "modified_on": _NOW.isoformat()}}
    eng.on_match(Soccer(), _EM, _KEV, bad_alt, _BOOKS, _NOW)
    assert eng.gateway.placed == []


def test_pull_window_cancels_everything(eng):
    eng.on_match(Soccer(), _EM, _KEV, _LADDER, _BOOKS, _NOW)
    near_kick = _KICK - datetime.timedelta(minutes=2)
    eng.on_match(Soccer(), _EM, _KEV, _LADDER, _BOOKS, near_kick)
    assert len(eng.gateway.cancelled) == 2 and eng.state.quotes_live() == 0


def test_fill_burst_triggers_cooloff(eng):
    eng.state.fill_times[1] = [_NOW - datetime.timedelta(seconds=i) for i in range(config.FILL_BURST_N + 1)]
    eng.on_match(Soccer(), _EM, _KEV, _LADDER, _BOOKS, _NOW)
    assert eng.gateway.placed == []
    assert eng.state.cooloff_until[1] > _NOW


def test_ledger_caps_size(eng, monkeypatch):
    monkeypatch.setattr(config, "MATCH_CAP_PCT", 0.02)      # $20 cap
    eng.on_match(Soccer(), _EM, _KEV, _LADDER, _BOOKS, _NOW)
    # yes quote at 40c: alone it could take floor(20/0.40)=50, but the no-side
    # resting quote is also assumed filled; every placement stays within cap
    for _t, s, p, n in eng.gateway.placed:
        assert n * (p / 100.0) <= 300.0 + 1e-9              # never above MAX_QUOTE_PCT
    from unabated_edge.maker import ledger
    assert ledger.worst_case(eng.state.exposure_fills(1)) >= -20.0 - 1e-9


def test_watchdog_pulls_on_stale_feed(eng):
    eng.on_match(Soccer(), _EM, _KEV, _LADDER, _BOOKS, _NOW)
    eng.note_success("soccer", _NOW)
    eng.watchdog(_NOW + datetime.timedelta(seconds=config.MAX_STALENESS_SEC + 5))
    assert eng.state.quotes_live() == 0


def test_sweep_pulls_unpaired_events(eng):
    eng.on_match(Soccer(), _EM, _KEV, _LADDER, _BOOKS, _NOW)
    eng.sweep("soccer", set(), _NOW)                        # event 1 not seen this tick
    assert eng.state.quotes_live() == 0


def test_daily_halt_pulls_everything(eng, monkeypatch):
    eng.on_match(Soccer(), _EM, _KEV, _LADDER, _BOOKS, _NOW)
    eng.state.roll_day(_NOW)
    eng.state.settled_pnl_today = -config.DAILY_LOSS_HALT_PCT * config.BANKROLL - 1
    eng.on_match(Soccer(), _EM, _KEV, _LADDER, _BOOKS, _NOW + datetime.timedelta(seconds=5))
    assert eng.state.quotes_live() == 0
    assert eng.stats()["halted"] is True


def test_fixed_contract_cap(eng):
    """_BOOK gives plenty of room; every placed quote must still be capped at
    MAKER_MAX_CONTRACTS (the tuition-run ceiling), regardless of budget room."""
    eng.on_match(Soccer(), _EM, _KEV, _LADDER, _BOOKS, _NOW)
    assert eng.gateway.placed                                   # sanity: quotes actually placed
    for _t, _s, _p, n in eng.gateway.placed:
        assert n <= config.MAKER_MAX_CONTRACTS


def test_hard_stop_pulls_all(eng):
    eng.state.fills = {1: [(2.5, "yes", 100, 0.90)]}
    eng._fair_by_event = {1: {2.5: 0.10}}   # unrealized ~= 100*(0.10-0.90) = -80 < -50
    eng.on_match(Soccer(), _EM, _KEV, _LADDER, _BOOKS, _NOW)
    assert eng.gateway.placed == []
    assert eng.stats()["halted"] is True


def test_anchor_stale_pulls_near_kickoff(eng):
    """Within ANCHOR_STALE_FARK_SEC of kickoff the sharp total should churn, so a
    10-min-old anchor means the feed is lagging -> pull."""
    near_kick = EventMeta(event_id=1, league_key="lg21",
                          start_utc=_NOW + datetime.timedelta(minutes=30),
                          home_id=10, away_id=20, home="TeamA", away="TeamB")
    stale = {2.5: {**_LADDER[2.5],
                  "modified_on": (_NOW - datetime.timedelta(minutes=10)).isoformat()}}
    eng.on_match(Soccer(), near_kick, _KEV, stale, _BOOKS, _NOW)
    assert eng.gateway.placed == []


def test_anchor_stale_tolerated_far_from_kickoff(eng):
    """Far from kickoff (_EM is 3h out, beyond the 2h FARK cutoff) the sharp
    total legitimately sits unchanged for hours -> quote against the latest
    number anyway; the poll-success watchdog is the dead-feed guard here."""
    stale = {2.5: {**_LADDER[2.5],
                  "modified_on": (_NOW - datetime.timedelta(minutes=90)).isoformat()}}
    eng.on_match(Soccer(), _EM, _KEV, stale, _BOOKS, _NOW)
    assert eng.gateway.placed


def test_anchor_no_timestamp_always_pulls(eng):
    """No parseable modifiedOn on any rung -> can't prove provenance -> pull,
    even far from kickoff."""
    no_ts = {2.5: {**_LADDER[2.5], "modified_on": None}}
    eng.on_match(Soccer(), _EM, _KEV, no_ts, _BOOKS, _NOW)
    assert eng.gateway.placed == []


def test_anchor_fresh_still_quotes(eng):
    fresh = {2.5: {**_LADDER[2.5], "modified_on": _NOW.isoformat()}}
    eng.on_match(Soccer(), _EM, _KEV, fresh, _BOOKS, _NOW)
    assert eng.gateway.placed


def test_no_crowd_side_skipped(eng):
    """no_bids empty -> yes_ask_from_book is None -> the YES side must be
    skipped (no_crowd), not quoted at some unconstrained price."""
    book = {"T-O25": {"yes_bids": [(0.30, 5.0)], "no_bids": []}}
    eng.on_match(Soccer(), _EM, _KEV, _LADDER, book, _NOW)
    assert not any(s == "yes" for _t, s, _p, _n in eng.gateway.placed)


def test_baseline_blocked(eng):
    eng.state.position_baseline = {"T-O25": 5.0}
    eng.on_match(Soccer(), _EM, _KEV, _LADDER, _BOOKS, _NOW)
    assert eng.gateway.placed == []


# ---------- touch-join (spec 2026-07-18) ----------

def _tj(eng, book, fair_cents, ticker="T-O25", side="yes", alt=False, opp_ask_cents=99):
    return eng._touch_join_price(ticker, side, fair_cents, alt, book, opp_ask_cents)

def test_touch_join_entry_within_band(eng):
    book = {"yes_bids": [(0.78, 500.0)], "no_bids": [(0.21, 500.0)]}
    assert _tj(eng, book, 80.0) == 78          # edge ~1.7c in [1, 5]

def test_touch_join_rejects_thin_and_suspicious_edges(eng):
    book = {"yes_bids": [(0.78, 500.0)], "no_bids": []}
    assert _tj(eng, book, 79.0) is None        # edge ~0.7c < 1c entry
    deep = {"yes_bids": [(0.60, 500.0)], "no_bids": []}
    assert _tj(eng, deep, 80.0) is None        # edge ~19c > MAX_MARGIN suspicion guard

def test_touch_join_alt_needs_bigger_edge(eng):
    book = {"yes_bids": [(0.78, 500.0)], "no_bids": []}
    assert _tj(eng, book, 79.5, alt=True) is None   # ~1.2c < 1.5c alt entry
    assert _tj(eng, book, 80.0, alt=True) == 78     # ~1.5c >= 1.5c

def test_touch_join_self_exclusion(eng):
    # our own 2-lot at 78c is the only order at the level -> not "the crowd"
    eng.state.on_place("T-O25", "yes", "o-x", 78, 2, mode="quote")
    book = {"yes_bids": [(0.78, 2.0), (0.63, 500.0)], "no_bids": []}
    assert _tj(eng, book, 80.0) is None        # real touch 63 -> edge 16.7c: suspicious

def test_touch_join_hysteresis_hold_and_exit(eng):
    eng.state.on_place("T-O25", "yes", "o-x", 78, 2, mode="touch_join")
    book = {"yes_bids": [(0.78, 2.0)], "no_bids": []}   # only us at the level
    assert _tj(eng, book, 79.0) == 78          # edge ~0.7c in [0.25, 5] -> hold
    assert _tj(eng, book, 78.2) is None        # edge ~-0.1c < 0.25c -> exit

def test_touch_join_hold_ignored_for_legacy_orders(eng):
    eng.state.on_place("T-O25", "yes", "o-x", 78, 2, mode="quote")
    book = {"yes_bids": [(0.78, 2.0)], "no_bids": []}
    assert _tj(eng, book, 79.0) is None        # hysteresis only holds touch-joins

def test_touch_join_never_crosses(eng):
    book = {"yes_bids": [(0.78, 500.0)], "no_bids": []}
    assert _tj(eng, book, 80.0, opp_ask_cents=78) is None   # touch == ask-0 -> cross
