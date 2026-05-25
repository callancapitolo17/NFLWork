"""Integration tests for symmetric side selection in _evaluate_quote.

We exercise the dry-run path so we don't need to mock the accept_quote /
position-reconciliation chain — side selection happens before the dry-run
short-circuit, so quote_log captures the chosen side and ev_yes_pct /
ev_no_pct just the same as a live accept would.
"""
from __future__ import annotations

from datetime import datetime, timezone
from unittest.mock import patch

import duckdb
import pytest

from kalshi_mlb_rfq import db, main


@pytest.fixture
def tmpdb(tmp_path, monkeypatch):
    p = tmp_path / "state.duckdb"
    monkeypatch.setattr(db, "DB_PATH", p)
    db.init_database()
    return p


def _seed_rfq(tmpdb, rfq_id: str = "rfq-1"):
    """Insert one RFQ + combo_cache row so _evaluate_quote has a target."""
    combo = "KX-COMBO-1"
    legs_json = "[]"
    leg_set_hash = "hash-1"
    game_id = "GAME-1"
    con = duckdb.connect(str(tmpdb))
    try:
        con.execute(
            "INSERT INTO live_rfqs (rfq_id, combo_market_ticker, leg_set_hash, "
            "game_id, status, submitted_at) VALUES (?, ?, ?, ?, ?, ?)",
            [rfq_id, combo, leg_set_hash, game_id, "open",
             datetime.now(timezone.utc)],
        )
        con.execute(
            "INSERT INTO combo_cache (leg_set_hash, collection_ticker, "
            "combo_market_ticker, combo_event_ticker, legs_json, game_id) "
            "VALUES (?, ?, ?, ?, ?, ?)",
            [leg_set_hash, "MVE-1", combo, "EVT-1", legs_json, game_id],
        )
    finally:
        con.close()
    return {"rfq_id": rfq_id, "combo": combo, "leg_set_hash": leg_set_hash}


def _read_log(tmpdb, quote_id: str) -> dict:
    con = duckdb.connect(str(tmpdb), read_only=True)
    try:
        row = con.execute(
            "SELECT decision, chosen_side, ev_yes_pct, ev_no_pct, "
            "post_fee_ev_pct, blended_fair_at_eval "
            "FROM quote_log WHERE quote_id=?",
            [quote_id]
        ).fetchone()
    finally:
        con.close()
    if row is None:
        return {}
    keys = ["decision", "chosen_side", "ev_yes_pct", "ev_no_pct",
            "post_fee_ev_pct", "blended_fair_at_eval"]
    return dict(zip(keys, row))


def _eval_with_stubs(tmpdb, quote: dict, fair: float):
    """Run _evaluate_quote in dry-run mode with all pre-EV gates stubbed-out."""
    with patch.object(main, "_fresh_blended_fair", return_value=fair), \
         patch.object(main, "_all_per_accept_gates_pass",
                      return_value=(True, "passed")):
        main._evaluate_quote(quote, dry_run=True)


def test_only_yes_eligible_selects_yes(tmpdb):
    _seed_rfq(tmpdb)
    # fair=0.60, yes_ask=0.45 → EV_buy_yes = 0.60 - 0.45 - fee ≈ +0.13 (well >5%)
    # fair=0.60 ⇒ fair_no=0.40, no_ask=0.55 → EV_buy_no = 0.40 - 0.55 - fee < 0
    quote = {"id": "q-yes-only", "rfq_id": "rfq-1",
             "yes_bid_dollars": 0.45, "no_bid_dollars": 0.55,
             "creator_id": "lp-A"}
    _eval_with_stubs(tmpdb, quote, fair=0.60)
    row = _read_log(tmpdb, "q-yes-only")
    assert row["decision"] == "declined_dry_run"
    assert row["chosen_side"] == "yes"
    assert row["ev_yes_pct"] > 0.05
    assert row["ev_no_pct"] < 0


def test_only_no_eligible_selects_no(tmpdb):
    _seed_rfq(tmpdb)
    # fair=0.40 ⇒ fair_no=0.60, yes_ask=0.55 → EV_buy_yes = 0.40-0.55-fee < 0
    # fair=0.40 ⇒ fair_no=0.60, no_ask=0.45 → EV_buy_no = 0.60-0.45-fee ≈ +0.13
    quote = {"id": "q-no-only", "rfq_id": "rfq-1",
             "yes_bid_dollars": 0.55, "no_bid_dollars": 0.45,
             "creator_id": "lp-A"}
    _eval_with_stubs(tmpdb, quote, fair=0.40)
    row = _read_log(tmpdb, "q-no-only")
    assert row["decision"] == "declined_dry_run"
    assert row["chosen_side"] == "no"
    assert row["ev_yes_pct"] < 0
    assert row["ev_no_pct"] > 0.05


def test_neither_eligible_declines_ev(tmpdb):
    _seed_rfq(tmpdb)
    # fair sits in the middle of the LP's spread → both sides -EV after fees
    quote = {"id": "q-neither", "rfq_id": "rfq-1",
             "yes_bid_dollars": 0.40, "no_bid_dollars": 0.40,
             "creator_id": "lp-A"}
    _eval_with_stubs(tmpdb, quote, fair=0.50)
    row = _read_log(tmpdb, "q-neither")
    assert row["decision"] == "declined_ev"
    assert row["chosen_side"] is None


def test_math_invariant_guard_fires_on_impossible_both_eligible(tmpdb):
    """Force the impossible case by hand-crafting a quote where the LP
    'pays more for both sides than they collect' (yes_bid + no_bid > 1).
    Real LPs never quote this; the test proves the defensive guard fires."""
    _seed_rfq(tmpdb)
    # yes_bid + no_bid = 1.30 > 1 (impossible IRL). yes_ask = no_ask = 0.35.
    # Against fair=0.50, both sides clear MIN_EV_PCT by a wide margin.
    quote = {"id": "q-impossible", "rfq_id": "rfq-1",
             "yes_bid_dollars": 0.65, "no_bid_dollars": 0.65,
             "creator_id": "lp-broken"}
    _eval_with_stubs(tmpdb, quote, fair=0.50)
    row = _read_log(tmpdb, "q-impossible")
    assert row["decision"] == "declined_math_invariant"
    assert row["chosen_side"] is None
    assert row["ev_yes_pct"] > 0
    assert row["ev_no_pct"] > 0


def test_hedge_diag_fires_when_opposite_position_held(tmpdb):
    """If we already hold YES on combo X and the bot is about to fill NO on
    the same combo, the quote_log row gets tagged hedge_added=True with the
    projected combined-net P&L on the matched contracts."""
    info = _seed_rfq(tmpdb)
    # Pretend we already hold 3 YES contracts at $0.42
    con = duckdb.connect(str(tmpdb))
    try:
        con.execute(
            "INSERT INTO positions (combo_market_ticker, side, game_id, "
            "net_contracts, weighted_price, legs_json, updated_at) VALUES "
            "(?, ?, ?, ?, ?, ?, ?)",
            [info["combo"], "yes", "GAME-1", 3.0, 0.42, "[]",
             datetime.now(timezone.utc)],
        )
    finally:
        con.close()

    # Quote that puts NO in the money: fair=0.40, no_ask=0.45 → NO is +EV.
    quote = {"id": "q-hedge", "rfq_id": "rfq-1",
             "yes_bid_dollars": 0.55, "no_bid_dollars": 0.45,
             "creator_id": "lp-A"}
    _eval_with_stubs(tmpdb, quote, fair=0.40)

    con = duckdb.connect(str(tmpdb), read_only=True)
    try:
        row = con.execute(
            "SELECT chosen_side, hedge_added, hedge_original_side, "
            "hedge_original_price, hedge_new_price, hedge_projected_net "
            "FROM quote_log WHERE quote_id='q-hedge'"
        ).fetchone()
    finally:
        con.close()
    keys = ["chosen_side", "hedge_added", "hedge_original_side",
            "hedge_original_price", "hedge_new_price", "hedge_projected_net"]
    r = dict(zip(keys, row))

    assert r["chosen_side"] == "no"
    assert r["hedge_added"] is True
    assert r["hedge_original_side"] == "yes"
    assert r["hedge_original_price"] == pytest.approx(0.42)
    assert r["hedge_new_price"] == pytest.approx(0.45)  # no_ask = 1 - yes_bid
    # Projected net on matched contracts (before fees): 1 - 0.42 - 0.45 = 0.13
    assert r["hedge_projected_net"] == pytest.approx(0.13, abs=1e-6)


def test_hedge_diag_absent_when_no_opposite_position(tmpdb):
    """No existing position → hedge_added should remain NULL on the log row."""
    _seed_rfq(tmpdb)
    quote = {"id": "q-no-hedge", "rfq_id": "rfq-1",
             "yes_bid_dollars": 0.45, "no_bid_dollars": 0.55,
             "creator_id": "lp-A"}
    _eval_with_stubs(tmpdb, quote, fair=0.60)
    con = duckdb.connect(str(tmpdb), read_only=True)
    try:
        hedge_added = con.execute(
            "SELECT hedge_added FROM quote_log WHERE quote_id='q-no-hedge'"
        ).fetchone()[0]
    finally:
        con.close()
    assert hedge_added is None


def test_chosen_side_uses_min_ev_pct_threshold(tmpdb):
    """Side just under MIN_EV_PCT shouldn't qualify even if it's positive."""
    from kalshi_mlb_rfq import config
    _seed_rfq(tmpdb)
    # yes_ask = 0.50, fair = 0.515 → EV_yes ≈ +0.015 (≈3% pct, below 5% default)
    # no_ask  = 0.50, fair_no = 0.485 → EV_no ≈ -0.015
    quote = {"id": "q-marginal", "rfq_id": "rfq-1",
             "yes_bid_dollars": 0.50, "no_bid_dollars": 0.50,
             "creator_id": "lp-A"}
    _eval_with_stubs(tmpdb, quote, fair=0.515)
    row = _read_log(tmpdb, "q-marginal")
    # Whichever way the EV pct lands, the result should be consistent with
    # config.MIN_EV_PCT. We just assert no false-positive accept.
    if row["decision"] == "declined_dry_run":
        assert row["chosen_side"] in ("yes", "no")
        assert max(row["ev_yes_pct"] or -1, row["ev_no_pct"] or -1) >= config.MIN_EV_PCT
    else:
        assert row["decision"] == "declined_ev"
