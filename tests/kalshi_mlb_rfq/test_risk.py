from datetime import datetime, timedelta, timezone

from kalshi_mlb_rfq import risk


# ---- Basic gates ----------------------------------------------------------

def test_staleness_gate_blocks_old_data():
    fresh = datetime.now(timezone.utc) - timedelta(seconds=300)  # 5 min
    old   = datetime.now(timezone.utc) - timedelta(seconds=900)  # 15 min
    assert risk.staleness_ok(fresh, max_age_sec=600)
    assert not risk.staleness_ok(old, max_age_sec=600)


def test_staleness_gate_handles_negative_age():
    future = datetime.now(timezone.utc) + timedelta(seconds=60)
    assert not risk.staleness_ok(future, max_age_sec=600)


def test_fair_bounds_gate():
    assert risk.fair_in_bounds(0.30, lower=0.05, upper=0.95)
    assert not risk.fair_in_bounds(0.03, lower=0.05, upper=0.95)
    assert not risk.fair_in_bounds(0.97, lower=0.05, upper=0.95)


def test_sanity_bound_gate():
    assert risk.sanity_bound_ok(quote_implied=0.27, blended_fair=0.30, max_deviation=0.15)
    assert not risk.sanity_bound_ok(quote_implied=0.10, blended_fair=0.30, max_deviation=0.15)


def test_tipoff_gate():
    now = datetime.now(timezone.utc)
    far_future_game = now + timedelta(minutes=60)
    near_game     = now + timedelta(minutes=3)
    assert risk.tipoff_ok(commence_time=far_future_game, cancel_min=5, now=now)
    assert not risk.tipoff_ok(commence_time=near_game, cancel_min=5, now=now)
    assert not risk.tipoff_ok(commence_time=None, cancel_min=5, now=now)


# ---- Exposure ------------------------------------------------------------

def test_per_game_cap_blocks_at_threshold():
    # bankroll=1000, MAX_GAME_EXPOSURE_PCT=0.10 → cap = $100/game
    today_fills = [
        {"game_id": "G1", "contracts": 100, "price_dollars": 0.30},
        {"game_id": "G1", "contracts": 200, "price_dollars": 0.20},
        {"game_id": "G2", "contracts": 100, "price_dollars": 0.50},
    ]
    assert risk.per_game_cap_ok(game_id="G1", today_fills=today_fills,
                                 bankroll=1000.0, max_pct=0.10)
    today_fills.append({"game_id": "G1", "contracts": 100, "price_dollars": 0.40})
    assert not risk.per_game_cap_ok(game_id="G1", today_fills=today_fills,
                                     bankroll=1000.0, max_pct=0.10)


def test_daily_cap_blocks_at_threshold():
    today_fills = [
        {"game_id": "G1", "contracts": 100, "price_dollars": 0.50},
        {"game_id": "G2", "contracts": 200, "price_dollars": 0.50},
    ]
    assert risk.daily_cap_ok(today_fills=today_fills, daily_cap_usd=200.0)
    today_fills.append({"game_id": "G3", "contracts": 200, "price_dollars": 0.50})
    assert not risk.daily_cap_ok(today_fills=today_fills, daily_cap_usd=200.0)


def test_cooldown_gate():
    now = datetime.now(timezone.utc)
    cooled = {"abc": now.replace(year=now.year + 1)}
    assert not risk.cooldown_ok(leg_set_hash="abc", cooldown_map=cooled, now=now)
    assert risk.cooldown_ok(leg_set_hash="other", cooldown_map=cooled, now=now)


def test_inverse_combo_guard():
    legs = [
        {"market_ticker": "X", "side": "yes"},
        {"market_ticker": "Y", "side": "no"},
    ]
    inverse_legs = [
        {"market_ticker": "X", "side": "no"},
        {"market_ticker": "Y", "side": "yes"},
    ]
    assert not risk.inverse_combo_ok(
        legs=legs, open_position_hashes={risk.inverse_leg_set_hash(inverse_legs)})


# ---- Advanced ------------------------------------------------------------

def test_line_move_detection():
    ref = {"spread": -1.5, "total": 8.5}
    assert risk.line_move_ok(ref_lines=ref, current_lines={"spread": -1.5, "total": 8.5},
                             threshold=0.5)
    assert not risk.line_move_ok(ref_lines=ref, current_lines={"spread": -2.0, "total": 8.5},
                                  threshold=0.5)
    assert risk.line_move_ok(ref_lines=ref, current_lines={"spread": -1.7, "total": 8.5},
                             threshold=0.5)


def test_kill_switch_present(tmp_path, monkeypatch):
    kill_path = tmp_path / ".kill"
    monkeypatch.setattr(risk, "KILL_FILE", kill_path)
    assert risk.kill_switch_ok()
    kill_path.touch()
    assert not risk.kill_switch_ok()


def test_fill_ratio_halt():
    log = [{"decision": "accepted"}] * 25 + [{"decision": "failed_quote_walked"}] * 25
    assert risk.fill_ratio_ok(quote_log_window=log, min_ratio=0.50)

    log = [{"decision": "accepted"}] * 20 + [{"decision": "failed_quote_walked"}] * 30
    assert not risk.fill_ratio_ok(quote_log_window=log, min_ratio=0.50)


def test_fill_ratio_ok_when_window_too_small():
    log = [{"decision": "accepted"}] * 5
    assert risk.fill_ratio_ok(quote_log_window=log, min_ratio=0.50, min_window=10)
