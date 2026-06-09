from datetime import datetime, timedelta, timezone
from kalshi_mlb_mm import risk

NOW = datetime(2026, 5, 27, 18, 0, tzinfo=timezone.utc)

def test_staleness_ok_true_when_fresh():
    assert risk.staleness_ok(NOW - timedelta(seconds=100), 600, now=NOW) is True

def test_staleness_ok_false_when_stale():
    assert risk.staleness_ok(NOW - timedelta(seconds=700), 600, now=NOW) is False

def test_staleness_ok_false_when_none():
    assert risk.staleness_ok(None, 600, now=NOW) is False

def test_staleness_ok_handles_naive_local_timestamp():
    # Regression: R writes naive LOCAL time into mlb_samples_meta.generated_at.
    # A fresh naive timestamp must read as FRESH, not be off by the UTC offset
    # (the old code forced naive->UTC, which made every sample look hours stale
    # in any non-UTC timezone and the bot never quoted).
    naive_recent = datetime.now() - timedelta(seconds=60)  # naive local, no tzinfo
    assert risk.staleness_ok(naive_recent, 600) is True
    naive_old = datetime.now() - timedelta(seconds=900)
    assert risk.staleness_ok(naive_old, 600) is False

def test_size_gate():
    assert risk.size_ok(5, max_contracts=5) is True
    assert risk.size_ok(6, max_contracts=5) is False
    assert risk.size_ok(0, max_contracts=5) is False

def test_book_move_circuit_breaker():
    assert risk.book_move_triggered(0.40, 0.45, threshold=0.03) is True
    assert risk.book_move_triggered(0.40, 0.41, threshold=0.03) is False

def test_last_look_ok_passes_when_still_ev_and_no_drift():
    assert risk.last_look_ok(side="yes", price=0.52, fee=0.005, current_fair=0.55,
                             prev_fair=0.55, drift_tol=0.02) is True

def test_last_look_voids_when_not_ev():
    assert risk.last_look_ok(side="yes", price=0.52, fee=0.005, current_fair=0.50,
                             prev_fair=0.55, drift_tol=0.02) is False

def test_last_look_voids_on_excess_drift():
    # still +EV (fair 0.60 > price+fee) but fair moved 0.05 > tol 0.02 -> void
    assert risk.last_look_ok(side="yes", price=0.52, fee=0.005, current_fair=0.60,
                             prev_fair=0.55, drift_tol=0.02) is False

def test_daily_cap():
    assert risk.daily_cap_ok([{"price": 4.0}, {"price": 5.0}], cap_usd=375.0) is True
    assert risk.daily_cap_ok([{"price": 380.0}], cap_usd=375.0) is False

def test_per_game_cap():
    fills = [{"game_id": "g1", "price": 40.0}, {"game_id": "g2", "price": 10.0}]
    assert risk.per_game_cap_ok("g1", fills, bankroll=500.0, pct=0.10) is True   # $40 < $50
    assert risk.per_game_cap_ok("g1", fills + [{"game_id": "g1", "price": 15.0}], bankroll=500.0, pct=0.10) is False  # $55 > $50
