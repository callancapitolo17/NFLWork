"""Risk gates for the maker bot. Pure functions, time injected for testability."""
from datetime import datetime, timedelta, timezone


def _now_matching(dt: datetime, now: datetime | None) -> datetime:
    """Return a `now` comparable to `dt`.

    CRITICAL: a NAIVE `dt` is treated as LOCAL time — the R answer-key pipeline
    writes naive local time into `mlb_samples_meta.generated_at`, and the Kalshi
    line tables likewise. Forcing naive→UTC (the old bug) made every sample look
    ~offset hours old and the bot never quoted. Aware `dt` is compared in UTC.
    An injected `now` (tests) is aligned to `dt`'s awareness. Mirrors the taker's
    `kalshi_mlb_rfq.risk.staleness_ok`.
    """
    if now is None:
        return datetime.now(timezone.utc) if dt.tzinfo is not None else datetime.now()
    if dt.tzinfo is not None and now.tzinfo is None:
        now = now.replace(tzinfo=timezone.utc)
    elif dt.tzinfo is None and now.tzinfo is not None:
        now = now.astimezone().replace(tzinfo=None)
    return now


def staleness_ok(generated_at: datetime | None, max_age_sec: int,
                 now: datetime | None = None) -> bool:
    if generated_at is None:
        return False
    now = _now_matching(generated_at, now)
    age = (now - generated_at).total_seconds()
    if age < 0:
        return False  # clock skew — fail safe
    return age <= max_age_sec


def tipoff_ok(commence_time: datetime | None, cancel_min: int,
              now: datetime | None = None) -> bool:
    if commence_time is None:
        return False
    now = _now_matching(commence_time, now)
    return now < commence_time - timedelta(minutes=cancel_min)


def size_ok(requested_contracts: int, max_contracts: int) -> bool:
    return 0 < requested_contracts <= max_contracts


def book_move_triggered(prev_fair: float, current_fair: float, threshold: float) -> bool:
    """True if the underlying book fair jumped by more than `threshold` between scrapes."""
    return abs(current_fair - prev_fair) > threshold


def last_look_ok(side: str, price: float, fee: float, current_fair: float,
                 prev_fair: float, drift_tol: float) -> bool:
    """Confirm only if the filled side is still +EV against current fair AND fair
    hasn't drifted against it past tolerance."""
    p = current_fair if side == "yes" else (1.0 - current_fair)
    if p <= price + fee:                      # no longer +EV
        return False
    if abs(current_fair - prev_fair) > drift_tol:
        return False
    return True


def daily_cap_ok(today_fills: list[dict], cap_usd: float) -> bool:
    used = sum(float(f.get("price", 0)) for f in today_fills)
    return used < cap_usd


def per_game_cap_ok(game_id: str, today_fills: list[dict], bankroll: float,
                    pct: float) -> bool:
    used = sum(float(f.get("price", 0)) for f in today_fills
               if f.get("game_id") == game_id)
    return used < bankroll * pct


def kill_switch_ok(kill_file) -> bool:
    return not kill_file.exists()
