"""Per-accept gates and risk sweep for the Kalshi MLB RFQ bot."""

import hashlib
from datetime import datetime, timedelta, timezone

from kalshi_mlb_rfq.config import KILL_FILE


# ---- Basic gates ----------------------------------------------------------

def staleness_ok(generated_at: datetime, max_age_sec: int) -> bool:
    """True if data is younger than max_age_sec. False on negative-age (clock skew)."""
    if generated_at.tzinfo is None:
        generated_at = generated_at.replace(tzinfo=timezone.utc)
    age = (datetime.now(timezone.utc) - generated_at).total_seconds()
    if age < 0:
        return False  # clock skew — fail safe
    return age <= max_age_sec


def fair_in_bounds(blended_fair: float, lower: float, upper: float) -> bool:
    return lower <= blended_fair <= upper


def sanity_bound_ok(quote_implied: float, blended_fair: float,
                    max_deviation: float) -> bool:
    return abs(quote_implied - blended_fair) <= max_deviation


def tipoff_ok(commence_time: datetime | None, cancel_min: int,
              now: datetime | None = None) -> bool:
    if commence_time is None:
        return False  # unknown tipoff = fail-safe refuse
    if commence_time.tzinfo is None:
        commence_time = commence_time.replace(tzinfo=timezone.utc)
    if now is None:
        now = datetime.now(timezone.utc)
    return (commence_time - now) > timedelta(minutes=cancel_min)


# ---- Exposure / cooldown / inverse ---------------------------------------

def per_game_cap_ok(game_id: str, today_fills: list[dict],
                     bankroll: float, max_pct: float) -> bool:
    cap = bankroll * max_pct
    spent = sum(f["contracts"] * f["price_dollars"]
                for f in today_fills if f["game_id"] == game_id)
    return spent < cap


def daily_cap_ok(today_fills: list[dict], daily_cap_usd: float) -> bool:
    spent = sum(f["contracts"] * f["price_dollars"] for f in today_fills)
    return spent < daily_cap_usd


def cooldown_ok(leg_set_hash: str, cooldown_map: dict, now: datetime | None = None) -> bool:
    if now is None:
        now = datetime.now(timezone.utc)
    cooled_until = cooldown_map.get(leg_set_hash)
    if cooled_until is None:
        return True
    if cooled_until.tzinfo is None:
        cooled_until = cooled_until.replace(tzinfo=timezone.utc)
    return cooled_until <= now


def inverse_leg_set_hash(legs: list[dict]) -> str:
    """Hash of leg-set after flipping every side. Used to detect held inverse positions."""
    flipped = [{"market_ticker": l["market_ticker"],
                "side": "no" if l["side"] == "yes" else "yes"} for l in legs]
    keys = sorted(f"{l['market_ticker']}|{l['side']}" for l in flipped)
    return hashlib.sha256("\n".join(keys).encode()).hexdigest()


def inverse_combo_ok(legs: list[dict], open_position_hashes: set[str]) -> bool:
    """True if no open position exists on the leg-set's inverse."""
    return inverse_leg_set_hash(legs) not in open_position_hashes


# ---- Advanced gates -------------------------------------------------------

def line_move_ok(ref_lines: dict, current_lines: dict, threshold: float) -> bool:
    """All lines moved by less than threshold? (None entries treated as no-data, blocks accept.)

    Uses >= for the threshold comparison so a move of exactly threshold trips the gate
    (conservative: matches CBB MM's behavior).
    """
    for key, ref in ref_lines.items():
        cur = current_lines.get(key)
        if cur is None or ref is None:
            return False
        if abs(cur - ref) >= threshold:
            return False
    return True


def kill_switch_ok() -> bool:
    return not KILL_FILE.exists()


def fill_ratio_ok(quote_log_window: list[dict], min_ratio: float,
                  min_window: int = 10) -> bool:
    """True if we have enough samples AND the accepted-fraction is >= min_ratio."""
    relevant = [r for r in quote_log_window
                if r["decision"] in ("accepted", "failed_quote_walked")]
    if len(relevant) < min_window:
        return True   # not enough data; don't halt
    accepted = sum(1 for r in relevant if r["decision"] == "accepted")
    return accepted / len(relevant) >= min_ratio
