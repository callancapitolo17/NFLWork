"""One RFQ -> one decision. Pure: the pricer call is injected (the leg table
for the first pass, the Kalshi-checked runner call for a confirm).

Phase 1 sends nothing, so the outcome is "would_quote" (every leg priced, a
parlay fair exists) or "decline" with a reason. The margin formula is #284;
everything it needs is logged here so margins can be replayed: per-leg fair,
its probBounds bracket, sigma (with c and tau), Kalshi's own price for the leg,
and hours to the earliest kickoff.

Parlay fair = product of leg fairs (legs are on different games, treated as
independent). fair_lo / fair_hi multiply the per-leg brackets.
"""
import math
import time
from datetime import datetime, timezone

from kalshi_parlay_quoter import scope
from kalshi_parlay_quoter.pricer import PricerBacklog, PricerError, PriceResult

LEG_FIELDS = ("ok", "reason", "detail", "priced_at_ms", "cache_age_ms", "fair", "fair_lo", "fair_hi", "bracketed",
              "sigma", "c", "tau", "line", "kickoff_ms", "fair_age_ms",
              "kalshi_status", "kalshi_yes_bid", "kalshi_yes_ask",
              "kalshi_no_ask")


def _num(v):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return f if math.isfinite(f) else None


def _prod(values):
    if not values or any(v is None for v in values):
        return None
    return math.prod(values)


def _parse_ts(value):
    if not value:
        return None
    try:
        return datetime.fromisoformat(str(value).replace("Z", "+00:00"))
    except ValueError:
        return None


def base_row(rfq: dict, seen_at: datetime) -> dict:
    return {
        "rfq_id": rfq.get("id"),
        "seen_at": seen_at,
        "created_ts": _parse_ts(rfq.get("created_ts")),
        "combo_ticker": rfq.get("market_ticker"),
        "creator_id": rfq.get("creator_id"),
        "contracts": _num(rfq.get("contracts_fp")),
        "target_cost_usd": _num(rfq.get("target_cost_dollars")),
        "n_legs": None, "leagues": None, "n_games": None, "same_game": None,
        "decision": "decline", "reason": None, "reason_detail": None,
        "parlay_fair": None, "parlay_fair_lo": None, "parlay_fair_hi": None,
        "max_leg_sigma": None, "legs_kalshi_ask_product": None,
        "earliest_kickoff_ms": None, "hours_to_kickoff": None,
        "pricer_ms": None, "priced_at_ms": None, "queue_ms": None,
        "confirmed": False,
    }


def decide(rfq: dict, legs: list[dict] | None, price_fn, *,
           seen_at: datetime | None = None,
           now_ms: int | None = None) -> tuple[dict, list[dict]]:
    """price_fn(rfq_id, [(ticker, side), ...]) -> PriceResult, may raise
    PricerError. Returns (decision row, leg rows)."""
    seen_at = seen_at or datetime.now(timezone.utc)
    now_ms = now_ms if now_ms is not None else int(time.time() * 1000)
    row = base_row(rfq, seen_at)
    verdict = scope.evaluate(legs)
    row.update(n_legs=len(verdict.legs), leagues=",".join(verdict.leagues) or None,
               n_games=verdict.n_games, same_game=verdict.same_game)
    leg_rows = [{
        "rfq_id": row["rfq_id"], "leg_idx": i, "market_ticker": li.market_ticker,
        "side": li.side, "series": li.series, "league": li.league,
        "market": li.market, "game_key": li.game_key,
        **{f: None for f in LEG_FIELDS},
    } for i, li in enumerate(verdict.legs)]

    if not verdict.in_scope:
        row.update(reason=verdict.reason, reason_detail=verdict.detail)
        return row, leg_rows

    try:
        res: PriceResult = price_fn(row["rfq_id"],
                                    [(l.market_ticker, l.side) for l in verdict.legs])
    except PricerBacklog as e:
        row.update(reason="pricer_backlog", reason_detail=str(e)[:500] or None)
        return row, leg_rows
    except PricerError as e:
        row.update(reason="pricer_error", reason_detail=str(e)[:500])
        return row, leg_rows

    row.update(pricer_ms=res.elapsed_ms, priced_at_ms=res.priced_at_ms)
    for lr, got in zip(leg_rows, res.legs):
        for f in LEG_FIELDS:
            lr[f] = got.get(f)
        lr["ok"] = bool(got.get("ok"))
        for f in ("fair", "fair_lo", "fair_hi", "sigma", "c", "tau",
                  "kalshi_yes_bid", "kalshi_yes_ask", "kalshi_no_ask",
                  "cache_age_ms"):
            lr[f] = _num(lr[f])
        kickoff = _num(lr["kickoff_ms"])
        if lr["ok"] and kickoff is not None and kickoff <= now_ms:
            # The runner's started gate ran when it priced the leg; a fair
            # held in the leg table must not outlive kickoff.
            lr.update(ok=False, reason="started", detail="started after it was priced")

    kickoffs = [k for k in (_num(lr["kickoff_ms"]) for lr in leg_rows) if k]
    if kickoffs:
        row["earliest_kickoff_ms"] = int(min(kickoffs))
        row["hours_to_kickoff"] = (min(kickoffs) - now_ms) / 3_600_000.0
    sigmas = [lr["sigma"] for lr in leg_rows if lr["sigma"] is not None]
    row["max_leg_sigma"] = max(sigmas) if sigmas else None
    asks = [lr["kalshi_yes_ask"] if lr["side"] == "yes" else lr["kalshi_no_ask"]
            for lr in leg_rows]
    row["legs_kalshi_ask_product"] = _prod(asks)

    failed = [lr for lr in leg_rows if not lr["ok"] or lr["fair"] is None]
    if failed:
        first = failed[0]
        row.update(reason=f"leg_{first['reason'] or 'no_fair'}",
                   reason_detail=";".join(
                       f"{lr['market_ticker']}:{lr['reason'] or 'no_fair'}"
                       for lr in failed)[:500])
        return row, leg_rows

    row["parlay_fair"] = _prod([lr["fair"] for lr in leg_rows])
    row["parlay_fair_lo"] = _prod([
        lr["fair_lo"] if lr["fair_lo"] is not None else lr["fair"] for lr in leg_rows])
    row["parlay_fair_hi"] = _prod([
        lr["fair_hi"] if lr["fair_hi"] is not None else lr["fair"] for lr in leg_rows])
    row["decision"] = "would_quote"
    return row, leg_rows
