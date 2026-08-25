"""One-sided Kalshi YRFI maker — daemon entry point.

Every CYCLE_SEC: refresh the open KXMLBRFI board (one Kalshi GET), poll our
fills, and per in-horizon game keep exactly one YES (YRFI) bid resting at
book-consensus fair minus MARGIN_CENTS, gated by dispersion/caps/deadlines.

Reads:  Kalshi REST (markets, portfolio), the 4 books' live SGP endpoints.
Writes: kalshi_rfi/kalshi_rfi.duckdb (snapshots/orders/fills/settlements —
        appends; fills upsert on trade_id). Live mode places/cancels real
        Kalshi limit orders. Shadow mode (default) writes the same rows,
        places nothing.

Run: python -m kalshi_rfi.main   (from repo root; see run.sh)
"""
import logging
import time
from datetime import datetime, timezone

from kalshi_common import auth_client

from kalshi_rfi import config, discovery, engine, state as state_mod, storage
from kalshi_rfi.fair import FairService
from kalshi_rfi.gateway import make_gateway
from kalshi_rfi.log_setup import setup_logging

log = logging.getLogger(__name__)


class _CachedFair:
    def __init__(self, result, mid_cents, now_mono):
        self.result = result          # FairResult | None (gate declined)
        self.mid_cents = mid_cents    # Kalshi mid when the fair was fetched
        self.fetched_mono = now_mono


def _mid_cents(game) -> float | None:
    if game.yes_bid_cents is None or game.yes_ask_cents is None:
        return None
    return (game.yes_bid_cents + game.yes_ask_cents) / 2.0


def _cancel(gateway, st, con, ticker: str, reason: str):
    resting = st.resting.get(ticker)
    if resting is None:
        return
    if gateway.cancel(resting["order_id"]):
        st.on_cancel(ticker)
        storage.log_order(con, ticker, resting["order_id"], "cancel",
                          resting["price_cents"], resting["count"],
                          "live" if gateway.is_live else "shadow", reason)
    else:
        # Keep local state: the order may still rest (retry next cycle) or
        # may have filled (poll_fills reconciles). Never place on top of it.
        log.warning("cancel failed for %s (%s) — holding state",
                    ticker, reason)


def _apply_decision(gateway, st, con, game, decision, now):
    mode = "live" if gateway.is_live else "shadow"
    resting = st.resting.get(game.ticker)
    if decision.action == "no_quote":
        _cancel(gateway, st, con, game.ticker, decision.reason)
        return
    if resting is not None and resting["price_cents"] == decision.price_cents:
        return                                       # hold — hysteresis
    if resting is not None:
        _cancel(gateway, st, con, game.ticker, "reprice")
        if game.ticker in st.resting:
            return                                   # cancel failed; retry next cycle
    coid = st.next_client_order_id(now)
    order_id = gateway.place_yes_bid(game.ticker, decision.price_cents,
                                     decision.count, coid)
    if order_id is None:
        storage.log_order(con, game.ticker, None, "place",
                          decision.price_cents, decision.count, mode,
                          "place_failed")
        return
    st.on_place(game.ticker, order_id, decision.price_cents, decision.count)
    storage.log_order(con, game.ticker, order_id, "place",
                      decision.price_cents, decision.count, mode,
                      decision.reason)
    log.info("QUOTE %s yes bid %dc x%d (%s)", game.ticker,
             decision.price_cents, decision.count, mode)


def run():
    setup_logging()
    if not config.KALSHI_API_KEY_ID or not config.KALSHI_PRIVATE_KEY_PATH:
        raise SystemExit("KALSHI_API_KEY_ID / KALSHI_PRIVATE_KEY_PATH not set")
    auth_client.configure(config.KALSHI_API_KEY_ID,
                          config.KALSHI_PRIVATE_KEY_PATH,
                          config.KALSHI_BASE_URL, config.PROJECT_ROOT)
    gateway = make_gateway(config.RFI_MODE, config.RFI_LIVE_ACK)
    if gateway is None:
        raise SystemExit("RFI_MODE=off — nothing to do")
    log.info("starting YRFI maker: mode=%s margin=%dc per_game=$%.0f "
             "daily=$%.0f books=%s", config.RFI_MODE, config.MARGIN_CENTS,
             config.PER_GAME_CAP_USD, config.DAILY_CAP_USD, config.BOOKS)

    con = storage.connect()
    st = state_mod.RfiState()
    fair_service = FairService()
    fair_cache: dict[str, _CachedFair] = {}
    if gateway.is_live:
        state_mod.sweep_orphan_orders(st, gateway)
    last_settle_poll = 0.0

    try:
        while True:
            if config.KILL_FILE.exists():
                log.warning("kill file present — exiting")
                break
            cycle_start = time.monotonic()
            now = datetime.now(timezone.utc)
            now_naive = now.replace(tzinfo=None)

            games = discovery.fetch_open_rfi_games()
            by_ticker = {g.ticker: g for g in games}
            for ticker in list(st.resting):
                if ticker not in by_ticker:
                    _cancel(gateway, st, con, ticker, "market_gone")

            if gateway.is_live:
                state_mod.poll_fills(st, con, now)
                if time.monotonic() - last_settle_poll >= config.SETTLEMENT_POLL_SEC:
                    state_mod.poll_settlements(st, con)
                    last_settle_poll = time.monotonic()

            for game in games:
                sec_to_start = (game.commence_utc - now_naive).total_seconds()
                if (sec_to_start > config.QUOTE_HORIZON_HOURS * 3600
                        or sec_to_start < -3600):
                    continue                # too far out / long since started
                if game.status not in ("active", "open"):
                    _cancel(gateway, st, con, game.ticker,
                            f"status_{game.status}")
                    continue

                cached = fair_cache.get(game.ticker)
                mid = _mid_cents(game)

                # Kalshi-side jump guard between book refreshes: the market
                # moving under our fair means the fair is stale (single-market
                # version of the MM's constituent_jump).
                if (cached is not None and mid is not None
                        and cached.mid_cents is not None
                        and abs(mid - cached.mid_cents) >= config.KALSHI_JUMP_CENTS):
                    log.info("kalshi_jump %s: mid %.1f -> %.1f — refetching",
                             game.ticker, cached.mid_cents, mid)
                    _cancel(gateway, st, con, game.ticker, "kalshi_jump")
                    cached = None

                past_pull = sec_to_start <= config.PULL_BEFORE_START_SEC
                ttl = (config.FAIR_REFRESH_NEAR_SEC
                       if sec_to_start <= config.NEAR_WINDOW_MIN * 60
                       else config.FAIR_REFRESH_FAR_SEC)
                if (not past_pull
                        and (cached is None
                             or time.monotonic() - cached.fetched_mono > ttl)):
                    result = fair_service.fetch(game)
                    cached = _CachedFair(result, mid, time.monotonic())
                    fair_cache[game.ticker] = cached

                fair_result = cached.result if cached is not None else None
                fair_yes = (fair_result.consensus_yes
                            if fair_result is not None else None)
                decision = engine.decide(
                    fair_yes, game.yes_ask_cents,
                    margin_cents=config.MARGIN_CENTS,
                    min_edge_cents=config.MIN_EDGE_CENTS,
                    per_game_cap_usd=config.PER_GAME_CAP_USD,
                    game_committed_usd=st.game_filled_cost_usd(game.ticker),
                    daily_remaining_usd=(config.DAILY_CAP_USD
                                         - st.daily_committed_usd(
                                             now, exclude_ticker=game.ticker)),
                    sec_to_start=sec_to_start,
                    pull_before_start_sec=config.PULL_BEFORE_START_SEC)
                _apply_decision(gateway, st, con, game, decision, now)
                storage.log_snapshot(con, game, fair_result, decision,
                                     sec_to_start)

            elapsed = time.monotonic() - cycle_start
            time.sleep(max(1.0, config.CYCLE_SEC - elapsed))
    finally:
        # Never leave unmanaged GTC orders resting after the daemon stops.
        for ticker in list(st.resting):
            _cancel(gateway, st, con, ticker, "shutdown")
        fair_service.close()
        con.close()
        log.info("YRFI maker stopped")


if __name__ == "__main__":
    run()
