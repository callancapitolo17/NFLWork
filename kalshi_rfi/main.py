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
import math
import time
from datetime import datetime, timezone

from kalshi_common import auth_client

from kalshi_rfi import config, discovery, engine, state as state_mod, storage
from kalshi_rfi.fair import FairService
from kalshi_rfi.gateway import NETWORK_ERRORS, make_gateway
from kalshi_rfi.log_setup import setup_logging

log = logging.getLogger(__name__)


# Floor between the cycle's fill poll and a pre-placement re-poll.
FILL_REPOLL_MIN_INTERVAL_SEC = 2.0


class _CachedFair:
    def __init__(self, result, ref_bid, ref_ask, now_mono):
        self.result = result          # FairResult | None (gate declined)
        self.ref_bid = ref_bid        # Kalshi touch when the fair was fetched
        self.ref_ask = ref_ask
        self.fetched_mono = now_mono


def _cancel(gateway, st, con, ticker: str, reason: str) -> bool:
    """Cancel EVERY order resting on one game. Returns True only if the game
    is left flat — the caller must never place on top of a cancel that
    failed. Normally there is exactly one order; more than one means state
    drifted and both must go."""
    ok = True
    for order_id, resting in st.resting_for(ticker):
        if gateway.cancel(order_id, ticker=ticker,
                          exchange_index=resting.get("exchange_index")):
            st.on_cancel(order_id)
            storage.log_order(con, ticker, order_id, "cancel",
                              resting["price_cents"], resting["count"],
                              "live" if gateway.is_live else "shadow", reason)
        else:
            # Keep local state: the order may still rest (retry next cycle)
            # or may have filled (poll_fills reconciles). Never place on
            # top of it.
            log.warning("cancel failed for %s (%s) — holding state",
                        ticker, reason)
            ok = False
    return ok


def _fit_to_caps(st, ticker: str, price_cents: int, count: int,
                 now, per_game_cap_usd: float, daily_cap_usd: float) -> int:
    """Hard pre-placement backstop: the largest order that still leaves the
    game and the day inside their caps, given everything we currently know
    is filled or resting. 0 means place nothing.

    This is the last line before the wire and it re-reads state rather than
    trusting the count the engine sized at the top of the cycle — that count
    came from a fill poll that may be a full cycle old, which is how one
    game reached 45 contracts against a $5 cap (bug 1)."""
    headroom = min(per_game_cap_usd - st.game_exposure_usd(ticker),
                   daily_cap_usd - st.daily_committed_usd(now))
    if headroom <= 0 or price_cents < 1:
        return 0
    return max(0, min(count, math.floor(headroom * 100.0 / price_cents)))


def _apply_decision(gateway, st, con, game, decision, now, *,
                    per_game_cap_usd, daily_cap_usd, resize=None):
    """Make the game's resting order match `decision`.

    resize(price_cents) -> int re-sizes against fill state refreshed RIGHT
    NOW, and runs before EVERY placement — not just replacements. It used to
    run only when an order was already resting, so a completely filled order
    (which poll_fills drops from state) took the fresh-quote path and sized
    off the stale top-of-cycle poll (bug 1).
    """
    mode = "live" if gateway.is_live else "shadow"
    resting = st.resting_for(game.ticker)
    if decision.action == "no_quote":
        _cancel(gateway, st, con, game.ticker, decision.reason)
        return
    count = decision.count
    if len(resting) == 1 and resting[0][1]["price_cents"] == decision.price_cents:
        if resting[0][1]["count"] <= count:
            return                          # hold — hysteresis (up-drift ok)
        # Down-drift means fills elsewhere consumed cap headroom: an
        # oversized resting order IS a cap breach, so replace despite the
        # unchanged price.
        reason = "downsize"
    else:
        reason = "reprice"
    if resting:
        if not _cancel(gateway, st, con, game.ticker, reason):
            return                          # cancel failed; retry next cycle
    if resize is not None:
        count = min(count, resize(decision.price_cents))
    count = _fit_to_caps(st, game.ticker, decision.price_cents, count, now,
                         per_game_cap_usd, daily_cap_usd)
    if count < 1:
        storage.log_order(con, game.ticker, None, "place", None, None, mode,
                          "caps_exhausted_pre_place")
        return
    coid = st.next_client_order_id(now)
    order_id = gateway.place_yes_bid(game.ticker, decision.price_cents,
                                     count, coid, game.exchange_index)
    if order_id is None:
        storage.log_order(con, game.ticker, None, "place",
                          decision.price_cents, count, mode,
                          "place_failed")
        return
    st.on_place(game.ticker, order_id, decision.price_cents, count,
                game.exchange_index)
    storage.log_order(con, game.ticker, order_id, "place",
                      decision.price_cents, count, mode, decision.reason)
    log.info("QUOTE %s yes bid %dc x%d (%s)", game.ticker,
             decision.price_cents, count, mode)


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
    if gateway.is_live:
        st.hydrate_from_db(con)
        state_mod.reconcile_resting_orders(st, gateway, con)
    fair_service = FairService()
    fair_cache: dict[str, _CachedFair] = {}
    last_settle_poll = 0.0

    try:
        while True:
            if config.KILL_FILE.exists():
                log.warning("kill file present — exiting")
                break
            cycle_start = time.monotonic()
            try:
                _run_cycle(gateway, st, con, fair_service, fair_cache)
                if (gateway.is_live and time.monotonic() - last_settle_poll
                        >= config.SETTLEMENT_POLL_SEC):
                    state_mod.poll_settlements(st, con)
                    last_settle_poll = time.monotonic()
            except NETWORK_ERRORS as e:
                # Transient network failure anywhere in the cycle: hold all
                # state and retry — a crash here would strand live GTC
                # orders with nothing watching them (finding 2).
                log.error("cycle network failure (state held): %s", e)
            elapsed = time.monotonic() - cycle_start
            time.sleep(max(1.0, config.CYCLE_SEC - elapsed))
    finally:
        # Never leave unmanaged GTC orders resting after the daemon stops —
        # and never let one failing cancel strand the rest (finding 2).
        for ticker in st.resting_tickers():
            try:
                _cancel(gateway, st, con, ticker, "shutdown")
            except Exception as e:
                log.error("shutdown cancel failed for %s: %s", ticker, e)
        fair_service.close()
        con.close()
        log.info("YRFI maker stopped")


def _run_cycle(gateway, st, con, fair_service,
               fair_cache: dict[str, _CachedFair]):
    now = datetime.now(timezone.utc)
    now_naive = now.replace(tzinfo=None)

    games = discovery.fetch_open_rfi_games()
    if games is None:
        return                # API failure, not an empty board: hold state
    by_ticker = {g.ticker: g for g in games}
    for ticker in st.resting_tickers():
        if ticker not in by_ticker:
            _cancel(gateway, st, con, ticker, "market_gone")
    for ticker in list(fair_cache):
        if ticker not in by_ticker:
            fair_cache.pop(ticker, None)

    if gateway.is_live:
        state_mod.poll_fills(st, con, now)
        # Local order tracking drifts (a misrouted cancel, a lost place
        # response); every cycle Kalshi's resting listing is the truth.
        state_mod.reconcile_resting_orders(st, gateway, con)

    for game in games:
        sec_to_start = (game.commence_utc - now_naive).total_seconds()
        if (sec_to_start > config.QUOTE_HORIZON_HOURS * 3600
                or sec_to_start < -3600):
            continue                # too far out / long since started
        if game.status not in ("active", "open"):
            _cancel(gateway, st, con, game.ticker, f"status_{game.status}")
            continue

        cached = fair_cache.get(game.ticker)
        resting = st.resting_for(game.ticker)

        # Kalshi-side jump guard between book refreshes: the market moving
        # under our fair means the fair is stale. Per-side comparison that
        # ignores our own resting bid — see engine.kalshi_book_moved.
        if cached is not None:
            moved = engine.kalshi_book_moved(
                cached.ref_bid, cached.ref_ask,
                game.yes_bid_cents, game.yes_ask_cents,
                resting[0][1]["price_cents"] if len(resting) == 1 else None,
                config.KALSHI_JUMP_CENTS)
            if moved is not None:
                log.info("kalshi_jump %s: %s — refetching", game.ticker,
                         moved)
                _cancel(gateway, st, con, game.ticker, f"kalshi_{moved}")
                cached = None

        past_pull = sec_to_start <= config.PULL_BEFORE_START_SEC
        ttl = (config.FAIR_REFRESH_NEAR_SEC
               if sec_to_start <= config.NEAR_WINDOW_MIN * 60
               else config.FAIR_REFRESH_FAR_SEC)
        if (not past_pull
                and (cached is None
                     or time.monotonic() - cached.fetched_mono > ttl)):
            result = fair_service.fetch(game)
            cached = _CachedFair(result, game.yes_bid_cents,
                                 game.yes_ask_cents, time.monotonic())
            fair_cache[game.ticker] = cached

        fair_result = cached.result if cached is not None else None
        fair_yes = (fair_result.consensus_yes
                    if fair_result is not None else None)

        def _resize(price_cents: int, _t=game.ticker) -> int:
            # Pre-placement resize: refresh fills first so a fill that landed
            # after the cycle's poll shrinks the order (finding 4 / bug 1).
            # Throttled: this now runs before EVERY placement, so a full
            # slate would otherwise poll fills once per quoted game. A poll
            # seconds old still closes the bug (the stale one was a full
            # cycle old).
            if (gateway.is_live and time.monotonic() - st.last_fill_poll_mono
                    > FILL_REPOLL_MIN_INTERVAL_SEC):
                state_mod.poll_fills(st, con, datetime.now(timezone.utc))
            return engine.size_contracts(
                price_cents, config.PER_GAME_CAP_USD,
                st.game_filled_cost_usd(_t),
                config.DAILY_CAP_USD
                - st.daily_committed_usd(now, exclude_ticker=_t))

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
        if (decision.action == "quote"
                and st.in_fill_cooldown(game.ticker,
                                        config.POST_FILL_COOLDOWN_SEC)):
            # Just filled: stand down (and pull what's left resting) until
            # the fill has settled into the caps. Re-quoting straight into a
            # fill is what turned one game into 5 fills in 60s (bug 1).
            decision = engine.QuoteDecision("no_quote", "post_fill_cooldown")
        _apply_decision(gateway, st, con, game, decision, now,
                        per_game_cap_usd=config.PER_GAME_CAP_USD,
                        daily_cap_usd=config.DAILY_CAP_USD,
                        resize=_resize)
        storage.log_snapshot(con, game, fair_result, decision, sec_to_start)


if __name__ == "__main__":
    run()
