"""Kalshi MLB MM (maker) bot — REST-polling daemon.

Composes the prior modules into three timed loops:
  - discovery: poll open RFQs, scope/price/quote in-scope spread×total combos
  - confirm:   last-look gate on accepted quotes before confirming the fill
  - risk:      kill-switch, book-freshness auto-pull, tipoff-cancel, drift sweep

Fair value is now pure book-consensus median — the model was removed in the v1
hardening pass (`fairs.py` docstring explains why). Book fairs come from the
maker's own market DB (its own sgp_runner cadence) and are filtered through
`_consensus_filter` (median + ±BOOK_CONSENSUS_BAND outlier rejection) so a
single rogue book can't anchor the quote.
"""
import argparse
import json
import os
import signal
import statistics
import threading
import time
import uuid
from datetime import datetime, timezone

import duckdb

from kalshi_common import auth_client, sgp_runner
from kalshi_common.ev_calc import maker_fee_per_contract
from kalshi_common.fair_value import devig_book
from kalshi_common.leg_types import (
    _parse_event_suffix,
    _spread_line_from_legs,
    _total_line_from_legs,
    _MLB_CODE_TO_TEAM,
)
from kalshi_mlb_mm import config, db, notify, pricing, risk, scope, fairs
from kalshi_mlb_mm.rfq_source import RestRFQSource
from kalshi_mlb_mm.quote_gateway import RestQuoteGateway

_running = threading.Event()
_running.set()
_SGP_ODDS = None         # pd.DataFrame
_PREV_BOOK_FAIR = {}     # combo_market_ticker -> last blended book fair (circuit breaker)
_SCOPE_CACHE = {}        # market_ticker -> (in_scope, game_id, legs)


def _signal_handler(_s, _f):
    _running.clear()


def _configure_auth():
    auth_client.configure(config.KALSHI_API_KEY_ID, config.KALSHI_PRIVATE_KEY_PATH,
                          config.KALSHI_BASE_URL, config.PROJECT_ROOT)


def _refresh_sgp():
    global _SGP_ODDS
    if not config.MARKET_DB.exists():
        return
    try:
        con = duckdb.connect(str(config.MARKET_DB), read_only=True)
    except duckdb.IOException:
        return
    try:
        tables = {t[0] for t in con.execute("SHOW TABLES").fetchall()}
        if "mlb_sgp_odds" not in tables:
            return
        _SGP_ODDS = con.execute(
            "SELECT game_id, combo, period, bookmaker, sgp_decimal, fetch_time, "
            "spread_line, total_line FROM mlb_sgp_odds WHERE period='FG' "
            "AND fetch_time > NOW() - INTERVAL (CAST(? AS BIGINT)) SECOND",
            [config.MAX_BOOK_STALENESS_SEC]).fetchdf()
    finally:
        con.close()


def _consensus_filter(book_fairs: dict[str, float]) -> dict[str, float]:
    """Book-consensus-band gate (v1 correlation defense).

    Algorithm:
      1. If fewer than MIN_AGREEING_BOOKS supplied → return {} (no quote).
      2. Compute median of all supplied book devigged fairs.
      3. Keep only books within ±BOOK_CONSENSUS_BAND of that median.
      4. If fewer than MIN_AGREEING_BOOKS survive → return {} (no quote).
      5. Otherwise return the surviving (agreeing) books.

    Mirrors the MLB answer-key dashboard's consensus-band logic. The caller
    then medians the surviving books to get the fair (see `fairs.blended_fair`).
    """
    if len(book_fairs) < config.MIN_AGREEING_BOOKS:
        return {}
    med = statistics.median(book_fairs.values())
    agreeing = {b: f for b, f in book_fairs.items()
                if abs(f - med) <= config.BOOK_CONSENSUS_BAND}
    if len(agreeing) < config.MIN_AGREEING_BOOKS:
        return {}
    return agreeing


# book fairs per (game, spread_line, total_line) — mirrors taker _load_book_fairs,
# but REQUIRES full 4-side devig (no fallback) per accepted-risk #6, then runs
# through _consensus_filter (v1 correlation defense, replaces MIN_BOOK_COUNT_FOR_BLEND).
def _book_fairs(game_id, spread_line, total_line):
    if _SGP_ODDS is None or _SGP_ODDS.empty:
        return {}
    rows = _SGP_ODDS[(_SGP_ODDS.game_id == game_id)
                     & (_SGP_ODDS.spread_line.astype(float).round(2) == round(spread_line, 2))
                     & (_SGP_ODDS.total_line.astype(float).round(2) == round(total_line, 2))]
    if rows.empty:
        return {}
    out = {}
    for book in rows.bookmaker.unique():
        sub = rows[rows.bookmaker == book]
        if len(sub) < 4:        # require full 4-side coverage (no crude fallback)
            continue
        f = devig_book(sub, combo="Home Spread + Over", vig_fallback=0.0)
        if f is not None:
            out[book] = f
    return _consensus_filter(out)


def _commence_time(game_id):
    # read from mlb_target_lines (written by sgp_runner) for tipoff gating
    if not config.MARKET_DB.exists():
        return None
    try:
        con = duckdb.connect(str(config.MARKET_DB), read_only=True)
    except duckdb.IOException:
        return None
    try:
        row = con.execute(
            "SELECT commence_time FROM mlb_target_lines WHERE game_id=? LIMIT 1",
            [game_id]).fetchone()
        return row[0] if row else None
    except duckdb.CatalogException:
        return None
    finally:
        con.close()


def _today_fills():
    start = datetime.now(timezone.utc).replace(hour=0, minute=0, second=0, microsecond=0)
    with db.connect(read_only=True) as con:
        return [{"game_id": g, "price": p * c}
                for g, p, c in con.execute(
                    "SELECT game_id, price, contracts FROM fills WHERE filled_at >= ?",
                    [start]).fetchall()]


def _resolve_game_and_lines(market_ticker, legs):
    """Resolve game_id + (spread_line, total_line) from decoded legs.

    - Parse the spread leg's event-ticker suffix -> (away_code, home_code).
    - Compute spread_line / total_line directly from the legs.
    - Map team codes -> canonical names, then look up game_id in
      MARKET_DB::mlb_target_lines by (home_team, away_team).

    Returns (game_id | None, spread_line, total_line).
    """
    spread_leg = next(
        (l for l in legs if str(l.get("market_ticker", "")).startswith("KXMLBSPREAD-")), None)
    if not spread_leg:
        return None, None, None
    # Event ticker self-encodes away/home: ...-{YYMMMDDHHMM}{Away}{Home}.
    suffix = str(spread_leg.get("event_ticker", "")).rsplit("-", 1)[-1]
    away, home = _parse_event_suffix(suffix)
    if not away or not home:
        return None, None, None
    spread_line = _spread_line_from_legs(legs)
    total_line = _total_line_from_legs(legs)
    home_name = _MLB_CODE_TO_TEAM.get(home)
    away_name = _MLB_CODE_TO_TEAM.get(away)
    if not home_name or not away_name or not config.MARKET_DB.exists():
        return None, spread_line, total_line
    try:
        con = duckdb.connect(str(config.MARKET_DB), read_only=True)
    except duckdb.IOException:
        return None, spread_line, total_line
    try:
        row = con.execute(
            "SELECT game_id FROM mlb_target_lines WHERE home_team=? AND away_team=? LIMIT 1",
            [home_name, away_name]).fetchone()
    except duckdb.CatalogException:
        row = None
    finally:
        con.close()
    return (row[0] if row else None), spread_line, total_line


def _get_position_contracts(market_ticker: str, timeout: int = 5) -> int | None:
    """Authoritative current position via /portfolio/positions. Returns a SIGNED
    int (positive = long YES, negative = long NO). None on API failure or timeout.
    Short default timeout so a stuck positions call can't block the confirm loop.
    Mirrors kalshi_mlb_rfq.rfq_client.get_position_contracts."""
    try:
        status, body, _ = auth_client.api(
            "GET", f"/portfolio/positions?ticker={market_ticker}&limit=10",
            timeout=timeout)
    except Exception:
        return None
    if status != 200 or not isinstance(body, dict):
        return None
    for p in body.get("market_positions") or []:
        if p.get("ticker") == market_ticker:
            fp = p.get("position_fp")
            if fp is not None:
                return int(float(fp))
            return int(p.get("position", 0))
    return 0


def _log_decision(decision, *, rfq_id=None, quote_id=None, ticker=None, game_id=None,
                  reason=None, model=None, book=None, blended=None, yb=None, nb=None):
    with db.connect() as con:
        con.execute(
            "INSERT INTO quote_decisions (decision_id, rfq_id, quote_id, "
            "combo_market_ticker, game_id, decision, reason, model_fair, book_fair, "
            "blended_fair, yes_bid, no_bid, observed_at) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
            [str(uuid.uuid4()), rfq_id, quote_id, ticker, game_id, decision, reason,
             model, book, blended, yb, nb, datetime.now(timezone.utc)])


def _discovery_tick(source, gateway, dry_run):
    # FIX M2: kill-switch — stop quoting immediately if the kill file exists
    if config.KILL_FILE.exists():
        return
    # Book-freshness gate: if our books are stale/missing we cannot price anything.
    # Replaces the old samples-staleness gate (model was removed in v1 hardening).
    if _SGP_ODDS is None or _SGP_ODDS.empty:
        return
    rfqs = source.poll()
    with db.connect(read_only=True) as con:
        open_count = con.execute(
            "SELECT COUNT(*) FROM live_quotes WHERE status='open'").fetchone()[0]
    # FIX I3: load today's fills once before the loop for cap checks
    fills_today = _today_fills()
    for rfq in rfqs:
        if open_count >= config.MAX_OPEN_QUOTES:
            break
        rid = rfq.get("id")
        ticker = rfq.get("market_ticker")
        if not rid or not ticker:
            continue
        with db.connect(read_only=True) as con:
            seen = con.execute("SELECT in_scope FROM seen_rfqs WHERE rfq_id=?", [rid]).fetchone()
        # scope (cache the market lookup verdict)
        if ticker in _SCOPE_CACHE:
            in_scope, game_id, legs = _SCOPE_CACHE[ticker]
        else:
            market = source.get_market(ticker)
            legs = scope.decode_legs(market) if market else None
            in_scope = bool(legs and scope.is_spread_total_2leg(legs))
            game_id = None
            _SCOPE_CACHE[ticker] = (in_scope, game_id, legs)
        if not in_scope:
            if not seen:
                _log_decision("skipped", rfq_id=rid, ticker=ticker, reason="out_of_scope")
                with db.connect() as con:
                    con.execute("INSERT OR REPLACE INTO seen_rfqs VALUES (?,?,?,?,?,?,?)",
                                [rid, ticker, False, None, None,
                                 datetime.now(timezone.utc), "out_of_scope"])
            continue
        # size gate
        if not risk.size_ok(int(rfq.get("contracts", 0) or 0), config.MAX_RFQ_CONTRACTS):
            _log_decision("skipped", rfq_id=rid, ticker=ticker, reason="size_gate")
            continue
        game_id, spread_line, total_line = _resolve_game_and_lines(ticker, legs)
        if game_id is None:
            _log_decision("skipped", rfq_id=rid, ticker=ticker, reason="no_game")
            continue
        # FIX I3: exposure cap gates (wired — were defined but never called)
        if not risk.daily_cap_ok(fills_today, config.daily_exposure_cap_usd()):
            _log_decision("skipped", rfq_id=rid, ticker=ticker, game_id=game_id,
                          reason="daily_cap")
            continue
        if not risk.per_game_cap_ok(game_id, fills_today, config.BANKROLL,
                                    config.MAX_GAME_EXPOSURE_PCT):
            _log_decision("skipped", rfq_id=rid, ticker=ticker, game_id=game_id,
                          reason="per_game_cap")
            continue
        # tipoff gate
        ct = _commence_time(game_id)
        if not risk.tipoff_ok(ct, config.TIPOFF_CANCEL_MIN):
            _log_decision("skipped", rfq_id=rid, ticker=ticker, game_id=game_id, reason="tipoff")
            continue
        book_fairs = _book_fairs(game_id, spread_line, total_line)
        book_med, blended = fairs.blended_fair(legs, game_id, book_fairs)
        if blended is None or not (config.MIN_FAIR_PROB <= blended <= config.MAX_FAIR_PROB):
            _log_decision("skipped", rfq_id=rid, ticker=ticker, game_id=game_id,
                          reason="no_fair", model=None, book=book_med, blended=blended)
            continue
        # circuit breaker bookkeeping
        prev = _PREV_BOOK_FAIR.get(ticker)
        _PREV_BOOK_FAIR[ticker] = book_med
        if (prev is not None and book_med is not None
                and risk.book_move_triggered(prev, book_med, config.BOOK_MOVE_CB_THRESHOLD)):
            # H1: circuit breaker now actually CANCELS open quotes on this combo
            # (previously it only blocked re-quotes — resting quotes were exposed).
            with db.connect(read_only=True) as con:
                opens = con.execute(
                    "SELECT quote_id FROM live_quotes WHERE combo_market_ticker=? AND status='open'",
                    [ticker]).fetchall()
            for (open_qid,) in opens:
                try:
                    gateway.cancel(open_qid)
                except Exception:
                    pass
                with db.connect() as con:
                    con.execute(
                        "UPDATE live_quotes SET status='cancelled', closed_at=? WHERE quote_id=?",
                        [datetime.now(timezone.utc), open_qid])
            _log_decision("circuit_breaker", rfq_id=rid, ticker=ticker, game_id=game_id,
                          reason=f"book_move_pulled_{len(opens)}_quotes")
            continue
        q = pricing.quote(blended, config.TARGET_ROI)
        if q is None:
            _log_decision("skipped", rfq_id=rid, ticker=ticker, game_id=game_id,
                          reason="unpriceable")
            continue
        if dry_run:
            _log_decision("dry_run_quote", rfq_id=rid, ticker=ticker, game_id=game_id,
                          model=None, book=book_med, blended=blended, yb=q.yes_bid, nb=q.no_bid)
            continue
        # FIX I1 + M1: dedup + hysteresis — skip re-quote if price unchanged, else
        # mark the old row 'replaced' so we never leave ghost-open rows.
        with db.connect(read_only=True) as con:
            existing = con.execute(
                "SELECT quote_id, yes_bid, no_bid FROM live_quotes "
                "WHERE rfq_id=? AND status='open'",
                [rid]).fetchone()
        if existing:
            eqid, eyb, enb = existing
            if (abs(q.yes_bid - eyb) < config.QUOTE_HYSTERESIS
                    and abs(q.no_bid - enb) < config.QUOTE_HYSTERESIS):
                continue  # price essentially unchanged — leave the resting quote in place
            # price moved beyond hysteresis → replace: mark old row closed, then submit new
            with db.connect() as con:
                con.execute(
                    "UPDATE live_quotes SET status='replaced', closed_at=? WHERE quote_id=?",
                    [datetime.now(timezone.utc), eqid])
        qid = gateway.submit_quote(rid, q.yes_bid, q.no_bid)
        if qid:
            with db.connect() as con:
                con.execute("INSERT INTO live_quotes VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
                            [qid, rid, ticker, game_id, q.yes_bid, q.no_bid, None, book_med,
                             blended, "open", datetime.now(timezone.utc), None])
                con.execute("INSERT OR REPLACE INTO seen_rfqs VALUES (?,?,?,?,?,?,?)",
                            [rid, ticker, True, game_id, json.dumps(legs),
                             datetime.now(timezone.utc), "quoted"])
            _log_decision("quoted", rfq_id=rid, quote_id=qid, ticker=ticker, game_id=game_id,
                          model=None, book=book_med, blended=blended, yb=q.yes_bid, nb=q.no_bid)
            open_count += 1


def _confirm_tick(gateway, dry_run):
    with db.connect(read_only=True) as con:
        # FIX I2: also select model_fair and book_fair so we can carry them into fills
        live = con.execute(
            "SELECT quote_id, rfq_id, combo_market_ticker, game_id, yes_bid, no_bid, "
            "blended_fair, model_fair, book_fair "
            "FROM live_quotes WHERE status='open'").fetchall()
    for qid, rid, ticker, game_id, yb, nb, prev_fair, model_fair_at_q, book_fair_at_q in live:
        status, body, _ = auth_client.api("GET", f"/communications/quotes/{qid}")
        q = body.get("quote") if isinstance(body, dict) else None
        st = (q or {}).get("status")
        if st == "accepted":
            # Which side do we end up HOLDING? The taker accepted one side of our
            # two-sided quote; we hold the opposite. If they took our yes_bid (they
            # sold us YES), we hold YES; the API marks the accepted_side as theirs.
            accepted_side = (q or {}).get("accepted_side")
            side_held = "no" if accepted_side == "yes" else "yes"
            price = (1 - nb) if side_held == "yes" else (1 - yb)  # ask we transact at
            fee = maker_fee_per_contract(price)
            with db.connect(read_only=True) as con:
                legs_row = con.execute(
                    "SELECT legs_json FROM seen_rfqs WHERE rfq_id=?", [rid]).fetchone()
            legs = json.loads(legs_row[0]) if legs_row and legs_row[0] else []
            # H2: last-look HARD-FAIL when we can't re-price.
            # Principle: "can't re-price ⇒ don't confirm." Previously a missing
            # legs row or empty book_fairs silently fell back to prev_fair, which
            # made the drift check a no-op for exactly the cases that matter most.
            if not legs:
                _log_decision("voided_no_legs", rfq_id=rid, quote_id=qid, ticker=ticker,
                              game_id=game_id)
                with db.connect() as con:
                    con.execute(
                        "UPDATE live_quotes SET status='voided', closed_at=? WHERE quote_id=?",
                        [datetime.now(timezone.utc), qid])
                continue
            spread_line = _spread_line_from_legs(legs)
            total_line = _total_line_from_legs(legs)
            book_fairs_now = _book_fairs(game_id, spread_line, total_line)
            if not book_fairs_now:
                _log_decision("voided_no_fresh_books", rfq_id=rid, quote_id=qid, ticker=ticker,
                              game_id=game_id)
                with db.connect() as con:
                    con.execute(
                        "UPDATE live_quotes SET status='voided', closed_at=? WHERE quote_id=?",
                        [datetime.now(timezone.utc), qid])
                continue
            _, cur_fair = fairs.blended_fair(legs, game_id, book_fairs_now)
            if cur_fair is None:
                _log_decision("voided_blend_failed", rfq_id=rid, quote_id=qid, ticker=ticker,
                              game_id=game_id)
                with db.connect() as con:
                    con.execute(
                        "UPDATE live_quotes SET status='voided', closed_at=? WHERE quote_id=?",
                        [datetime.now(timezone.utc), qid])
                continue
            if risk.last_look_ok(side_held, price, fee, cur_fair, prev_fair,
                                 config.FAIR_DRIFT_TOLERANCE):
                if not dry_run and gateway.confirm(qid):
                    contracts = int(float((q or {}).get("contracts", 1) or 1))
                    # H6 (refined for N1/N2): reconcile against Kalshi positions as
                    # ground truth, but compute THIS fill's delta from our prior
                    # recorded total — Kalshi reports AGGREGATE per ticker, not
                    # per-fill. Retry briefly to wait out eventual consistency
                    # (positions read right after confirm can return the pre-fill total).
                    with db.connect(read_only=True) as con:
                        prior_total = con.execute(
                            "SELECT COALESCE(SUM(CASE WHEN side_held='yes' THEN contracts "
                            "ELSE -contracts END), 0) FROM fills WHERE combo_market_ticker=?",
                            [ticker]).fetchone()[0]
                    expected_signed = contracts if side_held == "yes" else -contracts

                    actual = None
                    for attempt in range(3):
                        actual = _get_position_contracts(ticker)
                        if actual is not None and actual != prior_total:
                            break  # Kalshi has caught up
                        time.sleep(0.2 * (2 ** attempt))  # 200ms, 400ms, 800ms

                    if actual is None:
                        # API down or persistently lagging → trust expected as fallback
                        print(f"[position_reconcile_unavailable] ticker={ticker} "
                              f"using expected={expected_signed}", flush=True)
                        delta = expected_signed
                    else:
                        delta = actual - int(prior_total)  # Kalshi total - prior = this fill
                        if delta != expected_signed:
                            print(f"[position_mismatch] ticker={ticker} "
                                  f"expected={expected_signed} actual_delta={delta} "
                                  f"(kalshi_total={actual} prior={prior_total}) "
                                  f"— trusting Kalshi", flush=True)

                    # Apply the trusted delta to this fill's recorded side/size
                    if delta > 0:
                        side_held, contracts = "yes", delta
                    elif delta < 0:
                        side_held, contracts = "no", abs(delta)
                    # delta == 0: only reachable when API is down AND expected was 0;
                    # keep side_held/contracts as-is.
                    with db.connect() as con:
                        con.execute(
                            "INSERT INTO fills (fill_id, quote_id, rfq_id, combo_market_ticker, "
                            "game_id, side_held, contracts, price, fee, model_fair_at_quote, "
                            "book_fair_at_quote, blended_fair_at_quote, fair_at_confirm, "
                            "realized_pnl, filled_at) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                            # FIX I2: carry model_fair_at_q / book_fair_at_q from the live_quotes row
                            [str(uuid.uuid4()), qid, rid, ticker, game_id, side_held,
                             contracts, price, fee, model_fair_at_q, book_fair_at_q,
                             prev_fair, cur_fair, None,
                             datetime.now(timezone.utc)])
                        con.execute(
                            "UPDATE live_quotes SET status='filled', closed_at=? WHERE quote_id=?",
                            [datetime.now(timezone.utc), qid])
                    notify.fill(ticker, side_held, contracts, price)
                    _log_decision("confirmed", rfq_id=rid, quote_id=qid, ticker=ticker,
                                  game_id=game_id)
            else:
                _log_decision("voided_last_look", rfq_id=rid, quote_id=qid, ticker=ticker,
                              game_id=game_id)
                with db.connect() as con:
                    con.execute(
                        "UPDATE live_quotes SET status='voided', closed_at=? WHERE quote_id=?",
                        [datetime.now(timezone.utc), qid])
        elif st in ("cancelled", "expired", "executed"):
            with db.connect() as con:
                con.execute("UPDATE live_quotes SET status=?, closed_at=? WHERE quote_id=?",
                            [st, datetime.now(timezone.utc), qid])


def _risk_sweep_tick(gateway):
    if config.KILL_FILE.exists():
        notify.halt("kill_switch")
        return
    # Freshness / auto-pull: if our book odds are stale/missing, pull EVERY open
    # quote — we can no longer trust our fair value, so we must stop resting risk.
    # Replaces the old samples-staleness gate (model was removed in v1 hardening).
    books_stale = _SGP_ODDS is None or _SGP_ODDS.empty
    with db.connect(read_only=True) as con:
        # Pull the per-quote fields we need for the drift-since-quote check too.
        live = con.execute(
            "SELECT lq.quote_id, lq.game_id, lq.combo_market_ticker, lq.book_fair, "
            "       lq.rfq_id, sr.legs_json "
            "FROM live_quotes lq LEFT JOIN seen_rfqs sr ON lq.rfq_id = sr.rfq_id "
            "WHERE lq.status='open'").fetchall()
    for qid, game_id, ticker, book_fair_at_q, rid, legs_json in live:
        ct = _commence_time(game_id)
        cancel = False
        if books_stale or not risk.tipoff_ok(ct, config.TIPOFF_CANCEL_MIN):
            cancel = True
        else:
            # H1 (part 2): per-open-quote drift-since-quote sweep.
            # Recompute current book consensus for this combo; if it has drifted
            # more than BOOK_MOVE_CB_THRESHOLD from book_fair-at-quote, cancel.
            # Catches gradual drift that the per-tick (last vs current) circuit
            # breaker misses (e.g., 1¢ moves over several ticks adding to 4¢).
            if book_fair_at_q is not None and legs_json:
                try:
                    legs = json.loads(legs_json)
                    sl = _spread_line_from_legs(legs)
                    tl = _total_line_from_legs(legs)
                    bf_now = _book_fairs(game_id, sl, tl)
                    if bf_now:
                        cur_med = statistics.median(bf_now.values())
                        if risk.book_move_triggered(book_fair_at_q, cur_med,
                                                    config.BOOK_MOVE_CB_THRESHOLD):
                            cancel = True
                except Exception:
                    pass
        if cancel:
            try:
                gateway.cancel(qid)
            except Exception:
                pass
            with db.connect() as con:
                con.execute(
                    "UPDATE live_quotes SET status='cancelled', closed_at=? WHERE quote_id=?",
                    [datetime.now(timezone.utc), qid])


def main_loop(dry_run: bool):
    _configure_auth()
    db.init_database()
    sid = db.start_session(pid=os.getpid(), dry_run=dry_run)
    print(f"=== MM bot session {sid} dry_run={dry_run} ===", flush=True)
    source, gateway = RestRFQSource(), RestQuoteGateway()
    # synchronous warm-up: one SGP cycle
    try:
        sgp_runner.sgp_cycle(bot_market_db=str(config.MARKET_DB),
                             scraper_dir=str(config.MLB_SGP_DIR),
                             venv_python=str(config.MLB_SGP_DIR / "venv" / "bin" / "python"),
                             timeout_sec=config.SGP_SCRAPER_TIMEOUT_SEC)
    except Exception as e:
        print(f"  warmup sgp failed: {e}", flush=True)
    _refresh_sgp()
    last = {"disc": 0.0, "conf": 0.0, "risk": 0.0, "sgp": time.time()}
    try:
        while _running.is_set():
            now = time.time()
            if now - last["disc"] >= config.DISCOVERY_SEC:
                try:
                    _discovery_tick(source, gateway, dry_run)
                except Exception as e:
                    print(f"  disc err: {e}", flush=True)
                last["disc"] = now
            if now - last["conf"] >= config.CONFIRM_SEC:
                try:
                    _confirm_tick(gateway, dry_run)
                except Exception as e:
                    print(f"  conf err: {e}", flush=True)
                last["conf"] = now
            if now - last["risk"] >= config.RISK_SWEEP_SEC:
                try:
                    _risk_sweep_tick(gateway)
                except Exception as e:
                    print(f"  risk err: {e}", flush=True)
                last["risk"] = now
            if now - last["sgp"] >= config.SGP_REFRESH_SEC:
                try:
                    sgp_runner.sgp_cycle(
                        bot_market_db=str(config.MARKET_DB),
                        scraper_dir=str(config.MLB_SGP_DIR),
                        venv_python=str(config.MLB_SGP_DIR / "venv" / "bin" / "python"),
                        timeout_sec=config.SGP_SCRAPER_TIMEOUT_SEC)
                    _refresh_sgp()
                except Exception as e:
                    print(f"  sgp err: {e}", flush=True)
                last["sgp"] = now
            time.sleep(0.25)   # short sleep → responsive SIGTERM
    finally:
        with db.connect(read_only=True) as con:
            live = [r[0] for r in con.execute(
                "SELECT quote_id FROM live_quotes WHERE status='open'").fetchall()]
        for qid in live:
            try:
                gateway.cancel(qid)
                with db.connect() as con:
                    con.execute(
                        "UPDATE live_quotes SET status='cancelled', closed_at=? WHERE quote_id=?",
                        [datetime.now(timezone.utc), qid])
            except Exception:
                pass
        db.end_session(sid)
        print("=== shutdown complete ===", flush=True)


def cli():
    p = argparse.ArgumentParser()
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args()
    signal.signal(signal.SIGINT, _signal_handler)
    signal.signal(signal.SIGTERM, _signal_handler)
    main_loop(dry_run=args.dry_run)


if __name__ == "__main__":
    cli()
