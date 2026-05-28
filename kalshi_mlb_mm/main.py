"""Kalshi MLB MM (maker) bot — REST-polling daemon.

Composes the prior modules into three timed loops:
  - discovery: poll open RFQs, scope/price/quote in-scope spread×total combos
  - confirm:   last-look gate on accepted quotes before confirming the fill
  - risk:      kill-switch, samples-staleness auto-pull, tipoff-cancel

Fair value reuses the taker's exact math via kalshi_common: model samples come
from mlb_mm.duckdb::mlb_game_samples (read-only), book fairs from the maker's
own market DB (its own sgp_runner cadence), blended via fair_value.blend.
"""
import argparse
import json
import os
import signal
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
_SAMPLES = {}            # game_id -> df
_SAMPLES_GEN_AT = None   # datetime | None — mlb_samples_meta.generated_at (freshness gate)
_SGP_ODDS = None         # pd.DataFrame
_PREV_BOOK_FAIR = {}     # combo_market_ticker -> last blended book fair (circuit breaker)
_SCOPE_CACHE = {}        # market_ticker -> (in_scope, game_id, legs)


def _signal_handler(_s, _f):
    _running.clear()


def _configure_auth():
    auth_client.configure(config.KALSHI_API_KEY_ID, config.KALSHI_PRIVATE_KEY_PATH,
                          config.KALSHI_BASE_URL, config.PROJECT_ROOT)


def _refresh_samples():
    """Load model samples + their generation timestamp from the answer-key DB.

    Reads mlb_game_samples (keyed into _SAMPLES by game_id) and the latest
    mlb_samples_meta.generated_at (into _SAMPLES_GEN_AT, used by the freshness
    gate). Missing DB / tables leave the prior caches untouched; _SAMPLES_GEN_AT
    is set to None only when the meta row is genuinely unavailable.
    """
    global _SAMPLES, _SAMPLES_GEN_AT
    if not config.ANSWER_KEY_DB.exists():
        return
    try:
        con = duckdb.connect(str(config.ANSWER_KEY_DB), read_only=True)
    except duckdb.IOException:
        return
    try:
        sdf = con.execute(
            "SELECT game_id, sim_idx, home_margin, total_final_score, "
            "home_margin_f5, total_f5 FROM mlb_game_samples").fetchdf()
        _SAMPLES = {g: d.reset_index(drop=True) for g, d in sdf.groupby("game_id")}
        try:
            row = con.execute(
                "SELECT generated_at FROM mlb_samples_meta "
                "ORDER BY generated_at DESC LIMIT 1").fetchone()
            _SAMPLES_GEN_AT = row[0] if row else None
        except duckdb.CatalogException:
            _SAMPLES_GEN_AT = None
    except duckdb.CatalogException:
        # mlb_game_samples mid-replace or missing — keep prior caches.
        return
    finally:
        con.close()


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


# book fairs per (game, spread_line, total_line) — mirrors taker _load_book_fairs,
# but REQUIRES full 4-side devig (no fallback) per accepted-risk #6.
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
    return out if len(out) >= config.MIN_BOOK_COUNT_FOR_BLEND else {}


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
    # FIX I4: samples-staleness gate — don't quote anything when model is stale
    if not risk.staleness_ok(_SAMPLES_GEN_AT, config.MAX_PREDICTION_STALENESS_SEC):
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
        samples = _SAMPLES.get(game_id)
        # tipoff gate
        ct = _commence_time(game_id)
        if not risk.tipoff_ok(ct, config.TIPOFF_CANCEL_MIN):
            _log_decision("skipped", rfq_id=rid, ticker=ticker, game_id=game_id, reason="tipoff")
            continue
        book_fairs = _book_fairs(game_id, spread_line, total_line)
        model, book_med, blended = fairs.blended_fair(legs, game_id, samples, book_fairs)
        if blended is None or not (config.MIN_FAIR_PROB <= blended <= config.MAX_FAIR_PROB):
            _log_decision("skipped", rfq_id=rid, ticker=ticker, game_id=game_id,
                          reason="no_fair", model=model, book=book_med, blended=blended)
            continue
        # circuit breaker bookkeeping
        prev = _PREV_BOOK_FAIR.get(ticker)
        _PREV_BOOK_FAIR[ticker] = book_med
        if (prev is not None and book_med is not None
                and risk.book_move_triggered(prev, book_med, config.BOOK_MOVE_CB_THRESHOLD)):
            _log_decision("skipped", rfq_id=rid, ticker=ticker, game_id=game_id,
                          reason="circuit_breaker")
            continue
        q = pricing.quote(blended, config.TARGET_ROI)
        if q is None:
            _log_decision("skipped", rfq_id=rid, ticker=ticker, game_id=game_id,
                          reason="unpriceable")
            continue
        if dry_run:
            _log_decision("dry_run_quote", rfq_id=rid, ticker=ticker, game_id=game_id,
                          model=model, book=book_med, blended=blended, yb=q.yes_bid, nb=q.no_bid)
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
                            [qid, rid, ticker, game_id, q.yes_bid, q.no_bid, model, book_med,
                             blended, "open", datetime.now(timezone.utc), None])
                con.execute("INSERT OR REPLACE INTO seen_rfqs VALUES (?,?,?,?,?,?,?)",
                            [rid, ticker, True, game_id, json.dumps(legs),
                             datetime.now(timezone.utc), "quoted"])
            _log_decision("quoted", rfq_id=rid, quote_id=qid, ticker=ticker, game_id=game_id,
                          model=model, book=book_med, blended=blended, yb=q.yes_bid, nb=q.no_bid)
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
            samples = _SAMPLES.get(game_id)
            with db.connect(read_only=True) as con:
                legs_row = con.execute(
                    "SELECT legs_json FROM seen_rfqs WHERE rfq_id=?", [rid]).fetchone()
            legs = json.loads(legs_row[0]) if legs_row and legs_row[0] else []
            # FIX M3: recompute book_fairs at confirm time so cur_fair is real
            # (previously called with {} → blend returned None → drift was always 0)
            if legs:
                spread_line = _spread_line_from_legs(legs)
                total_line = _total_line_from_legs(legs)
                book_fairs_now = _book_fairs(game_id, spread_line, total_line)
                _, _, cur_fair = fairs.blended_fair(legs, game_id, samples, book_fairs_now)
            else:
                cur_fair = prev_fair
            cur_fair = cur_fair if cur_fair is not None else prev_fair
            if risk.last_look_ok(side_held, price, fee, cur_fair, prev_fair,
                                 config.FAIR_DRIFT_TOLERANCE):
                if not dry_run and gateway.confirm(qid):
                    contracts = int(float((q or {}).get("contracts", 1) or 1))
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
    # Freshness / auto-pull: if model samples have gone stale, pull EVERY open
    # quote — we can no longer trust our fair value, so we must stop resting risk.
    samples_stale = not risk.staleness_ok(_SAMPLES_GEN_AT, config.MAX_PREDICTION_STALENESS_SEC)
    with db.connect(read_only=True) as con:
        live = con.execute(
            "SELECT quote_id, game_id FROM live_quotes WHERE status='open'").fetchall()
    for qid, game_id in live:
        ct = _commence_time(game_id)
        if samples_stale or not risk.tipoff_ok(ct, config.TIPOFF_CANCEL_MIN):
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
    # synchronous warm-up: one SGP cycle + sample load
    try:
        sgp_runner.sgp_cycle(bot_market_db=str(config.MARKET_DB),
                             scraper_dir=str(config.MLB_SGP_DIR),
                             venv_python=str(config.MLB_SGP_DIR / "venv" / "bin" / "python"),
                             timeout_sec=config.SGP_SCRAPER_TIMEOUT_SEC)
    except Exception as e:
        print(f"  warmup sgp failed: {e}", flush=True)
    _refresh_sgp()
    _refresh_samples()
    last = {"disc": 0.0, "conf": 0.0, "risk": 0.0, "sgp": time.time(), "samp": time.time()}
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
            if now - last["samp"] >= config.SAMPLES_REFRESH_SEC:
                _refresh_samples()
                last["samp"] = now
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
