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
from datetime import datetime, timedelta, timezone

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
_VOID_HALT_ACTIVE = False  # N12: track prior void-rate halt state for edge-triggered notify


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
    # N8: treat unreconciled fills conservatively — use GREATEST(contracts,
    # MAX_RFQ_CONTRACTS) so a fill whose confirm-response shape is still
    # unverified can't under-count and let extra quotes through the cap gates
    # during the 30s reconcile window.
    start = datetime.now(timezone.utc).replace(hour=0, minute=0, second=0, microsecond=0)
    with db.connect(read_only=True) as con:
        rows = con.execute(
            "SELECT game_id, price, "
            "CASE WHEN reconciled THEN contracts "
            "     ELSE GREATEST(contracts, ?) END AS effective_contracts "
            "FROM fills WHERE filled_at >= ?",
            [config.MAX_RFQ_CONTRACTS, start]).fetchall()
    return [{"game_id": g, "price": p * c} for g, p, c in rows]


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


def _void_rate_halt_triggered() -> bool:
    """H4 (global): if voided/(voided+confirmed) > threshold over the recent
    window, halt new quotes. A high void rate means we're either getting picked
    off (book moves we cancel on) or chronically can't re-price — both are
    signals to step back instead of farming volume."""
    cutoff = datetime.now(timezone.utc) - timedelta(hours=config.VOID_RATE_WINDOW_HOURS)
    with db.connect(read_only=True) as con:
        row = con.execute(
            "SELECT "
            "  SUM(CASE WHEN decision LIKE 'voided_%' THEN 1 ELSE 0 END), "
            "  SUM(CASE WHEN decision='confirmed' THEN 1 ELSE 0 END) "
            "FROM quote_decisions WHERE observed_at >= ?",
            [cutoff]).fetchone()
    voided = int(row[0] or 0)
    confirmed = int(row[1] or 0)
    denom = voided + confirmed
    if denom == 0:
        return False
    return (voided / denom) > config.VOID_RATE_HALT_THRESHOLD


def _creator_halt_active(creator_id: str) -> bool:
    """H4 (per-creator): if this counterparty has already farmed us for >=N
    fills in the recent window, refuse to quote them again."""
    if not creator_id:
        return False
    cutoff = datetime.now(timezone.utc) - timedelta(hours=config.PER_CREATOR_WINDOW_HOURS)
    with db.connect(read_only=True) as con:
        count = con.execute(
            "SELECT COUNT(*) FROM fills f JOIN seen_rfqs s ON f.rfq_id = s.rfq_id "
            "WHERE s.creator_id = ? AND f.filled_at >= ?",
            [creator_id, cutoff]).fetchone()[0]
    return int(count or 0) >= config.PER_CREATOR_FILL_HALT


def _discovery_tick(source, gateway, dry_run):
    global _VOID_HALT_ACTIVE
    # FIX M2: kill-switch — stop quoting immediately if the kill file exists
    if config.KILL_FILE.exists():
        return
    # Book-freshness gate: if our books are stale/missing we cannot price anything.
    # Replaces the old samples-staleness gate (model was removed in v1 hardening).
    if _SGP_ODDS is None or _SGP_ODDS.empty:
        return
    # H4 + N12: global void-rate halt — step back if we're voiding too often.
    # N12: edge-triggered notification so operator gets a single ping on state
    # transitions (False→True = halt; True→False = resume) rather than silence.
    void_halt_now = _void_rate_halt_triggered()
    if void_halt_now and not _VOID_HALT_ACTIVE:
        # Transition False → True: halt started.
        cutoff = datetime.now(timezone.utc) - timedelta(hours=config.VOID_RATE_WINDOW_HOURS)
        with db.connect(read_only=True) as con:
            row = con.execute(
                "SELECT "
                "  SUM(CASE WHEN decision LIKE 'voided_%' THEN 1 ELSE 0 END), "
                "  SUM(CASE WHEN decision='confirmed' THEN 1 ELSE 0 END) "
                "FROM quote_decisions WHERE observed_at >= ?",
                [cutoff]).fetchone()
        voided = int(row[0] or 0)
        confirmed = int(row[1] or 0)
        denom = voided + confirmed
        rate = voided / denom if denom else 0.0
        notify.halt("void_rate",
                    detail=f"void_rate={rate:.0%} in last {config.VOID_RATE_WINDOW_HOURS}h")
    elif not void_halt_now and _VOID_HALT_ACTIVE:
        # Transition True → False: halt lifted.
        notify.resume("void_rate", detail="void_rate fell back below threshold")
    _VOID_HALT_ACTIVE = void_halt_now
    if void_halt_now:
        _log_decision("halted_high_void_rate")
        return
    rfqs = source.poll()
    with db.connect(read_only=True) as con:
        open_count = con.execute(
            "SELECT COUNT(*) FROM live_quotes WHERE status='open'").fetchone()[0]
    # FIX I3: load today's fills once before the loop for cap checks
    fills_today = _today_fills()
    now_utc = datetime.now(timezone.utc)
    for rfq in rfqs:
        if open_count >= config.MAX_OPEN_QUOTES:
            break
        rid = rfq.get("id")
        ticker = rfq.get("market_ticker")
        if not rid or not ticker:
            continue
        # Kalshi RFQ poll responses tag the creator as creator_user_id; fall
        # back to creator_id for forward-compatibility if the API name changes.
        creator_id = rfq.get("creator_user_id") or rfq.get("creator_id") or ""
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
                    con.execute(
                        "INSERT OR REPLACE INTO seen_rfqs "
                        "(rfq_id, market_ticker, in_scope, game_id, legs_json, "
                        "first_seen_at, last_decision, creator_id) "
                        "VALUES (?,?,?,?,?,?,?,?)",
                        [rid, ticker, False, None, None, now_utc,
                         "out_of_scope", creator_id])
            continue
        # H4 (per-creator): refuse to quote a counterparty that's already
        # farmed us this window. Check BEFORE size gate so we don't waste cycles.
        if _creator_halt_active(creator_id):
            _log_decision("skipped", rfq_id=rid, ticker=ticker, reason="skipped_creator_halt")
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
        # H8 + N7 + N8: per-combo exposure cap — block multi-RFQ concentration
        # on one combo. N7: also count outstanding open live_quotes' worst-case
        # exposure so a burst of RFQs on the same combo can't all be accepted
        # before any fills register. N8: treat unreconciled fills conservatively
        # (GREATEST(contracts, MAX_RFQ_CONTRACTS)) for the same reason.
        with db.connect(read_only=True) as con:
            combo_exp = con.execute(
                "SELECT COALESCE(SUM(price * "
                "CASE WHEN reconciled THEN contracts "
                "     ELSE GREATEST(contracts, ?) END), 0) "
                "FROM fills WHERE combo_market_ticker=?",
                [config.MAX_RFQ_CONTRACTS, ticker]).fetchone()[0]
            inflight_count = con.execute(
                "SELECT COUNT(*) FROM live_quotes "
                "WHERE combo_market_ticker=? AND status='open'",
                [ticker]).fetchone()[0]
        worst_inflight = float(inflight_count) * config.MAX_RFQ_CONTRACTS * 1.0
        this_quote_max = config.MAX_RFQ_CONTRACTS * 1.0
        if float(combo_exp or 0) + worst_inflight + this_quote_max > config.MAX_COMBO_EXPOSURE_USD:
            _log_decision("skipped", rfq_id=rid, ticker=ticker, game_id=game_id,
                          reason="per_combo_cap")
            continue
        # H9: per-combo cooldown — refuse new quotes for COMBO_COOLDOWN_SEC
        # after a recent fill on this combo (same-price re-pickoff defense).
        # Compare in SQL: DuckDB normalizes the AWARE bind param into the
        # stored naive-local frame for us, so this is tz-consistent end-to-end.
        with db.connect(read_only=True) as con:
            cd_row = con.execute(
                "SELECT 1 FROM combo_cooldown "
                "WHERE combo_market_ticker=? AND cooled_until > ?",
                [ticker, now_utc]).fetchone()
        if cd_row is not None:
            _log_decision("skipped", rfq_id=rid, ticker=ticker, game_id=game_id,
                          reason="in_cooldown")
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
                con.execute(
                    "INSERT OR REPLACE INTO seen_rfqs "
                    "(rfq_id, market_ticker, in_scope, game_id, legs_json, "
                    "first_seen_at, last_decision, creator_id) "
                    "VALUES (?,?,?,?,?,?,?,?)",
                    [rid, ticker, True, game_id, json.dumps(legs),
                     datetime.now(timezone.utc), "quoted", creator_id])
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
                    # N5: record the fill IMMEDIATELY with EXPECTED side/size
                    # (no positions retry in the hot path — that lived here
                    # before and could spend up to 15s/fill, blowing the 30s
                    # confirm window during accept bursts). reconciled=FALSE
                    # marks the row for the background sweep to verify against
                    # Kalshi positions and correct side/size if needed.
                    now_ts = datetime.now(timezone.utc)
                    with db.connect() as con:
                        con.execute(
                            "INSERT INTO fills (fill_id, quote_id, rfq_id, combo_market_ticker, "
                            "game_id, side_held, contracts, price, fee, model_fair_at_quote, "
                            "book_fair_at_quote, blended_fair_at_quote, fair_at_confirm, "
                            "realized_pnl, filled_at, reconciled) "
                            "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                            # FIX I2: carry model_fair_at_q / book_fair_at_q from the live_quotes row
                            [str(uuid.uuid4()), qid, rid, ticker, game_id, side_held,
                             contracts, price, fee, model_fair_at_q, book_fair_at_q,
                             prev_fair, cur_fair, None, now_ts, False])
                        con.execute(
                            "UPDATE live_quotes SET status='filled', closed_at=? WHERE quote_id=?",
                            [now_ts, qid])
                        # H9: arm the per-combo cooldown — block fresh quotes
                        # on the same combo for COMBO_COOLDOWN_SEC.
                        con.execute(
                            "INSERT OR REPLACE INTO combo_cooldown "
                            "(combo_market_ticker, cooled_until) VALUES (?, ?)",
                            [ticker, now_ts + timedelta(seconds=config.COMBO_COOLDOWN_SEC)])
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


def _reconcile_sweep_tick():
    """N5: background reconciliation sweep. Picks up fills written with
    reconciled=FALSE and verifies side/size against Kalshi `/portfolio/positions`.

    For each unreconciled fill on ticker T:
      - prior_total = signed sum of OTHER fills on T (excluding this one)
      - actual      = current Kalshi position on T (signed)
      - delta       = actual - prior_total  → THIS fill's true side/size
    If delta disagrees with what we recorded, UPDATE the fill row. Either way,
    set reconciled=TRUE so we don't re-check it.

    N10: phantom-fill detection — if delta==0 but expected was non-zero, the
    confirm succeeded but no position materialized (e.g. cancel race, internal
    Kalshi mismatch). Mark contracts=0 + reconciled=TRUE + log [phantom_fill].
    Caps are no longer inflated by confirms that didn't actually fill.

    N11: max-age fallback — if positions API is persistently down, a fill older
    than MAX_RECONCILE_AGE_SEC gets marked reconciled=TRUE with its recorded
    values + [reconcile_max_age_fallback] warning. Prevents indefinite
    cap-blocking on Kalshi outages.

    Lifted out of the confirm hot path (was up to 15s/fill via retries on
    eventual consistency). Sweep runs every RECONCILE_SWEEP_SEC.
    """
    now = datetime.now(timezone.utc)
    with db.connect(read_only=True) as con:
        rows = con.execute(
            "SELECT fill_id, combo_market_ticker, side_held, contracts, filled_at "
            "FROM fills WHERE reconciled = FALSE").fetchall()
    if not rows:
        return
    for fill_id, ticker, recorded_side, recorded_ct, filled_at in rows:
        # N11: compute fill age before the API call so we can apply the max-age
        # fallback if the positions API is persistently unreachable.
        # DuckDB TIMESTAMP (not TIMESTAMPTZ) returns naive datetimes stored in
        # local time. Use naive local `now` on both sides to stay consistent.
        if filled_at is not None:
            now_local = datetime.now()
            filled_at_naive = (filled_at.replace(tzinfo=None)
                               if filled_at.tzinfo is not None else filled_at)
            fill_age = max(0.0, (now_local - filled_at_naive).total_seconds())
        else:
            fill_age = 0.0
        actual = _get_position_contracts(ticker)
        if actual is None:
            if fill_age > config.MAX_RECONCILE_AGE_SEC:
                # N11: API persistently down — fall back to recorded values so
                # caps don't block new quotes indefinitely during an outage.
                with db.connect() as con:
                    con.execute("UPDATE fills SET reconciled=TRUE WHERE fill_id=?",
                                [fill_id])
                print(f"[reconcile_max_age_fallback] fill_id={fill_id} ticker={ticker} "
                      f"age={fill_age:.0f}s — positions API unreachable, marking "
                      f"reconciled with recorded values", flush=True)
            else:
                # API down — leave reconciled=FALSE; next sweep will retry.
                print(f"[position_reconcile_unavailable] ticker={ticker} "
                      f"fill_id={fill_id} — will retry next sweep", flush=True)
            continue
        # prior_total = signed sum of all OTHER fills on this ticker.
        with db.connect(read_only=True) as con:
            prior_total = con.execute(
                "SELECT COALESCE(SUM(CASE WHEN side_held='yes' THEN contracts "
                "ELSE -contracts END), 0) FROM fills "
                "WHERE combo_market_ticker=? AND fill_id<>?",
                [ticker, fill_id]).fetchone()[0]
        delta = int(actual) - int(prior_total)
        recorded_signed = (int(recorded_ct) if recorded_side == "yes"
                           else -int(recorded_ct))
        # N10: phantom-fill — confirm succeeded but no position materialized.
        if delta == 0 and recorded_signed != 0:
            with db.connect() as con:
                con.execute(
                    "UPDATE fills SET contracts=0, reconciled=TRUE WHERE fill_id=?",
                    [fill_id])
            print(f"[phantom_fill] fill_id={fill_id} ticker={ticker} "
                  f"expected={recorded_signed} actual_delta=0 (kalshi_total={actual} "
                  f"prior={prior_total}) — marking contracts=0", flush=True)
            continue
        if delta != recorded_signed:
            print(f"[position_mismatch] ticker={ticker} fill_id={fill_id} "
                  f"recorded={recorded_signed} actual_delta={delta} "
                  f"(kalshi_total={actual} prior={prior_total}) "
                  f"— trusting Kalshi", flush=True)
            if delta > 0:
                new_side, new_ct = "yes", delta
            elif delta < 0:
                new_side, new_ct = "no", abs(delta)
            else:
                # delta == 0 AND recorded_signed == 0: already correct, just mark done.
                new_side, new_ct = recorded_side, 0
            with db.connect() as con:
                con.execute(
                    "UPDATE fills SET side_held=?, contracts=?, reconciled=TRUE "
                    "WHERE fill_id=?",
                    [new_side, new_ct, fill_id])
        else:
            with db.connect() as con:
                con.execute("UPDATE fills SET reconciled=TRUE WHERE fill_id=?",
                            [fill_id])


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
    last = {"disc": 0.0, "conf": 0.0, "risk": 0.0, "reconcile": 0.0,
            "sgp": time.time()}
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
            if (not dry_run
                    and now - last["reconcile"] >= config.RECONCILE_SWEEP_SEC):
                try:
                    _reconcile_sweep_tick()
                except Exception as e:
                    print(f"  reconcile err: {e}", flush=True)
                last["reconcile"] = now
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
