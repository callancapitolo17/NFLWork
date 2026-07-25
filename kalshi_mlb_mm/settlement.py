"""Settlement sweep — populate fills.realized_pnl once markets settle (issue #12, audit M-1).

For every combo_market_ticker holding reconciled fills with realized_pnl IS
NULL, polls Kalshi `GET /markets/{ticker}`. When the market reports a terminal
yes/no result, in ONE transaction:
  - fills.realized_pnl = contracts * ((result == side_held ? 1-price : -price) - fee)
  - one `settlements` audit row per ticker, storing the RAW market payload —
    the verification record for Kalshi's (still unverified) settlement
    response shape (same pattern as the accept_observed research event).
Then emits a `settlement_recorded` research event per settled fill, and — on
the FIRST settled ticker ever — best-effort verifies the maker-fee assumption
(`kalshi_common.ev_calc.maker_fee_per_contract`, assumed 25% of taker) against
the /portfolio endpoints, logging LOUDLY on mismatch.

Inputs:  kalshi_mlb_mm.duckdb::fills (reconciled, unsettled rows); Kalshi API.
Outputs: UPDATEs fills.realized_pnl; INSERT OR REPLACE settlements;
         research events (separate research DB, buffered).
Safety:  never raises into the trading loop — per-ticker failures are logged
         and retried on the next sweep; only reconciled fills are settled so
         a later reconcile correction can never invalidate a written P&L.
"""
import json
import logging
from datetime import datetime, timezone

from kalshi_common import auth_client
from kalshi_mlb_mm import db, research

log = logging.getLogger("kalshi_mlb_mm.settlement")

# Kalshi market lifecycle vocabulary is UNVERIFIED for MVE combos (audit O-3).
# We key settlement on `result` in {yes,no} rather than on `status`, but a
# terminal-looking status with a non-yes/no result (voided market?) deserves a
# loud log so it isn't silently re-polled forever.
_TERMINAL_STATUSES = {"settled", "finalized", "determined"}

# The sweep shares the trading loop's single thread — a hung API call at the
# default 30s timeout, multiplied over several unsettled tickers, would starve
# the 2s confirm loop (whose accepts expire in ~30s). Same rationale as
# _get_position_contracts' short timeout in main.py.
API_TIMEOUT_SEC = 10


def realized_pnl_usd(side_held: str, contracts: float, price: float,
                     fee: float, result: str) -> float:
    """Signed settlement P&L in dollars for one fill.

    `price` and `fee` are per-contract dollars (as stored on the fills row);
    a binary pays $1 on a win, $0 on a loss:
      win:  contracts * (1 - price - fee)
      loss: contracts * (-price - fee)
    """
    won = (result == side_held)
    outcome_per_contract = (1.0 - price) if won else (-price)
    return contracts * (outcome_per_contract - fee)


def _unsettled_tickers() -> list[str]:
    """Tickers with at least one reconciled, not-yet-settled fill.

    reconciled=TRUE only: the reconcile sweep may still rewrite side/size of a
    fresh fill, and a P&L computed from pre-correction values would be wrong
    forever. Settlement is hours away; reconciliation is seconds — waiting
    costs nothing.
    """
    with db.connect(read_only=True) as con:
        rows = con.execute(
            "SELECT DISTINCT combo_market_ticker FROM fills "
            "WHERE realized_pnl IS NULL AND reconciled").fetchall()
    return [r[0] for r in rows]


def _any_settlements_recorded() -> bool:
    with db.connect(read_only=True) as con:
        n = con.execute("SELECT COUNT(*) FROM settlements").fetchone()[0]
    return int(n or 0) > 0


def _fetch_market(ticker: str) -> dict | None:
    """GET /markets/{ticker} → market dict, or None on any API failure."""
    status_code, body, _ = auth_client.api("GET", f"/markets/{ticker}",
                                           timeout=API_TIMEOUT_SEC)
    if status_code != 200 or not isinstance(body, dict):
        log.warning("[settlement_fetch_failed] ticker=%s status=%s — will "
                    "retry next sweep", ticker, status_code)
        return None
    market = body.get("market")
    if not isinstance(market, dict):
        log.warning("[settlement_bad_shape] ticker=%s body_keys=%s — will "
                    "retry next sweep", ticker, list(body.keys()))
        return None
    return market


def _parse_settled_time(market: dict):
    """Best-effort aware datetime from the market payload, else None.

    Field name is unverified (O-3) — the raw payload lands in `settlements`
    either way, so a miss here only costs the convenience column.
    """
    for key in ("settled_time", "settlement_time"):
        val = market.get(key)
        if not val:
            continue
        try:
            return datetime.fromisoformat(str(val).replace("Z", "+00:00"))
        except ValueError:
            continue
    return None


def _settle_ticker(ticker: str) -> bool:
    """Settle every eligible fill on one ticker. Returns True iff settled."""
    market = _fetch_market(ticker)
    if market is None:
        return False
    result = str(market.get("result") or "").strip().lower()
    if result not in ("yes", "no"):
        market_status = str(market.get("status") or "").strip().lower()
        if market_status in _TERMINAL_STATUSES:
            # Terminal but not yes/no (voided/scratched market?). P&L handling
            # for voids is out of scope until we see one — surface it loudly.
            log.warning("[settlement_unrecognized_result] ticker=%s status=%s "
                        "result=%r — fills left unsettled, inspect payload",
                        ticker, market_status, result)
        return False
    now_utc = datetime.now(timezone.utc)
    settled_at = _parse_settled_time(market)
    raw_payload = json.dumps(market, default=repr)
    settled_fills = []
    with db.connect() as con:
        con.execute("BEGIN TRANSACTION")
        rows = con.execute(
            "SELECT fill_id, side_held, contracts, price, fee FROM fills "
            "WHERE combo_market_ticker=? AND realized_pnl IS NULL AND reconciled",
            [ticker]).fetchall()
        for fill_id, side_held, contracts, price, fee in rows:
            pnl = realized_pnl_usd(side_held, contracts, price, fee, result)
            con.execute("UPDATE fills SET realized_pnl=? WHERE fill_id=?",
                        [pnl, fill_id])
            settled_fills.append(dict(fill_id=fill_id, side_held=side_held,
                                      contracts=contracts, price=price,
                                      fee=fee, realized_pnl=pnl))
        con.execute(
            "INSERT OR REPLACE INTO settlements (combo_market_ticker, result, "
            "settled_at, raw_payload, recorded_at) VALUES (?,?,?,?,?)",
            [ticker, result, settled_at, raw_payload, now_utc])
        con.execute("COMMIT")
    for f in settled_fills:
        research.emit("settlement_recorded", ticker=ticker,
                      payload=dict(result=result, **f))
        log.info("[settlement_recorded] ticker=%s fill_id=%s result=%s "
                 "side_held=%s contracts=%s realized_pnl=%.4f",
                 ticker, f["fill_id"], result, f["side_held"],
                 f["contracts"], f["realized_pnl"])
    return True


def _collect_fee_fields(obj, prefix: str = "", out: dict | None = None) -> dict:
    """Recursively harvest numeric values under any key containing 'fee'."""
    if out is None:
        out = {}
    if isinstance(obj, dict):
        for key, val in obj.items():
            path = f"{prefix}.{key}" if prefix else key
            if isinstance(val, (dict, list)):
                _collect_fee_fields(val, path, out)
            elif "fee" in key.lower():
                try:
                    out[path] = float(val)
                except (TypeError, ValueError):
                    pass
    elif isinstance(obj, list):
        for i, val in enumerate(obj):
            _collect_fee_fields(val, f"{prefix}[{i}]", out)
    return out


def _filter_entries_to_ticker(obj, ticker: str):
    """Copy of `obj` with list entries about OTHER tickers dropped.

    Any dict inside a list that carries a "ticker" key is kept only if it
    matches `ticker`. Dicts without a "ticker" key are kept (fail-open:
    better to over-collect than to silently verify against nothing).
    """
    if isinstance(obj, dict):
        return {k: _filter_entries_to_ticker(v, ticker) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_filter_entries_to_ticker(item, ticker) for item in obj
                if not (isinstance(item, dict) and "ticker" in item
                        and item.get("ticker") != ticker)]
    return obj


def _verify_fee_assumption(ticker: str) -> None:
    """One-time fee check (issue #12 item 3): does the fee Kalshi ACTUALLY
    charged match maker_fee_per_contract (assumed 25% of taker)? The
    assumption silently shapes every quoted price, so a mismatch is loud.

    Best-effort: the portfolio payload shapes are unverified, so we harvest
    any fee-named field, compare at dollar AND cent scales, capture the raw
    payloads in a `fee_verification` research event for manual inspection,
    and never raise.
    """
    try:
        with db.connect(read_only=True) as con:
            row = con.execute(
                "SELECT COALESCE(SUM(fee * contracts), 0), "
                "COALESCE(SUM(contracts), 0) FROM fills "
                "WHERE combo_market_ticker=?", [ticker]).fetchone()
        expected_fee_usd = float(row[0])
        total_contracts = float(row[1])
        payloads = {}
        for name, path in (
                ("portfolio_fills", f"/portfolio/fills?ticker={ticker}&limit=100"),
                ("portfolio_settlements", f"/portfolio/settlements?limit=100")):
            try:
                status_code, body, _ = auth_client.api("GET", path,
                                                       timeout=API_TIMEOUT_SEC)
                payloads[name] = dict(status=status_code, body=body)
            except Exception as exc:
                payloads[name] = dict(error=repr(exc))
        # Harvest fee fields only from entries about THIS ticker: the
        # /portfolio/settlements payload has no ticker filter and covers the
        # whole account, so other markets' fees would otherwise pollute the
        # sum and fire a false [FEE_MISMATCH]. The research event still gets
        # the unfiltered payloads.
        observed_fields = _collect_fee_fields(
            _filter_entries_to_ticker(payloads, ticker))
        verdict = "no_fee_field_found"
        if observed_fields:
            observed_total = sum(observed_fields.values())
            tolerance = max(0.01, 0.005 * total_contracts)
            # Unit of the fee field is unverified — accept dollars or cents.
            if any(abs(observed_total * scale - expected_fee_usd) <= tolerance
                   for scale in (1.0, 0.01)):
                verdict = "matched"
                log.info("[fee_verified] ticker=%s expected=$%.4f observed=%s "
                         "— maker_fee_per_contract assumption holds",
                         ticker, expected_fee_usd, observed_fields)
            else:
                verdict = "MISMATCH"
                log.warning("[FEE_MISMATCH] ticker=%s expected=$%.4f "
                            "observed=%s — maker_fee_per_contract (25%% of "
                            "taker) does NOT match what Kalshi charged; every "
                            "quoted price is built on this number. Fix "
                            "kalshi_common/ev_calc.py before trusting margins.",
                            ticker, expected_fee_usd, observed_fields)
        else:
            log.warning("[fee_unverified] ticker=%s: no fee field found in "
                        "portfolio payloads — inspect the fee_verification "
                        "research event; maker fee assumption remains "
                        "UNVERIFIED", ticker)
        research.emit("fee_verification", ticker=ticker,
                      payload=dict(expected_fee_usd=expected_fee_usd,
                                   total_contracts=total_contracts,
                                   observed_fee_fields=observed_fields,
                                   verdict=verdict, raw=payloads))
    except Exception as exc:  # verification is advisory — never break the sweep
        log.warning("[fee_verification_failed] ticker=%s: %s", ticker, exc)


def settlement_sweep_tick() -> None:
    """One sweep pass. Idempotent; API-failure tolerant; never raises."""
    try:
        tickers = _unsettled_tickers()
    except Exception as exc:
        log.warning("[settlement_sweep_error] listing unsettled fills: %s — "
                    "will retry next sweep", exc)
        return
    if not tickers:
        return
    try:
        had_settlements = _any_settlements_recorded()
    except Exception:
        had_settlements = True   # fail-safe: skip fee check rather than crash
    for ticker in tickers:
        try:
            if _settle_ticker(ticker) and not had_settlements:
                _verify_fee_assumption(ticker)
                had_settlements = True
        except Exception as exc:
            log.warning("[settlement_sweep_error] ticker=%s: %s — will retry "
                        "next sweep", ticker, exc)
