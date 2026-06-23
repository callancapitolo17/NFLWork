"""Read-only, lock-safe query helpers for the Kalshi MLB monitor.

Every function opens a short-lived ``read_only=True`` DuckDB handle, runs one
query, and closes it. If the live bot is mid-checkpoint and the file is briefly
locked, ``_read`` returns the ``LOCKED`` sentinel instead of raising — callbacks
treat that as "no new data, keep the last render".

Bot timestamps are naive local wall-clock (verified: the maker's newest decision
time matches the local clock). So time-window cutoffs are computed in Python with
``datetime.now()`` and passed as parameters — never compared against SQL ``now()``
(which is timezone-aware and would mis-cast).
"""

from __future__ import annotations

import subprocess
import time
from datetime import datetime, timedelta

import duckdb
import pandas as pd

from .bots import BOTS


class _Locked:
    """Sentinel returned when a read can't proceed (lock/IO error)."""

    def __bool__(self) -> bool:  # so `if df is LOCKED` reads naturally
        return False


LOCKED = _Locked()

WINDOWS = {"1h": timedelta(hours=1), "24h": timedelta(days=1),
           "7d": timedelta(days=7), "all": None}


def window_cutoff(window: str):
    """Return the naive-local datetime cutoff for a window key, or None for 'all'."""
    delta = WINDOWS.get(window)
    return None if delta is None else datetime.now() - delta


def _read(db_path: str, sql: str, params=None, retries: int = 3):
    """Run a read-only query → DataFrame, or LOCKED after retries.

    The live bot briefly holds an exclusive lock while DuckDB checkpoints, which
    makes a fresh read_only connect fail for a few milliseconds. We retry a few
    times with short backoff so those sub-second locks don't surface as empty
    data (which would make a *live* bot look dormant). Only a genuinely stuck
    lock returns LOCKED.
    """
    for attempt in range(retries):
        con = None
        try:
            con = duckdb.connect(db_path, read_only=True)
            try:
                con.execute("PRAGMA disable_progress_bar")  # keep logs clean
            except Exception:
                pass
            return con.execute(sql, params or []).fetchdf()
        except Exception:
            if attempt < retries - 1:
                time.sleep(0.06)  # ~60ms; checkpoint locks clear fast
                continue
            return LOCKED
        finally:
            if con is not None:
                try:
                    con.close()
                except Exception:
                    pass


def signature(bot_key: str) -> str | None:
    """Cheap fingerprint of a bot's state — used by the dashboard's poll guard to
    skip re-rendering (and the flash that comes with it) when nothing changed.

    Combines newest decision time + row counts. Counts over a few hundred k rows
    are fast in DuckDB. Returns None if the DB is locked so callers re-render.
    """
    bot = BOTS[bot_key]
    df = _read(bot["state_db"],
               f"SELECT max({bot['decision_ts']}) AS t, count(*) AS n "
               f"FROM {bot['decision_table']}")
    if df is LOCKED or not len(df):
        return None
    r = df.iloc[0]
    f = _read(bot["state_db"], f"SELECT count(*) AS n FROM {bot['fills_table']}")
    p = _read(bot["state_db"], f"SELECT count(*) AS n FROM {bot['positions_table']}")
    fn = None if (f is LOCKED or not len(f)) else int(f.iloc[0]["n"])
    pn = None if (p is LOCKED or not len(p)) else int(p.iloc[0]["n"])
    t = r["t"]
    nn = int(r["n"]) if pd.notna(r["n"]) else 0
    return f"{t}|{nn}|{fn}|{pn}"


def is_running(bot_key: str) -> bool:
    """True if the bot's process is currently running (best-effort pgrep)."""
    pat = BOTS[bot_key]["proc_match"]
    try:
        out = subprocess.run(["pgrep", "-f", pat], capture_output=True, text=True, timeout=3)
        return out.returncode == 0 and bool(out.stdout.strip())
    except Exception:
        return False


def _category_expr(bot) -> str:
    """SQL expression for the meaningful 'why' bucket for this bot.

    Maker: the reason matters (out_of_scope / no_fair / size_gate…); a quoted row
    has NULL reason so we fold it to its decision. Taker: the decision string is
    already the granular outcome (declined_ev, failed_quote_walked, …).
    """
    if bot["kind"] == "maker":
        return f"COALESCE({bot['reason_col']}, {bot['decision_col']})"
    return bot["decision_col"]


def _win_clause(ts_col: str, cutoff):
    """Return (sql_fragment, params) restricting ts_col to >= cutoff (or no-op)."""
    if cutoff is None:
        return "", []
    return f" AND {ts_col} >= ?", [cutoff]


# --------------------------------------------------------------------------- #
# Status / freshness
# --------------------------------------------------------------------------- #
def status(bot_key: str) -> dict:
    bot = BOTS[bot_key]
    out = {"running": is_running(bot_key), "last_activity": None,
           "session_start": None, "session_end": None, "dry_run": None}

    df = _read(bot["state_db"],
               f"SELECT max({bot['decision_ts']}) AS t FROM {bot['decision_table']}")
    if df is not LOCKED and len(df) and pd.notna(df.iloc[0]["t"]):
        out["last_activity"] = df.iloc[0]["t"].to_pydatetime()

    df = _read(bot["state_db"],
               "SELECT started_at, ended_at, dry_run FROM sessions "
               "ORDER BY started_at DESC LIMIT 1")
    if df is not LOCKED and len(df):
        r = df.iloc[0]
        out["session_start"] = r["started_at"].to_pydatetime() if pd.notna(r["started_at"]) else None
        out["session_end"] = r["ended_at"].to_pydatetime() if pd.notna(r["ended_at"]) else None
        out["dry_run"] = bool(r["dry_run"]) if pd.notna(r["dry_run"]) else None
    return out


# --------------------------------------------------------------------------- #
# KPIs + funnel
# --------------------------------------------------------------------------- #
def kpis(bot_key: str, cutoff) -> dict:
    bot = BOTS[bot_key]
    wc, wp = _win_clause(bot["rfq_ts"], cutoff)
    wcd, wpd = _win_clause(bot["decision_ts"], cutoff)
    wcf, wpf = _win_clause(bot["fill_ts"], cutoff)
    db = bot["state_db"]

    def scalar(path, sql, params):
        df = _read(path, sql, params)
        if df is LOCKED or not len(df):
            return None
        v = df.iloc[0, 0]
        return None if pd.isna(v) else v

    k = {}
    if bot["kind"] == "maker":
        k["rfqs"] = scalar(db, f"SELECT count(*) FROM {bot['rfq_table']} WHERE 1=1{wc}", wp)
        k["in_scope"] = scalar(db, f"SELECT count(*) FROM {bot['rfq_table']} WHERE in_scope{wc}", wp)
        k["quotes"] = scalar(
            db,
            f"SELECT count(DISTINCT rfq_id) FROM {bot['decision_table']} "
            f"WHERE {bot['decision_col']}='quoted'{wcd}", wpd)
    else:
        k["rfqs"] = scalar(db, f"SELECT count(*) FROM {bot['rfq_table']} WHERE 1=1{wc}", wp)
        k["in_scope"] = scalar(db, f"SELECT count(*) FROM {bot['decision_table']} WHERE 1=1{wcd}", wpd)
        k["quotes"] = k["in_scope"]  # taker "evaluated" == quotes considered

    k["fills"] = scalar(db, f"SELECT count(*) FROM {bot['fills_table']} WHERE 1=1{wcf}", wpf)
    k["staked"] = scalar(
        db, f"SELECT sum({bot['fill_price']}*contracts) FROM {bot['fills_table']} WHERE 1=1{wcf}", wpf)
    k["open_positions"] = scalar(db, f"SELECT count(*) FROM {bot['positions_table']}", [])
    # NaN-safe: DuckDB stores NaN (not NULL) in weighted_price, and SUM
    # propagates NaN, so treat NaN as 0 before summing.
    k["exposure"] = scalar(
        db,
        f"SELECT sum(abs(net_contracts)*CASE WHEN isnan(weighted_price) THEN 0 "
        f"ELSE weighted_price END) FROM {bot['positions_table']}", [])

    denom = k.get("quotes") or 0
    k["fill_rate"] = (k["fills"] / denom) if denom else None
    return k


def funnel(bot_key: str, cutoff) -> list:
    """Ordered [(stage_label, count), ...] for the conversion funnel."""
    bot = BOTS[bot_key]
    k = kpis(bot_key, cutoff)
    if bot["kind"] == "maker":
        return [("RFQs seen", k["rfqs"]), ("In scope", k["in_scope"]),
                ("Quoted", k["quotes"]), ("Filled", k["fills"])]
    return [("RFQs sent", k["rfqs"]), ("Quotes evaluated", k["in_scope"]),
            ("Filled (accepted)", k["fills"])]


# --------------------------------------------------------------------------- #
# Why-not-filled: decision breakdown + timeseries
# --------------------------------------------------------------------------- #
def decision_breakdown(bot_key: str, cutoff) -> "pd.DataFrame | _Locked":
    bot = BOTS[bot_key]
    cat = _category_expr(bot)
    wcd, wpd = _win_clause(bot["decision_ts"], cutoff)
    sql = (f"SELECT {cat} AS category, count(*) AS n "
           f"FROM {bot['decision_table']} WHERE 1=1{wcd} "
           f"GROUP BY 1 ORDER BY n DESC")
    return _read(bot["state_db"], sql, wpd)


def decision_timeseries(bot_key: str, cutoff, window: str) -> "pd.DataFrame | _Locked":
    bot = BOTS[bot_key]
    cat = _category_expr(bot)
    bucket = {"1h": "minute", "24h": "hour", "7d": "day", "all": "day"}[window]
    wcd, wpd = _win_clause(bot["decision_ts"], cutoff)
    sql = (f"SELECT date_trunc('{bucket}', {bot['decision_ts']}) AS bucket, "
           f"       {cat} AS category, count(*) AS n "
           f"FROM {bot['decision_table']} WHERE 1=1{wcd} "
           f"GROUP BY 1,2 ORDER BY 1")
    return _read(bot["state_db"], sql, wpd)


# --------------------------------------------------------------------------- #
# Fills & P&L
# --------------------------------------------------------------------------- #
def recent_fills(bot_key: str, limit: int = 200) -> "pd.DataFrame | _Locked":
    bot = BOTS[bot_key]
    pnl = f", {bot['fill_pnl']} AS realized_pnl" if bot["fill_pnl"] else ""
    sql = (f"SELECT {bot['fill_ts']} AS filled_at, game_id, combo_market_ticker, "
           f"       {bot['fill_side']} AS side, contracts, "
           f"       {bot['fill_price']} AS price, {bot['fill_fee']} AS fee, "
           f"       {bot['fill_fair_quote']} AS fair_at_quote{pnl} "
           f"FROM {bot['fills_table']} ORDER BY {bot['fill_ts']} DESC LIMIT {int(limit)}")
    return _read(bot["state_db"], sql)


def fills_cumulative(bot_key: str) -> "pd.DataFrame | _Locked":
    bot = BOTS[bot_key]
    sql = (f"SELECT {bot['fill_ts']} AS filled_at, "
           f"       {bot['fill_price']}*contracts AS stake "
           f"FROM {bot['fills_table']} ORDER BY {bot['fill_ts']}")
    df = _read(bot["state_db"], sql)
    if df is LOCKED or not len(df):
        return df
    df["cum_fills"] = range(1, len(df) + 1)
    df["cum_stake"] = df["stake"].cumsum()
    return df


# --------------------------------------------------------------------------- #
# Positions & open orders
# --------------------------------------------------------------------------- #
def positions(bot_key: str) -> "pd.DataFrame | _Locked":
    bot = BOTS[bot_key]
    legs = ", legs_json" if bot["has_position_legs"] else ""
    # NaN-safe exposure (see kpis): NaN weighted_price -> 0 so sorting/sum behave.
    sql = (f"SELECT combo_market_ticker, side, game_id, net_contracts, weighted_price, "
           f"       abs(net_contracts)*CASE WHEN isnan(weighted_price) THEN 0 "
           f"       ELSE weighted_price END AS exposure, updated_at{legs} "
           f"FROM {bot['positions_table']} ORDER BY exposure DESC")
    return _read(bot["state_db"], sql)


def open_orders(bot_key: str) -> "pd.DataFrame | _Locked":
    bot = BOTS[bot_key]
    if bot["kind"] == "maker":
        sql = ("SELECT submitted_at, combo_market_ticker, game_id, "
               "       yes_bid, no_bid, blended_fair, status, closed_at "
               "FROM live_quotes ORDER BY submitted_at DESC LIMIT 200")
    else:
        sql = ("SELECT submitted_at, combo_market_ticker, game_id, intended_side, "
               "       edge_at_submit, blended_fair_at_submit, status, cancellation_reason "
               "FROM live_rfqs WHERE status='open' ORDER BY submitted_at DESC LIMIT 200")
    return _read(bot["state_db"], sql)


def open_orders_summary(bot_key: str) -> "pd.DataFrame | _Locked":
    """Status counts for working orders (so we can flag orphaned 'open' rows)."""
    bot = BOTS[bot_key]
    if bot["kind"] == "maker":
        sql = "SELECT status, count(*) AS n FROM live_quotes GROUP BY 1 ORDER BY n DESC"
    else:
        sql = ("SELECT status, count(*) AS n FROM live_rfqs GROUP BY 1 ORDER BY n DESC")
    return _read(bot["state_db"], sql)


# --------------------------------------------------------------------------- #
# Adverse selection / fill quality
# --------------------------------------------------------------------------- #
def adverse(bot_key: str) -> "pd.DataFrame | _Locked":
    """Maker: per-fill fair drift (confirm - quote) vs the quoted margin.
    Taker: per-fill edge_at_fill proxy via blended_fair vs price.
    """
    bot = BOTS[bot_key]
    if bot["kind"] == "maker":
        sql = (f"SELECT {bot['fill_ts']} AS filled_at, combo_market_ticker, "
               f"       {bot['fill_fair_quote']} AS fair_at_quote, "
               f"       {bot['fill_fair_confirm']} AS fair_at_confirm, "
               f"       {bot['fill_fair_confirm']} - {bot['fill_fair_quote']} AS fair_drift, "
               f"       {bot['fill_price']} AS price, {bot['fill_pnl']} AS realized_pnl "
               f"FROM {bot['fills_table']} ORDER BY {bot['fill_ts']} DESC LIMIT 500")
    else:
        sql = (f"SELECT {bot['fill_ts']} AS filled_at, combo_market_ticker, "
               f"       {bot['fill_fair_quote']} AS fair_at_quote, "
               f"       {bot['fill_price']} AS price, expected_ev_dollars "
               f"FROM {bot['fills_table']} ORDER BY {bot['fill_ts']} DESC LIMIT 500")
    return _read(bot["state_db"], sql)


def taker_quote_quality(bot_key: str, cutoff) -> "pd.DataFrame | _Locked":
    """Taker-only: walk vs accept counts over time + how close walks were."""
    bot = BOTS[bot_key]
    wcd, wpd = _win_clause(bot["decision_ts"], cutoff)
    sql = (f"SELECT date_trunc('hour', {bot['decision_ts']}) AS bucket, "
           f"       {bot['decision_col']} AS decision, count(*) AS n "
           f"FROM {bot['decision_table']} "
           f"WHERE {bot['decision_col']} IN "
           f"      ('accepted','failed_quote_walked','halted_low_fill_ratio'){wcd} "
           f"GROUP BY 1,2 ORDER BY 1")
    return _read(bot["state_db"], sql, wpd)
