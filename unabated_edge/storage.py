import time
import random
import json
import logging
import datetime
import duckdb
from contextlib import contextmanager

from unabated_edge import config

log = logging.getLogger("unabated_edge")
_BUFFER = []


@contextmanager
def connect(path, read_only=False, retries=10):
    attempt = 0
    while True:
        try:
            con = duckdb.connect(str(path), read_only=read_only)
            break
        except duckdb.IOException:
            if attempt >= retries:
                raise
            time.sleep(0.05 * 2**attempt + random.random() * 0.05)
            attempt += 1
    try:
        yield con
    finally:
        con.close()


def init():
    with connect(config.MARKET_DB_PATH) as c:
        c.execute("""CREATE TABLE IF NOT EXISTS line_snapshots(
            ts TIMESTAMPTZ, sport VARCHAR, event_id BIGINT, market_source_id INTEGER,
            bet_type VARCHAR, side VARCHAR, price DOUBLE, points DOUBLE,
            modified_on VARCHAR)""")
        c.execute("""CREATE TABLE IF NOT EXISTS flagged_edges(
            ts TIMESTAMPTZ, sport VARCHAR, event_id BIGINT, market_ticker VARCHAR, outcome VARCHAR,
            fair_prob DOUBLE, yes_ask DOUBLE, ev_pct DOUBLE, kelly_contracts INTEGER, dry_run BOOLEAN)""")
        # Kalshi-side microstructure capture (maker design data): top-of-book
        # columns for easy querying + the full depth ladder as JSON. Pre-kickoff
        # only, same close semantics as line_snapshots.
        c.execute("""CREATE TABLE IF NOT EXISTS book_snapshots(
            ts TIMESTAMPTZ, sport VARCHAR, market_ticker VARCHAR, floor_strike DOUBLE,
            yes_bid DOUBLE, yes_bid_qty DOUBLE, no_bid DOUBLE, no_bid_qty DOUBLE,
            yes_ask DOUBLE, no_ask DOUBLE, volume DOUBLE, open_interest DOUBLE,
            depth JSON)""")
        # Executed-trades tape ("where flow trades" + taker direction). PK on
        # trade_id so overlapping poll windows dedup via INSERT OR IGNORE.
        # created_time kept as the raw ISO string (same convention as
        # line_snapshots.modified_on — ISO-8601 UTC sorts lexicographically).
        c.execute("""CREATE TABLE IF NOT EXISTS kalshi_trades(
            trade_id VARCHAR PRIMARY KEY, sport VARCHAR, market_ticker VARCHAR,
            created_time VARCHAR, yes_price DOUBLE, count DOUBLE, taker_side VARCHAR,
            captured_ts TIMESTAMPTZ)""")
    with connect(config.RESEARCH_DB_PATH) as c:
        c.execute("""CREATE TABLE IF NOT EXISTS research_events(
            ts TIMESTAMPTZ, event_type VARCHAR, sport VARCHAR, payload JSON)""")


def snapshot_lines(sport, rows):
    if not rows:
        return
    with connect(config.MARKET_DB_PATH) as c:
        c.executemany(
            "INSERT INTO line_snapshots VALUES (?,?,?,?,?,?,?,?,?)",
            [
                [r["ts"], sport, r["event_id"], r["market_source_id"],
                 r["bet_type"], r["side"], r["price"], r["points"], r.get("modified_on")]
                for r in rows
            ],
        )


def snapshot_books(sport, rows):
    if not rows:
        return
    with connect(config.MARKET_DB_PATH) as c:
        c.executemany(
            "INSERT INTO book_snapshots VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
            [
                [r["ts"], sport, r["market_ticker"], r["floor_strike"],
                 r["yes_bid"], r["yes_bid_qty"], r["no_bid"], r["no_bid_qty"],
                 r["yes_ask"], r["no_ask"], r["volume"], r["open_interest"],
                 json.dumps(r["depth"])]
                for r in rows
            ],
        )


def insert_trades(sport, rows, captured_ts):
    if not rows:
        return
    with connect(config.MARKET_DB_PATH) as c:
        c.executemany(
            "INSERT OR IGNORE INTO kalshi_trades VALUES (?,?,?,?,?,?,?,?)",
            [
                [r["trade_id"], sport, r["market_ticker"], r["created_time"],
                 r["yes_price"], r["count"], r["taker_side"], captured_ts]
                for r in rows
            ],
        )


def log_flagged(r):
    with connect(config.MARKET_DB_PATH) as c:
        c.execute(
            "INSERT INTO flagged_edges VALUES (?,?,?,?,?,?,?,?,?,?)",
            [r["ts"], r["sport"], r["event_id"], r["market_ticker"], r["outcome"],
             r["fair_prob"], r["yes_ask"], r["ev_pct"], r["kelly_contracts"], r["dry_run"]],
        )


def emit(event_type, sport, **payload):
    try:
        _BUFFER.append((
            datetime.datetime.now(datetime.timezone.utc),
            event_type,
            sport,
            json.dumps(payload),
        ))
    except Exception:
        pass


def flush():
    if not _BUFFER:
        return
    batch, _BUFFER[:] = list(_BUFFER), []
    try:
        with connect(config.RESEARCH_DB_PATH) as c:
            c.executemany("INSERT INTO research_events VALUES (?,?,?,?)", batch)
    except Exception:
        # Never raise into the trading loop, but don't drop silently — the
        # firehose is the dry-run's analytical record; losing it must be visible.
        log.warning("research flush failed, dropped %d rows", len(batch))
