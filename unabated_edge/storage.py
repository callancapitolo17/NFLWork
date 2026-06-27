import time
import random
import json
import datetime
import duckdb
from contextlib import contextmanager

from unabated_edge import config

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
            bet_type VARCHAR, side VARCHAR, price DOUBLE, points DOUBLE)""")
        c.execute("""CREATE TABLE IF NOT EXISTS flagged_edges(
            ts TIMESTAMPTZ, sport VARCHAR, event_id BIGINT, market_ticker VARCHAR, outcome VARCHAR,
            fair_prob DOUBLE, yes_ask DOUBLE, ev_pct DOUBLE, kelly_contracts INTEGER, dry_run BOOLEAN)""")
    with connect(config.RESEARCH_DB_PATH) as c:
        c.execute("""CREATE TABLE IF NOT EXISTS research_events(
            ts TIMESTAMPTZ, event_type VARCHAR, sport VARCHAR, payload JSON)""")


def snapshot_lines(sport, rows):
    if not rows:
        return
    with connect(config.MARKET_DB_PATH) as c:
        c.executemany(
            "INSERT INTO line_snapshots VALUES (?,?,?,?,?,?,?,?)",
            [
                [r["ts"], sport, r["event_id"], r["market_source_id"],
                 r["bet_type"], r["side"], r["price"], r["points"]]
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
        pass
