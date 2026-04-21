"""End-to-end test: Kalshi trade tape dedups by trade_id PK via INSERT OR IGNORE.

The kalshi_trades table has trade_id as PRIMARY KEY, so re-inserting the same
trade_id twice should be a no-op on the second pass. We exercise that path
directly so regressions in the schema (e.g. someone dropping the PK) surface
immediately.
"""
import pytest
from datetime import datetime
from nfl_draft.lib import db as db_module
from nfl_draft.lib import seed
from nfl_draft.scrapers._base import TradeRow


@pytest.fixture
def seeded(monkeypatch, tmp_path):
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "test.duckdb")
    db_module.init_schema()
    seed.run()


def test_trades_insert_or_ignore_dedup(seeded):
    """Second INSERT OR IGNORE with the same trade_id must be a no-op."""
    t = TradeRow(
        trade_id="abc123",
        ticker="KXNFLDRAFT1-T1",
        side="yes",
        price_cents=50,
        count=10,
        traded_at=datetime.now(),
        fetched_at=datetime.now(),
    )
    with db_module.write_connection() as con:
        for _ in range(2):
            con.execute(
                "INSERT OR IGNORE INTO kalshi_trades VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                [
                    t.trade_id, t.ticker, t.side, t.price_cents, t.count,
                    t.count * t.price_cents * 0.01, t.traded_at, t.fetched_at,
                ],
            )
    with db_module.read_connection() as con:
        count = con.execute(
            "SELECT COUNT(*) FROM kalshi_trades WHERE trade_id='abc123'"
        ).fetchone()[0]
    assert count == 1  # second insert was deduped
