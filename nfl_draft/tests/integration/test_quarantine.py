"""End-to-end test: unmapped player + unmapped market land in quarantine, NOT draft_odds."""
import pytest
from datetime import datetime
from nfl_draft.lib import db as db_module
from nfl_draft.lib import seed
from nfl_draft.lib.quarantine import write_or_quarantine
from nfl_draft.scrapers._base import OddsRow


@pytest.fixture
def seeded(monkeypatch, tmp_path):
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "test.duckdb")
    db_module.init_schema()
    seed.run()


def test_unmapped_player_lands_in_quarantine(seeded):
    rows = [OddsRow(
        book="draftkings", book_label="First QB", book_subject="Made Up Player",
        american_odds=200, fetched_at=datetime.now(),
    )]
    write_or_quarantine(rows)
    with db_module.read_connection() as con:
        odds_count = con.execute("SELECT COUNT(*) FROM draft_odds").fetchone()[0]
        quarantine_count = con.execute("SELECT COUNT(*) FROM draft_odds_unmapped").fetchone()[0]
    assert odds_count == 0
    assert quarantine_count == 1


def test_mapped_row_lands_in_draft_odds(seeded):
    from nfl_draft.lib.db import write_connection
    with write_connection() as con:
        # Insert a market row first to satisfy any FK
        con.execute("INSERT INTO draft_markets (market_id, market_type, subject_player, position) VALUES ('first_qb_cam-ward', 'first_at_position', 'Cam Ward', 'QB')")
        con.execute("INSERT INTO market_map VALUES ('draftkings', 'First QB', 'Cam Ward', 'first_qb_cam-ward')")
    rows = [OddsRow(
        book="draftkings", book_label="First QB", book_subject="Cam Ward",
        american_odds=200, fetched_at=datetime.now(),
    )]
    write_or_quarantine(rows)
    with db_module.read_connection() as con:
        odds_count = con.execute("SELECT COUNT(*) FROM draft_odds").fetchone()[0]
        quarantine_count = con.execute("SELECT COUNT(*) FROM draft_odds_unmapped").fetchone()[0]
    assert odds_count == 1
    assert quarantine_count == 0
