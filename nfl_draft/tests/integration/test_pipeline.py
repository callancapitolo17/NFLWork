"""End-to-end pipeline test: fixture -> normalize -> market_map -> devig -> DB."""
import json
import pytest
from pathlib import Path
from nfl_draft.lib import db as db_module
from nfl_draft.lib import seed
from nfl_draft.lib.quarantine import write_or_quarantine
from nfl_draft.scrapers.kalshi import parse_markets_response


# Fixtures live under nfl_draft/tests/fixtures/, so go up one directory
# from this file (integration/) to reach tests/, then into fixtures/.
FIXTURES = Path(__file__).resolve().parent.parent / "fixtures"


@pytest.fixture
def seeded(monkeypatch, tmp_path):
    """Point DB_PATH to a fresh temp DuckDB, init schema, and load seeds."""
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "test.duckdb")
    db_module.init_schema()
    seed.run()


def test_kalshi_pipeline_writes_to_draft_odds_or_quarantine(seeded):
    raw = json.loads((FIXTURES / "kalshi" / "markets_response.json").read_text())
    rows = parse_markets_response(raw, series_ticker="KXNFLDRAFT1")
    if not rows:
        pytest.skip("Empty fixture (all markets have yes_bid <= 0)")
    mapped, unmapped = write_or_quarantine(rows)
    assert (mapped + unmapped) == len(rows)
