"""Tests for each scraper's parse() function - pure logic, no network."""

import json
from pathlib import Path

from nfl_draft.scrapers.kalshi import parse_markets_response


FIXTURES = Path(__file__).resolve().parent.parent / "fixtures"


def test_kalshi_parse_markets_returns_oddsrow_list():
    raw = json.loads((FIXTURES / "kalshi" / "markets_response.json").read_text())
    rows = parse_markets_response(raw, series_ticker="KXNFLDRAFT1")
    assert isinstance(rows, list)
    if rows:
        row = rows[0]
        assert row.book == "kalshi"
        assert row.book_label.startswith("KXNFLDRAFT1")
