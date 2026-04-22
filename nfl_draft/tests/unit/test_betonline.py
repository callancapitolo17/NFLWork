"""Offline parser tests for scrapers/betonline.py.

Drives the parser off the committed fixture; never hits the network.
"""

import json
from pathlib import Path
import pytest

from nfl_draft.scrapers.betonline import parse_response

FIXTURE = Path(__file__).resolve().parent.parent / "fixtures" / "betonline" / "draft_markets.json"


@pytest.fixture(scope="module")
def rows():
    envelope = json.loads(FIXTURE.read_text())
    return parse_response(envelope)


def test_emits_pick_outright_rows(rows):
    pick = [r for r in rows if r.market_group == "pick_outright"]
    # At least the 9 Nth-Overall Pick markets × their contestants.
    assert len(pick) >= 100, f"expected >= 100 pick_outright rows, got {len(pick)}"
    # book_label is the market description (e.g. '10th Overall Pick').
    labels = {r.book_label for r in pick}
    assert "10th Overall Pick" in labels
    assert "2nd Overall Pick" in labels
    # Spot-check one: Caleb Downs at +350 for 10th Overall (from brainstorm capture).
    cd = [r for r in pick if r.book_label == "10th Overall Pick" and r.book_subject == "Caleb Downs"]
    assert cd and cd[0].american_odds == 350


def test_all_rows_have_betonline_book_and_valid_odds(rows):
    for r in rows:
        assert r.book == "betonline"
        assert isinstance(r.american_odds, int)
        assert r.american_odds != 0
        assert r.book_subject
        assert r.book_label
