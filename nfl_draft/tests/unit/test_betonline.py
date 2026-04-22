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


def test_emits_first_at_position_from_to_be_drafted_1st(rows):
    """to-be-drafted-1st has descriptions like 'First Wide Receiver Drafted'
    with player contestants. Maps to canonical first_at_position."""
    fap = [r for r in rows if r.market_group == "first_at_position"]
    assert len(fap) >= 10, f"expected >= 10 first_at_position rows, got {len(fap)}"
    labels = {r.book_label for r in fap}
    # The fixture has CB, OL, WR first-drafted markets at minimum.
    assert any("cornerback" in lbl.lower() for lbl in labels), labels
    assert any("wide receiver" in lbl.lower() for lbl in labels), labels


def test_emits_nth_at_position_2(rows):
    """to-be-drafted-2nd has '2nd X Selected' with player runners."""
    nth = [r for r in rows if r.market_group == "nth_at_position_2"]
    assert len(nth) >= 20, f"expected >= 20 nth_at_position_2 rows, got {len(nth)}"
    for r in nth:
        assert "2nd" in r.book_label.lower(), r.book_label


def test_emits_top_n_range_from_to_be_selected(rows):
    """to-be-selected has 'Drafted Top 5/10' and 'Drafted in Round 1' with
    player runners. Maps to top_N_range for N in {5, 10, 32}."""
    top_rows = [r for r in rows
                if r.market_group.startswith("top_") and r.market_group.endswith("_range")]
    assert len(top_rows) >= 30, f"expected >= 30 top_N_range rows, got {len(top_rows)}"
    groups = {r.market_group for r in top_rows}
    assert groups == {"top_5_range", "top_10_range", "top_32_range"}, (
        f"expected exactly top_{{5,10,32}}_range, got {groups}"
    )


def test_emits_mr_irrelevant_position(rows):
    """mr-irrelevant contestants are position words, not players."""
    mi = [r for r in rows if r.market_group == "mr_irrelevant_position"]
    assert len(mi) >= 5, f"expected >= 5 mr_irrelevant rows, got {len(mi)}"
    subjects = {r.book_subject.lower() for r in mi}
    assert any("receiver" in s or "lineman" in s or "back" in s or "linebacker" in s
               for s in subjects), subjects


def test_1st_round_props_emits_totals_props_with_lines(rows):
    """1st-round-props are totals-style ('Total X Drafted in 1st Round' with
    O/U + GroupLine) — no canonical join today, so we emit as props with
    line encoded in subject for future canonicalization."""
    totals = [r for r in rows if r.market_group == "prop_first_round_total_ou"]
    assert len(totals) >= 20, f"expected >= 20 1st-round-props rows, got {len(totals)}"
    for r in totals[:5]:
        assert r.book_subject.startswith(("Over ", "Under ")), (
            f"subject {r.book_subject!r} must start with 'Over '/'Under ' + line"
        )
