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


def test_emits_team_first_pick(rows):
    """team-to-draft has 'Team to Draft <Player>' with team-name runners.
    Market_group is ``team_first_pick`` (not ``team_drafts_player``) so
    ``outright_group_key`` routes the rows through the per-team mutex
    devig bucket.
    """
    tfp = [r for r in rows if r.market_group == "team_first_pick"]
    assert len(tfp) >= 200, f"expected >= 200 team_first_pick rows, got {len(tfp)}"
    for r in tfp[:5]:
        assert "team to draft" in r.book_label.lower()


def test_emits_team_first_pick_position(rows):
    """teams-1st-drafted-position has '<Team> 1st Drafted Player Position'
    with position-word runners."""
    tfp = [r for r in rows if r.market_group == "team_first_pick_position"]
    assert len(tfp) >= 250, f"expected >= 250 team_first_pick_position rows, got {len(tfp)}"
    for r in tfp[:5]:
        assert "1st drafted" in r.book_label.lower(), r.book_label


def test_emits_matchup_before_pairs(rows):
    """matchups has 'To Be Drafted First' pairs. Synthesized label must
    contain both player names so MARKET_MAP can re-pair them."""
    mb = [r for r in rows if r.market_group == "matchup_before"]
    assert len(mb) >= 8, f"expected >= 8 matchup_before rows, got {len(mb)}"
    for r in mb[:4]:
        assert " vs " in r.book_label, r.book_label
    # Sanity: every matchup pair has exactly 2 rows sharing the same label.
    from collections import Counter
    label_counts = Counter(r.book_label for r in mb)
    assert all(count == 2 for count in label_counts.values()), label_counts


def test_emits_draft_position_ou(rows):
    """draft-position encodes the line into subject as 'Over 9.5' / 'Under 9.5'
    so the MARKET_MAP builder can reconstruct (player, line, direction) from
    (label, subject) alone."""
    dpo = [r for r in rows if r.market_group == "draft_position_over_under"]
    assert len(dpo) >= 30, f"expected >= 30 draft_position_ou rows, got {len(dpo)}"
    for r in dpo[:5]:
        assert r.book_subject.startswith(("Over ", "Under ")), (
            f"subject {r.book_subject!r} must start with 'Over '/'Under ' + line"
        )
        assert "draft position" in r.book_label.lower()


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


def test_specials_bucket_emits_prop_rows(rows):
    """specials are heterogeneous one-offs; the generic prop fallback
    (slug-keyed) emits prop_specials rows with the raw Description as
    book_label — operator can audit and canonicalize later."""
    specials = [r for r in rows if r.market_group == "prop_specials"]
    assert len(specials) >= 10, f"expected >= 10 prop_specials rows, got {len(specials)}"
    # Some specials include 'How Many Trades Will Occur in Round 1 of Draft'.
    assert any("trades" in r.book_label.lower() for r in specials), (
        "expected at least one trade-count special in the fixture"
    )


def test_total_row_count_close_to_snapshot(rows):
    """Regression guard: fixture snapshot had 902 runners (captured 2026-04-21).
    Allow +/-15% drift for minor parser logic adjustments."""
    assert 750 <= len(rows) <= 1050, f"unexpected total row count: {len(rows)}"


def test_v1_market_groups_are_allowlisted(rows):
    """Every emitted market_group should be either a recognized v1 canonical
    OR a prop_* fallback. Anything else is a bug."""
    allowed_canonical = {
        "pick_outright", "first_at_position",
        "top_5_range", "top_10_range", "top_32_range",
        "nth_at_position_2", "nth_at_position_3",
        "mr_irrelevant_position",
        "team_first_pick", "team_first_pick_position",
        "matchup_before",
        "draft_position_over_under",
    }
    actual = {r.market_group for r in rows}
    props = {g for g in actual if g.startswith("prop_")}
    unexpected = actual - allowed_canonical - props
    assert not unexpected, (
        f"unexpected market_groups not in v1 allowlist: {unexpected}"
    )
