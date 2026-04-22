"""Fixture sanity checks — captured offline, parser subagent lives here."""

import json
from pathlib import Path

FIXTURES = Path(__file__).resolve().parent.parent / "fixtures"
BETONLINE_FIXTURE = FIXTURES / "betonline" / "draft_markets.json"


def test_betonline_fixture_present_and_complete():
    """Recon script should have captured all 11 NFL Draft buckets."""
    assert BETONLINE_FIXTURE.exists(), (
        f"Missing fixture at {BETONLINE_FIXTURE}. "
        "Run: python nfl_draft/scrapers/recon_betonline.py"
    )
    envelope = json.loads(BETONLINE_FIXTURE.read_text())
    assert envelope["book"] == "betonline"
    buckets = envelope["data"]
    assert isinstance(buckets, dict)

    # All 11 market buckets captured, each a full offering envelope.
    expected_slugs = {
        "1st-round", "1st-round-props", "draft-position", "matchups",
        "mr-irrelevant", "specials", "team-to-draft",
        "teams-1st-drafted-position", "to-be-drafted-1st",
        "to-be-drafted-2nd", "to-be-selected",
    }
    assert set(buckets.keys()) == expected_slugs, (
        f"fixture buckets != expected. got {sorted(buckets)}"
    )

    # Each bucket should have a ContestOfferings dict (not null/empty).
    for slug, payload in buckets.items():
        co = payload.get("ContestOfferings") or {}
        dg = co.get("DateGroup") or []
        assert dg, f"bucket {slug!r} has empty DateGroup"
