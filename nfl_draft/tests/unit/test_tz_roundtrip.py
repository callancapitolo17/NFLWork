from datetime import datetime
from nfl_draft.lib.db import local_to_utc_iso, utc_iso_to_local


def test_local_to_utc_iso_returns_string_with_offset():
    local = datetime(2026, 4, 23, 12, 0, 0)
    iso = local_to_utc_iso(local)
    assert isinstance(iso, str)
    assert iso.endswith("+00:00") or iso.endswith("Z")


def test_utc_iso_to_local_roundtrip():
    original_local = datetime(2026, 4, 23, 12, 0, 0)
    iso = local_to_utc_iso(original_local)
    roundtripped = utc_iso_to_local(iso)
    assert roundtripped == original_local


def test_utc_iso_to_local_handles_kalshi_format():
    # Kalshi returns 2026-04-23T16:00:00Z — should convert to local
    result = utc_iso_to_local("2026-04-23T16:00:00Z")
    assert result.tzinfo is None  # naive local
