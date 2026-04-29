from datetime import datetime, timezone

from kalshi_mlb_rfq import ticker_map


def test_event_suffix_round_trip():
    # Apr 28 2026 at 8:05 PM ET → 26APR282005
    commence = datetime(2026, 4, 29, 0, 5, tzinfo=timezone.utc)  # 8:05 PM ET
    suffix = ticker_map.format_event_suffix(commence, away_code="NYY", home_code="TEX")
    assert suffix == "26APR282005NYYTEX"


def test_spread_ticker_for_home_team():
    t = ticker_map.spread_ticker("26APR282005NYYTEX", team_code="TEX", line=-1.5)
    assert t == "KXMLBSPREAD-26APR282005NYYTEX-TEX2"


def test_spread_ticker_for_away_team_alt_line():
    t = ticker_map.spread_ticker("26APR282005NYYTEX", team_code="NYY", line=-2.5)
    assert t == "KXMLBSPREAD-26APR282005NYYTEX-NYY3"


def test_total_ticker():
    t = ticker_map.total_ticker("26APR282005NYYTEX", line=7.5)
    assert t == "KXMLBTOTAL-26APR282005NYYTEX-8"


def test_total_ticker_alt():
    t = ticker_map.total_ticker("26APR282005NYYTEX", line=10.5)
    assert t == "KXMLBTOTAL-26APR282005NYYTEX-11"
