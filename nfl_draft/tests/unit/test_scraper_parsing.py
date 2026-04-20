"""Tests for each scraper's parse() function - pure logic, no network."""

import json
from pathlib import Path

from nfl_draft.scrapers.kalshi import parse_markets_response
from nfl_draft.scrapers.draftkings import parse_response as dk_parse
from nfl_draft.scrapers.fanduel import parse_response as fd_parse
from nfl_draft.scrapers.bookmaker import parse_response as bm_parse
from nfl_draft.scrapers.wagerzon import parse_response as wz_parse
from nfl_draft.scrapers.hoop88 import parse_response as h88_parse


FIXTURES = Path(__file__).resolve().parent.parent / "fixtures"


def test_kalshi_parse_markets_returns_oddsrow_list():
    """Fresh fixture is a per-series dict; iterate each series through the
    parser and assert we get a non-trivial row count with non-zero bids."""
    raw = json.loads((FIXTURES / "kalshi" / "markets_response.json").read_text())
    series_responses = raw.get("series_responses") or {}
    assert series_responses, "kalshi fixture missing series_responses -- recapture needed"

    all_rows = []
    for ticker, resp in series_responses.items():
        all_rows.extend(parse_markets_response(resp, series_ticker=ticker))

    # Parser skips zero-bid markets; fresh fixture should still yield >= 50 rows
    # across all series. If this drops to 0, the fixture needs re-capture.
    assert len(all_rows) >= 50, f"expected >= 50 Kalshi rows, got {len(all_rows)}"
    assert all(r.book == "kalshi" for r in all_rows)
    assert all(r.book_label.startswith("KXNFL") for r in all_rows)
    assert all(r.book_subject for r in all_rows)


def test_kalshi_parse_trades():
    """Parser should turn Kalshi /markets/trades response into TradeRows.

    Live Kalshi v2 ships `yes_price_dollars` (string) + `count_fp` (string)
    + `taker_side` + `created_time`; we normalize to int cents + int count.
    """
    from nfl_draft.scrapers.kalshi import parse_trades_response
    raw = json.loads((FIXTURES / "kalshi" / "trades_response.json").read_text())
    rows = parse_trades_response(raw)
    assert isinstance(rows, list)
    if rows:
        t = rows[0]
        assert t.trade_id and t.ticker
        assert t.price_cents is not None
        assert t.count is not None
        assert 0 <= t.price_cents <= 100
        assert t.count >= 1
        assert t.side in ("yes", "no")


def test_dk_parse_returns_oddsrow_list():
    """DK fixture has 14 subcategories and ~795 selections; smoke-check
    that the parser hits every structured market_type and returns >= 500 rows."""
    raw = json.loads((FIXTURES / "draftkings" / "draft_markets.json").read_text())
    rows = dk_parse(raw)
    assert isinstance(rows, list)
    assert len(rows) >= 500, f"expected >= 500 DK rows, got {len(rows)}"
    assert all(r.book == "draftkings" for r in rows)
    assert all(r.american_odds != 0 for r in rows)  # 0 is never a valid US odd
    assert all(isinstance(r.american_odds, int) for r in rows)
    # Sanity: every OddsRow has a non-empty label + subject.
    for r in rows:
        assert r.book_label, f"empty book_label on {r!r}"
        assert r.book_subject, f"empty book_subject on {r!r}"
    # Coverage: we saw all 3 structured market groups.
    groups = {r.market_group for r in rows}
    assert "pick_outright" in groups
    assert "first_at_position" in groups
    assert any(g.startswith("top_") for g in groups)


def test_dk_parse_handles_unicode_minus():
    """DK ships American odds with a Unicode minus ('\u2212') rather than ASCII '-'.
    The parser must convert those to negative ints."""
    raw = json.loads((FIXTURES / "draftkings" / "draft_markets.json").read_text())
    rows = dk_parse(raw)
    negatives = [r for r in rows if r.american_odds < 0]
    assert negatives, "fixture has no negative odds after parsing - minus sign not handled"


def test_fd_parse_returns_oddsrow_list():
    """FD fixture has 16 draft markets under tab 391 - should yield >= 200
    runner-level rows covering picks 1-10, top-N, and first-at-position."""
    raw = json.loads((FIXTURES / "fanduel" / "draft_markets.json").read_text())
    rows = fd_parse(raw)
    assert isinstance(rows, list)
    assert len(rows) >= 200, f"expected >= 200 FD rows, got {len(rows)}"
    assert all(r.book == "fanduel" for r in rows)
    # FD already ships signed ints; negative odds should round-trip as-is.
    assert any(r.american_odds < 0 for r in rows)
    assert any(r.american_odds > 0 for r in rows)
    groups = {r.market_group for r in rows}
    assert "pick_outright" in groups
    assert "first_at_position" in groups
    assert any(g.startswith("top_") for g in groups)


def test_bm_parse_returns_oddsrow_list():
    """BM fixture has 18 draft games (10 overall picks + 8 position markets)
    and should yield >= 200 line-level rows."""
    raw = json.loads((FIXTURES / "bookmaker" / "draft_markets.json").read_text())
    rows = bm_parse(raw)
    assert isinstance(rows, list)
    assert len(rows) >= 200, f"expected >= 200 BM rows, got {len(rows)}"
    assert all(r.book == "bookmaker" for r in rows)
    groups = {r.market_group for r in rows}
    assert "pick_outright" in groups
    assert "first_at_position" in groups
    # Sanity: BM ships signed int strings; parser must yield real ints.
    for r in rows:
        assert isinstance(r.american_odds, int)
        assert r.american_odds != 0


def test_wz_parse_returns_oddsrow_list():
    """WZ fixture has 32 active draft leagues. Parser must handle all three
    shapes (single-game multi-runner, multi-game per-runner, embedded-in-htm)
    and yield >= 400 rows covering picks, top-N, first-at-position, and props."""
    raw = json.loads((FIXTURES / "wagerzon" / "draft_markets.json").read_text())
    rows = wz_parse(raw)
    assert isinstance(rows, list)
    assert len(rows) >= 400, f"expected >= 400 WZ rows, got {len(rows)}"
    assert all(r.book == "wagerzon" for r in rows)
    groups = {r.market_group for r in rows}
    assert "pick_outright" in groups
    assert "first_at_position" in groups
    assert any(g.startswith("top_") for g in groups)
    # WZ seeds every market with '*ALL BETS ACTION*' placeholder lines that
    # have empty odds -- parser must drop them.
    assert not any(r.book_subject.startswith("*ALL BETS ACTION*") for r in rows)


def test_kalshi_fetch_draft_odds_walks_cursor(monkeypatch):
    """fetch_draft_odds must paginate: if page 1 returns a cursor, page 2
    must be requested and its rows merged into the output. Regression guard
    for the bug where only the first 100 markets were captured (KXNFLDRAFTPICK
    has ~649 open markets — the other ~549 were silently dropped).
    """
    from nfl_draft.scrapers import kalshi as kalshi_mod

    # Mock discover_draft_series to return one fake series.
    monkeypatch.setattr(
        kalshi_mod.legacy_fetcher,
        "discover_draft_series",
        lambda: [{"series_ticker": "KXFAKE", "title": "Fake"}],
    )
    # Mock the legacy side-effect call (keeps this a pure unit test, no
    # network and no DB).
    monkeypatch.setattr(
        kalshi_mod.legacy_fetcher,
        "fetch_markets_for_series",
        lambda ticker: [],
    )

    # Build two mock pages. Page 1 has a cursor, page 2 does not.
    # Each page has markets with valid yes_bid so parse_markets_response
    # yields rows.
    def mk_market(ticker, bid):
        return {
            "ticker": ticker,
            "yes_bid": bid,
            "yes_sub_title": f"Candidate {ticker}",
        }

    page1 = {
        "markets": [mk_market(f"KXFAKE-P1-{i}", 50) for i in range(3)],
        "cursor": "CURSOR_P2",
    }
    page2 = {
        "markets": [mk_market(f"KXFAKE-P2-{i}", 50) for i in range(2)],
        "cursor": None,
    }

    call_log = []

    def fake_public_request(path):
        call_log.append(path)
        if "cursor=CURSOR_P2" in path:
            return page2
        return page1

    monkeypatch.setattr(kalshi_mod, "public_request", fake_public_request)

    rows = kalshi_mod.fetch_draft_odds()

    # 3 rows from page 1 + 2 rows from page 2
    assert len(rows) == 5, f"expected 5 rows across 2 pages, got {len(rows)}"
    assert len(call_log) == 2, f"expected 2 API calls (2 pages), got {len(call_log)}"
    assert "cursor=" not in call_log[0], "page 1 should have no cursor"
    assert "cursor=CURSOR_P2" in call_log[1], "page 2 should carry cursor"


def test_h88_parse_returns_oddsrow_list():
    """H88 fixture has 13 markets (9 pick_outright + 3 first_at_position +
    1 prop 'Mr Irrelevant Position') and ~155 total runners. Parser must
    emit real signed ints for every contestant and classify markets by
    propDescription label."""
    raw = json.loads((FIXTURES / "hoop88" / "draft_markets.json").read_text())
    rows = h88_parse(raw)
    assert isinstance(rows, list)
    assert len(rows) >= 100, f"expected >= 100 H88 rows, got {len(rows)}"
    assert all(r.book == "hoop88" for r in rows)
    groups = {r.market_group for r in rows}
    assert "pick_outright" in groups
    assert "first_at_position" in groups
    # Sanity: H88 ships signed ints; parser must yield real non-zero ints.
    for r in rows:
        assert isinstance(r.american_odds, int)
        assert r.american_odds != 0
        assert r.book_label, f"empty book_label on {r!r}"
        assert r.book_subject, f"empty book_subject on {r!r}"
        # H88 pads ContestantName with trailing spaces -- parser must strip.
        assert r.book_subject == r.book_subject.strip()
