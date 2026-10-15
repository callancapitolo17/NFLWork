"""Unit tests for scraper_betonline.py — reading BetOnline bets from the bets service.

The round trip pins the claim the sheet relies on: the report fixture, stored
the way the VM stores it (unabated_ticket/bets_service/sources/betonline.py),
then read back through /bets.json, gives the same sheet rows as parsing the
report directly. No network: requests.get is replaced.
"""
import json
import sys
from pathlib import Path

import pytest

import scraper_betonline
from scraper_betonline import parse_api_bets, report_rows_from_payload, scrape_betonline

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))
from unabated_ticket.bets_service.sources.betonline import normalize_betonline  # noqa: E402

FIXTURE = REPO_ROOT / "unabated_ticket" / "tests" / "fixtures" / "bets" / "betonline_history.json"


def load_report_rows() -> list:
    data = json.loads(FIXTURE.read_text())
    return data if isinstance(data, list) else data["Data"]


def payload_of(report_rows: list, ok: bool = True) -> dict:
    """A /bets.json body holding the rows as the VM stores them."""
    return {
        "sources": {"betonline": {"ok": ok, "error": None if ok else "token refresh failed",
                                  "fetchedAt": "2026-09-14T00:00:00Z", "count": len(report_rows)}},
        "bets": normalize_betonline(report_rows, None)
        + [{"id": "kalshi:x", "venue": "kalshi", "raw": {"nativeId": "x"}}],
    }


class FakeResponse:
    def __init__(self, body: dict, status_code: int = 200):
        self._body = body
        self.status_code = status_code
        self.text = json.dumps(body)

    def json(self) -> dict:
        return self._body


def test_round_trip_gives_the_same_sheet_rows_as_the_report():
    report_rows = load_report_rows()
    from_date = "2026-09-01"
    expected = parse_api_bets(sorted((row for row in report_rows if row["Date"][:10] >= from_date),
                                     key=lambda row: row["Date"]))
    actual = parse_api_bets(report_rows_from_payload(payload_of(report_rows), from_date))
    assert actual == expected
    assert expected, "fixture should hold settled bets after from_date"


def test_parlay_legs_collapse_to_one_wager():
    report_rows = [row for row in load_report_rows() if row["WagerType"] == "Same Game Parlay"]
    payload = payload_of(report_rows)
    assert len(payload["bets"]) > len(report_rows) + 1, "fixture parlay should be stored as legs"
    assert len(report_rows_from_payload(payload, "2000-01-01")) == len(report_rows)


def test_rows_placed_before_from_date_are_left_out():
    rows = report_rows_from_payload(payload_of(load_report_rows()), "2026-09-12")
    assert rows and all(row["Date"][:10] >= "2026-09-12" for row in rows)


def test_scrape_reads_the_service_with_a_window_covering_from_date(monkeypatch):
    calls = []

    def fake_get(url, params, timeout):
        calls.append((url, params))
        return FakeResponse(payload_of(load_report_rows()))

    monkeypatch.setattr(scraper_betonline.requests, "get", fake_get)
    bets = scrape_betonline(from_date="2026-09-01", service_url="https://vm.example/")
    assert calls[0][0] == "https://vm.example/bets.json"
    assert calls[0][1]["days"] >= 2
    assert all(bet["platform"] == "BetOnline" and bet["result"] for bet in bets)


def test_failed_vm_poll_fails_loudly(monkeypatch):
    monkeypatch.setattr(scraper_betonline.requests, "get",
                        lambda url, params, timeout: FakeResponse(payload_of(load_report_rows(), ok=False)))
    with pytest.raises(RuntimeError, match="last BetOnline poll failed"):
        scrape_betonline(days_back=7, service_url="https://vm.example")


def test_missing_source_and_http_error_fail_loudly(monkeypatch):
    monkeypatch.setattr(scraper_betonline.requests, "get",
                        lambda url, params, timeout: FakeResponse({"sources": {}, "bets": []}))
    with pytest.raises(RuntimeError, match="no 'betonline' source"):
        scrape_betonline(service_url="https://vm.example")
    monkeypatch.setattr(scraper_betonline.requests, "get",
                        lambda url, params, timeout: FakeResponse({}, status_code=403))
    with pytest.raises(RuntimeError, match="HTTP 403"):
        scrape_betonline(service_url="https://vm.example")
