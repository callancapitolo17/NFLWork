"""check_sources: per-venue verdicts without touching any venue or the network."""
import io

from unabated_ticket.bets_service import check_sources


class FakeSource:
    name = "fake"
    poll_sec = 60.0

    def __init__(self, records=None, error=None):
        self._records = records or []
        self._error = error

    def fetch(self):
        if self._error:
            raise self._error
        return self._records


def fixed_clock():
    ticks = iter([10.0, 12.5])
    return lambda: next(ticks)


def test_ok_reports_count_and_seconds_never_contents():
    result = check_sources.check_source("bfa", lambda: FakeSource(records=[{"id": "bfa:1", "stake": 50}] * 3),
                                        clock=fixed_clock())
    assert result == check_sources.CheckResult("bfa", "ok", "3 records in 2.5s")


def test_missing_credentials_is_not_a_failure():
    result = check_sources.check_source("wagerzon", lambda: None)
    assert result.status == "not configured"


def test_fetch_error_is_a_failure_with_its_type():
    result = check_sources.check_source("betonline", lambda: FakeSource(error=PermissionError("403 Forbidden")),
                                        clock=fixed_clock())
    assert result.status == "FAILED"
    assert result.detail == "PermissionError: 403 Forbidden after 2.5s"


def test_factory_that_raises_is_a_setup_failure():
    def broken():
        raise RuntimeError("KALSHI_API_KEY_ID / KALSHI_PRIVATE_KEY_PATH not set")
    result = check_sources.check_source("kalshi", broken)
    assert result.status == "FAILED"
    assert result.detail.startswith("setup: RuntimeError: KALSHI_API_KEY_ID")


def test_long_errors_are_truncated():
    result = check_sources.check_source("novig", lambda: FakeSource(error=ValueError("x" * 1000)),
                                        clock=fixed_clock())
    assert len(result.detail) < check_sources.MAX_ERROR_CHARS + 50


class FakeResponse(io.BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def test_unabated_feed_ok_and_block_page():
    def big(url, timeout):
        return FakeResponse(b"{" + b" " * 5000 + b"}")

    def tiny(url, timeout):
        return FakeResponse(b"<html>blocked</html>")

    assert check_sources.check_unabated_feed(big).status == "ok"
    blocked = check_sources.check_unabated_feed(tiny)
    assert blocked.status == "FAILED" and "got 20" in blocked.detail


def test_every_service_venue_is_checked(monkeypatch):
    monkeypatch.setattr(check_sources, "KalshiSource", lambda: None)
    monkeypatch.setattr(check_sources, "check_unabated_feed",
                        lambda: check_sources.CheckResult("unabated feed", "ok", "stub"))
    monkeypatch.setattr(check_sources.service, "OPTIONAL_SOURCE_FACTORIES",
                        tuple((venue, lambda: None) for venue, _ in check_sources.service.OPTIONAL_SOURCE_FACTORIES))
    venues = [result.venue for result in check_sources.run_checks()]
    assert venues == ["kalshi", "betonline", "novig", "bfa", "wagerzon", "polymarket_us", "unabated feed"]


def test_table_aligns_columns():
    table = check_sources.format_results([
        check_sources.CheckResult("kalshi", "ok", "4 records in 1.0s"),
        check_sources.CheckResult("polymarket_us", "FAILED", "boom"),
    ])
    lines = table.splitlines()
    assert lines[0].index("ok") == lines[1].index("FAILED")
