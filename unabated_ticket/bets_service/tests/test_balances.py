"""Balance by book: the per-venue balance parsers (fields verified on Cal's
account 2026-10-13), the in-memory BalanceBook, and GET /balances.json."""
import http.client
import json
import threading
from http.server import ThreadingHTTPServer

import pytest

from unabated_ticket.bets_service import service
from unabated_ticket.bets_service.balances import BalanceBook, read_balance_once, validate_reading
from unabated_ticket.bets_service.sources.kalshi import KalshiSource, balance_of_kalshi
from unabated_ticket.bets_service.sources.novig import NovigSource, balance_of_novig_summary
from unabated_ticket.bets_service.sources.polymarket_us import balance_of_polymarket_us
from unabated_ticket.bets_service.store import BetsStore


def test_kalshi_balance_is_cents_to_dollars():
    assert balance_of_kalshi({"balance": 826312, "portfolio_value": 109050}) == {
        "cash": 8263.12, "inBets": 1090.5, "available": 8263.12, "creditLimit": None}


def test_kalshi_balance_refuses_a_missing_field():
    with pytest.raises(RuntimeError, match="portfolio_value"):
        balance_of_kalshi({"balance": 100})


def test_kalshi_source_balance_raises_on_non_200():
    source = KalshiSource(api=lambda method, path: (401, {"error": "no"}, {}), configure_auth=False)
    with pytest.raises(RuntimeError, match="401"):
        source.balance()


def test_novig_balance_reads_cash_wallet_and_open_positions():
    body = {"wallets": {"cashBalance": 320.5, "coinBalance": 9000, "bonusBalance": 0},
            "openPositionsValue": 19857.25, "totalValue": 20177.75, "tradeCreditsValue": 0}
    assert balance_of_novig_summary(body) == {
        "cash": 320.5, "inBets": 19857.25, "available": 320.5, "creditLimit": None}


def test_novig_balance_refuses_a_missing_field():
    with pytest.raises(RuntimeError, match="openPositionsValue"):
        balance_of_novig_summary({"wallets": {"cashBalance": 1}})


def test_novig_source_balance_asks_for_the_cash_summary():
    class Token:
        token = "t"

    class Auth:
        def token(self):
            return Token()

    urls = []

    def http_get(url, bearer):
        urls.append(url)
        return {"wallets": {"cashBalance": 1}, "openPositionsValue": 2}

    source = NovigSource(auth=Auth(), http_get=http_get, poll_sec=1, base_url="https://x/nbx/v1")
    assert source.balance()["inBets"] == 2
    assert urls == ["https://x/nbx/v1/portfolio/summary?currency=CASH"]


def test_polymarket_us_balance_reads_the_one_entry():
    body = {"balances": [{"currentBalance": 1082.4, "buyingPower": 1080.0, "assetAvailable": 0,
                          "pendingCredit": 0, "displayedCash": 1082.4, "availableToWithdraw": 1082.4}]}
    assert balance_of_polymarket_us(body) == {
        "cash": 1082.4, "inBets": None, "available": 1080.0, "creditLimit": None}


def test_polymarket_us_balance_refuses_more_than_one_entry():
    entry = {"currentBalance": 1, "buyingPower": 1}
    with pytest.raises(RuntimeError, match="exactly one"):
        balance_of_polymarket_us({"balances": [entry, entry]})


def test_validate_reading_needs_cash():
    with pytest.raises(RuntimeError, match="no cash"):
        validate_reading({"cash": None, "inBets": 5})


class FakeSource:
    name = "fakebook"

    def __init__(self):
        self.answers = []

    def balance(self):
        answer = self.answers.pop(0)
        if isinstance(answer, Exception):
            raise answer
        return answer


def test_a_failed_read_keeps_the_last_numbers_and_says_why():
    book = BalanceBook()
    source = FakeSource()
    source.answers = [{"cash": 100, "inBets": 20, "available": 100}, RuntimeError("login refused")]
    read_balance_once(source, book)
    first = book.snapshot()["fakebook"]
    read_balance_once(source, book)
    held = book.snapshot()["fakebook"]
    assert held["cash"] == 100.0 and held["fetchedAt"] == first["fetchedAt"]
    assert held["error"] == "RuntimeError: login refused"


def test_a_venue_never_read_has_no_cash():
    book = BalanceBook()
    source = FakeSource()
    source.answers = [RuntimeError("login refused")]
    read_balance_once(source, book)
    assert book.snapshot()["fakebook"]["cash"] is None


def test_a_source_without_balance_is_skipped():
    book = BalanceBook()

    class NoBalance:
        name = "nobalance"

    read_balance_once(NoBalance(), book)
    assert book.snapshot() == {}


def test_balances_route(tmp_path):
    book = BalanceBook()
    book.record_success("kalshi", {"cash": 10.0, "inBets": 1.0, "available": 10.0, "creditLimit": None})
    store = BetsStore(tmp_path / "bets.duckdb", 30)
    handler = service.make_handler(store, 0.0, ["kalshi"], extra_hosts=(), balance_book=book)
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    try:
        port = server.server_address[1]
        connection = http.client.HTTPConnection("127.0.0.1", port, timeout=5)
        connection.request("GET", "/balances.json", headers={"Host": f"127.0.0.1:{port}"})
        response = connection.getresponse()
        body = json.loads(response.read())
    finally:
        server.shutdown()
        store.close()
    assert response.status == 200
    assert body["balances"]["kalshi"]["cash"] == 10.0
