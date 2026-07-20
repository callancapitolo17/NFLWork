from unabated_edge.venues import kalshi


def test_best_yes_ask_from_fp_dollars(monkeypatch):
    """Live API: orderbook_fp.no_dollars as string dollars (verified 2026-06-28)."""
    monkeypatch.setattr(kalshi.auth_client, "api",
        lambda m, p, *a, **k: (200, {"orderbook_fp": {"no_dollars": [["0.3300", "1"], ["0.3500", "2"]]}}, {}))
    assert kalshi.best_yes_ask("X") == 0.65    # 1 - 0.35 (highest no bid)


def test_best_yes_ask_legacy_cents_fallback(monkeypatch):
    monkeypatch.setattr(kalshi.auth_client, "api",
        lambda m, p, *a, **k: (200, {"orderbook": {"no": [[30, 100], [28, 50]]}}, {}))
    assert kalshi.best_yes_ask("X") == 0.70    # 1 - 0.30


def test_best_yes_ask_empty_book(monkeypatch):
    monkeypatch.setattr(kalshi.auth_client, "api",
        lambda m, p, *a, **k: (200, {"orderbook_fp": {"no_dollars": []}}, {}))
    assert kalshi.best_yes_ask("X") is None


def test_best_no_ask_from_fp_dollars(monkeypatch):
    """best_no_ask = 1 - best yes bid (from yes_dollars). Under = buy NO on Over rung."""
    monkeypatch.setattr(kalshi.auth_client, "api",
        lambda m, p, *a, **k: (200, {"orderbook_fp": {"yes_dollars": [["0.6400", "1"], ["0.6500", "2"]]}}, {}))
    # best yes bid = 0.65 → no_ask = 1 - 0.65 = 0.35
    assert kalshi.best_no_ask("X") == round(1 - 0.65, 2)


def test_get_book_normalizes_and_sorts_best_first(monkeypatch):
    monkeypatch.setattr(kalshi.auth_client, "api",
        lambda m, p, *a, **k: (200, {"orderbook_fp": {
            "yes_dollars": [["0.6400", "10"], ["0.6500", "2"]],
            "no_dollars": [["0.3300", "1"], ["0.3500", "73972.6"]]}}, {}))
    book = kalshi.get_book("X")
    assert book["yes_bids"] == [(0.65, 2.0), (0.64, 10.0)]      # best first, floats
    assert book["no_bids"] == [(0.35, 73972.6), (0.33, 1.0)]
    assert kalshi.yes_ask_from_book(book) == 0.65               # 1 - best no bid
    assert kalshi.no_ask_from_book(book) == 0.35                # 1 - best yes bid


def test_get_book_legacy_cents_fallback(monkeypatch):
    monkeypatch.setattr(kalshi.auth_client, "api",
        lambda m, p, *a, **k: (200, {"orderbook": {"yes": [[64, 10]], "no": [[30, 100], [28, 50]]}}, {}))
    book = kalshi.get_book("X")
    assert book["no_bids"] == [(0.30, 100.0), (0.28, 50.0)]
    assert kalshi.yes_ask_from_book(book) == 0.70


def test_ask_from_empty_book_side_is_none():
    book = {"yes_bids": [], "no_bids": [(0.30, 5.0)]}
    assert kalshi.no_ask_from_book(book) is None    # no yes bids → can't price NO
    assert kalshi.yes_ask_from_book(book) == 0.70
    assert kalshi.yes_ask_from_book(None) is None


def test_recent_trades_parses_dollars_and_cents_variants(monkeypatch):
    """Trades tape: *_dollars/count_fp strings and legacy integer cents both
    normalize; rows without a trade_id are dropped."""
    payload = {"trades": [
        {"trade_id": "t1", "ticker": "X", "created_time": "2026-07-10T17:00:00Z",
         "yes_price_dollars": "0.4200", "count_fp": "12.5", "taker_side": "yes"},
        {"trade_id": "t2", "ticker": "X", "created_time": "2026-07-10T17:00:01Z",
         "yes_price": 58, "count": 3, "taker_side": "no"},
        {"ticker": "X"},                                         # no trade_id → dropped
    ]}
    monkeypatch.setattr(kalshi.auth_client, "api", lambda m, p, *a, **k: (200, payload, {}))
    rows = kalshi.recent_trades("X", 60)
    assert [r["trade_id"] for r in rows] == ["t1", "t2"]
    assert rows[0]["yes_price"] == 0.42 and rows[0]["count"] == 12.5
    assert rows[1]["yes_price"] == 0.58 and rows[1]["count"] == 3.0


def test_recent_trades_failure_returns_empty(monkeypatch):
    monkeypatch.setattr(kalshi.auth_client, "api", lambda m, p, *a, **k: (500, None, {}))
    assert kalshi.recent_trades("X", 60) == []
