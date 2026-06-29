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
