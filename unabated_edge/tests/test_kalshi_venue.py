from unabated_edge.venues import kalshi


def test_best_yes_ask_from_no_bids(monkeypatch):
    monkeypatch.setattr(kalshi.auth_client, "api",
        lambda m,p,*a,**k: (200, {"orderbook":{"no":[[30,100],[28,50]]}}, {}))
    assert kalshi.best_yes_ask("X") == 0.70    # 1 - 0.30
