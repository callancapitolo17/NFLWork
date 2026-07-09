from unabated_edge import config

def test_core_defaults():
    assert config.SHARP_BOOK_PRICE_ID == 7
    assert {7, 6, 68}.issubset(set(config.ANCHOR_SOURCE_IDS))
    assert config.BANKROLL == 1000.0
    assert config.UNABATED_SNAPSHOT_URL.startswith("https://content.unabated.com")
    assert not hasattr(config, "WC_LEAGUE_KEY_PREFIX")   # sport specifics live in adapters
