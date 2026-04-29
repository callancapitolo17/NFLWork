import os
from pathlib import Path


def test_config_loads_from_env_example(monkeypatch):
    """Loading config from .env.example yields all expected knobs with correct types."""
    env_path = Path(__file__).parent.parent.parent / "kalshi_mlb_rfq" / ".env.example"
    for line in env_path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        k, v = line.split("=", 1)
        if v.strip():
            monkeypatch.setenv(k.strip(), v.strip())

    import importlib
    import kalshi_mlb_rfq.config as config_mod
    importlib.reload(config_mod)

    assert config_mod.KALSHI_BASE_URL == "https://api.elections.kalshi.com/trade-api/v2"
    assert config_mod.MVE_COLLECTION_TICKER == "KXMVECROSSCATEGORY-R"
    assert config_mod.BANKROLL == 1000.0
    assert config_mod.KELLY_FRACTION == 0.25
    assert config_mod.MIN_EV_PCT == 0.05
    assert config_mod.MAX_LIVE_RFQS == 80
    assert config_mod.MAX_PREDICTION_STALENESS_SEC == 600
    assert config_mod.RFQ_REFRESH_SEC == 30
    assert config_mod.MIN_FILL_RATIO == 0.50
    assert config_mod.MAX_GAME_EXPOSURE_PCT == 0.10
