"""Issue #50 item 5: Caesars WAF-token pre-warm.

``ensure_token`` refreshes only at expiry, so with a 240s TTL an RFQ that
lands just after expiry pays the ~1.5s node mint inline. The warming pass
instead calls ``ensure_fresh_token(min_remaining_sec=...)`` — re-mint in the
BACKGROUND whenever the remaining TTL would not cover the next warming gap,
so the mint cost never lands on the quote path.

No network, no node: ``_mint`` / ``_load_cache`` / ``_validate`` are stubbed;
the clock is monkeypatched at the module's ``time.time`` seam.
"""
import inspect

import pytest

from mlb_sgp import caesars_client
from mlb_sgp.caesars_client import TOKEN_TTL_SEC, CaesarsClient


@pytest.fixture
def client(monkeypatch, tmp_path):
    """A CaesarsClient with a controllable clock and no disk/wire access."""
    monkeypatch.setattr(caesars_client, "TOKEN_CACHE",
                        tmp_path / "token.json")
    clock = {"t": 1_000_000.0}
    monkeypatch.setattr(caesars_client.time, "time", lambda: clock["t"])
    c = CaesarsClient()
    c._clock = clock            # test handle, not production state
    minted = []

    def fake_mint():
        minted.append(clock["t"])
        c._token = "fresh-token"
        c._minted_at = clock["t"]
        return True
    monkeypatch.setattr(c, "_mint", fake_mint)
    monkeypatch.setattr(c, "_load_cache", lambda: False)
    c._mint_log = minted
    return c


def _age_token(c, age_sec):
    c._token = "aging-token"
    c._minted_at = c._clock["t"] - age_sec


def test_ensure_fresh_token_min_remaining_is_keyword_only():
    sig = inspect.signature(CaesarsClient.ensure_fresh_token)
    assert (sig.parameters["min_remaining_sec"].kind
            is inspect.Parameter.KEYWORD_ONLY)
    with pytest.raises(TypeError):
        sig.bind(object(), 120.0)


def test_fresh_enough_token_is_left_alone(client):
    _age_token(client, age_sec=30)             # 210s remaining of 240
    assert client.ensure_fresh_token(min_remaining_sec=120.0) is True
    assert client._mint_log == []              # no mint, no wire


def test_near_expiry_token_is_reminted_before_ensure_token_would(client):
    """The pre-warm case: still VALID (ensure_token would keep it), but not
    enough runway to cover the next warming gap — mint now, in background."""
    _age_token(client, age_sec=TOKEN_TTL_SEC - 60)     # 60s remaining
    assert client.ensure_token() is True               # old path: still fine
    assert client._mint_log == []
    assert client.ensure_fresh_token(min_remaining_sec=120.0) is True
    assert len(client._mint_log) == 1                  # pre-warm minted
    assert client._token == "fresh-token"


def test_expired_token_is_minted(client):
    _age_token(client, age_sec=TOKEN_TTL_SEC + 10)
    assert client.ensure_fresh_token(min_remaining_sec=120.0) is True
    assert len(client._mint_log) == 1


def test_failed_mint_returns_false_without_raising(client, monkeypatch):
    _age_token(client, age_sec=TOKEN_TTL_SEC - 10)
    monkeypatch.setattr(client, "_mint", lambda: False)
    assert client.ensure_fresh_token(min_remaining_sec=120.0) is False


def test_fresh_cache_from_another_process_avoids_a_mint(client, monkeypatch):
    """Two processes share the on-disk token cache; a sibling's fresh mint
    must satisfy the pre-warm without a duplicate node run."""
    _age_token(client, age_sec=TOKEN_TTL_SEC - 60)

    def fake_load_cache():
        client._token = "sibling-token"
        client._minted_at = client._clock["t"] - 5     # nearly new
        return True
    monkeypatch.setattr(client, "_load_cache", fake_load_cache)
    assert client.ensure_fresh_token(min_remaining_sec=120.0) is True
    assert client._mint_log == []
    assert client._token == "sibling-token"
