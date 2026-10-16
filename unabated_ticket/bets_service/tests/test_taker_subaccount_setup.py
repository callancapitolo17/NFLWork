"""The taker's one-time subaccount setup (taker_subaccount_setup.py) against a
scripted Kalshi: it makes or reuses the subaccount, funds it once, mints a key
locked to it with no transfer scope, saves the private key at 0600 without
ever printing it, and proves the lock the desk checks before arming. No network."""
import io
import stat
from datetime import datetime, timezone

import pytest

from unabated_ticket.bets_service import taker_subaccount_setup as setup

NOW = datetime(2026, 10, 16, 20, 0, tzinfo=timezone.utc)
NEW_KEY_ID = "new-key-0001"


def rsa_pem() -> str:
    from cryptography.hazmat.primitives import serialization
    from cryptography.hazmat.primitives.asymmetric import rsa
    key = rsa.generate_private_key(public_exponent=65537, key_size=2048)
    return key.private_bytes(serialization.Encoding.PEM, serialization.PrivateFormat.PKCS8,
                             serialization.NoEncryption()).decode()


PEM = rsa_pem()


class MainKalshi:
    """Cal's main key: balances per subaccount, and scripted answers to the
    create / upgrade / transfer / generate POSTs. Records every call."""

    def __init__(self, balances, create=(201,), upgrade=201, transfer=200, generate=201, generated=None):
        self.balances = dict(balances)
        self.create = list(create)
        self.upgrade = upgrade
        self.transfer = transfer
        self.generate = generate
        self.generated = generated if generated is not None else {"api_key_id": NEW_KEY_ID, "key_type": "rsa",
                                                                    "private_key": PEM}
        self.calls = []

    def request(self, method, path, body=None):
        self.calls.append((method, path, body))
        if (method, path) == ("GET", setup.SUBACCOUNT_BALANCES_PATH):
            return 200, {"subaccount_balances": [
                {"subaccount_number": number, "exchange_index": 0, "balance": f"{dollars:.4f}", "updated_ts": 1}
                for number, dollars in self.balances.items()]}
        if (method, path) == ("POST", setup.SUBACCOUNTS_PATH):
            status = self.create.pop(0)
            if status != 201:
                return status, {"error": {"code": "forbidden", "message": "advanced tier required"}}
            number = max(self.balances) + 1
            self.balances[number] = 0.0
            return 201, {"subaccount_number": number}
        if (method, path) == ("POST", setup.UPGRADE_PATH):
            return self.upgrade, {} if self.upgrade == 201 else {"error": "no API order in the last 100"}
        if (method, path) == ("POST", setup.TRANSFER_PATH):
            if self.transfer == 200:
                self.balances[body["to_subaccount"]] += body["amount_cents"] / 100
                self.balances[body["from_subaccount"]] -= body["amount_cents"] / 100
                return 200, {}
            return self.transfer, {"error": "insufficient balance"}
        if (method, path) == ("POST", setup.GENERATE_KEY_PATH):
            return (201, self.generated) if self.generate == 201 else (self.generate, {"error": "scope"})
        raise AssertionError(f"unexpected {method} {path}")

    def posts(self, path):
        return [body for method, called, body in self.calls if method == "POST" and called == path]


class BotKalshi:
    """A key's view of GET /portfolio/balance: locked = refused the primary."""

    def __init__(self, locked=True):
        self.locked = locked

    def request(self, method, path, body=None):
        if path.endswith("subaccount=0") and self.locked:
            return 403, {"error": "subaccount must match the key's restriction"}
        if path.endswith("subaccount=0"):
            return 200, {"balance": 1}
        return 200, {"balance": 50000}


def run(main, tmp_path, wanted=None, configured_key_id=None, configured_subaccount=None, locked=True, key_path=None):
    out = io.StringIO()
    made = []

    def make_bot_client(key_id, path):
        made.append((key_id, path))
        return BotKalshi(locked)
    key_path = key_path or tmp_path / "kalshi_taker.pem"
    code = setup.run_setup(main, make_bot_client, key_path, wanted, configured_key_id, configured_subaccount,
                           out, NOW, sleep=lambda _sec: None)
    return code, out.getvalue(), made


def test_a_fresh_account_gets_a_subaccount_500_and_a_locked_key_saved_unprinted(tmp_path):
    main = MainKalshi({0: 2000.0})
    code, output, made = run(main, tmp_path)
    assert code == 0, output
    assert main.balances == {0: 1500.0, 1: 500.0}
    [transfer] = main.posts(setup.TRANSFER_PATH)
    assert transfer == {"client_transfer_id": setup.funding_transfer_id(1), "from_subaccount": 0, "to_subaccount": 1,
                        "amount_cents": 50000}
    [generate] = main.posts(setup.GENERATE_KEY_PATH)
    assert generate == {"name": "unabated-taker-subaccount-1-20261016", "key_type": "rsa",
                        "scopes": ["read", "write::trade"], "subaccount": 1}
    key_file = tmp_path / "kalshi_taker.pem"
    assert key_file.read_text() == PEM
    assert stat.S_IMODE(key_file.stat().st_mode) == 0o600
    assert made == [(NEW_KEY_ID, key_file)]
    assert "KALSHI_TAKER_API_KEY_ID=new-key-0001\n" in output and "KALSHI_TAKER_SUBACCOUNT=1\n" in output
    assert "PRIVATE KEY" not in output and PEM.splitlines()[1] not in output


def test_the_funding_transfer_id_is_fixed_per_subaccount():
    assert setup.funding_transfer_id(1) == setup.funding_transfer_id(1)
    assert setup.funding_transfer_id(1) != setup.funding_transfer_id(2)


def test_a_refused_create_asks_for_the_advanced_tier_once_then_retries(tmp_path):
    main = MainKalshi({0: 2000.0}, create=(403, 201))
    code, output, _made = run(main, tmp_path)
    assert code == 0, output
    assert [call[1] for call in main.calls if call[0] == "POST"][:3] == [
        setup.SUBACCOUNTS_PATH, setup.UPGRADE_PATH, setup.SUBACCOUNTS_PATH]


def test_a_refused_tier_stops_before_any_money_moves(tmp_path):
    main = MainKalshi({0: 2000.0}, create=(403,), upgrade=403)
    code, output, _made = run(main, tmp_path)
    assert code == 1 and output.splitlines()[-1].startswith("STOPPED: the Advanced API tier was refused: HTTP 403")
    assert main.posts(setup.TRANSFER_PATH) == [] and main.posts(setup.GENERATE_KEY_PATH) == []


def test_the_one_existing_subaccount_is_reused_and_a_funded_one_gets_no_transfer(tmp_path):
    main = MainKalshi({0: 2000.0, 1: 480.0})
    code, output, _made = run(main, tmp_path)
    assert code == 0, output
    assert main.posts(setup.SUBACCOUNTS_PATH) == [] and main.posts(setup.TRANSFER_PATH) == []
    assert "already holds $480.00; no transfer" in output


def test_several_subaccounts_need_one_named_and_an_unlisted_name_stops(tmp_path):
    code, output, _made = run(MainKalshi({0: 10.0, 1: 0.0, 2: 0.0}), tmp_path)
    assert code == 1 and "name the bot's with --subaccount N" in output
    code, output, _made = run(MainKalshi({0: 10.0, 1: 0.0}), tmp_path, wanted=5)
    assert code == 1 and "Kalshi lists no subaccount 5" in output
    main = MainKalshi({0: 2000.0, 1: 0.0, 2: 0.0})
    code, output, _made = run(main, tmp_path, wanted=2)
    assert code == 0 and main.posts(setup.GENERATE_KEY_PATH)[0]["subaccount"] == 2


def test_a_refused_transfer_stops_before_a_key_is_made(tmp_path):
    main = MainKalshi({0: 100.0}, transfer=400)
    code, output, _made = run(main, tmp_path)
    assert code == 1 and "the $500 transfer to subaccount 1 was refused: HTTP 400" in output
    assert main.posts(setup.GENERATE_KEY_PATH) == []


def test_a_finished_setup_makes_no_second_key(tmp_path):
    key_file = tmp_path / "kalshi_taker.pem"
    key_file.write_text(PEM)
    main = MainKalshi({0: 1500.0, 1: 500.0})
    code, output, made = run(main, tmp_path, configured_key_id="key-live", configured_subaccount=1)
    assert code == 0 and "already locked to subaccount 1; no new key" in output
    assert main.posts(setup.GENERATE_KEY_PATH) == [] and made == [("key-live", key_file)]


def test_a_new_key_kalshi_did_not_lock_is_deleted_and_named_for_revoking(tmp_path):
    main = MainKalshi({0: 2000.0})
    code, output, _made = run(main, tmp_path, locked=False)
    assert code == 1
    assert "did not pass the lock check" in output and "Revoke new-key-0001 at kalshi.com" in output
    assert not (tmp_path / "kalshi_taker.pem").exists()


@pytest.mark.parametrize("generated", [{"private_key": PEM}, {"api_key_id": NEW_KEY_ID, "private_key": 7}])
def test_an_odd_generate_reply_never_echoes_the_private_key(tmp_path, generated):
    main = MainKalshi({0: 2000.0}, generated=generated)
    code, output, _made = run(main, tmp_path)
    assert code == 1 and "expected api_key_id and private_key" in output
    assert "PRIVATE KEY" not in output and PEM.splitlines()[1] not in output


def test_an_unusable_private_key_is_never_written(tmp_path):
    main = MainKalshi({0: 2000.0}, generated={"api_key_id": NEW_KEY_ID, "private_key": "not a key"})
    code, output, _made = run(main, tmp_path)
    assert code == 1 and "could not be saved" in output and "revoke new-key-0001 at kalshi.com" in output
    assert not (tmp_path / "kalshi_taker.pem").exists()
    assert "not a key" not in output
