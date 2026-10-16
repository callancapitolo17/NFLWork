"""The taker's one-time subaccount setup (taker_subaccount_setup.py) against a
scripted Kalshi: it puts the stop file in place, makes or reuses the
subaccount, funds it once and checks the money landed, mints a key locked to
it with no transfer scope, saves the private key at 0600 without ever
printing it, and proves the lock two ways (Kalshi's key listing, then the
desk's own balance check). No network."""
import io
import os
import stat
from datetime import datetime, timezone

import pytest

from unabated_ticket.bets_service import taker_subaccount_setup as setup

NOW = datetime(2026, 10, 16, 20, 0, tzinfo=timezone.utc)
NEW_KEY_ID = "new-key-0001"
LOCKED_RECORD = {"api_key_id": NEW_KEY_ID, "name": "bot", "scopes": ["read", "write::trade"], "subaccount": 1}


def rsa_pem() -> str:
    from cryptography.hazmat.primitives import serialization
    from cryptography.hazmat.primitives.asymmetric import rsa
    key = rsa.generate_private_key(public_exponent=65537, key_size=2048)
    return key.private_bytes(serialization.Encoding.PEM, serialization.PrivateFormat.PKCS8,
                             serialization.NoEncryption()).decode()


PEM = rsa_pem()


class MainKalshi:
    """Cal's main key: balances per subaccount, scripted answers to the create /
    upgrade / transfer / generate POSTs, and the key listing. Records every
    call. `listed` is what GET /api_keys carries (None: the minted key, locked)."""

    def __init__(self, balances, create=(201,), upgrade=201, transfer=200, generate=201, generated=None,
                 listed=None, transfer_lands=True, land_after_reads=0, create_error="forbidden"):
        self.balances = dict(balances)
        self.create = list(create)
        self.upgrade = upgrade
        self.transfer = transfer
        self.transfer_lands = transfer_lands
        # Balance reads after an accepted transfer that still show it missing (Kalshi's lag).
        self.land_after_reads = land_after_reads
        self.pending = {}
        self.create_error = create_error
        self.generate = generate
        self.generated = generated if generated is not None else {"api_key_id": NEW_KEY_ID, "key_type": "rsa",
                                                                    "private_key": PEM}
        self.listed = listed
        self.calls = []

    def request(self, method, path, body=None):
        self.calls.append((method, path, body))
        if (method, path) == ("GET", setup.SUBACCOUNT_BALANCES_PATH):
            if self.pending:
                if self.land_after_reads > 0:
                    self.land_after_reads -= 1
                else:
                    for number, dollars in self.pending.items():
                        self.balances[number] += dollars
                    self.pending = {}
            return 200, {"subaccount_balances": [
                {"subaccount_number": number, "exchange_index": 0, "balance": f"{dollars:.4f}", "updated_ts": 1}
                for number, dollars in self.balances.items()]}
        if (method, path) == ("GET", setup.API_KEYS_PATH):
            listed = self.listed if self.listed is not None else [dict(LOCKED_RECORD, subaccount=self.minted_in())]
            return 200, {"api_keys": listed}
        if (method, path) == ("POST", setup.SUBACCOUNTS_PATH):
            status = self.create.pop(0)
            if status != 201:
                return status, {"error": {"code": self.create_error, "message": "refused"}}
            number = max(self.balances) + 1
            self.balances[number] = 0.0
            return 201, {"subaccount_number": number}
        if (method, path) == ("POST", setup.UPGRADE_PATH):
            return self.upgrade, {} if self.upgrade == 201 else {"error": "no API order in the last 100"}
        if (method, path) == ("POST", setup.TRANSFER_PATH):
            if self.transfer == 200:
                if self.transfer_lands:
                    self.pending[body["to_subaccount"]] = body["amount_cents"] / 100
                    self.balances[body["from_subaccount"]] -= body["amount_cents"] / 100
                return 200, {}
            return self.transfer, {"error": "insufficient balance"}
        if (method, path) == ("POST", setup.GENERATE_KEY_PATH):
            return (201, self.generated) if self.generate == 201 else (self.generate, {"error": "scope"})
        raise AssertionError(f"unexpected {method} {path}")

    def posts(self, path):
        return [body for method, called, body in self.calls if method == "POST" and called == path]

    def minted_in(self):
        generated = self.posts(setup.GENERATE_KEY_PATH)
        return generated[-1]["subaccount"] if generated else 1


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


def run(main, tmp_path, wanted=None, configured_key_id=None, configured_subaccount=None, locked=True, key_path=None,
        stop_path=None):
    out = io.StringIO()
    made = []

    def make_bot_client(key_id, path):
        made.append((key_id, path))
        return BotKalshi(locked)
    key_path = key_path or tmp_path / "kalshi_taker.pem"
    stop_path = stop_path or tmp_path / "taker.stop"
    code = setup.run_setup(main, make_bot_client, key_path, stop_path, wanted, configured_key_id,
                           configured_subaccount, out, NOW, sleep=lambda _sec: None)
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
    # Kalshi's listing of the key was read on the main key before the bot key was tried.
    assert ("GET", setup.API_KEYS_PATH, None) in main.calls
    assert not (tmp_path / "kalshi_taker.pem.new").exists() and not (tmp_path / "kalshi_taker.pem.probe").exists()


def test_the_stop_file_is_put_in_place_first_and_left_alone_when_present(tmp_path):
    main = MainKalshi({0: 2000.0})
    code, output, _made = run(main, tmp_path)
    assert code == 0 and (tmp_path / "taker.stop").exists()
    assert output.startswith(f"Put the stop file {tmp_path / 'taker.stop'} in place")
    code, output, _made = run(MainKalshi({0: 1500.0, 1: 500.0}), tmp_path, configured_key_id=NEW_KEY_ID,
                              configured_subaccount=1)
    assert code == 0 and output.startswith(f"Stop file {tmp_path / 'taker.stop'} is in place")


def test_an_unwritable_key_directory_stops_before_kalshi_is_asked_anything(tmp_path):
    if os.geteuid() == 0:
        pytest.skip("root writes anywhere")
    locked_dir = tmp_path / "locked"
    locked_dir.mkdir()
    locked_dir.chmod(0o500)
    main = MainKalshi({0: 2000.0})
    try:
        code, output, _made = run(main, tmp_path, key_path=locked_dir / "kalshi_taker.pem")
    finally:
        locked_dir.chmod(0o700)
    assert code == 1 and "cannot write the bot key file" in output
    assert main.calls == []


def test_the_funding_transfer_id_is_fixed_per_subaccount():
    assert setup.funding_transfer_id(1) == setup.funding_transfer_id(1)
    assert setup.funding_transfer_id(1) != setup.funding_transfer_id(2)


def test_a_refused_create_asks_for_the_advanced_tier_once_then_retries(tmp_path):
    main = MainKalshi({0: 2000.0}, create=(403, 201))
    code, output, _made = run(main, tmp_path)
    assert code == 0, output
    assert [call[1] for call in main.calls if call[0] == "POST"][:3] == [
        setup.SUBACCOUNTS_PATH, setup.UPGRADE_PATH, setup.SUBACCOUNTS_PATH]
    assert main.posts(setup.UPGRADE_PATH) == [{}]


def test_a_refused_tier_stops_before_any_money_moves(tmp_path):
    main = MainKalshi({0: 2000.0}, create=(403,), upgrade=403)
    code, output, _made = run(main, tmp_path)
    assert code == 1 and output.splitlines()[-1].startswith("STOPPED: the Advanced API tier was refused: HTTP 403")
    assert main.posts(setup.TRANSFER_PATH) == [] and main.posts(setup.GENERATE_KEY_PATH) == []


def test_a_read_only_signing_key_stops_before_the_tier_is_asked_or_anything_changes(tmp_path):
    main = MainKalshi({0: 2000.0}, create=(403,), create_error="insufficient_scope:_write_required")
    code, output, _made = run(main, tmp_path)
    assert code == 1
    assert output.splitlines()[-1].startswith("STOPPED: the signing key (KALSHI_API_KEY_ID) is read-only")
    assert "deploy/README step 11" in output and "Advanced API tier" not in output
    assert (tmp_path / "taker.stop").exists()  # step 0 still ran, so the run says "at Kalshi"
    assert [call[1] for call in main.calls if call[0] == "POST"] == [setup.SUBACCOUNTS_PATH]


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
    assert code == 0, output
    assert main.posts(setup.GENERATE_KEY_PATH)[0]["subaccount"] == 2


def test_a_refused_transfer_stops_before_a_key_is_made(tmp_path):
    main = MainKalshi({0: 100.0}, transfer=400)
    code, output, _made = run(main, tmp_path)
    assert code == 1 and "the $500 transfer to subaccount 1 was refused: HTTP 400" in output
    assert main.posts(setup.GENERATE_KEY_PATH) == []


def test_a_transfer_kalshi_accepted_but_did_not_land_stops_before_a_key_is_made(tmp_path):
    main = MainKalshi({0: 2000.0}, transfer_lands=False)
    code, output, _made = run(main, tmp_path)
    assert code == 1 and "answered 200 to the $500 transfer but subaccount 1 holds $0.00" in output
    assert main.posts(setup.GENERATE_KEY_PATH) == []


def test_a_transfer_kalshi_shows_a_few_seconds_late_is_waited_for(tmp_path):
    main = MainKalshi({0: 2000.0, 1: 0.0}, land_after_reads=3)
    naps = []
    held = setup.fund_once(main, 1, io.StringIO(), sleep=naps.append)
    assert held == 500.0 and naps == [setup.FUNDING_LAND_RETRY_SEC] * 3
    assert len(main.posts(setup.TRANSFER_PATH)) == 1


def test_a_transfer_still_missing_after_the_wait_stops_without_a_second_transfer(tmp_path):
    main = MainKalshi({0: 2000.0, 1: 0.0}, land_after_reads=100)
    naps = []
    with pytest.raises(setup.SetupError, match=r"holds \$0.00 after 30 s"):
        setup.fund_once(main, 1, io.StringIO(), sleep=naps.append)
    assert sum(naps) == setup.FUNDING_LAND_WAIT_SEC and len(main.posts(setup.TRANSFER_PATH)) == 1


def test_a_finished_setup_makes_no_second_key(tmp_path):
    key_file = tmp_path / "kalshi_taker.pem"
    key_file.write_text(PEM)
    main = MainKalshi({0: 1500.0, 1: 500.0}, listed=[dict(LOCKED_RECORD, api_key_id="key-live")])
    code, output, made = run(main, tmp_path, configured_key_id="key-live", configured_subaccount=1)
    assert code == 0 and "already locked to subaccount 1; no new key" in output
    assert main.posts(setup.GENERATE_KEY_PATH) == [] and made == [("key-live", key_file)]


def test_a_configured_key_kalshi_no_longer_lists_is_replaced(tmp_path):
    key_file = tmp_path / "kalshi_taker.pem"
    key_file.write_text("old pem")
    main = MainKalshi({0: 1500.0, 1: 500.0}, listed=[LOCKED_RECORD])  # only the new key, never "key-gone"
    code, output, made = run(main, tmp_path, configured_key_id="key-gone", configured_subaccount=1)
    assert code == 0, output
    assert [body["subaccount"] for body in main.posts(setup.GENERATE_KEY_PATH)] == [1]
    assert key_file.read_text() == PEM and made == [(NEW_KEY_ID, key_file)]


@pytest.mark.parametrize("record, why", [
    (dict(LOCKED_RECORD, subaccount=None), "with subaccount None, not 1"),
    (dict(LOCKED_RECORD, subaccount=2), "with subaccount 2, not 1"),
    (dict(LOCKED_RECORD, scopes=["read", "write"]), "with scopes ['read', 'write']"),
    (dict(LOCKED_RECORD, scopes=["read", "write::trade", "write::transfer"]), "write::transfer"),
    (dict(LOCKED_RECORD, fcm_subtrader_id="fcm-1"), "bound to an FCM subtrader"),
    (dict(LOCKED_RECORD, api_key_id="someone-else"), "listing does not carry new-key-0001"),
])
def test_a_new_key_kalshi_lists_as_not_locked_is_deleted_and_named_for_revoking(tmp_path, record, why):
    main = MainKalshi({0: 2000.0}, listed=[record])
    code, output, made = run(main, tmp_path)
    assert code == 1 and why in output and "Revoke new-key-0001 at kalshi.com" in output
    assert not (tmp_path / "kalshi_taker.pem").exists()
    assert made == []  # the bot key is never tried on Kalshi once its listing is wrong


def test_a_new_key_the_balance_check_shows_unlocked_is_deleted_and_named_for_revoking(tmp_path):
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
