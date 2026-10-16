"""novig_subaccount_setup.py on a fake Novig: keygen prints only the public
half and stops the Novig desk first, open refuses a private key on stdin, finds
or opens the labelled subaccount, checks the trading key's record, and funds
under an idempotent transfer id polled to Applied."""
import io
import json
import os
import stat

import pytest

from unabated_ticket.bets_service import novig_subaccount_setup as setup
from unabated_ticket.bets_service.novig_orders import NovigOrderClient, public_key_pem_of

KEY_ID = "8f14e45f-ceea-467a-9b1c-3f2a51c8d7e0"


def _public_pem(tmp_path):
    key_path = tmp_path / "vm" / "novig_bot.pem"
    out, err = io.StringIO(), io.StringIO()
    assert setup.run_keygen("bot", key_path, tmp_path / "vm" / "novig_taker.stop", out, err) == 0
    return out.getvalue()


def _fake_novig(tmp_path, routes):
    """A management client whose transport answers (method, path) from `routes`
    (a list per route, consumed in order) and records every request."""
    management_pem = tmp_path / "mgmt.pem"
    setup.write_new_private_key(management_pem)
    sent = []

    def transport(method, url, headers, body):
        path = url.split("api.novig.com", 1)[1].split("?", 1)[0]
        sent.append((method, path, json.loads(body) if body else None))
        status, reply = routes[(method, path)].pop(0)
        return status, json.dumps(reply).encode()

    client = NovigOrderClient("mgmt-key", management_pem, "https://api.novig.com", transport=transport,
                              clock_ms=lambda: 1755000000000)
    return client, sent


def _key_record(public_pem, label="unabated-bot", scope="trading"):
    return {"keyId": KEY_ID, "scope": scope, "restriction": label, "fingerprint": setup.fingerprint_of(public_pem)}


def test_keygen_stops_the_desk_first_and_prints_only_the_public_half(tmp_path):
    key_path = tmp_path / "novig_bot.pem"
    stop_path = tmp_path / "novig_taker.stop"
    out, err = io.StringIO(), io.StringIO()
    assert setup.run_keygen("bot", key_path, stop_path, out, err) == 0
    assert stop_path.exists()
    assert out.getvalue().startswith("-----BEGIN PUBLIC KEY-----") and "PRIVATE" not in out.getvalue()
    assert stat.S_IMODE(os.stat(key_path).st_mode) == 0o600
    assert out.getvalue() == public_key_pem_of(key_path.read_bytes())


def test_keygen_keeps_an_existing_key_and_place_never_touches_the_stop_file(tmp_path):
    key_path = tmp_path / "novig_place.pem"
    stop_path = tmp_path / "novig_taker.stop"
    first, second = io.StringIO(), io.StringIO()
    assert setup.run_keygen("place", key_path, stop_path, first, io.StringIO()) == 0
    assert setup.run_keygen("place", key_path, stop_path, second, io.StringIO()) == 0
    assert first.getvalue() == second.getvalue()
    assert not stop_path.exists()


def test_keygen_refuses_an_unusable_key_file(tmp_path):
    key_path = tmp_path / "novig_bot.pem"
    key_path.write_text("not a key")
    out, err = io.StringIO(), io.StringIO()
    assert setup.run_keygen("bot", key_path, tmp_path / "stop", out, err) == 1
    assert out.getvalue() == "" and "STOPPED" in err.getvalue()
    assert key_path.read_text() == "not a key"


def test_open_refuses_a_private_key_on_stdin(tmp_path):
    client, sent = _fake_novig(tmp_path, {})
    out = io.StringIO()
    assert setup.run_open(client, "bot", "-----BEGIN PRIVATE KEY-----\nx\n", 500, "initial", out) == 1
    assert "PRIVATE" in out.getvalue() and sent == []


def test_open_opens_checks_and_funds(tmp_path):
    public_pem = _public_pem(tmp_path)
    routes = {
        ("GET", "/v3/account/subaccounts"): [(200, {"items": [{"keyId": "other", "label": "unabated-place"}]})],
        ("POST", "/v3/account/subaccounts"): [(201, {"keyId": KEY_ID, "label": "unabated-bot", "balance": "0"})],
        ("GET", f"/v3/keys/{KEY_ID}"): [(200, _key_record(public_pem))],
        ("POST", f"/v3/account/subaccounts/{KEY_ID}/transfer"): [(202, {"transferId": "t1", "status": "Requested"})],
        ("GET", f"/v3/account/subaccounts/{KEY_ID}/transfers/t1"): [
            (200, {"transferId": "t1", "status": "Requested"}),
            (200, {"transferId": "t1", "status": "Applied", "amount": "500.00000"})],
        ("GET", f"/v3/account/subaccounts/{KEY_ID}/balance"): [(200, {"keyId": KEY_ID, "balance": "500.00000"})],
    }
    client, sent = _fake_novig(tmp_path, routes)
    out = io.StringIO()
    assert setup.run_open(client, "bot", public_pem, 500, "initial", out, sleep=lambda _s: None) == 0
    opened = next(body for method, path, body in sent if (method, path) == ("POST", "/v3/account/subaccounts"))
    assert opened == {"label": "unabated-bot", "publicKey": public_pem, "algorithm": "Ed25519"}
    funded = next(body for method, path, body in sent if path.endswith("/transfer"))
    assert funded == {"direction": "fund", "amount": "500.00000", "clientTransferId": "unabated-bot-initial"}
    assert f"NOVIG_TAKER_KEY_ID={KEY_ID}" in out.getvalue() and "$500.00" in out.getvalue()


def test_open_reuses_the_labelled_subaccount_and_refuses_another_key(tmp_path):
    public_pem = _public_pem(tmp_path)
    other_pem = public_key_pem_of(_other_key(tmp_path))
    routes = {
        ("GET", "/v3/account/subaccounts"): [(200, [{"keyId": KEY_ID, "label": "unabated-bot", "balance": "500"}])],
        ("GET", f"/v3/keys/{KEY_ID}"): [(200, _key_record(other_pem))],
    }
    client, sent = _fake_novig(tmp_path, routes)
    out = io.StringIO()
    assert setup.run_open(client, "bot", public_pem, 500, "initial", out) == 1
    assert "not the VM's key" in out.getvalue()
    assert [method for method, _path, _body in sent] == ["GET", "GET"]  # nothing opened or moved


def _other_key(tmp_path):
    path = tmp_path / "other.pem"
    setup.write_new_private_key(path)
    return path.read_bytes()


@pytest.mark.parametrize("record_change, reason", [({"scope": "management"}, "scope"),
                                                   ({"restriction": "unabated-place"}, "restricted")])
def test_trading_key_problem_names_the_mismatch(tmp_path, record_change, reason):
    public_pem = _public_pem(tmp_path)
    record = {**_key_record(public_pem), **record_change}
    client, _sent = _fake_novig(tmp_path, {("GET", f"/v3/keys/{KEY_ID}"): [(200, record)]})
    assert reason in setup.trading_key_problem(client, KEY_ID, "unabated-bot", public_pem)


def test_a_rejected_transfer_stops_without_a_balance_line(tmp_path):
    public_pem = _public_pem(tmp_path)
    routes = {
        ("GET", "/v3/account/subaccounts"): [(200, {"items": [{"keyId": KEY_ID, "label": "unabated-bot"}]})],
        ("GET", f"/v3/keys/{KEY_ID}"): [(200, _key_record(public_pem))],
        ("POST", f"/v3/account/subaccounts/{KEY_ID}/transfer"): [
            (202, {"transferId": "t1", "status": "Rejected", "actualBalance": "12.00000"})],
    }
    client, _sent = _fake_novig(tmp_path, routes)
    out = io.StringIO()
    assert setup.run_open(client, "bot", public_pem, 500, "topup-1", out) == 1
    assert "rejected" in out.getvalue() and "NOVIG_TAKER_KEY_ID" not in out.getvalue()


def test_a_transfer_still_requested_after_polling_says_do_not_resend(tmp_path):
    public_pem = _public_pem(tmp_path)
    pending = (200, {"transferId": "t1", "status": "Requested"})
    routes = {
        ("GET", "/v3/account/subaccounts"): [(200, {"items": [{"keyId": KEY_ID, "label": "unabated-bot"}]})],
        ("GET", f"/v3/keys/{KEY_ID}"): [(200, _key_record(public_pem))],
        ("POST", f"/v3/account/subaccounts/{KEY_ID}/transfer"): [(202, {"transferId": "t1", "status": "Requested"})],
        ("GET", f"/v3/account/subaccounts/{KEY_ID}/transfers/t1"): [pending] * setup.TRANSFER_POLL_ATTEMPTS,
    }
    client, _sent = _fake_novig(tmp_path, routes)
    out = io.StringIO()
    assert setup.run_open(client, "bot", public_pem, 500, "initial", out, sleep=lambda _s: None) == 1
    assert "do not re-send" in out.getvalue()


def test_transfer_id_is_per_name_and_run():
    assert setup.transfer_id_of("place", "topup-2") == "unabated-place-topup-2"
    with pytest.raises(setup.SetupError):
        setup.transfer_id_of("bot", "x" * 80)
