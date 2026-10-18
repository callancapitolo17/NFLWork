"""novig_orders.py: NOVIG-V3 signing against Novig's own published test
vectors (tests/fixtures/bets/novig_signing_vectors.json, docs.novig.com
signing-vectors.json, 30 vectors and two test keypairs), the price grid,
1c-contract sizing, and the client's request shape on a fake transport."""
import base64
import json
from pathlib import Path

import pytest
from cryptography.hazmat.primitives import hashes, serialization
from cryptography.hazmat.primitives.asymmetric import ec

from unabated_ticket.bets_service import novig_orders
from unabated_ticket.bets_service.novig_orders import (NovigError, NovigOrderClient, NovigSigner, canonical_query,
                                                       cost_of, ioc_order_body, is_on_price_grid, order_totals,
                                                       price_string, qty_for_stake, string_to_sign)

VECTORS_PATH = Path(__file__).parents[2] / "tests" / "fixtures" / "bets" / "novig_signing_vectors.json"
VECTORS = json.loads(VECTORS_PATH.read_text())
KEYPAIRS = VECTORS["keypairs"]


def _vector_ids(algorithm):
    return [vector["id"] for vector in VECTORS["vectors"] if vector["algorithm"] == algorithm]


def _vector(vector_id):
    return next(vector for vector in VECTORS["vectors"] if vector["id"] == vector_id)


def test_the_fixture_is_novigs_thirty_vectors():
    assert len(VECTORS["vectors"]) == 30


@pytest.mark.parametrize("vector_id", [vector["id"] for vector in VECTORS["vectors"]])
def test_string_to_sign_matches_every_vector(vector_id):
    vector = _vector(vector_id)
    given = vector["input"]
    built = string_to_sign(given["timestamp"], given["method"], given["path"], given["query"], given["body"].encode())
    assert built == vector["string_to_sign"]


@pytest.mark.parametrize("vector_id", _vector_ids("ed25519"))
def test_ed25519_signature_matches_every_vector(vector_id):
    # Ed25519 is deterministic: the same key and string give Novig's exact signature.
    vector = _vector(vector_id)
    signer = NovigSigner(KEYPAIRS[vector["keypair_id"]]["private_key_pkcs8_pem"].encode())
    assert signer.algorithm == "Ed25519"
    assert signer.sign(vector["string_to_sign"]) == vector["signature"]


@pytest.mark.parametrize("vector_id", _vector_ids("ecdsa-p256"))
def test_p256_signature_verifies_as_der(vector_id):
    # ECDSA is randomised: check ours verifies under the published public key (DER r,s).
    vector = _vector(vector_id)
    pair = KEYPAIRS[vector["keypair_id"]]
    signer = NovigSigner(pair["private_key_pkcs8_pem"].encode())
    assert signer.algorithm == "P-256"
    signature = base64.b64decode(signer.sign(vector["string_to_sign"]))
    public_key = serialization.load_pem_public_key(pair["public_key_spki_pem"].encode())
    public_key.verify(signature, vector["string_to_sign"].encode(), ec.ECDSA(hashes.SHA256()))
    assert signature[0] == 0x30  # DER SEQUENCE, never raw r||s


def test_canonical_query_rules():
    assert canonical_query("") == ""
    assert canonical_query("b=2&a=1") == "a=1&b=2"
    assert canonical_query("q=foo+bar") == "q=foo%2Bbar"
    assert canonical_query("flag") == "flag="
    assert canonical_query("token=abc=") == "token=abc%3D"


def test_price_grid_bands():
    assert is_on_price_grid(0.001) and is_on_price_grid(0.05) and is_on_price_grid(0.999)
    assert is_on_price_grid(0.465) and is_on_price_grid(0.535) and is_on_price_grid(0.95)
    assert not is_on_price_grid(0.051) and not is_on_price_grid(0.467) and not is_on_price_grid(0.0)
    assert not is_on_price_grid(1.0) and not is_on_price_grid(0.4651)
    assert price_string(0.465) == "0.465"
    with pytest.raises(NovigError):
        price_string(0.467)


def test_one_minus_a_grid_price_is_on_the_grid():
    # The taker buys at 1 - the best opposite order; that must always be submittable.
    for thousandths in range(1, 1000):
        price = thousandths / 1000
        if is_on_price_grid(price):
            assert is_on_price_grid(round(1 - price, 3)), price


def test_contracts_pay_one_cent_not_a_dollar():
    # $10 at 0.465: 2150 contracts of 1c (a Kalshi-style $1 count would be 21).
    assert qty_for_stake(10, 0.465) == 2150
    assert cost_of(2150, 0.465) == pytest.approx(9.9975)
    assert qty_for_stake(0, 0.5) == 0
    assert qty_for_stake(0.004, 0.5) == 0


def test_ioc_order_body():
    body = ioc_order_body("3f2504e0-4f89-11d3-9a0c-0305e82c3301", 0.465, 2150, "0b3f5f20-9c1b-4d6e-8f3a-5c2d1e7b4a90")
    assert body == {"outcomeId": "3f2504e0-4f89-11d3-9a0c-0305e82c3301", "price": "0.465", "qty": 2150, "tif": "IOC",
                    "clientId": "0b3f5f20-9c1b-4d6e-8f3a-5c2d1e7b4a90"}
    with pytest.raises(NovigError):
        ioc_order_body("x", 0.465, 0, "y")


def test_order_totals_refuses_an_unreadable_fill():
    assert order_totals([{"fillId": "a", "qty": 100, "cost": "0.46500"}, {"fillId": "b", "qty": 50, "cost": "0.2", "fee": "0.01"}]) == (150, 0.665, 0.01)
    with pytest.raises(NovigError):
        order_totals([{"fillId": "a", "qty": 100}])


def _client(tmp_path, replies):
    key_path = tmp_path / "trading.pem"
    key_path.write_text(KEYPAIRS["ed25519-test-1"]["private_key_pkcs8_pem"])
    sent = []

    def transport(method, url, headers, body):
        sent.append({"method": method, "url": url, "headers": headers, "body": body})
        status, reply = replies.pop(0)
        return status, json.dumps(reply).encode()

    client = NovigOrderClient("8f14e45f-ceea-467a-9b1c-3f2a51c8d7e0", key_path, "https://api.novig.com/",
                              transport=transport, clock_ms=lambda: 1755000000000)
    return client, sent


def test_client_signs_the_exact_body_it_sends(tmp_path):
    client, sent = _client(tmp_path, [(201, {"orderId": "o1", "clientId": "c1"})])
    status, reply = client.place_ioc_order("3f2504e0-4f89-11d3-9a0c-0305e82c3301", 0.465, 2150, "c1")
    assert (status, reply["orderId"]) == (201, "o1")
    request = sent[0]
    assert request["method"] == "POST" and request["url"] == "https://api.novig.com/v3/orders"
    assert request["headers"]["Content-Type"] == "application/json"
    assert request["headers"]["Novig-Timestamp"] == "1755000000000"
    signer = NovigSigner(KEYPAIRS["ed25519-test-1"]["private_key_pkcs8_pem"].encode())
    expected = signer.sign(string_to_sign(1755000000000, "POST", "/v3/orders", "", request["body"]))
    assert request["headers"]["Novig-Signature"] == expected
    assert json.loads(request["body"])["tif"] == "IOC"


def test_paged_reads_follow_next_and_refuse_a_failed_page(tmp_path):
    client, sent = _client(tmp_path, [(200, {"items": [{"orderId": "a", "clientId": "x"}], "next": "p2"}),
                                      (200, {"items": [{"orderId": "b", "clientId": "c9"}]})])
    assert [order["orderId"] for order in client.orders("out-1", "FILLED")] == ["a", "b"]
    assert "after=p2" in sent[1]["url"] and "status=FILLED" in sent[1]["url"]
    client, _ = _client(tmp_path, [(500, {"code": "INTERNAL_ERROR"})])
    with pytest.raises(NovigError):
        client.fills(order_id="o1")


def test_find_order_by_client_id_looks_in_every_status(tmp_path):
    replies = [(200, {"items": []}), (200, {"items": []}), (200, {"items": [{"orderId": "o7", "clientId": "c7"}]})]
    client, sent = _client(tmp_path, replies)
    assert client.find_order_by_client_id("out-1", "c7")["orderId"] == "o7"
    assert len(sent) == 3
    client, _ = _client(tmp_path, [(200, {"items": []})] * len(novig_orders.ORDER_STATUSES))
    assert client.find_order_by_client_id("out-1", "c7") is None


def test_balance_reads_its_own_subaccount(tmp_path):
    client, sent = _client(tmp_path, [(200, {"keyId": "k", "balance": "500.00000"})])
    assert client.balance() == 500.0
    assert sent[0]["url"].endswith("/v3/account/subaccounts/8f14e45f-ceea-467a-9b1c-3f2a51c8d7e0/balance")


def test_market_if_listed_is_none_only_when_novig_says_market_not_found(tmp_path):
    market = {"marketId": "m1", "status": "OPEN"}
    client, sent = _client(tmp_path, [(200, market), (404, {"code": "MARKET_NOT_FOUND", "message": "Market not found"})])
    assert client.market_if_listed("m1") == market
    assert client.market_if_listed("m2") is None
    assert sent[1]["url"].endswith("/m2")
    client, _ = _client(tmp_path, [(404, {"code": "NOT_FOUND"})])
    with pytest.raises(NovigError):
        client.market_if_listed("m3")
    client, _ = _client(tmp_path, [(500, {"code": "INTERNAL_ERROR"})])
    with pytest.raises(NovigError):
        client.market_if_listed("m4")
