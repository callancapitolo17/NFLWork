"""The Edges settings the server runner reads (phone page plan step 1):
bets.duckdb::edge_settings (one row, NULL = the panel's default), the PUT
body's validation, and GET/PUT /settings.json over a real HTTP server with
the same Host and Content-Type guards as the other writes."""
import json
import threading
import urllib.request
from datetime import datetime, timezone
from http.server import ThreadingHTTPServer
from urllib.parse import urlparse

import pytest

from unabated_ticket.bets_service import service
from unabated_ticket.bets_service.store import EDGE_SETTINGS_FIELDS, BetsStore

NOTHING_HELD = dict.fromkeys(EDGE_SETTINGS_FIELDS)
SAVED_AT = datetime(2026, 10, 1, 12, 0, tzinfo=timezone.utc)


@pytest.fixture
def store(tmp_path):
    bets_store = BetsStore(tmp_path / "bets.duckdb", 7)
    yield bets_store
    bets_store.close()


@pytest.fixture
def http_server(store):
    server = ThreadingHTTPServer(("127.0.0.1", 0), service.make_handler(store, 0.0, ["kalshi"]))
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    yield f"http://127.0.0.1:{server.server_address[1]}"
    server.shutdown()
    server.server_close()


def request_json(method: str, url: str, payload: object = None, content_type: str = "application/json",
                 host: str | None = None) -> tuple[int, dict]:
    body = None if payload is None else (payload if isinstance(payload, bytes) else json.dumps(payload).encode())
    headers = {"Content-Type": content_type} if body is not None else {}
    if host:
        headers["Host"] = host
    request = urllib.request.Request(url, data=body, method=method, headers=headers)
    try:
        with urllib.request.urlopen(request, timeout=5) as response:
            return response.status, json.loads(response.read())
    except urllib.error.HTTPError as error:
        return error.code, json.loads(error.read())


# ---- store -----------------------------------------------------------------------------

def test_before_any_write_every_field_is_the_default(store):
    assert store.load_edge_settings() == {"settings": NOTHING_HELD, "updatedAt": None}


def test_settings_round_trip_lists_booleans_and_nulls(store):
    settings = {**NOTHING_HELD, "bankroll": 8000, "multiplier": 0.5, "leagues": [1, 2], "periods": [1, 2],
                "betTypes": [2, 3], "bookMode": "custom", "bookIds": [89, 105], "includeAlts": True, "sortBy": "stake"}
    store.save_edge_settings(settings, SAVED_AT)
    loaded = store.load_edge_settings()
    assert loaded["updatedAt"] == "2026-10-01T12:00:00Z"
    assert loaded["settings"] == settings
    # The one row is replaced whole, never a second row.
    store.save_edge_settings({**settings, "bankroll": None}, SAVED_AT)
    assert store.load_edge_settings()["settings"]["bankroll"] is None
    assert store._con.execute("SELECT count(*) FROM edge_settings").fetchone()[0] == 1


def test_saving_without_every_field_fails_loudly(store):
    with pytest.raises(ValueError, match=r"missing \['groupByMarket'\]"):
        store.save_edge_settings({field: None for field in EDGE_SETTINGS_FIELDS if field != "groupByMarket"}, SAVED_AT)


# ---- validation ------------------------------------------------------------------------

@pytest.mark.parametrize("update, error", [
    ({"bankrol": 1}, f"unknown settings field(s) ['bankrol']; expected some of {list(EDGE_SETTINGS_FIELDS)}"),
    ({"bankroll": 0}, "settings.bankroll must be a number above 0, got 0"),
    ({"bankroll": "8000"}, "settings.bankroll must be a number above 0, got '8000'"),
    ({"multiplier": 1.5}, "settings.multiplier must be above 0 and at most 1, got 1.5"),
    ({"maxLineAgeHours": 0}, "settings.maxLineAgeHours must be a number above 0, got 0"),
    ({"minEdgePct": -1}, "settings.minEdgePct must be a number, 0 or more, got -1"),
    ({"minLiquidityToWin": True}, "settings.minLiquidityToWin must be a number, 0 or more, got True"),
    ({"includeAlts": "yes"}, "settings.includeAlts must be true or false, got 'yes'"),
    ({"periods": []}, "settings.periods must name at least one id, got []"),
    ({"periods": [9]}, "settings.periods must hold only [1, 2, 3, 4, 5, 6, 7], got [9]"),
    ({"betTypes": [4]}, "settings.betTypes must hold only [1, 2, 3], got [4]"),
    ({"leagues": [1, True]}, "settings.leagues must be a list of whole ids, got [1, True]"),
    ({"sortBy": "profit"}, "settings.sortBy must be one of ['edge', 'stake', 'start', 'exposure'], got 'profit'"),
    ({"bookMode": "mine"}, "settings.bookMode must be one of ['default', 'all', 'custom'], got 'mine'"),
    ({"bookMode": "custom"}, "settings.bookIds must be a list exactly when bookMode is 'custom', got bookMode 'custom' with bookIds None"),
    ({"bookIds": [89]}, "settings.bookIds must be a list exactly when bookMode is 'custom', got bookMode None with bookIds [89]"),
])
def test_a_bad_update_names_the_first_problem(update, error):
    assert service.validate_settings_update({"settings": update}, NOTHING_HELD) == error


def test_an_update_keeps_what_it_leaves_out_and_null_resets_to_the_default():
    held = {**NOTHING_HELD, "bankroll": 8000.0, "bookMode": "custom", "bookIds": [89]}
    merged = service.validate_settings_update({"settings": {"minEdgePct": 2, "bankroll": None}}, held)
    assert merged == {**held, "minEdgePct": 2, "bankroll": None}
    # Leaving custom books drops the ticks in the same request.
    assert service.validate_settings_update({"settings": {"bookMode": "all", "bookIds": None}}, held)["bookMode"] == "all"
    assert service.validate_settings_update([], held) == "body must be an object with a `settings` object, got list"


# ---- HTTP ------------------------------------------------------------------------------

def test_http_get_put_settings(store, http_server):
    status, reply = request_json("GET", f"{http_server}/settings.json")
    assert (status, reply) == (200, {"settings": NOTHING_HELD, "updatedAt": None})
    status, reply = request_json("PUT", f"{http_server}/settings.json", {"settings": {"bankroll": 8000, "includeAlts": True}})
    assert status == 200 and reply["ok"] is True and reply["updatedAt"] is not None
    assert (reply["settings"]["bankroll"], reply["settings"]["includeAlts"], reply["settings"]["multiplier"]) == (8000, True, None)
    status, reply = request_json("PUT", f"{http_server}/settings.json",
                                 {"settings": {"bookMode": "custom", "bookIds": [89, 105]}}, "application/json; charset=utf-8")
    assert status == 200
    assert (reply["settings"]["bankroll"], reply["settings"]["bookIds"]) == (8000, [89, 105])
    assert request_json("GET", f"{http_server}/settings.json")[1]["settings"] == reply["settings"]


def test_http_settings_refuses_non_json_bad_values_other_paths_and_foreign_hosts(store, http_server):
    update = {"settings": {"bankroll": 8000}}
    # The CSRF guard: a web page's cross-origin write can only be form- or text-encoded.
    status, reply = request_json("PUT", f"{http_server}/settings.json", update, "text/plain")
    assert status == 415 and "application/json" in reply["error"]
    status, reply = request_json("PUT", f"{http_server}/settings.json", {"settings": {"multiplier": 2}})
    assert status == 400 and reply["error"].startswith("settings.multiplier must be")
    status, reply = request_json("PUT", f"{http_server}/settings.json", b"{not json")
    assert status == 400 and reply["error"].startswith("body is not JSON")
    status, _reply = request_json("PUT", f"{http_server}/bets.json", update)
    assert status == 404
    status, _reply = request_json("POST", f"{http_server}/settings.json", update)
    assert status == 404
    for method, payload in [("GET", None), ("PUT", update)]:
        status, reply = request_json(method, f"{http_server}/settings.json", payload, host="evil.example")
        assert status == 403 and "Host must be one of" in reply["error"], method
    # Nothing was written by any refused request.
    assert store.load_edge_settings() == {"settings": NOTHING_HELD, "updatedAt": None}
    port = urlparse(http_server).port
    assert request_json("GET", f"{http_server}/settings.json", host=f"localhost:{port}")[0] == 200
