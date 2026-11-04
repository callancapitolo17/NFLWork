"""ProphetX's order stream (prophetx_order_stream.py) against a scripted
Pusher-protocol socket and client. Frame shapes are built from ProphetX's
websocket docs (no live capture yet). No network."""
import base64
import json
import threading

import pytest

from unabated_ticket.bets_service.prophetx_order_stream import (ProphetXOrderStream, order_of_event, parse_frame,
                                                                pusher_url)

CHANNEL = "private-service=6-device_type=5-user=42"


def frame(event, data, channel=None):
    message = {"event": event, "data": json.dumps(data)}
    if channel:
        message["channel"] = channel
    return json.dumps(message)


def orders_frame(info):
    payload = base64.b64encode(json.dumps({"info": info}).encode()).decode()
    return frame("orders", {"timestamp": 1, "change_type": "orders", "op": "u", "payload": payload}, CHANNEL)


class FakeClock:
    def __init__(self):
        self.now = 0.0

    def __call__(self):
        return self.now


class FakeConnection:
    """Hands out scripted frames, one per recv, then raises (the socket closed)."""

    def __init__(self, frames, clock):
        self.frames = list(frames)
        self.sent = []
        self.clock = clock

    def send(self, text):
        self.sent.append(json.loads(text))

    def recv(self, timeout_sec):
        self.clock.now += 0.1
        if not self.frames:
            raise ConnectionError("closed")
        item = self.frames.pop(0)
        if isinstance(item, Exception):
            raise item
        return item

    def close(self):
        pass


class FakeClient:
    def __init__(self):
        self.registered = []

    def websocket_config(self):
        return {"app_id": "1", "key": "abc", "cluster": "mt1", "service": "6", "ws_host": "ws.prophetx.co"}

    def register_websocket(self, socket_id, service):
        self.registered.append((socket_id, service))
        return {"authenticated": {"auth": "sig", "user_data": "{\"id\":\"42\"}"},
                "authorized_channel": [{"channel_name": "public", "auth": None},
                                       {"channel_name": CHANNEL, "auth": "chan-sig"}]}


HANDSHAKE = [frame("pusher:connection_established", {"socket_id": "1.2", "activity_timeout": 120}),
             frame("pusher:signin_success", {}),
             frame("pusher_internal:subscription_succeeded", {}, CHANNEL)]


def stream_with(frames):
    clock = FakeClock()
    connection = FakeConnection(frames, clock)
    client = FakeClient()
    stream = ProphetXOrderStream(client, connect=lambda url, headers: connection, clock=clock)
    return stream, connection, client, clock


def test_url_uses_ws_host_else_the_pusher_cluster():
    assert pusher_url({"key": "k", "ws_host": "wss://ws.prophetx.co/"}).startswith("wss://ws.prophetx.co/app/k?protocol=7")
    assert pusher_url({"key": "k", "cluster": "mt1"}).startswith("wss://ws-mt1.pusher.com/app/k?")
    with pytest.raises(ValueError):
        pusher_url({"key": "k"})


def test_an_orders_event_reads_as_the_order_with_its_running_total():
    data = parse_frame(orders_frame({"id": "px-1", "filled_quantity": 5, "total_filled_quantity": 50,
                                     "matching_status": "fully_matched"}))[2]
    order = order_of_event(data)
    assert order["order_id"] == "px-1" and order["filled_quantity"] == 50
    assert order_of_event({"payload": base64.b64encode(b'{"info": {}}').decode()}) is None


def test_the_handshake_signs_in_subscribes_the_private_channel_and_answers_pings():
    stream, connection, client, _clock = stream_with(HANDSHAKE + [frame("pusher:ping", {})])
    stream._open()
    assert client.registered == [("1.2", "6")]
    assert connection.sent[0] == {"event": "pusher:signin", "data": {"auth": "sig", "user_data": "{\"id\":\"42\"}"}}
    assert connection.sent[1] == {"event": "pusher:subscribe", "data": {"channel": CHANNEL, "auth": "chan-sig"}}
    assert stream.is_live()
    with pytest.raises(ConnectionError):
        stream._pump(threading.Event())
    assert connection.sent[2] == {"event": "pusher:pong", "data": {}}


def test_place_gets_each_newer_pushed_state_and_none_once_the_stream_is_down():
    stream, _connection, _client, _clock = stream_with(HANDSHAKE + [
        orders_frame({"id": "px-1", "matching_status": "unmatched", "total_filled_quantity": 0}),
        orders_frame({"id": "px-1", "matching_status": "fully_matched", "total_filled_quantity": 50})])
    stream._open()
    with pytest.raises(ConnectionError):
        stream._pump(threading.Event())
    order, version = stream.wait_for_order("px-1", 0, 1.0)
    assert order["matching_status"] == "fully_matched" and version == 2
    stream._drop("closed")
    assert not stream.is_live()
    assert stream.wait_for_order("px-1", version, 1.0) == (None, version)
    assert stream.status()["ordersEvents"] == 2 and stream.status()["error"] == "closed"


def test_a_silent_socket_is_dropped():
    stream, connection, _client, clock = stream_with(HANDSHAKE)
    stream._open()
    connection.frames = [None] * 400
    with pytest.raises(TimeoutError):
        stream._pump(threading.Event())


def test_a_register_without_a_private_channel_fails_the_connect():
    stream, _connection, client, _clock = stream_with(HANDSHAKE)
    client.register_websocket = lambda socket_id, service: {"authenticated": {"auth": "sig"},
                                                            "authorized_channel": [{"channel_name": "public"}]}
    with pytest.raises(ValueError, match="no private user channel"):
        stream._open()


def test_the_websocket_endpoints_answer_flat_or_wrapped():
    # 2026-11-04 live: connection-config came back with no `data` wrapper.
    from unabated_ticket.bets_service.prophetx_client import data_or_flat
    flat = {"app_id": "1810913", "cluster": "mt1", "key": "k", "service": "pusher", "ws_host": "ws.prophetx.co"}
    assert data_or_flat(flat) == flat
    assert data_or_flat({"data": flat}) == flat
    assert pusher_url(flat).startswith("wss://ws.prophetx.co/app/k?protocol=7")
