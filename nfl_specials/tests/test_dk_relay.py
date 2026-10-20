import asyncio
import socket

import pytest

from nfl_specials import dk_relay


@pytest.mark.parametrize("head, expected", [
    (b"CONNECT sportsbook-nash.draftkings.com:443 HTTP/1.1\r\nHost: x\r\n\r\n",
     ("sportsbook-nash.draftkings.com", 443)),
    (b"CONNECT Sportsbook.DraftKings.com:443 HTTP/1.0\r\n\r\n", ("sportsbook.draftkings.com", 443)),
    (b"GET http://sportsbook.draftkings.com/ HTTP/1.1\r\n\r\n", None),
    (b"CONNECT sportsbook.draftkings.com HTTP/1.1\r\n\r\n", None),
    (b"CONNECT :443 HTTP/1.1\r\n\r\n", None),
    (b"CONNECT sportsbook.draftkings.com:443\r\n\r\n", None),
    ("CONNECT sportsbook.draftkings.com:\u00b2 HTTP/1.1\r\n\r\n".encode("latin-1"), None),
])
def test_connect_target_is_read_from_the_request_line(head, expected):
    assert dk_relay.parse_connect_target(head) == expected


@pytest.mark.parametrize("target, expected", [
    (("sportsbook-nash.draftkings.com", 443), True),
    (("sportsbook.draftkings.com", 443), True),
    (("sportsbook.draftkings.com", 80), False),
    (("example.com", 443), False),
    (None, False),
])
def test_only_draftkings_on_443_is_tunnelled(target, expected):
    assert dk_relay.is_allowed(target) is expected


async def _relay_round_trip(request_head: bytes, payload: bytes) -> tuple[bytes, bytes]:
    """Start an echo server and the relay on loopback, send `request_head`
    then `payload` through the relay, and return (status line, echoed bytes)."""
    async def echo(reader, writer):
        writer.write(await reader.read(1024))
        await writer.drain()
        writer.close()

    upstream = await asyncio.start_server(echo, "127.0.0.1", 0)
    relay = await asyncio.start_server(dk_relay.handle_client, "127.0.0.1", 0, limit=dk_relay.MAX_HEADER_BYTES)
    upstream_port = upstream.sockets[0].getsockname()[1]
    relay_port = relay.sockets[0].getsockname()[1]
    async with upstream, relay:
        reader, writer = await asyncio.open_connection("127.0.0.1", relay_port)
        writer.write(request_head.replace(b"{port}", str(upstream_port).encode()))
        await writer.drain()
        status_line = (await reader.readuntil(b"\r\n\r\n")).split(b"\r\n", 1)[0]
        echoed = b""
        if status_line.endswith(b"200 Connection established"):
            writer.write(payload)
            await writer.drain()
            echoed = await reader.read(1024)
        writer.close()
        return status_line, echoed


def test_an_allowed_target_gets_a_byte_for_byte_tunnel(monkeypatch):
    head = b"CONNECT 127.0.0.1:{port} HTTP/1.1\r\nHost: 127.0.0.1\r\n\r\n"
    monkeypatch.setattr(dk_relay, "is_allowed", lambda target: target is not None and target[0] == "127.0.0.1")
    assert asyncio.run(_relay_round_trip(head, b"tls bytes")) == (b"HTTP/1.1 200 Connection established", b"tls bytes")


def test_any_other_target_is_refused():
    head = b"CONNECT example.com:443 HTTP/1.1\r\nHost: example.com\r\n\r\n"
    assert asyncio.run(_relay_round_trip(head, b"never sent")) == (b"HTTP/1.1 403 Forbidden", b"")


def test_an_unreachable_target_answers_502(monkeypatch):
    with socket.socket() as probe:
        probe.bind(("127.0.0.1", 0))
        closed_port = probe.getsockname()[1]
    monkeypatch.setattr(dk_relay, "is_allowed", lambda target: target == ("127.0.0.1", closed_port))
    head = f"CONNECT 127.0.0.1:{closed_port} HTTP/1.1\r\n\r\n".encode()
    assert asyncio.run(_relay_round_trip(head, b"never sent")) == (b"HTTP/1.1 502 Bad Gateway", b"")
