"""HTTPS proxy on the Mac that lets the pricer on the Oracle VM reach DraftKings.

DraftKings answers "Access Denied" to data-center IPs, so the VM sends its
DraftKings requests here over Tailscale and they leave from this Mac's home IP.
It is a CONNECT tunnel: TLS runs end to end between the VM's client and
DraftKings, so this process sees only the target host, never a request.

    python3 -m nfl_specials.dk_relay     (launchd: nfl_specials/com.nflwork.dk-relay.plist)

Listens on this Mac's Tailscale IPv4 (`tailscale ip -4`, waiting until
Tailscale is up) at config.DK_RELAY_PORT, so only devices on the tailnet can
reach it, and tunnels only to dk_book.DK_HOSTS on 443; any other target gets
403. No DB, no files: logs go to stderr (start, waiting for Tailscale, refused
targets, upstream failures).
"""
from __future__ import annotations

import asyncio
import logging
import subprocess
import time

from nfl_specials import config
from nfl_specials.dk_book import DK_HOSTS

log = logging.getLogger("nfl_specials.dk_relay")

HTTPS_PORT = 443
# A CONNECT head is one short line plus a Host header.
MAX_HEADER_BYTES = 8192
HEADER_TIMEOUT_SECONDS = 10
UPSTREAM_CONNECT_TIMEOUT_SECONDS = 10
CHUNK_BYTES = 65536
TAILSCALE_RETRY_SECONDS = 30


def parse_connect_target(head: bytes) -> tuple[str, int] | None:
    """(host, port) from a "CONNECT host:port HTTP/1.1" request head, else None."""
    request_line = head.split(b"\r\n", 1)[0].decode("latin-1")
    parts = request_line.split(" ")
    if len(parts) != 3 or parts[0] != "CONNECT" or not parts[2].startswith("HTTP/1."):
        return None
    host, separator, port = parts[1].rpartition(":")
    if not separator or not host or not (port.isascii() and port.isdigit()):
        return None
    return host.lower(), int(port)


def is_allowed(target: tuple[str, int] | None) -> bool:
    return target is not None and target[0] in DK_HOSTS and target[1] == HTTPS_PORT


async def _pipe(reader: asyncio.StreamReader, writer: asyncio.StreamWriter) -> None:
    try:
        while chunk := await reader.read(CHUNK_BYTES):
            writer.write(chunk)
            await writer.drain()
    except OSError:
        pass
    finally:
        writer.close()


async def _answer_and_close(writer: asyncio.StreamWriter, status_line: bytes) -> None:
    try:
        writer.write(status_line + b"\r\nContent-Length: 0\r\n\r\n")
        await writer.drain()
    except OSError:
        pass
    finally:
        writer.close()


async def handle_client(reader: asyncio.StreamReader, writer: asyncio.StreamWriter) -> None:
    peer = writer.get_extra_info("peername")
    try:
        head = await asyncio.wait_for(reader.readuntil(b"\r\n\r\n"), HEADER_TIMEOUT_SECONDS)
    except (TimeoutError, asyncio.IncompleteReadError, asyncio.LimitOverrunError, OSError):
        writer.close()
        return
    target = parse_connect_target(head)
    if not is_allowed(target):
        log.warning("refused %s from %s", target or "a request that is not CONNECT", peer)
        await _answer_and_close(writer, b"HTTP/1.1 403 Forbidden")
        return
    try:
        upstream_reader, upstream_writer = await asyncio.wait_for(
            asyncio.open_connection(*target), UPSTREAM_CONNECT_TIMEOUT_SECONDS)
    except (TimeoutError, OSError) as exc:
        log.warning("could not reach %s:%d for %s: %s", target[0], target[1], peer, exc)
        await _answer_and_close(writer, b"HTTP/1.1 502 Bad Gateway")
        return
    try:
        writer.write(b"HTTP/1.1 200 Connection established\r\n\r\n")
        await writer.drain()
    except OSError:
        upstream_writer.close()
        writer.close()
        return
    await asyncio.gather(_pipe(reader, upstream_writer), _pipe(upstream_reader, writer))


def tailscale_ipv4() -> str:
    """This Mac's Tailscale IPv4. Raises when Tailscale is not up."""
    result = subprocess.run(["tailscale", "ip", "-4"], capture_output=True, text=True, timeout=10, check=False)
    lines = result.stdout.split()
    if result.returncode != 0 or not lines or not lines[0].startswith("100."):
        raise RuntimeError(f"expected a 100.x address from `tailscale ip -4`, got exit {result.returncode}: "
                           f"{result.stdout.strip()!r} {result.stderr.strip()!r}")
    return lines[0]


async def serve(host: str, port: int) -> None:
    server = await asyncio.start_server(handle_client, host, port, limit=MAX_HEADER_BYTES)
    log.info("dk_relay listening on %s:%d, tunnelling to %s on %d",
             host, port, ", ".join(sorted(DK_HOSTS)), HTTPS_PORT)
    async with server:
        await server.serve_forever()


def wait_for_tailscale_ipv4() -> str:
    """tailscale_ipv4(), retried until Tailscale is up; logs the first failure only."""
    waiting_logged = False
    while True:
        try:
            return tailscale_ipv4()
        except (RuntimeError, OSError, subprocess.TimeoutExpired) as exc:
            if not waiting_logged:
                log.warning("waiting for Tailscale: %s", exc)
                waiting_logged = True
            time.sleep(TAILSCALE_RETRY_SECONDS)


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    asyncio.run(serve(wait_for_tailscale_ipv4(), config.DK_RELAY_PORT))


if __name__ == "__main__":
    main()
