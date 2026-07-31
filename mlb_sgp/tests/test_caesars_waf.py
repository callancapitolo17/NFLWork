"""Regression tests for issue #41 — Caesars minted no AWS-WAF token, silently.

Caesars went to zero rows with zero log mentions. The ticket assumed AWS had
rotated challenge.js. It had not: the bundle fetched fine (200, ~654 KB) and
evaluated fine. The defect was OURS.

THE ROOT CAUSE. ``caesars_waf_node.js`` used to pre-stub
``window.AwsWafIntegration`` before eval'ing challenge.js, so the SDK's
auto-init would not crash on a not-yet-assigned global. AWS's bundle decides
whether to install itself by asking whether that global already exists:

    _0x5120cd = void 0x0 !== window['AwsWafIntegration']    // "already integrated?"

With our stub present the SDK concluded it was already integrated and never
installed the real object. ``getToken()`` then returned the STUB's empty
string, the broker reported "empty token", and the whole book was skipped.
Our own defensive shim was the poison.

Removing the stub re-exposes the benign auto-init race the stub was papering
over — the SDK reads ``window.AwsWafIntegration.checkForceRefresh`` while the
global is still undefined, producing an unhandled rejection. That rejection is
harmless (we drive ``forceRefreshToken`` ourselves), so it must NOT abort the
mint. The real gate is, and stays, "did a token come out".

Verified live 2026-07-30: mint 0.68s, token len 250, tabs feed 200 / 178 KB
JSON, ``price_sgps`` 8 rows, ``price_on_demand`` 2.84s.

No network: the node tests serve a stub challenge.js from a local HTTP server
that reproduces AWS's already-integrated guard, and the broker tests stub the
subprocess call entirely.
"""
from __future__ import annotations

import inspect
import json
import shutil
import subprocess
import threading
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path

import pytest

from mlb_sgp import caesars_waf

_HARNESS = Path(caesars_waf.__file__).resolve().parent / "caesars_waf_node.js"
_NODE = shutil.which("node")
needs_node = pytest.mark.skipif(_NODE is None, reason="node not on PATH")


# --------------------------------------------------------------------------- #
# A local stand-in for the AWS WAF edge issuer                                 #
# --------------------------------------------------------------------------- #

# Reproduces the two behaviours that matter, and nothing else:
#   1. it refuses to install itself if window.AwsWafIntegration already exists
#      (the exact guard that made our pre-stub fatal), and
#   2. it fires the benign auto-init rejection that the pre-stub existed to
#      suppress, so a fix that merely deletes the stub without tolerating the
#      rejection still fails here.
STUB_CHALLENGE_JS = """
(function () {
  // (2) the auto-init race: the SDK reads the global before it is assigned.
  Promise.resolve().then(function () {
    var captured = window.__never_assigned__;
    return captured.checkForceRefresh();      // -> unhandled rejection
  });
  // (1) AWS's real guard — bail out if something already claimed the name.
  if (typeof window.AwsWafIntegration !== 'undefined') { return; }
  window.AwsWafIntegration = {
    forceRefreshToken: function () { return Promise.resolve(); },
    getToken: function () { return Promise.resolve('STUB-WAF-TOKEN'); },
    hasToken: function () { return true; },
    checkForceRefresh: function () { return Promise.resolve(false); },
    saveReferrer: function () {},
  };
})();
"""

# A bundle that installs nothing at all — models a genuinely broken/rotated
# challenge.js, which must still produce a clean failure rather than a hang.
EMPTY_CHALLENGE_JS = "(function () { /* installs nothing */ })();"

# Installs an integration whose forceRefreshToken never settles — models a
# wedged SDK. Ignoring stray rejections means nothing else would stop this.
HANGING_CHALLENGE_JS = """
window.AwsWafIntegration = {
  forceRefreshToken: function () { return new Promise(function () {}); },
  getToken: function () { return Promise.resolve('NEVER'); },
};
"""


class _IssuerHandler(BaseHTTPRequestHandler):
    body = STUB_CHALLENGE_JS

    def do_GET(self):
        payload = self.body.encode()
        self.send_response(200)
        self.send_header("Content-Type", "text/javascript")
        self.send_header("Content-Length", str(len(payload)))
        self.end_headers()
        self.wfile.write(payload)

    def log_message(self, *args):
        pass            # keep pytest output clean


class _Issuer:
    """Serves one canned challenge.js over plain HTTP on a loopback port."""

    def __init__(self, body: str):
        handler = type("_H", (_IssuerHandler,), {"body": body})
        self.server = HTTPServer(("127.0.0.1", 0), handler)
        self.thread = threading.Thread(target=self.server.serve_forever,
                                       daemon=True)

    def __enter__(self) -> str:
        self.thread.start()
        host, port = self.server.server_address
        return f"http://{host}:{port}"

    def __exit__(self, *exc):
        self.server.shutdown()
        self.server.server_close()


def _run_harness(issuer_base: str, timeout: int = 30, env: dict | None = None) -> dict:
    import os
    proc = subprocess.run([_NODE, str(_HARNESS), issuer_base, "testkey", "UA/1.0"],
                          capture_output=True, text=True, timeout=timeout,
                          env={**os.environ, **(env or {})})
    line = (proc.stdout or "").strip().splitlines()
    assert line, f"harness produced no stdout; stderr={proc.stderr[:400]!r}"
    return json.loads(line[-1])


# --------------------------------------------------------------------------- #
# 1. The root cause: the harness must not defeat AWS's already-integrated guard #
# --------------------------------------------------------------------------- #

@needs_node
def test_mint_installs_the_real_integration_and_returns_a_token():
    """The regression. A bundle that bails when ``AwsWafIntegration`` already
    exists must still yield a token — i.e. we must not pre-stub that global."""
    with _Issuer(STUB_CHALLENGE_JS) as base:
        res = _run_harness(base)
    assert res.get("ok") is True, (
        f"mint failed against a stub issuer: {res} — the harness is most likely "
        "pre-defining window.AwsWafIntegration again (issue #41)")
    assert res.get("token") == "STUB-WAF-TOKEN"


@needs_node
def test_mint_survives_the_sdk_auto_init_unhandled_rejection():
    """Deleting the pre-stub re-exposes the SDK's own auto-init race. That
    rejection is benign and must not abort a mint that otherwise succeeds."""
    with _Issuer(STUB_CHALLENGE_JS) as base:
        res = _run_harness(base)
    assert res.get("ok") is True, (
        "the SDK's auto-init rejection aborted an otherwise-good mint")


@needs_node
def test_mint_fails_cleanly_when_the_bundle_installs_nothing():
    """A genuinely rotated/broken challenge.js must fail fast and loudly, not
    hang or emit a bogus token."""
    with _Issuer(EMPTY_CHALLENGE_JS) as base:
        res = _run_harness(base)
    assert res.get("ok") is False
    assert not res.get("token")
    assert "integration" in (res.get("error") or "")


@needs_node
def test_a_wedged_sdk_hits_the_harness_deadline_instead_of_hanging():
    """Stray rejections are ignored now, so nothing else would stop a wedged
    SDK — and the mint can sit on the on-demand quote path. The harness's own
    deadline must end it well before the Python subprocess timeout."""
    with _Issuer(HANGING_CHALLENGE_JS) as base:
        res = _run_harness(base, timeout=20,
                           env={"CZR_WAF_DEADLINE_MS": "1500"})
    assert res.get("ok") is False
    assert "deadline" in (res.get("error") or "")


def test_harness_does_not_predefine_the_aws_waf_integration():
    """Source-level guard. Re-adding the pre-stub silently re-breaks Caesars
    on the real bundle, where no test issuer can catch it."""
    src = _HARNESS.read_text()
    offenders = [
        f"{n}: {line.strip()}"
        for n, line in enumerate(src.splitlines(), 1)
        if "AwsWafIntegration" in line
        and "=" in line.split("AwsWafIntegration")[-1][:3]
        and not line.lstrip().startswith("//")
    ]
    assert not offenders, (
        "caesars_waf_node.js assigns window.AwsWafIntegration before eval'ing "
        f"challenge.js; AWS's bundle then declines to install itself: {offenders}")


# --------------------------------------------------------------------------- #
# 2. The broker: failures must be LOUD and carry node's stderr                 #
# --------------------------------------------------------------------------- #

class _Proc:
    def __init__(self, stdout="", stderr="", returncode=0):
        self.stdout, self.stderr, self.returncode = stdout, stderr, returncode


def _stub_run(monkeypatch, proc):
    calls = []

    def fake_run(cmd, **kwargs):
        calls.append((cmd, kwargs))
        if isinstance(proc, Exception):
            raise proc
        return proc

    monkeypatch.setattr(caesars_waf.subprocess, "run", fake_run)
    monkeypatch.setattr(caesars_waf.shutil, "which", lambda _: "/usr/bin/node")
    return calls


def test_mint_success_returns_token_and_a_fresh_device_id(monkeypatch):
    _stub_run(monkeypatch, _Proc(stdout='{"ok":true,"token":"T123"}\n'))
    token, device = caesars_waf.mint_token_browser_free()
    assert token == "T123"
    assert len(device) == 36 and device.count("-") == 4     # a uuid4


def test_missing_node_binary_logs_error(monkeypatch, caplog):
    monkeypatch.setattr(caesars_waf.shutil, "which", lambda _: None)
    with caplog.at_level("ERROR"):
        assert caesars_waf.mint_token_browser_free() == ("", "")
    assert any(r.levelname == "ERROR" for r in caplog.records)
    assert "node" in caplog.text


def test_nonzero_exit_logs_error_including_node_stderr(monkeypatch, caplog):
    _stub_run(monkeypatch, _Proc(stdout="", stderr="SyntaxError: boom",
                                 returncode=1))
    with caplog.at_level("ERROR"):
        assert caesars_waf.mint_token_browser_free() == ("", "")
    assert any(r.levelname == "ERROR" for r in caplog.records), (
        "a dead mint means zero Caesars rows — that is an ERROR, not a warning")
    assert "SyntaxError: boom" in caplog.text, (
        "node's stderr is the only diagnostic the operator gets (issue #41.3)")


def test_unparseable_node_output_logs_error_with_stderr(monkeypatch, caplog):
    _stub_run(monkeypatch, _Proc(stdout="not json at all",
                                 stderr="node: bad flag", returncode=0))
    with caplog.at_level("ERROR"):
        assert caesars_waf.mint_token_browser_free() == ("", "")
    assert any(r.levelname == "ERROR" for r in caplog.records)
    assert "node: bad flag" in caplog.text


def test_reported_mint_failure_logs_error_with_stderr(monkeypatch, caplog):
    """``{"ok":false}`` is the harness's own failure channel."""
    _stub_run(monkeypatch, _Proc(stdout='{"ok":false,"error":"empty token"}',
                                 stderr="[waf] challenge bailed", returncode=1))
    with caplog.at_level("ERROR"):
        assert caesars_waf.mint_token_browser_free() == ("", "")
    assert "empty token" in caplog.text
    assert "[waf] challenge bailed" in caplog.text


def test_timeout_logs_error(monkeypatch, caplog):
    _stub_run(monkeypatch,
              subprocess.TimeoutExpired(cmd="node", timeout=60))
    with caplog.at_level("ERROR"):
        assert caesars_waf.mint_token_browser_free() == ("", "")
    assert any(r.levelname == "ERROR" for r in caplog.records)


def test_empty_stdout_does_not_raise(monkeypatch, caplog):
    """``proc.stdout.strip().splitlines()[-1]`` used to IndexError on an empty
    stdout when returncode was 0."""
    _stub_run(monkeypatch, _Proc(stdout="", stderr="", returncode=0))
    with caplog.at_level("ERROR"):
        assert caesars_waf.mint_token_browser_free() == ("", "")


# --------------------------------------------------------------------------- #
# 3. Mint duration, for #50's token pre-warm                                   #
# --------------------------------------------------------------------------- #

def test_mint_reports_its_duration(monkeypatch):
    _stub_run(monkeypatch, _Proc(stdout='{"ok":true,"token":"T"}'))
    seen: list[float] = []
    caesars_waf.mint_token_browser_free(on_duration=seen.append)
    assert len(seen) == 1 and seen[0] >= 0.0


def test_mint_reports_duration_on_failure_too(monkeypatch):
    """#50 needs the cost of a FAILED mint too — that is the latency an RFQ
    actually pays when Caesars is down."""
    _stub_run(monkeypatch, _Proc(stdout='{"ok":false,"error":"x"}',
                                 returncode=1))
    seen: list[float] = []
    caesars_waf.mint_token_browser_free(on_duration=seen.append)
    assert len(seen) == 1


def test_on_duration_is_keyword_only():
    sig = inspect.signature(caesars_waf.mint_token_browser_free)
    with pytest.raises(TypeError):
        sig.bind("nj", "issuer", "key", "ua", None, 60, print)
    sig.bind(on_duration=print)


def test_a_broken_duration_callback_never_breaks_the_mint(monkeypatch):
    _stub_run(monkeypatch, _Proc(stdout='{"ok":true,"token":"T"}'))

    def boom(_):
        raise RuntimeError("measurement blew up")

    assert caesars_waf.mint_token_browser_free(on_duration=boom)[0] == "T"


# --------------------------------------------------------------------------- #
# 4. Token lifecycle invariants the ticket asks to pin                         #
# --------------------------------------------------------------------------- #

def test_client_exposes_the_last_mint_duration(monkeypatch):
    """#50 reads this off the live client to size the token pre-warm."""
    from mlb_sgp.caesars_client import CaesarsClient

    client = CaesarsClient.__new__(CaesarsClient)      # no session/network
    client.state = "nj"
    client.last_mint_sec = 0.0
    client._token = client._device = ""
    client._cookies = {}
    client._minted_at = 0.0
    def fake_mint(*, on_duration, **kwargs):
        on_duration(1.25)
        return "", ""

    monkeypatch.setattr(caesars_waf, "mint_token_browser_free", fake_mint)
    monkeypatch.setattr(CaesarsClient, "_validate", lambda self: False)
    client._mint()
    assert client.last_mint_sec == 1.25


def test_token_ttl_stays_short():
    """The WAF token outlives its TTL server-side, but a long TTL means a
    cycle spends its budget on 403 retries. 240s is the tuned value."""
    from mlb_sgp.caesars_client import TOKEN_TTL_SEC
    assert TOKEN_TTL_SEC == 240


def test_a_cached_token_is_validated_before_use():
    """A cached-but-rejected token must not be handed out: ``_load_cache``
    validates before returning True."""
    from mlb_sgp import caesars_client
    src = inspect.getsource(caesars_client.CaesarsClient._load_cache)
    assert "_validate()" in src


def test_auth_failure_raises_rather_than_yielding_an_empty_slate():
    """Issue #33's contract, re-pinned here: a book that cannot authenticate
    must raise, not return [] — an empty list reads as 'no markets today' and
    resets the failure counter."""
    from mlb_sgp import caesars_client
    src = inspect.getsource(caesars_client.CaesarsClient.list_events)
    assert "BookTransportError" in src and "auth" in src
