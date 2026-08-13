#!/usr/bin/env python3
"""
Browser-free AWS WAF token mint for Caesars (api.americanwagering.com).

Replaces a headless-Playwright mint with a far lighter node-subprocess that
runs AWS WAF's REAL challenge.js (no Chromium, no Playwright — just `node` on
PATH + curl_cffi for validation). Used by caesars_client.CaesarsClient._mint.

WHY NOT pure-Python PoW:
  Caesars' WebACL always issues the `NetworkBandwidth` challenge type
  (verified 12/12 /inputs calls). That is NOT a hashable proof-of-work — the
  SDK uploads timed multipart blobs to `/mp_verify` and the server measures
  bandwidth. Reproducing that blind is brittle. Instead we run the SDK's own
  code under node with minimal DOM shims; node supplies fetch / crypto.subtle
  / Blob / FormData / performance natively, and we shim only what it probes
  (document.scripts, navigator/screen fingerprint, FileReader for the blob).

FLOW (caesars_waf_node.js):
  1. fetch challenge.js from the AWS WAF edge issuer host
  2. eval it under node -> window.AwsWafIntegration
  3. forceRefreshToken()  -> SDK does GET /inputs + multipart POST /mp_verify
  4. getToken()           -> the aws-waf-token
  5. emit {"ok":true,"token":...} on stdout

The device id is a fresh UUID — CloudFront gates on the WAF token, not the
device id (the SPA just needs *some* uuid in x-unique-device-id).

DO NOT pre-define window.AwsWafIntegration in the harness (issue #41). AWS's
bundle installs itself only when that global is absent, so a placeholder makes
the SDK skip its own assignment; getToken() then returns the placeholder's
empty string and Caesars silently drops to zero rows. That is exactly how this
book died between 2026-06-24 and 2026-07-30. The auto-init rejection the
placeholder used to hide is now logged and ignored instead.

VERIFIED 2026-06-16: token_len=182, tabs feed -> 200 / 210 KB JSON. GREEN.
RE-VERIFIED 2026-07-30 after the issue #41 fix: mint 0.68s, token_len=250
(the token format rotated too), tabs feed -> 200 / 178 KB JSON. GREEN.
"""
from __future__ import annotations

import json
import logging
import shutil
import subprocess
import time
import uuid
from collections.abc import Callable
from pathlib import Path

logger = logging.getLogger(__name__)

_THIS_DIR = Path(__file__).resolve().parent
_NODE_HARNESS = _THIS_DIR / "caesars_waf_node.js"

DEFAULT_ISSUER = "4ad3fec456d9.edge.sdk.awswaf.com"
DEFAULT_KEY = "4ad3fec456d9"
REAL_UA = ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
           "(KHTML, like Gecko) Chrome/149.0.0.0 Safari/537.36")


def mint_token_browser_free(state: str = "nj",
                            issuer: str = DEFAULT_ISSUER,
                            issuer_key: str = DEFAULT_KEY,
                            ua: str = REAL_UA,
                            node_bin: str | None = None,
                            timeout: int = 60,
                            *,
                            on_duration: Callable[[float], None] | None = None
                            ) -> tuple[str, str]:
    """Mint an aws-waf-token without a browser. Returns (token, device_id);
    ("", "") on failure (caller then produces no rows, same as Playwright path).
    `state` is accepted for signature parity with CaesarsClient._mint.

    Every failure path logs at ERROR and includes node's stderr (issue #41): a
    mint that cannot complete means the whole book produces zero rows, which is
    indistinguishable from an empty slate unless it is loud. `on_duration`
    (keyword-only) receives the wall-clock seconds the mint cost, on success
    AND on failure — a failed mint is latency an RFQ still pays, which is what
    #50's token pre-warm needs to measure.
    """
    started = time.monotonic()

    def _done(token: str, device: str) -> tuple[str, str]:
        if on_duration is not None:
            try:
                on_duration(time.monotonic() - started)
            except Exception:       # measurement must never break the mint
                logger.debug("caesars-waf: on_duration callback raised",
                             exc_info=True)
        return token, device

    node = node_bin or shutil.which("node")
    if not node:
        # ERROR, not DEBUG: with no node there is no Caesars token and so no
        # Caesars data at all — an environment fault, not a slate condition.
        logger.error("caesars-waf: node not found on PATH")
        return _done("", "")
    if not _NODE_HARNESS.exists():
        logger.error("caesars-waf: missing node harness %s", _NODE_HARNESS)
        return _done("", "")
    try:
        proc = subprocess.run(
            [node, str(_NODE_HARNESS), issuer, issuer_key, ua],
            capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        logger.error("caesars-waf: node mint timed out after %ss", timeout)
        return _done("", "")

    stdout = proc.stdout or ""
    stderr = (proc.stderr or "")[:500]
    lines = stdout.strip().splitlines()
    if not lines:
        logger.error("caesars-waf: node mint produced no output "
                     "(exit=%s) stderr=%r", proc.returncode, stderr)
        return _done("", "")
    try:
        res = json.loads(lines[-1])
    except Exception:
        logger.error("caesars-waf: unparseable node output (exit=%s): %r "
                     "stderr=%r", proc.returncode, stdout[:200], stderr)
        return _done("", "")
    if not res.get("ok") or not res.get("token"):
        logger.error("caesars-waf: mint failed: %s (exit=%s) stderr=%r",
                     res.get("error"), proc.returncode, stderr)
        return _done("", "")
    return _done(res["token"], str(uuid.uuid4()))


# --- self-test --------------------------------------------------------------- #
def _validate(token: str, device: str, state: str = "nj") -> tuple[int, int, bool]:
    from curl_cffi import requests
    API = "https://api.americanwagering.com"
    tab = (f"{API}/regions/us/locations/{state}/brands/czr/sb/"
           f"v4/sports/baseball/tabs/SCHEDULE%7CGames%20%E2%9A%BE")
    h = {"User-Agent": REAL_UA, "Accept": "application/json, text/plain, */*",
         "Origin": "https://sportsbook.caesars.com",
         "Referer": "https://sportsbook.caesars.com/",
         "X-Aws-Waf-Token": token, "x-platform": "cordova-desktop",
         "x-app-version": "7.49.0", "x-unique-device-id": device}
    r = requests.get(tab, headers=h, cookies={"aws-waf-token": token},
                     impersonate="chrome", timeout=25)
    return r.status_code, len(r.text), r.text.strip().startswith("{")


if __name__ == "__main__":
    import time
    t0 = time.time()
    token, device = mint_token_browser_free()
    dt = time.time() - t0
    print(f"token_len={len(token)} device={device} mint={dt:.2f}s")
    if not token:
        print("MINT FAILED"); raise SystemExit(1)
    st, ln, is_json = _validate(token, device)
    print(f"tabs validate: status={st} bytes={ln} json={is_json}")
    print("VERDICT:", "GREEN — browser-free token VALIDATED"
          if (st == 200 and is_json)
          else "token minted but NOT validated (WAF reject or IP throttle)")
