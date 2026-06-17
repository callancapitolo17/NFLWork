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

VERIFIED 2026-06-16: token_len=182, tabs feed -> 200 / 210 KB JSON. GREEN.
"""
from __future__ import annotations

import json
import shutil
import subprocess
import uuid
from pathlib import Path

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
                            verbose: bool = False) -> tuple[str, str]:
    """Mint an aws-waf-token without a browser. Returns (token, device_id);
    ("", "") on failure (caller then produces no rows, same as Playwright path).
    `state` is accepted for signature parity with CaesarsClient._mint."""
    node = node_bin or shutil.which("node")
    if not node:
        if verbose:
            print("  [czr-waf] node not found on PATH")
        return "", ""
    if not _NODE_HARNESS.exists():
        if verbose:
            print(f"  [czr-waf] missing harness {_NODE_HARNESS}")
        return "", ""
    try:
        proc = subprocess.run(
            [node, str(_NODE_HARNESS), issuer, issuer_key, ua],
            capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        if verbose:
            print("  [czr-waf] node mint timed out")
        return "", ""
    line = (proc.stdout or "").strip().splitlines()[-1] if proc.stdout.strip() else ""
    try:
        res = json.loads(line)
    except Exception:
        if verbose:
            print(f"  [czr-waf] bad node output: {proc.stdout[:200]!r} "
                  f"stderr={proc.stderr[:200]!r}")
        return "", ""
    if not res.get("ok") or not res.get("token"):
        if verbose:
            print(f"  [czr-waf] mint failed: {res.get('error')}")
        return "", ""
    return res["token"], str(uuid.uuid4())


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
    token, device = mint_token_browser_free(verbose=True)
    dt = time.time() - t0
    print(f"token_len={len(token)} device={device} mint={dt:.2f}s")
    if not token:
        print("MINT FAILED"); raise SystemExit(1)
    st, ln, is_json = _validate(token, device)
    print(f"tabs validate: status={st} bytes={ln} json={is_json}")
    print("VERDICT:", "GREEN — browser-free token VALIDATED"
          if (st == 200 and is_json)
          else "token minted but NOT validated (WAF reject or IP throttle)")
