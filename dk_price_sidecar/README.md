# DraftKings price sidecar

A small local service that owns a real Chrome and makes DraftKings' SGP
price call (`calculateBets`) for local tools — today `nfl_specials/` (the
Wagerzon NFL trifecta/superfecta pricer). Exists because of issue #102:
since ~2026-08-20 Akamai denies that endpoint to anything that is not a
real, non-headless browser making one call at a time. Nothing at the HTTP
layer fixes it (see `mlb_sgp/README.md` § DraftKings price host for the
controlled runs), so the browser has to be real — and this process keeps it
out of the callers.

```
nfl_specials ──POST /price──▶ dk_price_sidecar ──in-page fetch──▶ DK
 (plain HTTP, loopback)         (minimized Chrome, serial)
```

The MLB bots do not use it: their wiring was reverted on 2026-09-02 by owner
decision (DK stays off same-game there).

## Run

```bash
dk_price_sidecar/run.sh
```

Run it from your own terminal. It opens a real Chrome window (minimized),
which a sandboxed process cannot draw — launched from a sandbox the page
loads but every price call hangs. Needs a Python with Playwright and a
Google Chrome install; `run.sh` defaults to the Framework 3.12 interpreter
and `DK_SIDECAR_PYTHON` overrides it.

## Interface

| | |
|---|---|
| `POST /price` `{"selections": ["<dk id>", ...]}` | DK's own status: **200** with `true_odds` (correlated decimal for the full set) / `422` declined / `403` blocked. `true_odds: null` on a 200 means DK returned singles only or a `combinabilityRestrictions` entry — a decline, not a price. |
| `GET /health` | `{ok, calls, last_status, relaunches, uptime_sec}` |
| `502` | the in-page fetch threw (no response reached the network layer) |
| `503` | the browser is gone; the NEXT call relaunches it |

## Config

| env | default | meaning |
|---|---|---|
| `DK_SIDECAR_PORT` | `8095` | loopback only, never bound off-machine |
| `DK_SIDECAR_MIN_INTERVAL_SEC` | `1.0` | floor between DK calls. Measured: 1.5s → 100 % pass, 4 concurrent → 5 %. Untested below 1.0 |
| `DK_SIDECAR_PAGE_RELOAD_SEC` | `600` | reload the sportsbook page to keep Akamai cookies fresh |
| `DK_SIDECAR_PROFILE_DIR` | `~/.dk_price_sidecar/profile` | persistent Chrome profile |
| `DK_SIDECAR_MINIMIZE` | `1` | minimize the window via CDP (verified 3/3 priced while minimized) |

## Why it is shaped this way

- **Single-threaded `HTTPServer` on purpose.** Concurrent requests queue in
  the socket backlog and reach DK one at a time. The no-burst rule is
  enforced where callers cannot bypass it — issue #93's concurrent partition
  cells are what drew Akamai's rule in the first place.
- **Not headless.** Headless Chrome is exactly what the rule detects (0/3,
  twice). The window is minimized instead.
- **Not Linux Chromium either.** Measured 2026-10-02: Playwright's bundled
  Chromium, headed under Xvfb in a Linux container, gets the preflight (200)
  and then Akamai's denial on the POST (no CORS header → `Failed to fetch`).
  Only a real Google Chrome has been seen to pass.
- **Warm-up absorbs the first-call throw.** The first fetch on a fresh page
  throws deterministically; `launch()` spends that one so no caller pays.
- **Relaunch on the next call, not in-line.** A Playwright failure closes
  the browser and answers 503; the following request relaunches.

## Measured (2026-09-02, pregame Yankees @ Angels, home −1.5 + O7.5)

`true_odds=7.75` 4/4, p50 1.40s (the 1.0s interval floor; first call
0.59s), 0 transport errors, 0 relaunches.

## Operational notes

- The machine must not sleep during a slate (macOS maintenance sleep
  freezes threads).
- `x-client-version` / `x-client-widget-version` are DK's betslip build
  numbers, hardcoded from the 2026-09-02 capture. DK rolls them; if prices
  stop while `/health` is green, re-capture them from DK's own request.
