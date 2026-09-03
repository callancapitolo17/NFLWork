# DraftKings price sidecar

A small local service that owns a real Chrome and makes DraftKings' SGP
price call (`calculateBets`) on behalf of the bots. Exists because of issue
#102: since ~2026-08-20 Akamai denies that endpoint to anything that is not
a real, non-headless browser making one call at a time. Nothing at the HTTP
layer fixes it (see `mlb_sgp/README.md` § DraftKings price host for the
controlled runs), so the browser has to be real — and this process keeps it
out of the bots.

```
maker / taker ──POST /price──▶ dk_price_sidecar ──in-page fetch──▶ DK
   (plain HTTP, loopback)        (minimized Chrome, serial)
```

## Run

```bash
dk_price_sidecar/run.sh
```

Needs a Python with Playwright and a Google Chrome install. The repo's
default `python3` has neither; `run.sh` defaults to the Framework 3.12
interpreter and `DK_SIDECAR_PYTHON` overrides it.

Then flip the bot's transport (config only, default `http` = the old path):

```bash
export DK_PRICE_TRANSPORT=sidecar        # DK_SIDECAR_URL defaults to http://127.0.0.1:8095
```

Rollback is unsetting `DK_PRICE_TRANSPORT`. The five other books are
untouched either way.

## Interface

| | |
|---|---|
| `POST /price` `{"selections": ["<dk id>", ...]}` | DK's own status: **200** with `true_odds` (correlated decimal for the full set) / `422` declined / `403` blocked. `true_odds: null` on a 200 means DK returned singles only or a `combinabilityRestrictions` entry — a decline, not a price. |
| `GET /health` | `{ok, calls, last_status, relaunches, uptime_sec}` |
| `502` | the in-page fetch threw (no response reached the network layer) |
| `503` | the browser is gone; the NEXT call relaunches it |

The bot maps these exactly as it maps curl_cffi's: float / None / counted
transport error (`mlb_sgp/draftkings.py::_price_via_sidecar`).

## Config

| env | default | meaning |
|---|---|---|
| `DK_SIDECAR_PORT` | `8095` | loopback only, never bound off-machine |
| `DK_SIDECAR_MIN_INTERVAL_SEC` | `1.0` | floor between DK calls. Measured: 1.5s → 100 % pass, 4 concurrent → 5 %. Untested below 1.0 |
| `DK_SIDECAR_PAGE_RELOAD_SEC` | `600` | reload the sportsbook page to keep Akamai cookies fresh |
| `DK_SIDECAR_PROFILE_DIR` | `~/.dk_price_sidecar/profile` | persistent Chrome profile |
| `DK_SIDECAR_MINIMIZE` | `1` | minimize the window via CDP (verified 3/3 priced while minimized) |

Bot side: `DK_PRICE_TRANSPORT=http|sidecar`, `DK_SIDECAR_URL`,
`DK_SIDECAR_TIMEOUT_SEC` (10 — a 4-cell partition queued behind another
flight legitimately takes several seconds).

## Why it is shaped this way

- **Single-threaded `HTTPServer` on purpose.** Concurrent bot requests queue
  in the socket backlog and reach DK one at a time. The no-burst rule is
  enforced where the bots cannot bypass it — issue #93's concurrent
  partition cells are what drew Akamai's rule in the first place.
- **Not headless.** Headless Chrome is exactly what the rule detects (0/3,
  twice). The window is minimized instead; "headed" is a process mode, not
  a screen requirement. On a Linux host, run it under `Xvfb`.
- **Warm-up absorbs the first-call throw.** The first fetch on a fresh page
  throws deterministically; `launch()` spends that one so no bot request
  pays for it.
- **Relaunch on the next call, not in-line.** A Playwright failure closes
  the browser and answers 503; the following request relaunches. The bot
  counts the 503 as a transport error (#102's first fix), so a crashed
  sidecar reads as "DK down", never as "DK's parser broke".

## Measured (2026-09-02, pregame Yankees @ Angels, home −1.5 + O7.5)

Through the maker's real `price_selection_set` with
`DK_PRICE_TRANSPORT=sidecar`: `true_odds=7.75` 4/4, p50 1.40s (the 1.0s
interval floor; first call 0.59s), 0 transport errors, 0 relaunches. A
4-cell partition is therefore ~5.6s serial — inside the 8–10s flight budget.

## Operational notes

- The machine must not sleep during a slate (see the MM overnight-sleep
  memory: macOS maintenance sleep freezes threads).
- `x-client-version` / `x-client-widget-version` are DK's betslip build
  numbers, hardcoded from the 2026-09-02 capture. DK rolls them; if prices
  stop while `/health` is green, re-capture them from DK's own request.
- Not started by any bot. Run it yourself, before the maker, like the
  monitor.
