# MLB Bets Tab — PR C (editable risk + WZ verified quote) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the bet card's hero-strip `Risk $XX` value click-to-edit. On commit, fire a Wagerzon `ConfirmWagerHelper` preflight to get the *actual* `Win` and *actual* current odds, then swap the local-math "to win" for the verified number. Surfaces line drift and other WZ rejections (integer-only amount, balance, line pulled) before the user clicks Place. Also fixes a latent bug where the existing placement flow conflates `actual_size` and `kelly_bet` into one `data-size` attribute.

**Architecture:** New Python module `wagerzon_odds/single_pricer.py` (mirror of `parlay_pricer.get_parlay_price` but with `WT="0"` for singles + `RiskWin="2"` for the no-balance-check preview). New Flask endpoint `POST /api/wz-quote-single` in `mlb_dashboard_server.py` that wraps `single_pricer.get_single_price`. Client side: `mlb_dashboard.R` gains click-to-edit on the Risk cell, splits the placement button's `data-size` (editable override) from a new `data-model-size` (immutable Kelly result), and a small JS coordinator that fires `/api/wz-quote-single` on commit and swaps the To Win value (with a verified-✓ badge for WZ rows or a clear error pill for drift / rejections).

**Tech Stack:** Python 3 + Flask (existing); requests + pytest with MagicMock for the new pricer unit test (matches `test_parlay_placer.py` pattern); plain DOM JavaScript with `contentEditable` for click-to-edit.

**Worktree:** Before Task 1, create a fresh worktree off `main` named `mlb-bets-tab-pr-c-editable-risk` (e.g. via `EnterWorktree({name: "mlb-bets-tab-pr-c-editable-risk"})` or `git worktree add .claude/worktrees/mlb-bets-tab-pr-c-editable-risk -b worktree-mlb-bets-tab-pr-c-editable-risk main`). All paths assume that worktree as the working directory. **PR A and PR B should already be merged to `main`** before this PR's manual verification step (Task 8).

---

## File Structure

```
wagerzon_odds/
├── single_pricer.py                          (CREATE: ConfirmWagerHelper preview)
├── tests/
│   └── test_single_pricer.py                 (CREATE: pytest with mocked session)
└── CLAUDE.md                                 (modify: Quick map entry)
Answer Keys/MLB Dashboard/
├── mlb_dashboard_server.py                   (modify: add /api/wz-quote-single)
├── mlb_dashboard.R                           (modify: 4 sections — see below)
└── README.md                                 (modify: feature list)
```

Each file's responsibility:

- **`wagerzon_odds/single_pricer.py`** (new) — single-bet preview helper. Mirrors `parlay_pricer.get_parlay_price` but with `WT="0"` (single) and a single-leg `sel`. Calls `ConfirmWagerHelper` with `RiskWin="2"` (skip balance validation). Returns `{win, current_wz_odds, error_msg, error_msg_key}`. ~70 lines including docstring.
- **`wagerzon_odds/tests/test_single_pricer.py`** (new) — pytest unit test with `requests.Session` mocked via `unittest.mock.MagicMock`. Mirrors `test_parlay_placer.py` patterns. Covers: happy path (returns `Win` from response), `ErrorMsg` path (rejection), HTML response (auth_error), network exception (network_error).
- **`mlb_dashboard_server.py`** — new endpoint `POST /api/wz-quote-single`. Body: `{bet_hash, amount, account}`. Looks up `idgm/play/line/odds` via the existing `_resolve_wagerzon_play_idgm` helper, calls `single_pricer.get_single_price`, returns parsed JSON.
- **`mlb_dashboard.R`** — four sections:
  1. The data-attribute builder (around line 1521-1529) splits `data-size` (editable) from a new `data-model-size` (immutable Kelly).
  2. `placeBet` JS (line 3790) reads `data.modelSize` for `kelly_bet` and `data.size` for `actual_size`.
  3. `render_hero_strip` makes the Risk value a click-to-edit `<span class="risk-value" contenteditable="false">` with a hidden `↶ snap-back` button and a hidden inline error pill, plus a verified-`✓` badge slot on To Win.
  4. New page-level `<script>` block with the click-to-edit coordinator, fire-quote handler, response-swap logic, snap-back handler, and the `MLB_QUOTE_SINGLE_URL` constant.
- **`wagerzon_odds/CLAUDE.md`** — add `single_pricer.py` to the "Quick map" list, mirroring the `parlay_pricer.py` entry.
- **`MLB Dashboard/README.md`** — add the editable-risk + verified-quote flow to the placing-bets section.

---

## Task 1 — `wagerzon_odds/single_pricer.py` (TDD)

**Files:**
- Create: `wagerzon_odds/single_pricer.py`
- Create: `wagerzon_odds/tests/test_single_pricer.py`

- [ ] **Step 1: Confirm the existing test pattern**

Look at `wagerzon_odds/tests/test_parlay_placer.py` for the MagicMock session pattern. We'll mirror it. (Reading-only step — no changes.)

- [ ] **Step 2: Write the failing tests**

Create `wagerzon_odds/tests/test_single_pricer.py`:

```python
"""Unit tests for wagerzon_odds.single_pricer.get_single_price.

Mirrors test_parlay_placer.py pattern: pass a MagicMock session into
get_single_price so we never hit the live WZ API.
"""

import json
from unittest.mock import MagicMock

import pytest

from wagerzon_odds.single_pricer import get_single_price


def _mock_response(json_body=None, status=200, content_type="application/json", raise_exc=None):
    """Build a MagicMock requests.Response with the given JSON body."""
    if raise_exc is not None:
        sess = MagicMock()
        sess.post.side_effect = raise_exc
        return sess
    resp = MagicMock()
    resp.status_code = status
    resp.headers = {"content-type": content_type}
    resp.json.return_value = json_body or {}
    sess = MagicMock()
    sess.post.return_value = resp
    return sess


def _wz_confirm_body(win, odds, error_key=None, error_msg=None):
    """Build a fake ConfirmWagerHelper response shape for a single."""
    body = {
        "result": {
            "details": [
                {
                    "Amount": 25.0,
                    "Risk": 25.0,
                    "Win": win,
                    "details": [{"Odds": odds, "IsMLine": True}],
                }
            ],
            "ErrorMsgKey": error_key or "",
            "ErrorMsg":    error_msg or "",
        }
    }
    return body


@pytest.fixture
def bet():
    return {
        "idgm": 5632938,
        "play": 5,            # home ML per WZ play codes
        "line": 0.0,
        "american_odds": -140,
        "amount": 25.0,
        "pitcher": 0,
    }


def test_happy_path_returns_win_and_current_odds(bet):
    sess = _mock_response(_wz_confirm_body(win=17.86, odds=-140))
    out = get_single_price(sess, bet, amount=25)
    assert out["win"] == pytest.approx(17.86)
    assert out["current_wz_odds"] == -140
    assert out["error_msg"] == ""
    assert out["error_msg_key"] == ""


def test_error_msg_surfaced(bet):
    sess = _mock_response(_wz_confirm_body(
        win=0, odds=0,
        error_key="MAXMONEYLINEEXCEEDED",
        error_msg="Maximum money line risk exceeded",
    ))
    out = get_single_price(sess, bet, amount=10000)
    assert out["error_msg_key"] == "MAXMONEYLINEEXCEEDED"
    assert out["error_msg"] == "Maximum money line risk exceeded"
    assert out["win"] is None


def test_html_response_returns_auth_error(bet):
    sess = _mock_response(json_body={}, content_type="text/html")
    out = get_single_price(sess, bet, amount=25)
    assert out["error_msg_key"] == "session_expired"


def test_network_exception_returns_network_error(bet):
    import requests
    sess = _mock_response(raise_exc=requests.RequestException("connreset"))
    out = get_single_price(sess, bet, amount=25)
    assert out["error_msg_key"] == "network_error"
    assert "connreset" in out["error_msg"]
```

- [ ] **Step 3: Run tests; confirm they fail**

```bash
cd "/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-bets-tab-pr-c-editable-risk" && \
  python -m pytest wagerzon_odds/tests/test_single_pricer.py -v
```

Expected: ImportError or ModuleNotFoundError on `from wagerzon_odds.single_pricer import get_single_price`.

- [ ] **Step 4: Implement `wagerzon_odds/single_pricer.py`**

Create `wagerzon_odds/single_pricer.py`:

```python
"""Wagerzon single-bet preview / pricer.

Mirror of `wagerzon_odds.parlay_pricer.get_parlay_price` for single
straight bets. Calls ConfirmWagerHelper with `RiskWin="2"` so WZ
returns a price quote without balance validation. Used by the MLB
Dashboard's editable-Risk feature to:

  1. Verify the actual `Win` (to-win) WZ would credit at any user-typed
     amount, instead of relying on local American-odds math.
  2. Detect line drift before the user clicks Place — the response
     includes WZ's *current* odds for the leg.
  3. Surface other WZ rejections (integer-only amounts, MAXRISK, line
     pulled, etc.) directly via `ErrorMsg` / `ErrorMsgKey`.

Does NOT place the bet. For actual placement, see `single_placer.py`.
"""

import json as _json

import requests

from wagerzon_odds.config import WAGERZON_BASE_URL
from wagerzon_odds.single_placer import build_sel_for_single

CONFIRM_URL = f"{WAGERZON_BASE_URL}/wager/ConfirmWagerHelper.aspx"


def _network_error(msg: str) -> dict:
    return {
        "win": None,
        "current_wz_odds": None,
        "error_msg": msg,
        "error_msg_key": "network_error",
    }


def _auth_error(msg: str = "Wagerzon returned HTML at preflight (session expired)") -> dict:
    return {
        "win": None,
        "current_wz_odds": None,
        "error_msg": msg,
        "error_msg_key": "session_expired",
    }


def get_single_price(session: requests.Session, bet: dict, amount: float) -> dict:
    """Get WZ's current price for a single bet at the given amount.

    Args:
        session: An authenticated requests.Session for WZ. In tests, pass a
            MagicMock — the function only exercises `.post`, `.json()`,
            `.headers`, `.status_code`.
        bet: Dict with keys `idgm, play, line, american_odds, pitcher`
            (pitcher optional, defaults to 0).
        amount: Risk amount to query at. Passed verbatim to WZ — if WZ
            rejects (e.g. integer-only rule), the rejection surfaces in
            the returned `error_msg_key`.

    Returns:
        Dict with `{win, current_wz_odds, error_msg, error_msg_key}`.
        `win` is the to-win amount WZ would credit; `current_wz_odds` is
        the leg's current American odds at WZ. Both None on failure.
    """
    # ConfirmWagerHelper payload mirrors single_placer._build_confirm_payload
    # but with RiskWin="2" (preview, no balance check) instead of "0".
    amount_str = str(int(amount)) if amount == int(amount) else str(amount)
    detail_data = [
        {
            "Amount": amount_str,
            "RiskWin": "2",
            "TeaserPointsPurchased": 0,
            "IdGame": bet["idgm"],
            "Play": bet["play"],
            "Pitcher": bet.get("pitcher", 0),
            "Points": {
                "BuyPoints": 0,
                "BuyPointsDesc": "",
                "LineDesc": "",
                "selected": True,
            },
        }
    ]
    payload = {
        "IDWT": "0",
        "WT": "0",
        "amountType": "0",
        "open": "0",
        "sameAmount": "false",
        "sameAmountNumber": amount_str,
        "useFreePlayAmount": "false",
        "sel": build_sel_for_single(bet),
        "detailData": _json.dumps(detail_data),
    }

    try:
        resp = session.post(
            CONFIRM_URL,
            data=payload,
            timeout=15,
            headers={"Accept": "application/json"},
        )
    except requests.RequestException as e:
        return _network_error(f"{type(e).__name__}: {e}")

    if "json" not in resp.headers.get("content-type", ""):
        return _auth_error()

    try:
        body = resp.json()
    except ValueError as e:
        return _network_error(f"json decode failed: {e}")

    result = body.get("result") or {}
    err_key = result.get("ErrorMsgKey") or result.get("ErrorCode") or ""
    err_msg = result.get("ErrorMsg") or result.get("ErrorMessage") or ""
    if err_key:
        return {
            "win": None,
            "current_wz_odds": None,
            "error_msg": err_msg or err_key,
            "error_msg_key": err_key,
        }

    details = result.get("details") or []
    if not details:
        return {
            "win": None,
            "current_wz_odds": None,
            "error_msg": "Wagerzon returned empty details (line pulled?)",
            "error_msg_key": "empty_details",
        }
    outer = details[0]
    legs = outer.get("details") or []
    win = outer.get("Win")
    odds_now = legs[0].get("Odds") if legs else None
    return {
        "win": win,
        "current_wz_odds": int(odds_now) if odds_now is not None else None,
        "error_msg": "",
        "error_msg_key": "",
    }
```

- [ ] **Step 5: Re-run tests; confirm green**

```bash
cd "/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-bets-tab-pr-c-editable-risk" && \
  python -m pytest wagerzon_odds/tests/test_single_pricer.py -v
```

Expected: 4 passed.

- [ ] **Step 6: Commit**

```bash
cd "/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-bets-tab-pr-c-editable-risk" && \
  git add wagerzon_odds/single_pricer.py wagerzon_odds/tests/test_single_pricer.py && \
  git commit -m "$(cat <<'EOF'
feat(wagerzon): add single_pricer.get_single_price for ConfirmWagerHelper preview

Mirror of parlay_pricer.get_parlay_price but for single straight bets
(WT=0). Uses RiskWin=2 to skip balance validation and return WZ's
quote at any amount — Win, current odds, ErrorMsg/Key. Used by the MLB
Dashboard's editable-Risk feature to verify to-win and detect line
drift before placement. Pure function; tests cover happy path, error
surface, HTML/auth response, and network exception.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 2 — `/api/wz-quote-single` Flask endpoint

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard_server.py`

- [ ] **Step 1: Add the endpoint**

Open `Answer Keys/MLB Dashboard/mlb_dashboard_server.py`. At the top of the file, add the import (near the other `from wagerzon_odds import ...` lines — search for `single_placer` import to find them):

```python
from wagerzon_odds import single_pricer
```

Then add a new endpoint. A natural location is immediately after the `/api/place-bet` route (search for `def place_bet():` to find it, scroll to the end of that function). Insert:

```python
@app.route("/api/wz-quote-single", methods=["POST"])
def wz_quote_single():
    """Preview a Wagerzon single bet at a user-typed amount.

    Calls ConfirmWagerHelper with RiskWin=2 so WZ returns the current
    quote without doing a balance check. Used by the bets-tab editable-
    Risk feature to swap the local-math To Win for the actual Win the
    user would receive, and to detect line drift / amount-rule
    rejections before the user clicks Place.

    Body:    {bet_hash, amount, account}
    Returns: {win, current_wz_odds, error_msg, error_msg_key}
             plus echo of {bet_hash, amount} so the client can correlate
             responses to their originating field.
    """
    data = request.json or {}
    bet_hash = data.get("bet_hash")
    amount   = data.get("amount")
    account  = data.get("account")

    if not bet_hash:
        return jsonify({"error_msg_key": "bad_request",
                        "error_msg": "bet_hash required"}), 400
    if amount is None:
        return jsonify({"error_msg_key": "bad_request",
                        "error_msg": "amount required"}), 400
    try:
        amount_f = float(amount)
    except (TypeError, ValueError):
        return jsonify({"error_msg_key": "bad_request",
                        "error_msg": "amount must be numeric"}), 400
    if not account:
        return jsonify({"error_msg_key": "bad_request",
                        "error_msg": "account required"}), 400

    # The dashboard caller passes the same bet payload it would send to
    # /api/place-bet (idgm + play already resolved client-side via the
    # placement modal flow). If those keys are missing, fall back to the
    # idgm/play resolver used by /api/place-bet.
    if "idgm" not in data or "play" not in data:
        wz_play = _resolve_wagerzon_play_idgm(data)
        if wz_play is None:
            return jsonify({"error_msg_key": "not_found",
                            "error_msg": "Could not find this bet in wagerzon_odds"}), 400
        data["idgm"] = wz_play["idgm"]
        data["play"] = wz_play["play"]

    bet_for_pricer = {
        "idgm":          data["idgm"],
        "play":          data["play"],
        "line":          data.get("line", 0.0) or 0.0,
        "american_odds": data.get("american_odds"),
        "pitcher":       data.get("pitcher", 0),
    }

    wz_account = get_account(account)
    sess = wagerzon_auth.get_session(wz_account)
    result = single_pricer.get_single_price(sess, bet_for_pricer, amount=amount_f)
    result["bet_hash"] = bet_hash
    result["amount"]   = amount_f
    return jsonify(result)
```

> If `get_account` and `wagerzon_auth` aren't already imported at the top of the file, add them — search for `from wagerzon_odds.wagerzon_accounts import` to find existing imports.

- [ ] **Step 2: Smoke-test the endpoint locally**

This requires a live WZ account configured in `bet_logger/.env`. Skip if unavailable; do the integration verification in Task 8 instead.

If you have a live account:
```bash
# Pick a bet_hash from the dashboard for a known WZ-pickable bet, then:
curl -s -X POST http://localhost:8083/api/wz-quote-single \
  -H 'Content-Type: application/json' \
  -d '{"bet_hash":"<hash>","amount":1,"account":"primary"}' | jq .
```
Expected: a JSON object with non-null `win` and `current_wz_odds`. If you get `error_msg_key: "session_expired"`, refresh WZ auth via the existing dashboard reauth flow.

- [ ] **Step 3: Commit**

```bash
cd "/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-bets-tab-pr-c-editable-risk" && \
  git add "Answer Keys/MLB Dashboard/mlb_dashboard_server.py" && \
  git commit -m "$(cat <<'EOF'
feat(server): /api/wz-quote-single — WZ ConfirmWagerHelper preview endpoint

Wraps wagerzon_odds.single_pricer.get_single_price for the bets-tab
editable-Risk flow. Resolves idgm/play via the same _resolve helper
used by /api/place-bet. Returns the parsed pricer result plus echoes
of bet_hash/amount so the client can correlate responses.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3 — Split `data-size` / `data-model-size` in placement attributes

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R` — the `data_attrs` builder and `placeBet` JS

**Why now (before the editable UI).** The fix is independent and tiny. Doing it first keeps each subsequent commit focused on UI without bundling backend semantics changes.

- [ ] **Step 1: Update the `data_attrs` builder in `create_bets_table`**

Open `Answer Keys/MLB Dashboard/mlb_dashboard.R`. Locate the `data_attrs` sprintf around lines 1521-1529:

```r
    data_attrs <- sprintf(
      'data-hash="%s" data-game-id="%s" data-home="%s" data-away="%s" data-time="%s" data-market="%s" data-bet-on="%s" data-line="%s" data-prob="%s" data-ev="%s" data-size="%s" data-odds="%s" data-book="%s" data-actual="%s" data-fill-status="%s"',
      row$bet_hash, row$id, row$home_team, row$away_team,
      as.character(row$pt_start_time), row$market, row$bet_on,
      ifelse(is.na(row$line), "", row$line),
      row$prob, row$ev, row$bet_size, row$odds, row$bookmaker_key,
      ifelse(is.na(placed_actual), "", placed_actual),
      row$fill_status
    )
```

Add a `data-model-size` attribute alongside `data-size`. Both start as `row$bet_size` — the JS-side click-to-edit will mutate `data-size` while leaving `data-model-size` untouched:

```r
    data_attrs <- sprintf(
      'data-hash="%s" data-game-id="%s" data-home="%s" data-away="%s" data-time="%s" data-market="%s" data-bet-on="%s" data-line="%s" data-prob="%s" data-ev="%s" data-size="%s" data-model-size="%s" data-odds="%s" data-book="%s" data-actual="%s" data-fill-status="%s"',
      row$bet_hash, row$id, row$home_team, row$away_team,
      as.character(row$pt_start_time), row$market, row$bet_on,
      ifelse(is.na(row$line), "", row$line),
      row$prob, row$ev, row$bet_size, row$bet_size, row$odds, row$bookmaker_key,
      ifelse(is.na(placed_actual), "", placed_actual),
      row$fill_status
    )
```

- [ ] **Step 2: Update `placeBet` JS to read both attributes**

In the same file, locate `function placeBet(btn)` around line 3790. Update the body construction:

```js
          var body = {
            bet_hash:         data.hash,
            bookmaker_key:    book,
            account:          account,
            bet_on:           data.betOn,
            line:             (data.line === '' || data.line === undefined) ? null : parseFloat(data.line),
            market:           data.market,
            american_odds:    parseInt(data.odds, 10),
            actual_size:      parseFloat(data.size),
            kelly_bet:        parseFloat(data.modelSize || data.size),
            wz_odds_at_place: parseInt(data.odds, 10),
            game_id:          data.gameId,
            home_team:        data.home,
            away_team:        data.away,
            game_time:        data.time,
            model_prob:       data.prob ? parseFloat(data.prob) : 0.0,
            model_ev:         data.ev   ? parseFloat(data.ev)   : 0.0
          };
```

The only line that changed is `kelly_bet: parseFloat(data.modelSize || data.size)`. The fallback to `data.size` keeps the legacy fallback table working (it doesn't set `data-model-size`).

- [ ] **Step 3: Commit**

```bash
cd "/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-bets-tab-pr-c-editable-risk" && \
  git add "Answer Keys/MLB Dashboard/mlb_dashboard.R" && \
  git commit -m "$(cat <<'EOF'
fix(bets-tab): split data-size (editable) from data-model-size (Kelly)

The placeBet JS sent both actual_size and kelly_bet from the same
data-size attribute, so when editable-Risk lands and the user
overrides the size, placed_bets.recommended_size would also become
the override — losing the model's recommended number. Add a sibling
data-model-size attribute that's never mutated; placeBet now reads
data.modelSize for kelly_bet (with data.size fallback for the legacy
fallback table that doesn't set the new attribute).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 4 — Editable Risk in the hero strip

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R` — `render_hero_strip` (around line 1205-1230)

- [ ] **Step 1: Update `render_hero_strip` to emit click-to-edit Risk + slots for snap-back, error, verified-badge**

Open `Answer Keys/MLB Dashboard/mlb_dashboard.R`. Locate `render_hero_strip` around line 1205. Replace the function body with:

```r
render_hero_strip <- function(pick_book, pick_odds, fair_odds,
                               ev_pct, risk_dollars, towin_dollars,
                               action_html) {
  ev_str   <- sprintf("+%.1f%%", ev_pct)
  risk_str <- sprintf("$%.0f", risk_dollars)
  win_str  <- sprintf("$%.0f", towin_dollars)
  fair_str <- .format_odds_signed(fair_odds)
  odds_str <- .format_odds_signed(pick_odds)

  # Risk and To Win get extra wrappers so the bets-tab JS coordinator
  # can attach click-to-edit + the WZ-verified-quote swap. data-model-risk
  # carries the original Kelly value so the snap-back button can revert.
  sprintf(
    '<div class="hero">
       <div class="pick">
         <span class="book">%s</span>
         <span class="odds">%s</span>
       </div>
       <div class="divider"></div>
       <div class="stat"><span class="lbl">Fair</span><span class="val fair">%s</span></div>
       <div class="stat"><span class="lbl">EV</span><span class="val ev">%s</span></div>
       <div class="stat risk-stat" data-model-risk="%.0f" data-american-odds="%d">
         <span class="lbl">Risk</span>
         <span class="val risk risk-value" tabindex="0" title="click to edit">%s</span>
         <button type="button" class="risk-reset" title="reset to model size">&#8634;</button>
         <span class="risk-error" hidden></span>
       </div>
       <div class="stat"><span class="lbl">To Win</span><span class="val win towin-value">%s</span><span class="towin-status" hidden></span></div>
       <div class="actions">%s</div>
     </div>',
    htmltools::htmlEscape(pick_book),
    odds_str, fair_str, ev_str,
    risk_dollars, as.integer(pick_odds),
    risk_str, win_str, action_html
  )
}
```

Notable changes from the existing version:
- The Risk `<span>` becomes `tabindex="0"` so it's keyboard-focusable; gains classes `risk-value` (for JS targeting) and a `title` tooltip.
- A `<button class="risk-reset">` (the `↶` snap-back) sits next to it, hidden by default via CSS.
- A `<span class="risk-error">` slot for inline drift / rejection messages.
- The To Win `<span>` gains `towin-value` class for JS targeting.
- A `<span class="towin-status">` slot for the verified-`✓` badge or spinner.
- The wrapping `<div class="stat">` becomes `risk-stat` and carries `data-model-risk` (original Kelly) and `data-american-odds` (so JS can compute optimistic local to-win when a non-WZ row is edited).

- [ ] **Step 2: Add CSS for editable Risk + snap-back + error + verified badge**

In the same file, find the `<style>` block where `.hero`, `.pick`, `.stat`, etc. are defined. Add:

```css
.hero .risk-stat { position: relative; }
.hero .risk-value {
  cursor: text;
  padding: 1px 6px;
  margin: -1px -6px;
  border-radius: 5px;
  border: 1px dashed transparent;
  display: inline-block;
  min-width: 48px;
  transition: all 0.15s ease;
}
.hero .risk-value:hover:not(.editing):not(.overridden) {
  border-color: #3a4658;
  background: rgba(255, 255, 255, 0.02);
}
.hero .risk-value.editing {
  border-style: solid;
  border-color: #3fb950;
  background: #0d1117;
  box-shadow: 0 0 0 3px rgba(63, 185, 80, 0.15);
  outline: none;
  cursor: text;
}
.hero .risk-value.overridden { color: #d29922; }

.hero .risk-reset {
  display: none;
  margin-left: 6px;
  width: 18px; height: 18px;
  padding: 0;
  border-radius: 50%;
  border: 1px solid #3a4658;
  background: #161b22;
  color: #7d8590;
  cursor: pointer;
  font-size: 11px;
  line-height: 1;
  vertical-align: middle;
}
.hero .risk-stat.is-overridden .risk-reset { display: inline-flex; align-items: center; justify-content: center; }
.hero .risk-reset:hover {
  color: #3fb950; border-color: #3fb950; background: rgba(63, 185, 80, 0.10);
}

.hero .risk-error {
  display: none;
  margin-left: 8px;
  padding: 2px 8px;
  border-radius: 5px;
  background: rgba(248, 81, 73, 0.12);
  color: #f85149;
  font-size: 11px;
  font-family: 'SF Mono', monospace;
}
.hero .risk-error[hidden] { display: none; }
.hero .risk-stat.has-error .risk-error { display: inline-block; }

.hero .towin-status { margin-left: 6px; font-size: 11px; }
.hero .towin-status.verified { color: #3fb950; }
.hero .towin-status.spinner  { color: #7d8590; }
.hero .towin-status[hidden]  { display: none; }
```

- [ ] **Step 3: Commit (markup + CSS only; JS coordinator follows)**

```bash
cd "/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-bets-tab-pr-c-editable-risk" && \
  git add "Answer Keys/MLB Dashboard/mlb_dashboard.R" && \
  git commit -m "$(cat <<'EOF'
feat(bets-tab): hero strip — editable Risk markup + snap-back/error slots

Risk value becomes a focusable .risk-value <span> with click-to-edit
affordance via CSS dashed border on hover. Sibling slots: risk-reset
(↶ snap-back, hidden until override), risk-error (inline error pill,
hidden until set), towin-status (verified-✓ badge or spinner). Wrapper
.stat.risk-stat carries data-model-risk (original Kelly) and
data-american-odds for the JS coordinator coming next.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 5 — JS coordinator: click-to-edit, fire WZ quote, swap to-win

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R` — append a `<script>` block to the bets-tab page template

- [ ] **Step 1: Append the coordinator script to the bets-tab page**

Find where the bets-tab page template adds inline scripts (search for the `<script>` block near `function placeBet(btn)` — that's the bets-tab page). Add this script alongside (in the same `<script>` tag or a sibling one):

```html
<script>
  (function setupEditableRisk() {
    var WZ_BOOK = 'wagerzon';

    function fmtMoney(x) {
      if (x == null || !isFinite(x)) return '$0';
      return '$' + Math.round(x).toLocaleString();
    }

    function localTowin(risk, americanOdds) {
      if (!isFinite(risk) || risk <= 0 || !isFinite(americanOdds) || americanOdds === 0) return 0;
      return americanOdds > 0 ? risk * americanOdds / 100 : risk * 100 / Math.abs(americanOdds);
    }

    function placeBtnFor(card) {
      // Place button is in .actions; only one per card.
      return card.querySelector('.actions .btn-place, .actions [onclick*="placeBet"]');
    }

    function bookFor(card) {
      var btn = placeBtnFor(card);
      return btn ? (btn.dataset.book || '').toLowerCase() : '';
    }

    function setOverride(card, riskValue, isOverride) {
      var stat    = card.querySelector('.risk-stat');
      var valueEl = stat.querySelector('.risk-value');
      var towinEl = card.querySelector('.towin-value');
      var amerOdds = parseInt(stat.dataset.americanOdds, 10);

      valueEl.textContent = fmtMoney(riskValue);
      stat.classList.toggle('is-overridden', isOverride);
      valueEl.classList.toggle('overridden', isOverride);

      // Optimistic local to-win update; verified-quote will overwrite if WZ.
      towinEl.textContent = fmtMoney(localTowin(riskValue, amerOdds));

      // Sync the placement button's data-size so the existing placeBet flow uses the override.
      var btn = placeBtnFor(card);
      if (btn) btn.dataset.size = String(riskValue);
    }

    function clearError(card) {
      var stat   = card.querySelector('.risk-stat');
      var errEl  = stat.querySelector('.risk-error');
      stat.classList.remove('has-error');
      errEl.hidden = true;
      errEl.textContent = '';
    }
    function showError(card, msg) {
      var stat  = card.querySelector('.risk-stat');
      var errEl = stat.querySelector('.risk-error');
      errEl.textContent = msg;
      errEl.hidden = false;
      stat.classList.add('has-error');
    }

    function setTowinStatus(card, kind, text) {
      // kind: 'verified' | 'spinner' | null (clear)
      var status = card.querySelector('.towin-status');
      if (!status) return;
      status.classList.remove('verified', 'spinner');
      if (!kind) {
        status.hidden = true;
        status.textContent = '';
        return;
      }
      status.classList.add(kind);
      status.textContent = text || '';
      status.hidden = false;
    }

    function startEdit(valueEl) {
      if (valueEl.classList.contains('editing')) return;
      var current = Number(valueEl.textContent.replace(/[^0-9.]/g, ''));
      valueEl.classList.add('editing');
      valueEl.contentEditable = 'true';
      valueEl.textContent = isFinite(current) ? String(current) : '';
      var sel = window.getSelection();
      var range = document.createRange();
      range.selectNodeContents(valueEl);
      sel.removeAllRanges(); sel.addRange(range);
      valueEl.focus();
    }

    function commitEdit(valueEl) {
      if (!valueEl.classList.contains('editing')) return;
      var card = valueEl.closest('.bet-card-v8');
      var stat = card.querySelector('.risk-stat');
      var modelRisk = Number(stat.dataset.modelRisk);
      var raw = valueEl.textContent.replace(/[^0-9.]/g, '');
      var amount = Number(raw);
      if (!isFinite(amount) || amount <= 0) amount = modelRisk;
      valueEl.classList.remove('editing');
      valueEl.contentEditable = 'false';
      var isOverride = Math.abs(amount - modelRisk) > 0.005;
      setOverride(card, amount, isOverride);
      clearError(card);

      if (bookFor(card) === WZ_BOOK) {
        verifyWithWz(card, amount);
      } else {
        // Non-WZ: local math is already shown; clear any verified badge.
        setTowinStatus(card, null);
      }
    }

    function verifyWithWz(card, amount) {
      var btn = placeBtnFor(card);
      if (!btn) return;
      var account = (window.WZ_SELECTED_ACCOUNT || null);
      if (!account) {
        showError(card, 'pick a WZ account first');
        return;
      }
      setTowinStatus(card, 'spinner', 'verifying...');
      var body = {
        bet_hash:      btn.dataset.hash,
        amount:        amount,
        account:       account,
        // Pass the same fields /api/place-bet sends so the resolver
        // can find idgm/play if not already known.
        bet_on:        btn.dataset.betOn,
        line:          (btn.dataset.line === '' || btn.dataset.line === undefined)
                          ? null : parseFloat(btn.dataset.line),
        market:        btn.dataset.market,
        american_odds: parseInt(btn.dataset.odds, 10),
        game_id:       btn.dataset.gameId,
        home_team:     btn.dataset.home,
        away_team:     btn.dataset.away,
      };
      fetch('/api/wz-quote-single', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify(body)
      })
        .then(function(r) { return r.json(); })
        .then(function(j) {
          var towinEl = card.querySelector('.towin-value');
          if (j.error_msg_key) {
            showError(card, j.error_msg || j.error_msg_key);
            setTowinStatus(card, null);
            return;
          }
          var expectedOdds = parseInt(btn.dataset.odds, 10);
          if (j.current_wz_odds && j.current_wz_odds !== expectedOdds) {
            showError(card, 'WZ now ' + (j.current_wz_odds > 0 ? '+' : '') +
                            j.current_wz_odds + ' (was ' + (expectedOdds > 0 ? '+' : '') +
                            expectedOdds + ')');
          }
          if (j.win != null) {
            towinEl.textContent = fmtMoney(j.win);
            setTowinStatus(card, 'verified', '✓ wz');
          } else {
            setTowinStatus(card, null);
          }
        })
        .catch(function(e) {
          showError(card, 'verify failed: ' + e.message);
          setTowinStatus(card, null);
        });
    }

    document.addEventListener('click', function(ev) {
      var valueEl = ev.target.closest && ev.target.closest('.risk-value');
      if (valueEl) {
        startEdit(valueEl);
        return;
      }
      var resetBtn = ev.target.closest && ev.target.closest('.risk-reset');
      if (resetBtn) {
        var card = resetBtn.closest('.bet-card-v8');
        var stat = card.querySelector('.risk-stat');
        var modelRisk = Number(stat.dataset.modelRisk);
        setOverride(card, modelRisk, false);
        clearError(card);
        setTowinStatus(card, null);
      }
    });
    document.addEventListener('focusout', function(ev) {
      var valueEl = ev.target.closest && ev.target.closest('.risk-value');
      if (valueEl && valueEl.classList.contains('editing')) commitEdit(valueEl);
    });
    document.addEventListener('keydown', function(ev) {
      var valueEl = ev.target.closest && ev.target.closest('.risk-value');
      if (!valueEl || !valueEl.classList.contains('editing')) return;
      if (ev.key === 'Enter')   { ev.preventDefault(); valueEl.blur(); }
      if (ev.key === 'Escape')  {
        ev.preventDefault();
        var card = valueEl.closest('.bet-card-v8');
        var stat = card.querySelector('.risk-stat');
        var modelRisk = Number(stat.dataset.modelRisk);
        valueEl.classList.remove('editing');
        valueEl.contentEditable = 'false';
        setOverride(card, modelRisk, false);
        clearError(card);
      }
    });
  })();
</script>
```

- [ ] **Step 2: Manually verify the click-to-edit + verify flow on a WZ bet**

Reload the dashboard. Pick a card whose pick book is WZ. Click `Risk $XX`:
1. The value becomes editable (green dashed border).
2. Type `1` and press Enter.
3. The value renders amber; a `↶` button appears next to it.
4. To Win briefly shows `verifying...`, then either swaps to a real number with `✓ wz` next to it OR shows an inline error pill (e.g. `WZ now -125 (was +120)`).
5. Click `↶` → snaps back to the model number, error and verified-badge clear.

If WZ is unreachable (auth expired), you'll see `verify failed: ...` — re-auth via the dashboard's existing flow and try again.

For a **non-WZ card** (Hoop88, BFA): editing should still work, To Win still updates from local math, but no `verifying...` and no `✓ wz`.

- [ ] **Step 3: Commit**

```bash
cd "/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-bets-tab-pr-c-editable-risk" && \
  git add "Answer Keys/MLB Dashboard/mlb_dashboard.R" && \
  git commit -m "$(cat <<'EOF'
feat(bets-tab): editable Risk + WZ verified to-win + drift surface

Click-to-edit on the hero-strip Risk value. On commit (Enter/blur):
local to-win updates immediately from American-odds math; for WZ
rows, a parallel /api/wz-quote-single fires to fetch the actual Win
and current odds. Real Win replaces local math with a ✓ wz badge.
Drift / rejection / amount-rule errors surface inline next to Risk.
↶ snap-back reverts to the model's recommended size. Refresh wipes
overrides (per design — no persistence). Non-WZ rows show local math
only.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 6 — Documentation updates

**Files:**
- Modify: `wagerzon_odds/CLAUDE.md`
- Modify: `Answer Keys/MLB Dashboard/README.md`
- Modify: `Answer Keys/CLAUDE.md`

- [ ] **Step 1: Update `wagerzon_odds/CLAUDE.md`**

Open `wagerzon_odds/CLAUDE.md`. In the "Quick map" section, add a bullet between `parlay_pricer.py` and `parlay_placer.py`:

```markdown
- **Pricing singles:** `single_pricer.py::get_single_price()` calls
  `ConfirmWagerHelper` with `WT=0` (single) and `RiskWin=2` (skip
  balance check). Returns `{win, current_wz_odds, error_msg, error_msg_key}`.
  Used by the MLB Dashboard's editable-Risk feature to verify to-win
  and detect line drift before placement.
```

- [ ] **Step 2: Update `Answer Keys/MLB Dashboard/README.md`**

Add to the placing-bets section:

```markdown
### Editable Risk + WZ verified to-win (PR C, 2026-05)

The hero-strip `Risk $XX` value on every bet card is click-to-edit.
Type any number and press Enter — for **Wagerzon** picks, the dashboard
fires a `ConfirmWagerHelper` preflight (`POST /api/wz-quote-single`)
and replaces the local-math "to win" with the actual `Win` WZ would
credit, plus a `✓ wz` verified badge. Inline error pill surfaces line
drift (`WZ now -125 (was +120)`), amount-rule rejections (integer-only,
min stake), or any other `ErrorMsg` WZ returns.

Non-WZ books (Hoop88, BFA, BetOnlineAG) still let you edit Risk, but
"to win" stays as local American-odds math (no preview API exposed
by those books).

The `↶` button next to an overridden Risk snaps back to the
Kelly-recommended size. Page refresh wipes overrides — they don't
persist across reloads.

Note: `placed_bets.recommended_size` continues to record the model's
Kelly number even when you override (`actual_size`), so CLV analysis
can later ask "when did I override and how did I do?"
```

- [ ] **Step 3: Update `Answer Keys/CLAUDE.md`**

In the "MLB Dashboard — Odds screen + WZ single-bet placer" section, add:

```markdown
### Editable Risk + WZ verified-quote (PR C, 2026-05)

`POST /api/wz-quote-single` (in `mlb_dashboard_server.py`) wraps
`wagerzon_odds.single_pricer.get_single_price()` for client-side
preview-on-edit. Body: `{bet_hash, amount, account}`; returns
`{win, current_wz_odds, error_msg, error_msg_key}`. The bets-tab JS
coordinator fires this on every Risk-edit commit for WZ picks.

Placement plumbing now keeps `data-size` (editable override → `actual_size`)
separate from `data-model-size` (immutable Kelly → `kelly_bet` /
`recommended_size`). The two attributes are siblings on every Place
button; legacy callers without `data-model-size` fall back to `data-size`
as a default to keep `create_bets_table_legacy` working.
```

- [ ] **Step 4: Commit**

```bash
cd "/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-bets-tab-pr-c-editable-risk" && \
  git add wagerzon_odds/CLAUDE.md "Answer Keys/MLB Dashboard/README.md" "Answer Keys/CLAUDE.md" && \
  git commit -m "$(cat <<'EOF'
docs: editable Risk + WZ verified-quote feature

Notes the new /api/wz-quote-single endpoint, single_pricer helper,
data-size / data-model-size split, and user-facing flow across the
three relevant CLAUDE.md / README.md files.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 7 — End-to-end placement smoke test (placement plumbing intact)

**Files:** none (verification only).

**Why this matters.** The `data-size` / `data-model-size` split (Task 3) and the new editable-Risk flow share the placement-button data attributes. We must confirm that an actual `$1` placement still completes successfully end-to-end.

- [ ] **Step 1: Re-run pytest for placement-adjacent tests**

```bash
cd "/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-bets-tab-pr-c-editable-risk" && \
  python -m pytest wagerzon_odds/test_parlay_placer.py wagerzon_odds/tests/test_single_pricer.py -v
```

Expected: all green. (No new tests for `single_placer.py` are needed; we didn't change it.)

- [ ] **Step 2: Place a $1 bet end-to-end**

In the dashboard:
1. Pick any WZ card.
2. Edit Risk → `1`, press Enter.
3. Verify the To Win swaps to a real number with `✓ wz`.
4. Click Place Bet.
5. Confirm the bet appears in WZ's open bets list (manual check via WZ web UI) AND in the dashboard's "Placed Bets" section.
6. Inspect `placed_bets` in `mlb_dashboard.duckdb`:
   ```bash
   python3 -c "
   import duckdb
   con = duckdb.connect('Answer Keys/MLB Dashboard/mlb_dashboard.duckdb', read_only=True)
   df = con.execute('''
     SELECT bet_hash, status, recommended_size, actual_size, ticket_number
     FROM placed_bets
     ORDER BY id DESC LIMIT 1
   ''').fetchdf()
   print(df.to_string())
   "
   ```
   Expected: `actual_size = 1.0`, `recommended_size = <model's Kelly number>` (NOT 1.0). The two columns differ — that's the bug fix from Task 3 working.

If `actual_size == recommended_size == 1.0`, the data-attribute split didn't reach `placeBet`. Re-check Task 3.

If WZ rejects placement, that's an unrelated issue (auth, balance, line) — surface the error and decide whether to retry or back out the test.

---

## Task 8 — Final verification + handoff

**Files:** none (verification only).

- [ ] **Step 1: Confirm both unit-test suites pass**

```bash
cd "/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-bets-tab-pr-c-editable-risk" && \
  python -m pytest wagerzon_odds/tests/test_single_pricer.py -v
```

- [ ] **Step 2: Sanity-check the full bets-tab flow**

Reload the dashboard. Confirm:
1. **Editable Risk** works for WZ and non-WZ cards (basic typing + Enter).
2. **Snap-back** restores the model number and clears verified badge / error.
3. **Drift detection** — find a card whose WZ odds have moved since the dashboard's last refresh; commit an edit; the verified-quote returns the new odds and surfaces the drift.
4. **Amount rejection** — try a non-integer like `1.5` (if WZ enforces integer-only). Inline error should show WZ's exact message.
5. **Refresh wipes overrides** — set an override on a card, hit dashboard Refresh; the override is gone, Risk shows the Kelly number.

- [ ] **Step 3: Summarize the branch state**

```bash
cd "/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-bets-tab-pr-c-editable-risk" && \
  git log --oneline main..HEAD && \
  git diff --stat main..HEAD
```

Expected: ~7 commits ahead of main. Diff stat shows the new `single_pricer.py`, the new test file, `mlb_dashboard_server.py`, `mlb_dashboard.R`, and the three doc files.

Hand off to the user with: "PR C is ready. Editable Risk + WZ verified-quote working end-to-end; placement plumbing now records actual vs recommended size separately. Ready to merge to main when you give the word."

---

## Out-of-band notes

- **Do NOT merge to main without explicit user approval.** Project policy.
- **Don't change `single_placer.py`.** We add a new `single_pricer.py` instead. Keeping the placement code stable (it already has production traffic) is intentional.
- **Auth dependency.** The endpoint requires a live WZ session. If the user's session expired between page load and edit, `verifyWithWz` surfaces `session_expired` cleanly and the user can re-auth via the existing flow without losing their override.
- **Per-keystroke firing was rejected** during brainstorming — only fire on commit (Enter / blur), not on `input`. Otherwise typing `25` becomes 1 call for `2`, 1 for `25`.
- **The verified-quote pattern is WZ-only** because no other supported book exposes a public preflight endpoint. If H88 / BFA / BetOnlineAG ever add one, we can extend.
- **`actual_size` vs `recommended_size` data integrity** — the two-attribute split is the core fix that makes future "did my overrides beat the model?" analysis possible. Don't merge a future change that re-conflates them without an explicit replacement plan.
