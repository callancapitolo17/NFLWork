# WC Maker Touch-Join Pricing Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the WC totals maker quote at the crowd's best bid (one-sided) whenever the devigged anchor fair clears that price by ≥1c after maker fees, with queue-preserving hysteresis, instead of always resting fair−margin behind the flow.

**Architecture:** A new `_touch_join_price` helper in `MakerEngine` computes a self-excluded crowd touch and applies entry/exit edge thresholds; `_desired` takes the better of the join price and the legacy fair−margin price; the resting-order state gains a `mode` field so held touch-joins are recognized across ticks and fills are attributable (`touch_join` vs `quote` in `maker_quotes`).

**Tech Stack:** Python 3.14, pytest, existing `unabated_edge.maker` package, `kalshi_common.ev_calc.maker_fee_per_contract`.

## Global Constraints

- Spec: `docs/superpowers/specs/2026-07-18-wc-touchjoin-pricing-design.md`. Deviation approved in-session: joins are additionally bounded above by the existing `MAX_MARGIN_CENTS` (5c) suspicion guard — entry requires `threshold ≤ edge ≤ 5c`; hysteresis holds require `0.25c ≤ edge ≤ 5c`.
- Touch-join is default behavior — **no enable flag**. Dials only: `TOUCH_JOIN_MIN_EDGE_CENTS=1.0`, `TOUCH_JOIN_ALT_MIN_EDGE_CENTS=1.5`, `TOUCH_JOIN_EXIT_EDGE_CENTS=0.25`.
- All work on branch `worktree-wc-maker-depth-capture` in its existing worktree (`/Users/callancapitolo/NFLWork/.claude/worktrees/wc-maker-depth-capture`). Never switch branches; never touch main. The live maker (pid 55049) runs from this cwd — code edits are inert until the coordinated restart in Task 4.
- Sizing, caps, ledger, pull triggers, gateway, shadow mode: unchanged.
- Book dict shape (from `venues/kalshi.py::get_book`): `{"yes_bids": [(price_dollars, qty), ...], "no_bids": [...]}`, each sorted best (highest) first.
- Run tests from the worktree root: `cd /Users/callancapitolo/NFLWork/.claude/worktrees/wc-maker-depth-capture && python3 -m pytest unabated_edge/tests/ -q`.

---

### Task 1: Config dials + `_touch_join_price` helper

**Files:**
- Modify: `unabated_edge/config.py` (maker block, after `ANCHOR_STALE_FARK_SEC`, ~line 71)
- Modify: `unabated_edge/maker/state.py:76-77` (`on_place` gains `mode`)
- Modify: `unabated_edge/maker/engine.py` (new method after `_margin_cents`, ~line 85)
- Test: `unabated_edge/tests/test_maker_engine.py` (append)

**Interfaces:**
- Consumes: `state.resting_for(ticker, side)` → `{"order_id","price_cents","count","mode"}` or None; `maker_fee_per_contract(p)` (dollars/contract).
- Produces: `MakerEngine._touch_join_price(ticker, side, fair_cents, alt, book, opp_ask_cents) -> int | None` (join price in cents, or None when touch-join does not apply); `MakerState.on_place(ticker, side, order_id, price_cents, count, mode="quote")`.

- [ ] **Step 1: Write the failing tests** (append to `unabated_edge/tests/test_maker_engine.py`):

```python
# ---------- touch-join (spec 2026-07-18) ----------

def _tj(eng, book, fair_cents, ticker="T-O25", side="yes", alt=False, opp_ask_cents=99):
    return eng._touch_join_price(ticker, side, fair_cents, alt, book, opp_ask_cents)

def test_touch_join_entry_within_band(eng):
    book = {"yes_bids": [(0.78, 500.0)], "no_bids": [(0.21, 500.0)]}
    assert _tj(eng, book, 80.0) == 78          # edge ~1.7c in [1, 5]

def test_touch_join_rejects_thin_and_suspicious_edges(eng):
    book = {"yes_bids": [(0.78, 500.0)], "no_bids": []}
    assert _tj(eng, book, 79.0) is None        # edge ~0.7c < 1c entry
    deep = {"yes_bids": [(0.60, 500.0)], "no_bids": []}
    assert _tj(eng, deep, 80.0) is None        # edge ~19c > MAX_MARGIN suspicion guard

def test_touch_join_alt_needs_bigger_edge(eng):
    book = {"yes_bids": [(0.78, 500.0)], "no_bids": []}
    assert _tj(eng, book, 79.5, alt=True) is None   # ~1.2c < 1.5c alt entry
    assert _tj(eng, book, 79.9, alt=True) == 78     # ~1.6c >= 1.5c

def test_touch_join_self_exclusion(eng):
    # our own 2-lot at 78c is the only order at the level -> not "the crowd"
    eng.state.on_place("T-O25", "yes", "o-x", 78, 2, mode="quote")
    book = {"yes_bids": [(0.78, 2.0), (0.63, 500.0)], "no_bids": []}
    assert _tj(eng, book, 80.0) is None        # real touch 63 -> edge 16.7c: suspicious

def test_touch_join_hysteresis_hold_and_exit(eng):
    eng.state.on_place("T-O25", "yes", "o-x", 78, 2, mode="touch_join")
    book = {"yes_bids": [(0.78, 2.0)], "no_bids": []}   # only us at the level
    assert _tj(eng, book, 79.0) == 78          # edge ~0.7c in [0.25, 5] -> hold
    assert _tj(eng, book, 78.2) is None        # edge ~-0.1c < 0.25c -> exit

def test_touch_join_hold_ignored_for_legacy_orders(eng):
    eng.state.on_place("T-O25", "yes", "o-x", 78, 2, mode="quote")
    book = {"yes_bids": [(0.78, 2.0)], "no_bids": []}
    assert _tj(eng, book, 79.0) is None        # hysteresis only holds touch-joins

def test_touch_join_never_crosses(eng):
    book = {"yes_bids": [(0.78, 500.0)], "no_bids": []}
    assert _tj(eng, book, 80.0, opp_ask_cents=78) is None   # touch == ask-0 -> cross
```

- [ ] **Step 2: Run to verify failure**

Run: `python3 -m pytest unabated_edge/tests/test_maker_engine.py -q -k touch_join`
Expected: errors — `AttributeError: 'MakerEngine' object has no attribute '_touch_join_price'` / `TypeError: on_place() got an unexpected keyword argument 'mode'`.

- [ ] **Step 3: Implement**

`unabated_edge/config.py`, in the maker block after `ANCHOR_STALE_FARK_SEC`:

```python
TOUCH_JOIN_MIN_EDGE_CENTS = float(_get("TOUCH_JOIN_MIN_EDGE_CENTS", "1.0"))       # net edge (fair - touch - maker fee) to join the crowd's best bid
TOUCH_JOIN_ALT_MIN_EDGE_CENTS = float(_get("TOUCH_JOIN_ALT_MIN_EDGE_CENTS", "1.5"))  # alt rungs: less-trusted fair demands more edge
TOUCH_JOIN_EXIT_EDGE_CENTS = float(_get("TOUCH_JOIN_EXIT_EDGE_CENTS", "0.25"))    # hold a resting touch-join (queue position) until edge decays below this
```

`unabated_edge/maker/state.py:76-77`:

```python
    def on_place(self, ticker, side, order_id, price_cents, count, mode="quote"):
        self.resting[(ticker, side)] = {"order_id": order_id, "price_cents": price_cents,
                                        "count": count, "mode": mode}
```

`unabated_edge/maker/engine.py`, new method after `_margin_cents`:

```python
    def _touch_join_price(self, ticker, side, fair_cents, alt, book, opp_ask_cents):
        """Crowd-touch join price (cents) or None. Entry: net edge within
        [threshold, MAX_MARGIN_CENTS] — the upper bound is the same too-good-
        to-be-true guard as the legacy crowd_tighter skip. Hysteresis: an
        already-resting touch-join holds (queue position is capital) until its
        edge decays below TOUCH_JOIN_EXIT_EDGE_CENTS or turns suspicious. The
        touch is self-excluded: our own resting qty never counts as crowd."""
        bids = book["yes_bids"] if side == "yes" else book["no_bids"]
        cur = self.state.resting_for(ticker, side)

        def edge(price_cents):
            return fair_cents - price_cents - maker_fee_per_contract(price_cents / 100.0) * 100

        touch = None
        for p, q in bids:                       # best-first
            pc = int(round(p * 100))
            if cur and cur["price_cents"] == pc:
                q -= cur["count"]
            if q > 0.5:
                touch = pc
                break
        thr = config.TOUCH_JOIN_ALT_MIN_EDGE_CENTS if alt else config.TOUCH_JOIN_MIN_EDGE_CENTS
        if touch is not None and touch <= opp_ask_cents - 1 and \
                thr <= edge(touch) <= config.MAX_MARGIN_CENTS:
            return touch
        if cur and cur.get("mode") == "touch_join" and cur["price_cents"] <= opp_ask_cents - 1 and \
                config.TOUCH_JOIN_EXIT_EDGE_CENTS <= edge(cur["price_cents"]) <= config.MAX_MARGIN_CENTS:
            return cur["price_cents"]
        return None
```

- [ ] **Step 4: Run to verify pass**

Run: `python3 -m pytest unabated_edge/tests/test_maker_engine.py -q -k touch_join`
Expected: 7 passed. Then full file: `python3 -m pytest unabated_edge/tests/test_maker_engine.py -q` — all pass (existing tests untouched by the helper).

- [ ] **Step 5: Commit**

```bash
git add unabated_edge/config.py unabated_edge/maker/state.py unabated_edge/maker/engine.py unabated_edge/tests/test_maker_engine.py
git commit -m "feat(unabated-edge): touch-join price helper + dials (spec 2026-07-18)"
```

---

### Task 2: Wire touch-join into `_desired` / `_sync` with mode attribution

**Files:**
- Modify: `unabated_edge/maker/engine.py::_desired` (lines 93-121) and `::_sync` (lines 125-156)
- Test: `unabated_edge/tests/test_maker_engine.py` (append)

**Interfaces:**
- Consumes: `_touch_join_price` from Task 1.
- Produces: `_desired` returns `((price_cents, count, mode), None)` or `(None, skip_reason)` where `mode` ∈ `{"quote", "touch_join"}`; `_sync` logs the success row's `reason` as that mode and passes it to `state.on_place`.

- [ ] **Step 1: Write the failing tests** (append):

```python
_JBOOK = {"T-O25": {"yes_bids": [(0.78, 500.0)], "no_bids": [(0.21, 500.0)]}}  # yes ask 79c
_JLADDER = {2.5: {"p_over": 0.80, "book": 7, "alt": False, "overround": 1.05,
                  "modified_on": _NOW.isoformat()}}

def test_desired_joins_touch_when_edge_in_band(eng):
    eng.on_match(Soccer(), _EM, _KEV, _JLADDER, _JBOOK, _NOW)
    yes = [(p, n) for _t, s, p, n in eng.gateway.placed if s == "yes"]
    assert yes and yes[0][0] == 78              # joined the touch, not legacy 77
    assert eng.state.resting_for("T-O25", "yes")["mode"] == "touch_join"

def test_desired_falls_back_to_legacy_when_no_join(eng):
    eng.on_match(Soccer(), _EM, _KEV, _LADDER, _BOOKS, _NOW)   # wide crowd fixture
    assert eng.state.resting_for("T-O25", "yes")["mode"] == "quote"

def test_join_holds_through_fair_wiggle_no_churn(eng):
    eng.on_match(Soccer(), _EM, _KEV, _JLADDER, _JBOOK, _NOW)
    n_placed, n_cancelled = len(eng.gateway.placed), len(eng.gateway.cancelled)
    wiggle = {2.5: {**_JLADDER[2.5], "p_over": 0.79}}           # entry fails, exit holds
    eng.on_match(Soccer(), _EM, _KEV, wiggle, _JBOOK, _NOW + datetime.timedelta(seconds=5))
    assert len(eng.gateway.placed) == n_placed and len(eng.gateway.cancelled) == n_cancelled

def test_join_exits_to_legacy_when_edge_gone(eng):
    eng.on_match(Soccer(), _EM, _KEV, _JLADDER, _JBOOK, _NOW)
    dead = {2.5: {**_JLADDER[2.5], "p_over": 0.782}}            # edge ~-0.1c
    eng.on_match(Soccer(), _EM, _KEV, dead, _JBOOK, _NOW + datetime.timedelta(seconds=5))
    assert eng.gateway.cancelled                                 # touch-join replaced
    yes_prices = [p for _t, s, p, _n in eng.gateway.placed if s == "yes"]
    assert yes_prices[-1] == 75                                  # legacy: floor(78.2) - 3c margin
    assert eng.state.resting_for("T-O25", "yes")["mode"] == "quote"
```

- [ ] **Step 2: Run to verify failure**

Run: `python3 -m pytest unabated_edge/tests/test_maker_engine.py -q -k "desired_joins or falls_back or holds_through or exits_to_legacy"`
Expected: FAIL — yes quote rests at 77 (legacy) and `mode` KeyError/mismatch.

- [ ] **Step 3: Implement**

`_desired` (replace lines 100-108, keep everything else):

```python
        opp_ask = yes_ask_from_book(book) if side == "yes" else no_ask_from_book(book)
        if opp_ask is None:
            return None, "no_crowd"
        opp_ask_cents = int(round(opp_ask * 100))
        legacy = math.floor(fair_cents + 1e-9) - self._margin_cents(fair_cents, alt)
        legacy = min(legacy, opp_ask_cents - 1)               # never cross
        join = self._touch_join_price(ticker, side, fair_cents, alt, book, opp_ask_cents)
        if join is not None and join > legacy:
            price, mode = join, "touch_join"
        else:
            price, mode = legacy, "quote"
            if price < fair_cents - config.MAX_MARGIN_CENTS - 1e-9:
                return None, "crowd_tighter"
        if price < 1:
            return None, "price_floor"
```

and the return becomes `return (price, n, mode), None`.

`_sync` (unpack + attribution; replace lines 139-156):

```python
        price, count, mode = desired
        if cur and cur["price_cents"] == price:
            return                                      # hold: queue position is capital
        action = "replace" if cur else "rest"
        if cur:
            if not self.gateway.cancel(cur["order_id"]):
                return                                  # keep local state; retry next tick
            self.state.on_cancel(ticker, side)
        coid = f"uemk-{ticker}-{side}-{price}-{int(now.timestamp())}"
        oid = self.gateway.place(ticker, side, price, count, coid)
        margin = round(fair - price / 100.0, 4)
        if oid is None:
            store.log_quote(now, sport, eid, ticker, side, "skip", price / 100.0, count,
                            fair, margin, alt, "place_failed", None)
            return
        self.state.on_place(ticker, side, oid, price, count, mode)
        store.log_quote(now, sport, eid, ticker, side, action, price / 100.0, count,
                        fair, margin, alt, mode, oid)
```

Note the hold-if-same-price branch is what makes hysteresis churn-free: `_touch_join_price` returns the current price, `_desired` returns it as desired, `_sync` no-ops.

- [ ] **Step 4: Run the full suite**

Run: `python3 -m pytest unabated_edge/tests/ -q`
Expected: all pass (118 pre-existing + 11 new). If any pre-existing `_desired`-shape test fails on the 3-tuple, fix the *test's unpacking only if it asserts on the tuple shape directly* — behavior assertions must pass unchanged.

- [ ] **Step 5: Commit**

```bash
git add unabated_edge/maker/engine.py unabated_edge/tests/test_maker_engine.py
git commit -m "feat(unabated-edge): quote at crowd touch when anchor edge >= 1c (touch_join attribution)"
```

---

### Task 3: README + docs

**Files:**
- Modify: `unabated_edge/README.md` (maker quoting-formula section, ~lines 86-92)

- [ ] **Step 1: Add a "Touch-join pricing" subsection** after the quoting formula: the join condition (`fair − touch − maker fee ∈ [1c, MAX_MARGIN_CENTS]`, alt 1.5c), self-exclusion, hysteresis (hold ≥0.25c), the three dials, `touch_join` vs `quote` attribution in `maker_quotes`, and a pointer to spec `2026-07-18-wc-touchjoin-pricing-design.md` + issues #1/#2 for the deferred breaker/markout work.

- [ ] **Step 2: Commit**

```bash
git add unabated_edge/README.md
git commit -m "docs(unabated-edge): touch-join pricing subsection"
```

---

### Task 4: Rollout to the live maker (final: ESP-ARG, Jul 19)

- [ ] **Step 1: Full suite green** — `python3 -m pytest unabated_edge/tests/ -q`, zero failures.
- [ ] **Step 2: FLAG IN CHAT before touching the live process** (user rule: announce kills). Then clean stop: `kill <pid of unabated_edge.runner>` (SIGTERM; escalate to `kill -9` only if it ignores SIGTERM >60s — startup phantom-cleanup makes this safe).
- [ ] **Step 3: Relaunch from this worktree cwd** with the same inline env as the current live launch (per memory: `MAKER_MODE=live MAKER_LIVE_ACK=1 MAKER_MAX_CONTRACTS=100000000 HARD_STOP_DOLLARS=100000000 nohup python3 -u -m unabated_edge.runner ...`), confirming the exact command from the running process's env before killing it (`ps eww <pid>`).
- [ ] **Step 4: Verify within ~2 min**: heartbeat line resumes in `bot.log`; `maker_quotes` (fresh read-only copy of `unabated_edge_maker.duckdb`) shows `touch_join` rests on the NO side of ESPARG-2/-3 at the live touch (cross-check price vs `book_snapshots` latest), and legacy `quote` rows elsewhere.
- [ ] **Step 5: Report** tickers/sides/prices/sizes rested + worst-case total to the user; update memory file `unabated_edge_wc_maker.md`.

---

## Version control

Branch `worktree-wc-maker-depth-capture` (existing worktree, live maker's cwd). One commit per task as above; no merge to main tonight; no push without user approval.

## Documentation

Task 3 (README) ships in the same branch before rollout; spec + this plan already committed; memory update in Task 4 Step 5.
