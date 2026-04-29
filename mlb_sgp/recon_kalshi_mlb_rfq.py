"""
Kalshi MLB RFQ recon — answers open questions O-1..O-6 from the design spec.

Pure read-mostly probe. The only side effect is creating up to 2 RFQs with
target_cost_dollars=0.01, both DELETE'd in a finally block before exit.
We never accept a quote.

Run from anywhere: python3 recon_kalshi_mlb_rfq.py
Requires: KALSHI_API_KEY_ID + KALSHI_PRIVATE_KEY_PATH in /Users/callancapitolo/NFLWork/kalshi_mm/.env
"""

import json
import os
import sys
import time
import urllib.request
import urllib.error
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path("/Users/callancapitolo/NFLWork")
KALSHI_BASE_URL = "https://api.elections.kalshi.com/trade-api/v2"

# Load creds from main repo's kalshi_mm/.env (gitignored, not in worktree)
ENV_PATH = REPO_ROOT / "kalshi_mm" / ".env"
for line in ENV_PATH.read_text().splitlines():
    line = line.strip()
    if not line or line.startswith("#") or "=" not in line:
        continue
    k, v = line.split("=", 1)
    os.environ.setdefault(k.strip(), v.strip().strip('"').strip("'"))

# auth.sign_request lives in main repo's kalshi_draft/auth.py
sys.path.insert(0, str(REPO_ROOT / "kalshi_draft"))
from auth import sign_request  # noqa: E402

API_KEY = os.environ["KALSHI_API_KEY_ID"]
PRIVATE_KEY = os.environ["KALSHI_PRIVATE_KEY_PATH"]


def api(method: str, path: str, body=None):
    """Authenticated Kalshi API call. Returns (status, body_dict_or_text, headers)."""
    ts = str(int(datetime.now(timezone.utc).timestamp() * 1000))
    full_path = f"/trade-api/v2{path.split('?')[0]}"  # signature uses path only
    sig = sign_request(PRIVATE_KEY, ts, method, full_path)

    url = f"{KALSHI_BASE_URL}{path}"
    data = json.dumps(body).encode() if body else None
    req = urllib.request.Request(url, data=data, method=method)
    req.add_header("KALSHI-ACCESS-KEY", API_KEY)
    req.add_header("KALSHI-ACCESS-SIGNATURE", sig)
    req.add_header("KALSHI-ACCESS-TIMESTAMP", ts)
    req.add_header("Content-Type", "application/json")

    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            text = r.read().decode()
            try:
                return r.status, json.loads(text), dict(r.headers)
            except json.JSONDecodeError:
                return r.status, text, dict(r.headers)
    except urllib.error.HTTPError as e:
        body_text = e.read().decode()
        try:
            return e.code, json.loads(body_text), dict(e.headers)
        except json.JSONDecodeError:
            return e.code, body_text, dict(e.headers)


def banner(msg):
    print()
    print("=" * 76)
    print(msg)
    print("=" * 76)


# Track RFQs we create so we can DELETE them at exit no matter what
RFQS_TO_CLEANUP = []


def cleanup():
    if not RFQS_TO_CLEANUP:
        return
    banner("CLEANUP — deleting any RFQs we created")
    for rfq_id in RFQS_TO_CLEANUP:
        status, body, _ = api("DELETE", f"/communications/rfqs/{rfq_id}")
        print(f"  DELETE rfq {rfq_id} → {status}")


def main():
    # ──────────────────────────────────────────────────────────────────────
    banner("PHASE 0 — auth sanity")
    status, body, _ = api("GET", "/exchange/status")
    print(f"  /exchange/status → {status}: {body if status != 200 else 'OK'}")
    if status != 200:
        print("  auth failed; aborting"); return

    # ──────────────────────────────────────────────────────────────────────
    banner("PHASE A — confirm cross-category MVE collections include MLB legs")
    # We saw earlier that KXMVECROSSCATEGORY-R and KXMVESPORTSMULTIGAMEEXTENDED-R
    # both have MLB events as eligible legs. Re-verify and pick one.
    status, body, _ = api("GET", "/multivariate_event_collections?limit=200")
    cols = body.get("multivariate_contracts", [])
    mlb_collections = []
    for c in cols:
        ae = c.get("associated_event_tickers", []) or []
        mlb_legs = [t for t in ae if "MLB" in t.upper()]
        if mlb_legs:
            mlb_collections.append((c.get("collection_ticker"), c.get("series_ticker"), len(ae), len(mlb_legs)))
    print("  MVE collections that include MLB events as legs:")
    for ct, st, n_total, n_mlb in mlb_collections:
        print(f"    {ct} (series={st})  total_events={n_total}  mlb_events={n_mlb}")

    if not mlb_collections:
        print("  no MLB-eligible MVE collections found; aborting"); return

    # Pick the first cross-category collection
    chosen_collection = mlb_collections[0][0]
    print(f"\n  → using collection: {chosen_collection}")

    # ──────────────────────────────────────────────────────────────────────
    banner("PHASE B — pick an MLB game with first pitch ≥4 hours from now")
    # Fetch open MLB game events; pick the latest first-pitch we can find
    status, body, _ = api("GET", "/events?series_ticker=KXMLBGAME&status=open&limit=50")
    events = body.get("events", []) if status == 200 else []
    print(f"  open KXMLBGAME events: {len(events)}")

    # Try each in order; we want one that is in the cross-category collection's
    # associated_event_tickers list. The associated_event_tickers list shows
    # only event_tickers (not market_tickers), so cross-reference.
    eligible_event_tickers = set(
        t for c in cols if c.get("collection_ticker") == chosen_collection
        for t in (c.get("associated_event_tickers") or [])
    )

    chosen_event = None
    for ev in events:
        et = ev.get("event_ticker")
        if et in eligible_event_tickers:
            chosen_event = ev
            break

    if not chosen_event:
        print("  no MLB event ticker found in chosen collection's eligible list")
        print(f"  first 5 open events: {[e.get('event_ticker') for e in events[:5]]}")
        print(f"  first 5 eligible MLB legs in collection: {[t for t in list(eligible_event_tickers)[:5] if 'KXMLBGAME' in t]}")
        return
    print(f"  → chosen event: {chosen_event.get('event_ticker')}  ({chosen_event.get('title')})")

    EVENT_TICKER = chosen_event["event_ticker"]

    # ──────────────────────────────────────────────────────────────────────
    banner("PHASE C — inspect leg market structures (O-1, O-2)")
    # Resolve the corresponding KXMLBSPREAD, KXMLBTOTAL, KXMLBRFI tickers
    # for this game. Conventions to test:
    #   suffix on KXMLBGAME: "26APR291940AZMIL" (date+time+away+home)
    #   spread/total/rfi event tickers should be KXMLBSPREAD-{same suffix} etc.
    suffix = EVENT_TICKER.replace("KXMLBGAME-", "")
    related = {
        "KXMLBGAME":   f"KXMLBGAME-{suffix}",
        "KXMLBSPREAD": f"KXMLBSPREAD-{suffix}",
        "KXMLBTOTAL":  f"KXMLBTOTAL-{suffix}",
        "KXMLBRFI":    f"KXMLBRFI-{suffix}",
    }
    market_summaries = {}
    for label, ev_ticker in related.items():
        status, body, _ = api("GET", f"/markets?event_ticker={ev_ticker}&limit=50")
        if status != 200:
            print(f"  {label}: GET event {ev_ticker} → {status}: {body}")
            continue
        markets = body.get("markets", []) if isinstance(body, dict) else []
        market_summaries[label] = markets
        print(f"\n  {label}  event={ev_ticker}  markets={len(markets)}")
        for m in markets[:6]:
            print(f"    ticker={m.get('ticker')}")
            print(f"      yes_sub_title={m.get('yes_sub_title')!r}  no_sub_title={m.get('no_sub_title')!r}")
            print(f"      title={m.get('title')!r}")
            print(f"      strike_type={m.get('strike_type')}  custom_strike={m.get('custom_strike')}")
            print(f"      yes_bid/ask=${m.get('yes_bid_dollars')}/${m.get('yes_ask_dollars')}  "
                  f"no_bid/ask=${m.get('no_bid_dollars')}/${m.get('no_ask_dollars')}")
        if len(markets) > 6:
            print(f"    ... ({len(markets) - 6} more)")

    # ──────────────────────────────────────────────────────────────────────
    banner("PHASE D — mint a 2-leg combo ticker via MVE lookup (O-5)")
    # Pick: home spread YES + total over YES.
    # We need 2 legs from the chosen collection's eligible event list.
    spread_markets = market_summaries.get("KXMLBSPREAD", [])
    total_markets = market_summaries.get("KXMLBTOTAL", [])
    if not spread_markets or not total_markets:
        print("  missing spread or total markets; can't mint combo"); return

    # Default leg: first spread market YES + first total market YES.
    # We'll try whatever the first markets are, with side='yes'. If there are
    # multiple spread/total markets per event (e.g., per team or per line),
    # we'll see how the lookup endpoint reacts.
    leg1 = spread_markets[0]
    leg2 = total_markets[0]
    selected = [
        {"market_ticker": leg1["ticker"], "event_ticker": leg1["event_ticker"], "side": "yes"},
        {"market_ticker": leg2["ticker"], "event_ticker": leg2["event_ticker"], "side": "yes"},
    ]
    print(f"  selected legs:")
    for s in selected:
        print(f"    {s}")

    # Try the documented lookup path. From the docs the path is:
    # POST /multivariate_event_collections/{collection_ticker}
    # (recon-confirmed 2026-04-27: no /lookup suffix — docs were misleading)
    lookup_path = f"/multivariate_event_collections/{chosen_collection}"
    status, body, _ = api("POST", lookup_path, body={"selected_markets": selected})
    print(f"\n  POST {lookup_path}")
    print(f"    status={status}")
    print(f"    body={json.dumps(body, indent=2) if isinstance(body, dict) else body}")

    if status != 200 or not isinstance(body, dict) or not body.get("market_ticker"):
        print("\n  lookup did not return a market_ticker; aborting recon")
        return

    combo_ticker = body["market_ticker"]
    combo_event = body.get("event_ticker")
    print(f"\n  → minted combo: market_ticker={combo_ticker}  event_ticker={combo_event}")

    # ──────────────────────────────────────────────────────────────────────
    banner("PHASE E — submit RFQ for the combo (target_cost=$0.01)")
    rfq_body = {
        "market_ticker": combo_ticker,
        "rest_remainder": False,
        "target_cost_dollars": "0.01",
        "replace_existing": True,
    }
    status, body, hdrs = api("POST", "/communications/rfqs", body=rfq_body)
    print(f"  POST /communications/rfqs → {status}")
    print(f"    body={json.dumps(body, indent=2) if isinstance(body, dict) else body}")
    if status != 200 or not isinstance(body, dict) or not body.get("id"):
        print("  RFQ create failed; aborting"); return
    rfq_id = body["id"]
    RFQS_TO_CLEANUP.append(rfq_id)
    print(f"  → rfq_id={rfq_id}")

    # ──────────────────────────────────────────────────────────────────────
    banner("PHASE F — poll for quotes (~20s @ 2s cadence) (O-4)")
    quote_history = []
    poll_start = time.time()
    while time.time() - poll_start < 20:
        status, body, _ = api("GET", f"/communications/quotes?rfq_id={rfq_id}")
        quotes = body.get("quotes", []) if isinstance(body, dict) else []
        elapsed = round(time.time() - poll_start, 1)
        if quotes:
            for q in quotes:
                key = (q.get("id"), q.get("status"))
                if key not in [(qq.get("id"), qq.get("status")) for qq in quote_history]:
                    quote_history.append(q)
                    print(f"  t+{elapsed}s  quote {q.get('id')} status={q.get('status')} "
                          f"yes_bid=${q.get('yes_bid_dollars')} no_bid=${q.get('no_bid_dollars')}")
        time.sleep(2)
    if not quote_history:
        print("  no quotes received in 20s window")
    else:
        print(f"\n  total distinct quotes: {len({q.get('id') for q in quote_history})}")

    # ──────────────────────────────────────────────────────────────────────
    banner("PHASE G — test replace_existing semantics (O-3)")
    # Submit a 2nd RFQ on the SAME market_ticker with replace_existing=True.
    # Expectation: prior RFQ should be auto-cancelled.
    status2, body2, _ = api("POST", "/communications/rfqs", body=rfq_body)
    print(f"  POST /communications/rfqs (replace_existing=True) → {status2}")
    print(f"    body={json.dumps(body2, indent=2) if isinstance(body2, dict) else body2}")
    if isinstance(body2, dict) and body2.get("id"):
        new_rfq_id = body2["id"]
        RFQS_TO_CLEANUP.append(new_rfq_id)
        # Check status of the original
        status3, body3, _ = api("GET", f"/communications/rfqs/{rfq_id}")
        print(f"\n  GET /communications/rfqs/{rfq_id} (original) → {status3}")
        if isinstance(body3, dict):
            print(f"    status={body3.get('rfq', {}).get('status')}  "
                  f"cancellation_reason={body3.get('rfq', {}).get('cancellation_reason')}")


if __name__ == "__main__":
    try:
        main()
    finally:
        cleanup()
