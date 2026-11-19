"""Kalshi NFL recon for the SGP compare Kalshi adapter.

Read-mostly. The only side effect is ONE RFQ (~$100, two legs) created on a
combo, polled for quotes for ~8s, then DELETEd in a finally block. It never
accepts a quote.

Run on the Mac from the repo root:
    python3 sgp_compare/tools/recon_kalshi.py
Credentials: kalshi_mlb_rfq/.env (KALSHI_API_KEY_ID, KALSHI_PRIVATE_KEY_PATH,
KALSHI_USER_ID), the same as the MLB taker.

Writes raw JSON to /tmp/kalshi_sgp_recon/ with every
user / creator id replaced by "REDACTED".
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from kalshi_mlb_rfq import config
from kalshi_mlb_rfq.auth_client import api

COLLECTION = "KXMVECROSSCATEGORY-R"
OUT = Path("/tmp/kalshi_sgp_recon")   # raw; the tests use trimmed copies (tests/fixtures/kalshi_*.json)
RFQ_DOLLARS = 100.0
QUOTE_WINDOW_S = 8.0


def scrub(obj):
    if isinstance(obj, dict):
        return {k: ("REDACTED" if ("user_id" in k or "creator_id" in k) and v else scrub(v))
                for k, v in obj.items()}
    if isinstance(obj, list):
        return [scrub(x) for x in obj]
    return obj


def save(name: str, obj) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{name}.json").write_text(json.dumps(scrub(obj), indent=1, sort_keys=True))
    print(f"  wrote {name}.json")


def get(path: str):
    status, body, _ = api("GET", path)
    if status != 200:
        print(f"  GET {path} -> {status}: {str(body)[:300]}")
        return None
    return body


def main() -> None:
    print("== auth")
    if get("/exchange/status") is None:
        sys.exit("auth failed")

    print("== collection", COLLECTION)
    col = get(f"/multivariate_event_collections/{COLLECTION}")
    if col is None:
        sys.exit("no collection")
    c = col.get("multivariate_contract", col)
    nfl_assoc = [a for a in (c.get("associated_events") or []) if "KXNFL" in json.dumps(a)]
    nfl_tickers = [t for t in (c.get("associated_event_tickers") or []) if t.startswith("KXNFL")]
    trimmed = {k: v for k, v in c.items() if k not in ("associated_events", "associated_event_tickers")}
    trimmed["associated_events_nfl"] = nfl_assoc
    trimmed["associated_event_tickers_nfl"] = nfl_tickers
    trimmed["associated_events_total"] = len(c.get("associated_events") or [])
    save("collection_nfl", trimmed)
    nfl_series = sorted({t.split("-")[0] for t in nfl_tickers} |
                        {str(a.get("ticker", a.get("event_ticker", ""))).split("-")[0] for a in nfl_assoc})
    print("  NFL series in collection:", nfl_series)

    print("== NFL games")
    games = get("/events?series_ticker=KXNFLGAME&status=open&limit=100&with_nested_markets=true") or {}
    events = games.get("events") or []
    print(f"  open KXNFLGAME events: {len(events)}")
    save("games_kxnflgame", {"events": events[:4], "cursor": games.get("cursor")})
    if not events:
        sys.exit("no open NFL games")
    # Soonest game that is in the collection
    elig = set(nfl_tickers) | {a.get("event_ticker") or a.get("ticker") for a in nfl_assoc}
    events.sort(key=lambda e: e.get("strike_date") or e.get("expected_expiration_time") or "")
    game = next((e for e in events if e["event_ticker"] in elig), events[0])
    key = game["event_ticker"].split("-", 1)[1]
    print("  chosen game:", game["event_ticker"], game.get("title"))

    print("== series fee info")
    for s in ("KXNFLGAME", "KXNFLSPREAD", "KXNFLTOTAL"):
        b = get(f"/series/{s}")
        if b:
            save(f"series_{s}", b)

    print("== every event for this game")
    game_events = sorted({t for t in elig if t and t.endswith(key)})
    if not game_events:
        game_events = [f"{s}-{key}" for s in nfl_series if s]
    print(f"  {len(game_events)} events:", game_events)
    markets_by_event: dict[str, list[dict]] = {}
    for et in game_events:
        b = get(f"/events/{et}?with_nested_markets=true")
        if not b:
            continue
        ev = b.get("event", b)
        mk = ev.get("markets") or b.get("markets") or []
        markets_by_event[et] = mk
        save(f"event_{et}", b)
        print(f"   {et}: {len(mk)} markets")

    def near_half(et: str):
        mk = [m for m in markets_by_event.get(et, []) if m.get("yes_ask_dollars")]
        mk.sort(key=lambda m: abs(float(m["yes_ask_dollars"]) - 0.5))
        return mk[0] if mk else None

    spread = near_half(f"KXNFLSPREAD-{key}")
    total = near_half(f"KXNFLTOTAL-{key}")
    print("== orderbooks")
    for m in [x for x in (spread, total) if x]:
        b = get(f"/markets/{m['ticker']}/orderbook")
        if b:
            save(f"orderbook_{m['ticker']}", b)

    if not (spread and total):
        print("  no spread/total pair; skipping RFQ test")
        return

    print("== ONE test RFQ (never accepted)")
    selected = [
        {"market_ticker": spread["ticker"], "event_ticker": f"KXNFLSPREAD-{key}", "side": "yes"},
        {"market_ticker": total["ticker"], "event_ticker": f"KXNFLTOTAL-{key}", "side": "yes"},
    ]
    status, mint, _ = api("POST", f"/multivariate_event_collections/{COLLECTION}",
                          body={"selected_markets": selected})
    save("mint", {"request": selected, "status": status, "response": mint})
    if status != 200:
        print("  mint failed", status, mint)
        return
    combo = mint["market_ticker"]
    print("  combo:", combo)
    cm = get(f"/markets/{combo}")
    if cm:
        save("combo_market", cm)

    rfq_id = None
    try:
        t0 = time.time()
        status, resp, _ = api("POST", "/communications/rfqs", body={
            "market_ticker": combo, "rest_remainder": False, "replace_existing": True,
            "target_cost_dollars": f"{RFQ_DOLLARS:.2f}"})
        save("rfq_create", {"status": status, "response": resp})
        if status not in (200, 201):
            print("  create_rfq failed", status, resp)
            return
        rfq_id = resp["id"]
        polls = []
        user = config.KALSHI_USER_ID
        q = f"/communications/quotes?rfq_id={rfq_id}" + (f"&rfq_creator_user_id={user}" if user else "")
        while time.time() - t0 < QUOTE_WINDOW_S:
            time.sleep(0.5)
            b = get(q) or {}
            quotes = b.get("quotes") or []
            polls.append({"t": round(time.time() - t0, 2), "quotes": quotes})
            print(f"  t={polls[-1]['t']}s quotes={len(quotes)}")
        save("rfq_quote_polls", polls)
        b = get(f"/communications/rfqs/{rfq_id}")
        if b:
            save("rfq_get", b)
    finally:
        if rfq_id:
            status, _body, _ = api("DELETE", f"/communications/rfqs/{rfq_id}")
            print(f"  DELETE rfq {rfq_id} -> {status}")


if __name__ == "__main__":
    main()
