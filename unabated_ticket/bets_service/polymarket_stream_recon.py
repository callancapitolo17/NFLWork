"""READ-ONLY recon of Polymarket US's markets websocket (2026-11-05).

Run inside the bets container (its config holds POLYMARKET_US_KEY_ID /
_SECRET_KEY, and it has websocket-client):  git show origin/claude/project-thread-d1vreh:unabated_ticket/bets_service/polymarket_stream_recon.py | docker compose exec -T bets python -
Opens ONE websocket to wss://api.polymarket.us/v1/ws/markets signed like the
REST calls, subscribes order books for a few market slugs, prints what comes
back for ~60 s, then closes. Sends no orders, writes nothing, prints no key.
"""
import base64
import json
import time

import websocket
from unabated_ticket.bets_service import config
from unabated_ticket.bets_service.polymarket_orders import signer_of_secret_key

URL = "wss://api.polymarket.us/v1/ws/markets"
PATH = "/v1/ws/markets"
SLUGS_FIRST = ["tsc-cfb-duke-ncst-2026-11-07-total-48pt5", "tsc-cfb-okst-kanst-2026-11-07-total-58pt5"]
SLUG_LATER = "tsc-cfb-pennst-wash-2026-11-07-total-38pt5"
RUN_SEC = 60
FRAMES_PRINTED_PER_KIND = 3

# The key comes from the bets service's own .env files, the way the service reads it.
key_id = config.POLYMARKET_US_KEY_ID
if not key_id or not config.POLYMARKET_US_SECRET_KEY:
    raise SystemExit("no POLYMARKET_US_KEY_ID / POLYMARKET_US_SECRET_KEY in the bets service config")
signer = signer_of_secret_key(config.POLYMARKET_US_SECRET_KEY)
stamp = str(int(time.time() * 1000))
signature = base64.b64encode(signer.sign(f"{stamp}GET{PATH}".encode())).decode()
headers = [f"X-PM-Access-Key: {key_id}", f"X-PM-Timestamp: {stamp}", f"X-PM-Signature: {signature}"]

started = time.time()
try:
    ws = websocket.create_connection(URL, header=headers, timeout=10)
except Exception as error:
    print("CONNECT FAILED:", type(error).__name__, str(error)[:400])
    raise SystemExit(1)
print(f"connected in {time.time() - started:.2f}s, key {key_id[:8]}")
ws.send(json.dumps({"subscribe": {"requestId": "recon-md-1", "subscriptionType": "SUBSCRIPTION_TYPE_MARKET_DATA", "marketSlugs": SLUGS_FIRST}}))
ws.send(json.dumps({"subscribe": {"requestId": "recon-lite-1", "subscriptionType": "SUBSCRIPTION_TYPE_MARKET_DATA_LITE", "marketSlugs": SLUGS_FIRST[:1]}}))

counts, printed, per_slug, last_heartbeat, added_later = {}, {}, {}, None, False
ws.settimeout(2)
while time.time() - started < RUN_SEC:
    if not added_later and time.time() - started > 20:
        ws.send(json.dumps({"subscribe": {"requestId": "recon-md-2", "subscriptionType": "SUBSCRIPTION_TYPE_MARKET_DATA", "marketSlugs": [SLUG_LATER]}}))
        added_later = True
        print(f"[{time.time() - started:5.1f}s] sent second subscription recon-md-2 for {SLUG_LATER}")
    try:
        raw = ws.recv()
    except websocket.WebSocketTimeoutException:
        continue
    except Exception as error:
        print(f"[{time.time() - started:5.1f}s] RECV ENDED:", type(error).__name__, str(error)[:300])
        break
    try:
        message = json.loads(raw)
    except ValueError:
        print("NON-JSON:", str(raw)[:300])
        continue
    kind = next((k for k in ("heartbeat", "error", "marketData", "marketDataLite", "trade") if k in message), "other:" + ",".join(sorted(message))[:80])
    counts[kind] = counts.get(kind, 0) + 1
    at = time.time() - started
    if kind == "heartbeat":
        if last_heartbeat is not None and counts[kind] <= 3:
            print(f"[{at:5.1f}s] heartbeat gap {at - last_heartbeat:.1f}s: {json.dumps(message)[:200]}")
        last_heartbeat = at
        continue
    if kind == "marketData":
        data = message["marketData"]
        slug = data.get("marketSlug")
        per_slug[slug] = per_slug.get(slug, 0) + 1
        summary = {"requestId": message.get("requestId"), "slug": slug, "state": data.get("state"),
                   "bids": len(data.get("bids") or []), "offers": len(data.get("offers") or []),
                   "topBid": (data.get("bids") or [None])[0], "topOffer": (data.get("offers") or [None])[0],
                   "transactTime": data.get("transactTime"), "keys": sorted(data)}
        if printed.get(kind, 0) < FRAMES_PRINTED_PER_KIND or per_slug[slug] == 1:
            print(f"[{at:5.1f}s] marketData {json.dumps(summary)[:700]}")
            if printed.get(kind, 0) == 0:
                print("   FULL FIRST FRAME:", json.dumps(message)[:2500])
            printed[kind] = printed.get(kind, 0) + 1
        continue
    if printed.get(kind, 0) < FRAMES_PRINTED_PER_KIND:
        print(f"[{at:5.1f}s] {kind}: {json.dumps(message)[:800]}")
        printed[kind] = printed.get(kind, 0) + 1
ws.close()
print("counts by kind:", counts)
print("marketData frames per slug:", per_slug)
