import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from kalshi_common import auth_client
from unabated_edge import config

def init():
    auth_client.configure(api_key_id=config.KALSHI_API_KEY_ID,
        private_key_path=config.KALSHI_PRIVATE_KEY_PATH, base_url=config.KALSHI_BASE_URL,
        project_root=str(Path(__file__).resolve().parent.parent.parent))

def list_events(series_ticker: str) -> list[dict]:
    status, body, _ = auth_client.api("GET", f"/events?series_ticker={series_ticker}&status=open&with_nested_markets=true")
    return body.get("events", []) if isinstance(body, dict) else []

def best_yes_ask(market_ticker: str) -> float | None:
    status, body, _ = auth_client.api("GET", f"/markets/{market_ticker}/orderbook")
    if status != 200 or not isinstance(body, dict): return None
    no = (body.get("orderbook") or {}).get("no") or []
    if not no: return None
    return round((100 - max(l[0] for l in no)) / 100.0, 2)
