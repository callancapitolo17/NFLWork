"""Read-only look at how Novig, Polymarket US, ProphetX and Unabated's free feed
show an NFL game WHILE IT IS BEING PLAYED, so live Place can be built on real
shapes instead of guesses. Run it during a live game (a break is best):

    cd ~/NFLWork/unabated_ticket/deploy
    docker compose exec -T bets python -m unabated_ticket.bets_service.live_recon

What it answers, per venue:
  Unabated   whether the free feed's live rows (lg1:pt*:live) carry each venue's
             line with its venue ids (Kalshi sourceKey, Novig / ProphetX sourceData).
  Novig      the event and market statuses of an in-play NFL game, and its
             market `fee` object (pregame Place refuses fee.charged ALWAYS).
  ProphetX   the in-play event status (pregame is "not_started") and whether
             its markets and offers stay up.
  Polymarket whether in-play games stay in the public event list, and their
             market state and book.

Inputs:  the bets service's own config (.env): the bot's Novig key
         (NOVIG_TAKER_*, catalog GETs only, through novig_catalog_reads'
         allowlist) and PROPHETX_ACCESS_KEY / _SECRET_KEY. Polymarket US and
         Unabated are public reads.
Outputs: one JSON file next to bets.duckdb (live_recon_<UTC stamp>.json) with
         raw samples, cut to MAX_SAMPLE_CHARS each, and a short summary on stdout.
Side effects: GETs only, no order, no DuckDB write. At most 33 Novig catalog reads, at most
         1 + MAX_GAMES ProphetX reads on Cal's key (it shares that key's 1 req/s
         with the running service, so a ProphetX Place clicked in the same few
         seconds may see a 429), about 4 Polymarket US and 1 Unabated read.
         Each venue fails on its own and is recorded under "error".
"""
import argparse
import json
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path

import requests

from unabated_ticket.bets_service import config
from unabated_ticket.bets_service.novig_catalog_reads import NovigCatalogReads, validate_read_request
from unabated_ticket.bets_service.novig_orders import NovigOrderClient
from unabated_ticket.bets_service.prophetx_client import ProphetXClient

NFL_LEAGUE_ID = 1
UNABATED_SNAPSHOT_URL = "https://content.unabated.com/markets/v2/league/{league_id}/odds.json"
POLYMARKET_GATEWAY_URL = "https://gateway.polymarket.us"
NOVIG_LEAGUE = "NFL"
# One page, the runner's limit (runner.js NOVIG_MARKETS_PAGE_LIMIT); no eventStatus filter, so every status shows.
NOVIG_MARKETS_PAGE_LIMIT = "5000"
# Event reads spend the bot key's read bucket (16/s) the runner shares.
MAX_NOVIG_EVENT_READS = 30
# Unabated book ids of the venues live Place would cover (extension/place.js).
VENUE_BOOK_IDS = {"kalshi": 105, "novig": 89, "polymarket_us": 111, "prophetx": 66}
# A game that started this long ago or less counts as possibly in play.
IN_PLAY_WINDOW = timedelta(hours=4, minutes=30)
MAX_GAMES = 2
MAX_SAMPLE_CHARS = 4000
HTTP_TIMEOUT_SEC = 20
# Extra gap between this script's own ProphetX reads, on top of the client's pacer.
PROPHETX_EXTRA_GAP_SEC = 2.0


def sample(value: object) -> object:
    """value as is when small, else its JSON cut to MAX_SAMPLE_CHARS (a string)."""
    text = json.dumps(value, default=str)
    return value if len(text) <= MAX_SAMPLE_CHARS else text[:MAX_SAMPLE_CHARS] + "...(cut)"


def scalars_of(record: dict) -> dict:
    """The top-level non-container fields of a record: where status words live."""
    return {key: value for key, value in record.items() if not isinstance(value, (dict, list))}


def parse_time(text: object) -> datetime | None:
    if not isinstance(text, str) or not text:
        return None
    try:
        parsed = datetime.fromisoformat(text.replace("Z", "+00:00"))
    except ValueError:
        return None
    return parsed if parsed.tzinfo else parsed.replace(tzinfo=timezone.utc)


def started_recently(start: datetime | None, now: datetime) -> bool:
    return start is not None and now - IN_PLAY_WINDOW <= start <= now


def get_json(url: str, params: dict | None = None) -> object:
    reply = requests.get(url, params=params, headers={"Accept": "application/json"}, timeout=HTTP_TIMEOUT_SEC)
    if reply.status_code != 200:
        raise RuntimeError(f"GET {url}: HTTP {reply.status_code} {reply.text[:200]}")
    return reply.json()


# ---- Unabated -------------------------------------------------------------------

def recon_unabated() -> dict:
    """Live rows of the NFL snapshot: per game, each venue's line fields."""
    snapshot = get_json(UNABATED_SNAPSHOT_URL.format(league_id=NFL_LEAGUE_ID))
    odds = snapshot.get("odds") if isinstance(snapshot, dict) else None
    if not isinstance(odds, dict):
        raise RuntimeError("Unabated snapshot: no `odds` map")
    games: dict[str, dict] = {}
    for league_key, rows in odds.items():
        if ":live" not in league_key or not isinstance(rows, list):
            continue
        for row in rows:
            if not isinstance(row, dict):
                continue
            game = games.setdefault(str(row.get("eventId")), {
                "eventName": row.get("eventName"), "eventStart": row.get("eventStart"), "rowKeys": [], "venues": {}})
            game["rowKeys"].append(row.get("key"))
            for side_key, books in (row.get("sides") or {}).items():
                for venue, book_id in VENUE_BOOK_IDS.items():
                    line = (books or {}).get(f"ms{book_id}")
                    if not isinstance(line, dict):
                        continue
                    found = game["venues"].setdefault(venue, {"lines": 0, "withSourceKey": 0, "withSourceData": 0, "examples": []})
                    found["lines"] += 1
                    found["withSourceKey"] += 1 if line.get("sourceKey") else 0
                    found["withSourceData"] += 1 if line.get("sourceData") else 0
                    if len(found["examples"]) < 2:
                        found["examples"].append({"rowKey": row.get("key"), "side": side_key, "line": sample(line)})
    return {"liveLeagueKeys": [key for key in odds if ":live" in key], "games": games}


# ---- Novig ----------------------------------------------------------------------

def novig_read(reads: NovigCatalogReads, path: str, query: dict | None = None) -> dict:
    checked = validate_read_request({"path": path, "query": query or {}})
    if isinstance(checked, str):
        raise RuntimeError(f"novig read refused: {checked}")
    return reads.read(*checked)


def recon_novig(now: datetime) -> dict:
    """NFL catalog markets of every status, then one in-play game's event, a market and its book."""
    if not config.NOVIG_TAKER_KEY_ID or not config.NOVIG_TAKER_PRIVATE_KEY_PATH:
        raise RuntimeError("no bot Novig key configured (NOVIG_TAKER_KEY_ID / NOVIG_TAKER_PRIVATE_KEY_PATH)")
    reads = NovigCatalogReads(NovigOrderClient(config.NOVIG_TAKER_KEY_ID, config.NOVIG_TAKER_PRIVATE_KEY_PATH,
                                               config.NOVIG_API_URL))
    listing = novig_read(reads, "/v3/catalog/markets", {"league": NOVIG_LEAGUE, "limit": NOVIG_MARKETS_PAGE_LIMIT})
    body = listing["body"]
    markets = body.get("items") if isinstance(body, dict) else None
    result: dict = {"listStatus": listing["status"], "listShape": sample(scalars_of(body)) if isinstance(body, dict) else type(body).__name__}
    if not isinstance(markets, list):
        result["listSample"] = sample(body)
        return result
    result["marketCount"] = len(markets)
    result["marketStatuses"] = sorted({str(market.get("status")) for market in markets if isinstance(market, dict)})
    result["feeCharged"] = sorted({str((market.get("fee") or {}).get("charged")) for market in markets if isinstance(market, dict)})
    event_ids = []
    for market in markets:
        if isinstance(market, dict) and market.get("eventId") not in event_ids:
            event_ids.append(market.get("eventId"))
    result["events"] = []
    result["moreListPages"] = bool(isinstance(body, dict) and body.get("next"))
    for event_id in event_ids[:MAX_NOVIG_EVENT_READS]:
        event_read = novig_read(reads, f"/v3/catalog/events/{event_id}")
        event = event_read["body"] if isinstance(event_read["body"], dict) else {}
        start = parse_time(event.get("startTime") or event.get("scheduledStart") or event.get("start"))
        entry = {"eventId": event_id, "readStatus": event_read["status"], "fields": sample(scalars_of(event))}
        result["events"].append(entry)
        if not started_recently(start, now) and str(event.get("status")) == "OPEN_PREGAME":
            continue
        game_markets = [market for market in markets if isinstance(market, dict) and market.get("eventId") == event_id]
        entry["inPlay"] = True
        entry["marketSamples"] = [sample(market) for market in game_markets[:3]]
        if game_markets:
            market_id = game_markets[0].get("marketId")
            entry["marketRead"] = sample(novig_read(reads, f"/v3/catalog/markets/{market_id}"))
            entry["bookRead"] = sample(novig_read(reads, f"/v3/catalog/markets/{market_id}/book"))
        if sum(1 for other in result["events"] if other.get("inPlay")) >= MAX_GAMES:
            break
    return result


# ---- ProphetX -------------------------------------------------------------------

def recon_prophetx(now: datetime) -> dict:
    """Every listed ProphetX game's status words, then the markets of up to MAX_GAMES in-play NFL games."""
    if not (config.PROPHETX_ACCESS_KEY and config.PROPHETX_SECRET_KEY):
        raise RuntimeError("no ProphetX key configured (PROPHETX_ACCESS_KEY / PROPHETX_SECRET_KEY)")
    client = ProphetXClient(config.PROPHETX_ACCESS_KEY, config.PROPHETX_SECRET_KEY, config.PROPHETX_BASE_URL)
    events = client.sport_events()
    result: dict = {"eventCount": len(events), "statuses": sorted({str(event.get("status")) for event in events}),
                    "inPlay": []}
    in_play = [event for event in events
               if started_recently(parse_time(event.get("scheduled")), now) or event.get("status") not in (None, "not_started")]
    result["inPlayFields"] = [sample(scalars_of(event)) for event in in_play[:10]]
    nfl_first = sorted(in_play, key=lambda event: "nfl" not in json.dumps(event, default=str).lower())
    for event in nfl_first[:MAX_GAMES]:
        time.sleep(PROPHETX_EXTRA_GAP_SEC)
        event_id = event.get("event_id")
        try:
            markets = client.markets(int(event_id))
        except Exception as error:  # recorded, the next game still runs
            result["inPlay"].append({"eventId": event_id, "error": str(error)[:300]})
            continue
        result["inPlay"].append({
            "eventId": event_id, "name": event.get("name"), "status": event.get("status"), "marketCount": len(markets),
            "marketStatuses": sorted({str(market.get("status")) for market in markets}),
            "marketSamples": [sample(market) for market in markets[:2]]})
    return result


# ---- Polymarket US --------------------------------------------------------------

def recon_polymarket(now: datetime) -> dict:
    """Public events that started within IN_PLAY_WINDOW, with each one's market states and one book."""
    window = {"limit": "100", "active": "true",
              "startTimeMin": (now - IN_PLAY_WINDOW).strftime("%Y-%m-%dT%H:%M:%SZ"),
              "startTimeMax": now.strftime("%Y-%m-%dT%H:%M:%SZ")}
    body = get_json(f"{POLYMARKET_GATEWAY_URL}/v1/events", window)
    events = body.get("events") if isinstance(body, dict) else None
    if not isinstance(events, list):
        raise RuntimeError(f"Polymarket US events: no events list in {str(body)[:200]}")
    result: dict = {"startedEventCount": len(events), "events": []}
    football = [event for event in events if "football" in json.dumps(event.get("tags") or event.get("sport") or event, default=str).lower()]
    for event in (football or events)[:MAX_GAMES]:
        markets = [market for market in event.get("markets") or [] if isinstance(market, dict)]
        entry = {"slug": event.get("slug"), "fields": sample(scalars_of(event)),
                 "marketStates": sorted({f"{market.get('sportsMarketType')}:{market.get('state') or market.get('status')}" for market in markets}),
                 "marketSample": sample(markets[0]) if markets else None}
        if markets and markets[0].get("slug"):
            entry["book"] = sample(get_json(f"{POLYMARKET_GATEWAY_URL}/v1/markets/{markets[0]['slug']}/book"))
        result["events"].append(entry)
    return result


# ---- run ------------------------------------------------------------------------

def run_venue(name: str, recon, *args) -> dict:
    try:
        return recon(*args)
    except Exception as error:  # one venue failing must not lose the others
        return {"error": f"{type(error).__name__}: {str(error)[:500]}"}


def summary_lines(report: dict) -> list[str]:
    lines = [f"live recon at {report['at']}"]
    unabated = report["unabated"]
    if "error" in unabated:
        lines.append(f"unabated: ERROR {unabated['error']}")
    else:
        lines.append(f"unabated: {len(unabated['games'])} live games, keys {unabated['liveLeagueKeys']}")
        for event_id, game in unabated["games"].items():
            venues = ", ".join(f"{venue} {found['lines']} lines ({found['withSourceKey']} sourceKey, {found['withSourceData']} sourceData)"
                               for venue, found in game["venues"].items()) or "no venue lines"
            lines.append(f"  {event_id} {game['eventName']}: {venues}")
    novig = report["novig"]
    lines.append(f"novig: ERROR {novig['error']}" if "error" in novig else
                 f"novig: {novig.get('marketCount')} markets, statuses {novig.get('marketStatuses')}, fee.charged {novig.get('feeCharged')}, "
                 f"in-play events {sum(1 for event in novig.get('events', []) if event.get('inPlay'))}")
    prophetx = report["prophetx"]
    lines.append(f"prophetx: ERROR {prophetx['error']}" if "error" in prophetx else
                 f"prophetx: {prophetx['eventCount']} events, statuses {prophetx['statuses']}, "
                 + "; ".join(f"{game.get('name')} {game.get('status')} markets {game.get('marketStatuses')}" for game in prophetx["inPlay"]))
    polymarket = report["polymarket_us"]
    lines.append(f"polymarket_us: ERROR {polymarket['error']}" if "error" in polymarket else
                 f"polymarket_us: {polymarket['startedEventCount']} started events listed; "
                 + "; ".join(f"{event['slug']} {event['marketStates']}" for event in polymarket["events"]))
    return lines


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out-dir", default=str(config.DB_PATH.parent), help="where the JSON report goes (default: next to bets.duckdb)")
    args = parser.parse_args()
    now = datetime.now(timezone.utc)
    report = {"at": now.isoformat(),
              "unabated": run_venue("unabated", recon_unabated),
              "novig": run_venue("novig", recon_novig, now),
              "prophetx": run_venue("prophetx", recon_prophetx, now),
              "polymarket_us": run_venue("polymarket_us", recon_polymarket, now)}
    out_path = Path(args.out_dir) / f"live_recon_{now.strftime('%Y%m%dT%H%M%SZ')}.json"
    out_path.write_text(json.dumps(report, indent=1, default=str))
    print("\n".join(summary_lines(report)))
    print(f"report: {out_path}")


if __name__ == "__main__":
    main()
