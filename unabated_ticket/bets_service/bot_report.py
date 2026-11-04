"""GET /bot.json: what the Bet Tracker's Bot tab shows about the two automated
takers (the Kalshi bot, server/taker.js + taker_orders.py, and the Novig bot,
novig_taker_desk.py).

Inputs: the taker log (taker.duckdb, read only through TakerStore), whether
each order desk is armed (and why not), the stop files, and the server
runner's /taker.json (its last pass: when it finished and what, if anything,
skipped it). Output: one JSON dict (`bot_payload`). Side effects: none; one GET
of the runner's /taker.json per call (no proxy, the service's runner timeout).

Status per bot, in order: "off" (the service has no armed desk for it: no key,
or the key was not proven locked to the bot's subaccount), "stopped" (a stop
file is present), else "running". The tracker adds the runner's view: a pass
that is late or skipped for that venue reads as a warning on a running bot.
"""
import json
from datetime import datetime
from pathlib import Path
from typing import Callable

from unabated_ticket.bets_service import taker_orders, novig_taker_desk
from unabated_ticket.bets_service.taker_store import TakerStore

RECENT_ORDERS_LIMIT = 25
STATUS_OFF = "off"
STATUS_STOPPED = "stopped"
STATUS_RUNNING = "running"


def desk_status(armed: bool, off_reason: str, stop_paths: tuple[Path, ...]) -> dict:
    """{status, detail} for one desk: off with the service's reason, stopped
    naming the stop file found, else running."""
    if not armed:
        return {"status": STATUS_OFF, "detail": off_reason}
    for stop_path in stop_paths:
        if Path(stop_path).exists():
            return {"status": STATUS_STOPPED, "detail": f"stop file {stop_path} is present"}
    return {"status": STATUS_RUNNING, "detail": None}


def runner_pass(fetch_taker_json: Callable[[], tuple[int, bytes]]) -> dict:
    """The runner's last taker pass as {finishedAt, kalshiSkipped, novigSkipped,
    rows, novigKellyBankroll, kalshiDesk, novigDesk, error}: error names why it
    could not be read (the rest null)."""
    empty = {"finishedAt": None, "kalshiSkipped": None, "novigSkipped": None, "rows": None,
             "novigKellyBankroll": None, "kalshiDesk": None, "novigDesk": None, "error": None}
    status, body = fetch_taker_json()
    try:
        parsed = json.loads(body)
    except ValueError:
        return {**empty, "error": f"the runner's /taker.json answered HTTP {status} with no JSON"}
    if status != 200 or not isinstance(parsed, dict):
        problem = parsed.get("error") if isinstance(parsed, dict) else None
        return {**empty, "error": problem or f"the runner's /taker.json answered HTTP {status}"}
    last_pass = parsed.get("pass") if isinstance(parsed.get("pass"), dict) else {}
    rules = parsed.get("rules") if isinstance(parsed.get("rules"), dict) else {}
    return {**empty, "finishedAt": last_pass.get("finishedAt"), "kalshiSkipped": last_pass.get("skipped"),
            "novigSkipped": last_pass.get("novigSkipped"), "rows": last_pass.get("rows"),
            "novigKellyBankroll": rules.get("novigKellyBankroll"),
            "kalshiDesk": parsed.get("desk"), "novigDesk": parsed.get("novigDesk")}


def bot_payload(store: TakerStore | None, kalshi_desk: dict, novig_desk: dict, runner: dict,
                window_start: datetime | None, window_end: datetime | None, generated_at: str) -> dict:
    """The /bot.json body. `kalshi_desk` / `novig_desk` are desk_status dicts,
    `runner` is runner_pass's; the activity covers [window_start, window_end)."""
    if store is None:
        return {"generatedAt": generated_at, "error": "the taker log is not open on this service",
                "bots": [], "recentOrders": [], "runner": runner}
    activity = store.bot_activity(window_start, window_end)
    open_orders = store.bot_open_orders()
    bots = [
        {"venue": "kalshi", "label": "Kalshi bot", **kalshi_desk,
         "capDollars": taker_orders.MAX_OPEN_COST, "openCost": store.open_exposure(),
         "kellyBankroll": taker_orders.KELLY_BANKROLL, "skipped": runner["kalshiSkipped"],
         "runnerDesk": runner["kalshiDesk"], **_iso_times(open_orders["kalshi"]), "activity": activity["kalshi"]},
        {"venue": "novig", "label": "Novig bot", **novig_desk,
         "capDollars": novig_taker_desk.OPEN_COST_CEILING, "openCost": store.novig_open_exposure(),
         "kellyBankroll": runner["novigKellyBankroll"], "skipped": runner["novigSkipped"],
         "runnerDesk": runner["novigDesk"], **_iso_times(open_orders["novig"]), "activity": activity["novig"]},
    ]
    recent = [_iso_times(order) for order in store.recent_bot_orders(RECENT_ORDERS_LIMIT)]
    return {"generatedAt": generated_at, "error": None,
            "window": {"start": _iso_or_none(window_start), "end": _iso_or_none(window_end)},
            "runner": {"finishedAt": runner["finishedAt"], "rows": runner["rows"], "error": runner["error"]},
            "bots": bots, "recentOrders": recent}


def parse_window(query: dict[str, list[str]]) -> tuple[datetime | None, datetime | None] | str:
    """(start, end) from ?start=<ISO instant>&end=<ISO instant> (each optional,
    timezone required), or what was wrong."""
    bounds = []
    for name in ("start", "end"):
        values = query.get(name)
        if not values:
            bounds.append(None)
            continue
        try:
            moment = datetime.fromisoformat(values[0].replace("Z", "+00:00"))
        except ValueError:
            return f"{name} must be an ISO instant like 2026-11-04T08:00:00Z, got {values[0]!r}"
        if moment.tzinfo is None:
            return f"{name} must carry a timezone, got {values[0]!r}"
        bounds.append(moment)
    start, end = bounds
    if start is not None and end is not None and end <= start:
        return f"end {end.isoformat()} must be after start {start.isoformat()}"
    return start, end


def _iso_or_none(moment: datetime | None) -> str | None:
    return None if moment is None else moment.isoformat().replace("+00:00", "Z")


def _iso_times(row: dict) -> dict:
    return {key: _iso_or_none(value) if isinstance(value, datetime) else value for key, value in row.items()}
