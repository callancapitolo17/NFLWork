"""ProphetX's live prices for the Edges list's live price overlay (Cal,
2026-11-05, option A): the server runner's live price pass names the board
games it lists ProphetX lines on (POST /prophetx_live.json), and this reader
reads each game's markets off ProphetX on a key pair of its OWN
(PROPHETX_LIVE_ACCESS_KEY / _SECRET_KEY), so the 1 request a second ProphetX
allows per key stays whole for the bet source and Place on PROPHETX_*.

The ask never waits: it marks the games wanted and answers each game's last
read. A thread reads one wanted game at a time, the one read longest ago
first, at most one ProphetX call every READ_GAP_SEC (the old VM address was
banned by ProphetX's CloudFront after a scraper read it hard, so the reader
stays well under the key's limit; the client backs off any refusal itself).
A game no longer asked for in WANTED_FOR_SEC is dropped.

Per game: ProphetX's game is found as Place finds it (prophetx_place.
event_of_line: Unabated's strike ids, else league, start and teams) and kept
EVENT_REUSE_SEC; its teams are checked as Place checks them (a total survives
home and away swapped; any other mismatch leaves the game unread); then ONE
GET /v4/mm/get_markets, every side read by prophetx_markets.lines_of_event and
priced at its best offer (best_offer_of; ProphetX's fee is unpublished, so not
counted, as Place does).

Reads only: no order, no DB. Network: ProphetX's Trading API on the live key
pair, and the runner's GET /team_keys.json for the team check.
"""
import logging
import threading
import time
from collections.abc import Callable
from datetime import datetime, timezone

from unabated_ticket.bets_service.prophetx_client import ProphetXBlocked, ProphetXClient, ProphetXError
from unabated_ticket.bets_service.prophetx_markets import best_offer_of, lines_of_event
from unabated_ticket.bets_service.prophetx_place import (
    PROPHETX_EVENT_NOT_STARTED, REVERSED_TEAMS_WORDS, TeamKeysReader, TeamKeysUnread, event_of_line, teams_mismatch,
)

log = logging.getLogger(__name__)

# Seconds between the reader's ProphetX calls: half the key's 1-a-second limit.
READ_GAP_SEC = 2.0
# A game the runner stopped asking for this long ago is no longer read (its pass runs every 30 s).
WANTED_FOR_SEC = 90.0
# How long a game's ProphetX event (found off 1-2 calls) is reused before it is found again.
EVENT_REUSE_SEC = 10 * 60.0
# How long the reader idles when nothing is wanted or ProphetX is refusing it.
IDLE_SEC = 1.0
# The most games one ask may name.
MAX_GAMES_PER_ASK = 300


def validate_prophetx_live_request(body: object) -> list[dict] | str:
    """A POST /prophetx_live.json body's games, or what was expected and found.
    {games: [{gameKey, eventStart, strikeIds: [..], line: {league, awayTeamId,
    homeTeamId, awayTeam, homeTeam}}]}"""
    if not isinstance(body, dict) or not isinstance(body.get("games"), list):
        return f"body must be {{games: [..]}}, got {type(body).__name__}"
    games = body["games"]
    if len(games) > MAX_GAMES_PER_ASK:
        return f"at most {MAX_GAMES_PER_ASK} games per ask, got {len(games)}"
    for game in games:
        if not isinstance(game, dict) or not isinstance(game.get("gameKey"), str) or not game["gameKey"]:
            return f"each game needs a non-empty gameKey, got {str(game)[:120]}"
        if not isinstance(game.get("line"), dict) or not isinstance(game.get("strikeIds"), list):
            return f"game {game['gameKey']} needs a line object and a strikeIds list"
    return games


def priced_lines(lines: list[dict]) -> list[dict]:
    """Each side as the runner matches it: {betType, period, side, points,
    strikeId, price (American, None when nothing rests), depth (dollars)}."""
    priced = []
    for line in lines:
        offer = best_offer_of(line["ladder"])
        priced.append({"betType": line["betType"], "period": line["period"], "side": line["side"],
                       "points": line["points"], "strikeId": line["strike_id"],
                       "price": offer["price"] if offer else None, "depth": offer["depth"] if offer else 0.0})
    return priced


class ProphetXLiveBooks:
    """The wanted games' last ProphetX reads, and the thread that reads them."""

    def __init__(self, client: ProphetXClient, team_keys_of: TeamKeysReader,
                 clock: Callable[[], float] = time.time, read_gap_sec: float = READ_GAP_SEC):
        self._client = client
        self._team_keys_of = team_keys_of
        self._clock = clock
        self._read_gap_sec = read_gap_sec
        self._lock = threading.Lock()
        self._wanted: dict[str, tuple[dict, float]] = {}
        self._reads: dict[str, dict] = {}
        self._events: dict[str, tuple[dict, float]] = {}
        self._last_error: str | None = None
        self._reads_done = 0

    def want(self, games: list[dict]) -> dict:
        """Marks the games wanted and answers {reader: status, games:
        {gameKey: read | None (not read yet)}}. A read is {readAt (ISO),
        prophetxEventId, totalsOnly, lines: priced_lines} or {readAt, error}."""
        with self._lock:
            now = self._clock()
            for game in games:
                self._wanted[game["gameKey"]] = (game, now)
            return {"reader": self._status_locked(),
                    "games": {game["gameKey"]: self._reads.get(game["gameKey"]) for game in games}}

    def status(self) -> dict:
        with self._lock:
            return self._status_locked()

    def _status_locked(self) -> dict:
        blocked = self._client.blocked_until
        return {"wanted": len(self._wanted), "read": len(self._reads), "readsDone": self._reads_done,
                "lastError": self._last_error, "accessKey": self._client.access_key[:8],
                "blockedUntil": datetime.fromtimestamp(blocked, timezone.utc).isoformat() if blocked else None}

    def _next_game(self) -> dict | None:
        """The wanted game read longest ago (never read first), dropping games no longer asked for."""
        with self._lock:
            now = self._clock()
            for key in [key for key, (_game, asked_at) in self._wanted.items() if now - asked_at > WANTED_FOR_SEC]:
                self._wanted.pop(key, None)
                self._reads.pop(key, None)
                self._events.pop(key, None)
            if not self._wanted:
                return None
            def read_at(key: str) -> float:
                read = self._reads.get(key)
                return read["readAtSec"] if read else -1.0
            key = min(self._wanted, key=read_at)
            return self._wanted[key][0]

    def run(self, stop: threading.Event) -> None:
        """Reads wanted games until `stop` is set."""
        while not stop.is_set():
            game = self._next_game()
            if game is None or self._client.blocked_until:
                stop.wait(IDLE_SEC)
                continue
            read = self.read_game(game)
            with self._lock:
                if game["gameKey"] in self._wanted:
                    self._reads[game["gameKey"]] = read
                self._reads_done += 1
                self._last_error = read.get("error") or self._last_error
            stop.wait(self._read_gap_sec)

    def read_game(self, game: dict) -> dict:
        """One game's read: {readAt, readAtSec, prophetxEventId, totalsOnly, lines} or {readAt, readAtSec, error}."""
        key = game["gameKey"]
        try:
            event = self._event_of(game)
            if isinstance(event, str):
                return self._error(event)
            if event.get("status") not in (None, PROPHETX_EVENT_NOT_STARTED):
                return self._error(f"ProphetX game {event.get('event_id')} is {event.get('status')}")
            mismatch = teams_mismatch(event, game["line"], self._team_keys_of)
            totals_only = bool(mismatch and mismatch.endswith(REVERSED_TEAMS_WORDS))
            if mismatch and not totals_only:
                return self._error(mismatch)
            lines = lines_of_event(event, self._client.markets(event["event_id"]))
        except ProphetXBlocked as error:
            return self._error(f"ProphetX is refusing the live key for now: {error}")
        except TeamKeysUnread as error:
            return self._error(f"ProphetX's teams could not be checked against the board's ({error})")
        except (OSError, ProphetXError) as error:
            with self._lock:
                self._events.pop(key, None)
            return self._error(f"ProphetX could not be read ({error})")
        if isinstance(lines, str):
            return self._error(lines)
        at = self._clock()
        return {"readAt": _iso(at), "readAtSec": at, "prophetxEventId": event["event_id"],
                "totalsOnly": totals_only, "lines": priced_lines(lines)}

    def _event_of(self, game: dict) -> dict | str:
        key = game["gameKey"]
        with self._lock:
            held = self._events.get(key)
        if held and self._clock() - held[1] <= EVENT_REUSE_SEC:
            return held[0]
        event = event_of_line(self._client, self._team_keys_of,
                              {"strikeIds": game.get("strikeIds"), "eventStart": game.get("eventStart")}, game["line"])
        if isinstance(event, dict):
            with self._lock:
                self._events[key] = (event, self._clock())
        return event

    def _error(self, why: str) -> dict:
        at = self._clock()
        return {"readAt": _iso(at), "readAtSec": at, "error": why.rstrip(".")}


def _iso(unix_sec: float) -> str:
    return datetime.fromtimestamp(unix_sec, timezone.utc).isoformat()


def reader_if_configured(access_key: str, secret_key: str, base_url: str,
                         team_keys_of: TeamKeysReader) -> ProphetXLiveBooks | None:
    """The reader on the live key pair, or None (logged) when it is unset."""
    if not (access_key and secret_key):
        log.info("prophetx live: no PROPHETX_LIVE_ACCESS_KEY / PROPHETX_LIVE_SECRET_KEY, ProphetX lines stay at Unabated's price")
        return None
    log.info("prophetx live: reading ProphetX's live prices on access key %s", access_key[:8])
    return ProphetXLiveBooks(ProphetXClient(access_key, secret_key, base_url), team_keys_of)
