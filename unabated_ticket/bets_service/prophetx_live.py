"""ProphetX's live prices off ProphetX's market websocket, for the Edges list's
live price overlay and the ProphetX Place desk (Cal, 2026-11-06, option B).

Why a stream: ProphetX allows about 1 REST request a second, and that limit
held across our two key pairs (2026-11-05 23:49Z, Place's get_markets got a
429 while this reader polled one game every 2 s on the live key from the same
VM). So the reader no longer polls: it reads each wanted game ONCE with
get_markets and then follows ProphetX's pushed `market_selections` frames for
it, so REST is used only to find a game, seed its book and register the socket.

The ask never waits: the server runner's live price pass names the board games
it lists ProphetX lines on (POST /prophetx_live.json); `want` marks them wanted
and answers each game's current book. A game no longer asked for in
WANTED_FOR_SEC is dropped (unsubscribed). Place reads `streamed_game`: the
game's event and markets while the stream holds a current book of it, else
None and Place reads ProphetX itself.

ProphetX facts this rests on (docs.prophetx.co "WebSocket Events" and
"Integrate with Websockets", and the read-only recon on the VM 2026-11-06
00:44Z, live key defff84a, 3 games, 1,104 frames in 120 s):
  * Pusher transport as the order stream (prophetx_order_stream.py). POST
    /v4/mm/websocket {socket_id, service, subscriptions: [{type: "event",
    ids: ["<event_id>", ...]}]} answers authorized_channel entries, one per
    game `private-broadcast-service=6-device_type=5-event=<id>` with its own
    `auth` and scope {event_id}; channel_limit 500 per socket, counting the
    union of every channel ever registered on it (a new socket starts clean).
    The registration names the complete wanted set each time.
  * `v2.market_selections` frames {change_type: market_selections, op: u,
    payload: base64 JSON {sport_event_id, market_id, info: {id, name,
    selections: [[side offers], [side offers]], strike, type, status,
    sequence_number}}}. One frame is the WHOLE top-of-book (up to 10 offers a
    side) of ONE strike pair, so it replaces that pair's `selections`.
  * No snapshot on subscribe: a frame comes only when a strike's book
    changes, so a quiet strike is known only from the seed read.
  * health_check every 5 s on the private user channel; silent 30 s = dead.

What it holds, in memory only: per wanted game its ProphetX event and the
markets list (get_markets' shape) with frames applied. A reconnect drops every
book (frames may have been missed) and each game is seeded again.

Side effects: one outbound websocket on the live key pair; REST on the same
pair (connection-config + register per connect or wanted-set change, event
lookups, one get_markets per game seed). Reads only: no order, no DB. The
runner's GET /team_keys.json is read for the team check.
"""
import base64
import copy
import functools
import json
import logging
import threading
import time
from collections.abc import Callable
from datetime import datetime, timezone

from unabated_ticket.bets_service.novig_book_stream import Connect, Connection, websocket_client_connect
from unabated_ticket.bets_service.prophetx_client import (
    CFTC_TERMINOLOGY_HEADER, WEBSOCKET_REGISTER_PATH, ProphetXBlocked, ProphetXClient, ProphetXError, data_or_flat,
)
from unabated_ticket.bets_service.prophetx_markets import best_offer_of, lines_of_event, sub_markets_of
from unabated_ticket.bets_service.prophetx_order_stream import (
    HANDSHAKE_TIMEOUT_SEC, HEALTH_EVENT, PRIVATE_USER_CHANNEL_MARK, RECONNECT_MAX_SEC, RECONNECT_MIN_SEC,
    RECV_POLL_SEC, SILENT_RECONNECT_SEC, STABLE_CONNECTION_SEC, parse_frame, pusher_url,
)
from unabated_ticket.bets_service.prophetx_place import (
    PROPHETX_EVENT_NOT_STARTED, REVERSED_TEAMS_WORDS, TeamKeysReader, TeamKeysUnread, event_of_line, teams_mismatch,
)

log = logging.getLogger(__name__)

# Seconds between the seeder's ProphetX calls: half the 1-a-second limit, the rest left to Place and the bet source.
READ_GAP_SEC = 2.0
# A game the runner stopped asking for this long ago is no longer followed (its pass runs every 30 s).
WANTED_FOR_SEC = 90.0
# How long a game's ProphetX event (status, teams) is reused before it is found again.
EVENT_REUSE_SEC = 10 * 60.0
# A seeded book is read again this often: a strike ProphetX adds later is in no frame we can place.
SEED_REFRESH_SEC = 10 * 60.0
# A frame for a strike the book does not hold asks for a fresh seed, but not sooner than this after the last one.
NEW_STRIKE_RESEED_SEC = 60.0
# How long the seeder idles when nothing is due or ProphetX is refusing the key.
IDLE_SEC = 1.0
# The most games one ask may name.
MAX_GAMES_PER_ASK = 300
# A new registration waits this long after the last (the wanted set changes in bursts every 30 s pass).
REGISTER_GAP_SEC = 5.0
# A socket whose channel count gets this close to its limit is replaced (a new socket starts at 0).
CHANNEL_HEADROOM = 50
# Frames held per game between its subscribe and its seed read landing.
MAX_PENDING_FRAMES = 2000
MARKET_SELECTIONS_CHANGE = "market_selections"
EVENT_CHANNEL_MARK = "-event="
ACTIVE_STATUS = "active"
STREAM_DOWN_WORDS = "ProphetX's price stream is down"


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


def frame_payload(data: object) -> dict | None:
    """A market_selections frame's decoded payload {sport_event_id, market_id, info}, or None."""
    if not isinstance(data, dict) or data.get("change_type") != MARKET_SELECTIONS_CHANGE:
        return None
    payload = data.get("payload")
    if isinstance(payload, str):
        payload = json.loads(base64.b64decode(payload))
    return payload if isinstance(payload, dict) and isinstance(payload.get("info"), dict) else None


def _side_strike_ids(selections: object) -> list[str | None]:
    """Each side's one strike id (None for a side naming none or two)."""
    ids = []
    for side in selections if isinstance(selections, list) else []:
        found = {offer.get("strike_id") for offer in side if isinstance(offer, dict) and offer.get("strike_id")} \
            if isinstance(side, list) else set()
        ids.append(found.pop() if len(found) == 1 else None)
    return ids


def _empty_side(old_side: list) -> list[dict]:
    """A side with nothing resting that still names its strike (so the pair keeps reading as a line)."""
    first = next((offer for offer in old_side if isinstance(offer, dict)), {})
    return [{**first, "price": None, "quantity": 0, "value": 0}]


def _sequence_of(info: dict) -> int:
    sequence = info.get("sequence_number")
    return sequence if isinstance(sequence, int) and not isinstance(sequence, bool) else -1


def _newest_offer_at(pair: dict) -> int:
    """The newest updated_at among a pair's offers (-1 when none names one)."""
    stamps = [offer.get("updated_at") for side in pair.get("selections") or [] if isinstance(side, list)
              for offer in side if isinstance(offer, dict)]
    return max((stamp for stamp in stamps if isinstance(stamp, int) and not isinstance(stamp, bool)), default=-1)


class GameBook:
    """One game's markets (get_markets' shape) with frames applied, keyed by strike id."""

    def __init__(self, markets: list[dict], seeded_at: float, epoch: int, frames_before: list[dict] = ()):
        """`frames_before`: frames that arrived between the subscribe and this
        seed; each is applied only when newer than every offer the seed shows
        on its pair (sequence_number against the offers' updated_at, both
        nanoseconds in the recon), so the seed is never overwritten by a frame
        it already includes."""
        self.markets = copy.deepcopy(markets)
        self.seeded_at = seeded_at
        self.epoch = epoch
        self.wants_reseed = False
        self._pair_of_strike: dict[str, dict] = {}
        self._sequence_of_pair: dict[int, int] = {}
        for market in self.markets:
            for sub in sub_markets_of(market):
                for strike_id in _side_strike_ids(sub.get("selections")):
                    if strike_id:
                        self._pair_of_strike[strike_id] = sub
        for info in frames_before:
            pair = self._pair_of(info)
            if pair is not None and _sequence_of(info) > _newest_offer_at(pair):
                self.apply(info)

    def _pair_of(self, info: dict) -> dict | None:
        frame_ids = _side_strike_ids([side for side in info.get("selections") or [] if isinstance(side, list)])
        return next((self._pair_of_strike[strike_id] for strike_id in frame_ids if strike_id in self._pair_of_strike), None)

    def apply(self, info: dict) -> bool:
        """Replaces one strike pair's selections with a frame's (a frame older
        than the last applied to that pair is skipped). False when the book
        holds none of the frame's strikes (a strike added since the seed)."""
        frame_sides = [side for side in info.get("selections") or [] if isinstance(side, list)]
        frame_ids = _side_strike_ids(frame_sides)
        pair = self._pair_of(info)
        if pair is None:
            return False
        sequence = _sequence_of(info)
        if sequence >= 0:
            if sequence <= self._sequence_of_pair.get(id(pair), -1):
                return True  # an older frame arriving late
            self._sequence_of_pair[id(pair)] = sequence
        frame_side_of = {strike_id: side for strike_id, side in zip(frame_ids, frame_sides) if strike_id}
        active = info.get("status") in (None, ACTIVE_STATUS)
        old_sides = pair.get("selections") if isinstance(pair.get("selections"), list) else []
        new_sides = []
        for strike_id, old_side in zip(_side_strike_ids(old_sides), old_sides):
            side = frame_side_of.get(strike_id) if active else None
            # A side the frame leaves out (or a placeholder with no strike) has nothing resting.
            new_sides.append(copy.deepcopy(side) if side else _empty_side(old_side if isinstance(old_side, list) else []))
        pair["selections"] = new_sides
        return True


class ProphetXLiveBooks:
    """The wanted games' streamed books, the socket thread and the seeder."""

    def __init__(self, client: ProphetXClient, team_keys_of: TeamKeysReader,
                 connect: Connect = functools.partial(websocket_client_connect, peer="ProphetX"),
                 clock: Callable[[], float] = time.time, read_gap_sec: float = READ_GAP_SEC):
        self._client = client
        self._team_keys_of = team_keys_of
        self._connect = connect
        self._clock = clock
        self._read_gap_sec = read_gap_sec
        self._lock = threading.Lock()
        self._wanted: dict[str, tuple[dict, float]] = {}
        # gameKey -> {event, foundAt, totalsOnly} or {error, foundAt}
        self._games: dict[str, dict] = {}
        self._books: dict[int, GameBook] = {}
        # event_id -> frames since its subscribe, until its seed lands
        self._pending: dict[int, list[dict]] = {}
        # Socket state; an epoch counts connects, so a book seeded on an older socket is never served.
        self._epoch = 0
        self._live = False
        self._connected_at: float | None = None
        self._last_frame_at: float | None = None
        self._subscribed: set[int] = set()
        self._registered_at = -1e18
        self._channel_count = 0
        self._channel_limit: int | None = None
        self._frames = 0
        self._unknown_strike_frames = 0
        self._seeds_done = 0
        self._last_error: str | None = None

    # ---- what the routes read ----

    def want(self, games: list[dict]) -> dict:
        """Marks the games wanted and answers {reader: status, games:
        {gameKey: read | None (not read yet)}}. A read is {readAt (ISO),
        prophetxEventId, totalsOnly, lines: priced_lines} or {readAt, error}."""
        with self._lock:
            now = self._clock()
            for game in games:
                self._wanted[game["gameKey"]] = (game, now)
            held = {game["gameKey"]: self._held_locked(game["gameKey"]) for game in games}
            status = self._status_locked()
        return {"reader": status, "games": {key: self._read_of(*item) if item else None for key, item in held.items()}}

    def streamed_game(self, game_key: str) -> tuple[dict, list[dict]] | None:
        """(event, markets) for Place while the stream holds a current book of
        the game, else None (Place then reads ProphetX itself)."""
        with self._lock:
            held = self._held_locked(game_key)
        if not held or isinstance(held[0], str):
            return None
        return held[0]["event"], held[1]

    def status(self) -> dict:
        with self._lock:
            return self._status_locked()

    def _stream_live_locked(self) -> bool:
        return (self._live and self._last_frame_at is not None
                and self._clock() - self._last_frame_at <= SILENT_RECONNECT_SEC)

    def _status_locked(self) -> dict:
        blocked = self._client.blocked_until
        now = self._clock()
        return {"wanted": len(self._wanted), "subscribed": len(self._subscribed),
                "seeded": sum(1 for book in self._books.values() if book.epoch == self._epoch),
                "streamLive": self._stream_live_locked(), "frames": self._frames,
                "unknownStrikeFrames": self._unknown_strike_frames, "seedsDone": self._seeds_done,
                "channelCount": self._channel_count, "channelLimit": self._channel_limit,
                "connectedForSec": None if self._connected_at is None else round(now - self._connected_at),
                "lastFrameAgeMs": None if self._last_frame_at is None else int((now - self._last_frame_at) * 1000),
                "lastError": self._last_error, "accessKey": self._client.access_key[:8],
                "blockedUntil": datetime.fromtimestamp(blocked, timezone.utc).isoformat() if blocked else None}

    def _held_locked(self, game_key: str) -> tuple[dict, list[dict]] | tuple[str, None] | None:
        """(game, markets copy) when servable, (why, None) for a game that reads as an error, None when not read yet."""
        game = self._games.get(game_key)
        if game is None:
            return None
        if "error" in game:
            return game["error"], None
        book = self._books.get(game["event"]["event_id"])
        if book is None or book.epoch != self._epoch:
            return None
        if not self._stream_live_locked():
            return STREAM_DOWN_WORDS, None
        return game, copy.deepcopy(book.markets)

    def _read_of(self, game_or_error: dict | str, markets: list[dict] | None) -> dict:
        at = self._clock()
        if isinstance(game_or_error, str):
            return {"readAt": _iso(at), "error": game_or_error.rstrip(".")}
        event = game_or_error["event"]
        start = _unix_of(event.get("scheduled"))
        if start is not None and start <= at:
            return {"readAt": _iso(at), "error": f"ProphetX game {event.get('event_id')} has started"}
        lines = lines_of_event(event, markets)
        if isinstance(lines, str):
            return {"readAt": _iso(at), "error": lines}
        # The stream is live and the book current, so it is as of now.
        return {"readAt": _iso(at), "prophetxEventId": event["event_id"], "totalsOnly": game_or_error["totalsOnly"],
                "lines": priced_lines(lines)}

    # ---- the seeder (REST, one call at a time) ----

    def run(self, stop: threading.Event) -> None:
        """Runs the socket thread and seeds wanted games until `stop` is set."""
        socket_thread = threading.Thread(target=self.run_socket, args=(stop,), name="prophetx-live-socket", daemon=True)
        socket_thread.start()
        while not stop.is_set():
            if self._client.blocked_until or not self.seed_one():
                stop.wait(IDLE_SEC)
                continue
            stop.wait(self._read_gap_sec)
        socket_thread.join(timeout=5)

    def seed_one(self) -> bool:
        """Does the next REST job (find a game's event, or seed a subscribed
        game's book); False when nothing is due."""
        job = self._next_job()
        if job is None:
            return False
        kind, game = job
        if kind == "event":
            self._find_event(game)
        else:
            self._seed(game)
        return True

    def _next_job(self) -> tuple[str, dict] | None:
        with self._lock:
            now = self._clock()
            self._drop_unwanted_locked(now)
            for game, _asked_at in self._wanted.values():
                held = self._games.get(game["gameKey"])
                if held is None or now - held["foundAt"] > EVENT_REUSE_SEC:
                    return "event", game
            if not self._stream_live_locked():
                return None
            due = []
            for game, _asked_at in self._wanted.values():
                held = self._games[game["gameKey"]]
                if "error" in held or held["event"]["event_id"] not in self._subscribed:
                    continue
                book = self._books.get(held["event"]["event_id"])
                if book is None or book.epoch != self._epoch:
                    due.append((-1.0, game))
                elif now - book.seeded_at > SEED_REFRESH_SEC or (
                        book.wants_reseed and now - book.seeded_at > NEW_STRIKE_RESEED_SEC):
                    due.append((book.seeded_at, game))
            return ("seed", min(due, key=lambda pair: pair[0])[1]) if due else None

    def _drop_unwanted_locked(self, now: float) -> None:
        for key in [key for key, (_game, asked_at) in self._wanted.items() if now - asked_at > WANTED_FOR_SEC]:
            self._wanted.pop(key, None)
            held = self._games.pop(key, None)
            if held and "event" in held:
                self._books.pop(held["event"]["event_id"], None)
                self._pending.pop(held["event"]["event_id"], None)

    def _find_event(self, game: dict) -> None:
        key = game["gameKey"]
        found_at = self._clock()
        try:
            event = event_of_line(self._client, self._team_keys_of,
                                  {"strikeIds": game.get("strikeIds"), "eventStart": game.get("eventStart")}, game["line"])
            held = self._checked(event, game) if isinstance(event, dict) else {"error": event}
        except ProphetXBlocked as error:
            held = {"error": f"ProphetX is refusing the live key for now: {error}"}
        except TeamKeysUnread as error:
            held = {"error": f"ProphetX's teams could not be checked against the board's ({error})"}
        except (OSError, ProphetXError) as error:
            held = {"error": f"ProphetX could not be read ({error})"}
        held["foundAt"] = found_at
        with self._lock:
            if key not in self._wanted:
                return
            old = self._games.get(key)
            if old and "event" in old and ("event" not in held or held["event"]["event_id"] != old["event"]["event_id"]):
                self._books.pop(old["event"]["event_id"], None)
            self._games[key] = held
            if "error" in held:
                self._last_error = held["error"]

    def _checked(self, event: dict, game: dict) -> dict:
        """{event, totalsOnly} for a game that reads as the board's, else {error}."""
        if _int_or_none(event.get("event_id")) is None:
            return {"error": f"ProphetX game has no numeric event_id: {str(event.get('event_id'))[:40]}"}
        event = {**event, "event_id": int(event["event_id"])}
        if event.get("status") not in (None, PROPHETX_EVENT_NOT_STARTED):
            return {"error": f"ProphetX game {event.get('event_id')} is {event.get('status')}"}
        mismatch = teams_mismatch(event, game["line"], self._team_keys_of)
        totals_only = bool(mismatch and mismatch.endswith(REVERSED_TEAMS_WORDS))
        if mismatch and not totals_only:
            return {"error": mismatch}
        return {"event": event, "totalsOnly": totals_only}

    def _seed(self, game: dict) -> None:
        with self._lock:
            held = self._games.get(game["gameKey"])
            epoch = self._epoch
        if not held or "event" not in held:
            return
        event_id = held["event"]["event_id"]
        try:
            markets = self._client.markets(event_id)
        except (OSError, ProphetXError) as error:
            with self._lock:
                self._last_error = f"ProphetX game {event_id} could not be seeded ({error})"
                # The event is found again first: a failed read may mean the game moved or ended.
                self._games.pop(game["gameKey"], None)
            return
        with self._lock:
            self._seeds_done += 1
            if epoch == self._epoch and event_id in self._subscribed:
                self._books[event_id] = GameBook(markets, self._clock(), epoch, self._pending.pop(event_id, []))

    # ---- the socket thread ----

    def run_socket(self, stop: threading.Event) -> None:
        """Connect, sign in, follow the wanted games' channels; reconnect with
        backoff whenever the socket drops or falls silent."""
        backoff = RECONNECT_MIN_SEC
        while not stop.is_set():
            opened_at = None
            connection = None
            try:
                connection, socket_id, service = self._open()
                opened_at = self._clock()
                self._pump(stop, connection, socket_id, service)
            except Exception as error:  # a dead thread would leave every ProphetX line unread until a restart
                self._drop(connection, f"{type(error).__name__}: {error}")
            if opened_at is not None and self._clock() - opened_at > STABLE_CONNECTION_SEC:
                backoff = RECONNECT_MIN_SEC
            if stop.wait(backoff):
                break
            backoff = min(backoff * 2, RECONNECT_MAX_SEC)
        self._drop(None, "stopped")

    def _open(self) -> tuple[Connection, str, str | None]:
        config = self._client.websocket_config()
        connection = self._connect(pusher_url(config), {})
        with self._lock:
            self._epoch += 1
            self._books.clear()
            self._pending.clear()
            self._subscribed.clear()
            self._channel_count = 0
            self._registered_at = -1e18
            self._connected_at = self._last_frame_at = self._clock()
        established = self._await_event(connection, "pusher:connection_established")
        socket_id = established.get("socket_id") if isinstance(established, dict) else None
        if not socket_id:
            raise ValueError(f"connection_established without socket_id: {str(established)[:200]}")
        service = config.get("service")
        registered = self._register(str(socket_id), service, [])
        authenticated = registered.get("authenticated")
        if not isinstance(authenticated, dict) or not authenticated.get("auth"):
            raise ValueError(f"websocket register: no authenticated block in {str(registered)[:300]}")
        connection.send(json.dumps({"event": "pusher:signin", "data": {
            "auth": authenticated["auth"], "user_data": authenticated.get("user_data")}}))
        self._await_event(connection, "pusher:signin_success")
        # The private channel carries health_check every 5 s, the sign of life on a quiet night.
        private = next((channel for channel in registered.get("authorized_channel") or []
                        if isinstance(channel, dict) and PRIVATE_USER_CHANNEL_MARK in str(channel.get("channel_name"))), None)
        if private is None:
            raise ValueError("websocket register named no private user channel")
        connection.send(json.dumps({"event": "pusher:subscribe", "data": {
            "channel": private["channel_name"], "auth": private.get("auth")}}))
        self._await_event(connection, "pusher_internal:subscription_succeeded")
        with self._lock:
            self._live = True
        log.info("prophetx live: stream live (key %s)", self._client.access_key[:8])
        return connection, str(socket_id), service

    def _register(self, socket_id: str, service: str | None, event_ids: list[int]) -> dict:
        """POST /v4/mm/websocket naming the complete wanted set."""
        body: dict = {"socket_id": socket_id, "subscriptions": [{"type": "event", "ids": [str(i) for i in event_ids]}]
                      if event_ids else []}
        if service:
            body["service"] = service
        status, reply = self._client.request("POST", WEBSOCKET_REGISTER_PATH, body=body, extra_headers=CFTC_TERMINOLOGY_HEADER)
        if status != 200:
            raise ProphetXError(f"websocket register: HTTP {status} {str(reply)[:300]}")
        data = data_or_flat(reply)
        if not isinstance(data, dict):
            raise ProphetXError(f"websocket register: expected an object, got {str(reply)[:300]}")
        with self._lock:
            self._registered_at = self._clock()
            if isinstance(data.get("channel_count"), int):
                self._channel_count = data["channel_count"]
            if isinstance(data.get("channel_limit"), int):
                self._channel_limit = data["channel_limit"]
        if data.get("rejected"):
            log.warning("prophetx live: register rejected %s", str(data["rejected"])[:300])
        return data

    def _await_event(self, connection: Connection, wanted: str) -> object:
        """Read frames until `wanted` (answering pings), within HANDSHAKE_TIMEOUT_SEC."""
        deadline = self._clock() + HANDSHAKE_TIMEOUT_SEC
        while self._clock() < deadline:
            text = connection.recv(RECV_POLL_SEC)
            if not text:
                continue
            event, channel, data = parse_frame(text)
            self._heard()
            if event == wanted:
                return data
            if event == "pusher:error" or event.endswith("subscription_error"):
                raise ConnectionError(f"ProphetX's socket answered {event}: {str(data)[:300]}")
            self._handle(connection, event, channel, data)
        raise TimeoutError(f"no {wanted} within {HANDSHAKE_TIMEOUT_SEC} s")

    def _pump(self, stop: threading.Event, connection: Connection, socket_id: str, service: str | None) -> None:
        while not stop.is_set():
            self._follow_wanted(connection, socket_id, service)
            text = connection.recv(RECV_POLL_SEC)
            if text:
                event, channel, data = parse_frame(text)
                self._heard()
                if event == "pusher:error" or event.endswith("subscription_error"):
                    raise ConnectionError(f"ProphetX's socket answered {event}: {str(data)[:300]}")
                self._handle(connection, event, channel, data)
            with self._lock:
                silent_sec = self._clock() - self._last_frame_at
            if silent_sec > SILENT_RECONNECT_SEC:
                raise TimeoutError(f"ProphetX's socket sent nothing for {round(silent_sec)} s")

    def _follow_wanted(self, connection: Connection, socket_id: str, service: str | None) -> None:
        """Subscribes newly wanted games (one registration naming the whole set)
        and unsubscribes games no longer wanted."""
        with self._lock:
            wanted_ids = {held["event"]["event_id"] for held in self._games.values() if "event" in held}
            to_add = wanted_ids - self._subscribed
            to_drop = self._subscribed - wanted_ids
            register_due = bool(to_add) and self._clock() - self._registered_at >= REGISTER_GAP_SEC
            for event_id in to_drop:
                self._subscribed.discard(event_id)
                self._books.pop(event_id, None)
                self._pending.pop(event_id, None)
        for event_id in to_drop:
            # Fire-and-forget per the docs; the channel still counts toward the limit until a new socket.
            connection.send(json.dumps({"event": "pusher:unsubscribe", "data": {"channel": _event_channel(event_id)}}))
        if not register_due:
            return
        if self._channel_limit is not None and self._channel_count + len(to_add) > self._channel_limit - CHANNEL_HEADROOM:
            raise ConnectionError(f"{self._channel_count} of {self._channel_limit} channels used: a new socket starts clean")
        registered = self._register(socket_id, service, sorted(wanted_ids))
        for channel in registered.get("authorized_channel") or []:
            if not isinstance(channel, dict):
                continue
            scope = channel.get("scope") if isinstance(channel.get("scope"), dict) else {}
            event_id = _int_or_none(scope.get("event_id"))
            if event_id not in to_add or EVENT_CHANNEL_MARK not in str(channel.get("channel_name")):
                continue
            connection.send(json.dumps({"event": "pusher:subscribe", "data": {
                "channel": channel["channel_name"], "auth": channel.get("auth")}}))
            with self._lock:
                # A book is seeded only after the subscribe is sent, so no frame after the seed is missed.
                self._subscribed.add(event_id)

    def _heard(self) -> None:
        with self._lock:
            self._last_frame_at = self._clock()

    def _handle(self, connection: Connection, event: str, channel: str | None, data: object) -> None:
        if event == "pusher:ping":
            connection.send(json.dumps({"event": "pusher:pong", "data": {}}))
            return
        if event == HEALTH_EVENT or not event.endswith(MARKET_SELECTIONS_CHANGE):
            return
        try:
            payload = frame_payload(data)
        except (ValueError, TypeError) as error:
            log.warning("prophetx live: unreadable market_selections frame (%s): %s", error, str(data)[:300])
            return
        if payload is None:
            return
        with self._lock:
            self._frames += 1
            event_id = _int_or_none(payload.get("sport_event_id"))
            book = self._books.get(event_id)
            if book is None or book.epoch != self._epoch:
                if event_id in self._subscribed:
                    # Held for the seed: it may land after this frame but show the book before it.
                    pending = self._pending.setdefault(event_id, [])
                    if len(pending) < MAX_PENDING_FRAMES:
                        pending.append(payload["info"])
                return
            if not book.apply(payload["info"]):
                self._unknown_strike_frames += 1
                book.wants_reseed = True

    def _drop(self, connection: Connection | None, reason: str) -> None:
        with self._lock:
            was_live = self._live
            repeated = reason == self._last_error
            self._live = False
            self._connected_at = None
            self._books.clear()
            self._pending.clear()
            self._subscribed.clear()
            if reason != "stopped":
                self._last_error = reason
        if connection is not None:
            connection.close()
        if reason == "stopped":
            return
        if was_live:
            log.warning("prophetx live: stream down: %s", reason)
        elif not repeated:
            log.error("prophetx live: cannot connect (logged once until it changes): %s", reason)


def _int_or_none(value: object) -> int | None:
    if isinstance(value, bool):
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _event_channel(event_id: int) -> str:
    return f"private-broadcast-service=6-device_type=5-event={event_id}"


def _iso(unix_sec: float) -> str:
    return datetime.fromtimestamp(unix_sec, timezone.utc).isoformat()


def _unix_of(value: object) -> float | None:
    try:
        return datetime.fromisoformat(str(value).replace("Z", "+00:00")).timestamp()
    except ValueError:
        return None


def reader_if_configured(access_key: str, secret_key: str, base_url: str,
                         team_keys_of: TeamKeysReader) -> ProphetXLiveBooks | None:
    """The reader on the live key pair, or None (logged) when it is unset."""
    if not (access_key and secret_key):
        log.info("prophetx live: no PROPHETX_LIVE_ACCESS_KEY / PROPHETX_LIVE_SECRET_KEY, ProphetX lines stay at Unabated's price")
        return None
    log.info("prophetx live: streaming ProphetX's live prices on access key %s", access_key[:8])
    return ProphetXLiveBooks(ProphetXClient(access_key, secret_key, base_url), team_keys_of)
