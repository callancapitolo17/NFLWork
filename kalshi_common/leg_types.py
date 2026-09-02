"""Kalshi MLB leg-typing helpers shared between taker and maker bots.

Converts raw Kalshi market-ticker dicts to typed fair_value.SpreadLeg /
fair_value.TotalLeg instances, and extracts canonical spread / total line
values from a legs list.
"""
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from zoneinfo import ZoneInfo

from kalshi_common import fair_value

# The 4-cell devig families stored in mlb_sgp_odds. Each family partitions the
# outcome space, so devig_book n-way-devigs the four cells together.
SPREAD_TOTAL_FAMILY = (
    "Home Spread + Over", "Home Spread + Under",
    "Away Spread + Over", "Away Spread + Under",
)
ML_TOTAL_FAMILY = (
    "Home ML + Over", "Home ML + Under",
    "Away ML + Over", "Away ML + Under",
)

# 3-letter Kalshi team code → mlb_parlay_lines.home_team / away_team canonical name.
# Kalshi uses 3-letter codes; mlb_parlay_lines stores Odds-API canonical names.
_MLB_CODE_TO_TEAM = {
    "ARI": "Arizona Diamondbacks", "ATL": "Atlanta Braves", "BAL": "Baltimore Orioles",
    "BOS": "Boston Red Sox", "CHC": "Chicago Cubs", "CWS": "Chicago White Sox",
    "CIN": "Cincinnati Reds", "CLE": "Cleveland Guardians", "COL": "Colorado Rockies",
    "DET": "Detroit Tigers", "HOU": "Houston Astros", "KC": "Kansas City Royals",
    "LAA": "Los Angeles Angels", "LAD": "Los Angeles Dodgers", "MIA": "Miami Marlins",
    "MIL": "Milwaukee Brewers", "MIN": "Minnesota Twins", "NYM": "New York Mets",
    "NYY": "New York Yankees", "OAK": "Athletics", "ATH": "Athletics",
    "AZ": "Arizona Diamondbacks", "PHI": "Philadelphia Phillies",
    "PIT": "Pittsburgh Pirates", "SD": "San Diego Padres", "SF": "San Francisco Giants",
    "SEA": "Seattle Mariners", "STL": "St. Louis Cardinals", "TB": "Tampa Bay Rays",
    "TEX": "Texas Rangers", "TOR": "Toronto Blue Jays",
    "WAS": "Washington Nationals", "WSH": "Washington Nationals",
}


_ET = ZoneInfo("America/New_York")
_SUFFIX_MONTHS = {"JAN": 1, "FEB": 2, "MAR": 3, "APR": 4, "MAY": 5, "JUN": 6,
                  "JUL": 7, "AUG": 8, "SEP": 9, "OCT": 10, "NOV": 11, "DEC": 12}


def parse_suffix_start_utc(suffix: str) -> datetime | None:
    """KXMLB* event-suffix YYMMMDDHHMM prefix (US/Eastern) -> naive-UTC, or None.

    The suffix is the ONLY reliable first-pitch source on a Kalshi market:
    ``close_time`` is first pitch + 72h (see the kalshi_close_time_not_start
    note). Naive UTC matches the ``GameRef.commence_time`` convention the
    per-book ``match_events`` helpers bucket on.

    Lives here rather than in a bot package because both ``kalshi_rfi`` and the
    maker's leg surface key games on it; ``kalshi_rfi.discovery`` re-exports it.
    """
    if len(suffix) < 11:
        return None
    try:
        year = 2000 + int(suffix[0:2])
        month = _SUFFIX_MONTHS[suffix[2:5].upper()]
        day = int(suffix[5:7])
        hour = int(suffix[7:9])
        minute = int(suffix[9:11])
        local = datetime(year, month, day, hour, minute, tzinfo=_ET)
    except (KeyError, ValueError):
        return None
    return local.astimezone(timezone.utc).replace(tzinfo=None)


# Doubleheader marker: Kalshi appends G1/G2 to BOTH games' event suffixes
# (live 2026-09-01: KXMLBGAME-26SEP041410DETCLEG1 and ...1915DETCLEG2). Team
# codes are letters only, so a trailing "G<digits>" can never be part of one.
_GAME_NUMBER_RE = re.compile(r"G(\d+)$")


def split_game_number(suffix: str) -> tuple[str, int | None]:
    """A KXMLB* event suffix -> (suffix without the G-marker, game number).

    Returns (suffix, None) for the ordinary single-game grammar. The game
    number is NOT part of the team block, but it IS part of the game's
    identity — ``legset.game_id_of`` keeps the whole suffix, so the two games
    of a doubleheader stay distinct keys everywhere downstream.
    """
    match = _GAME_NUMBER_RE.search(suffix)
    if not match:
        return suffix, None
    return suffix[:match.start()], int(match.group(1))


def game_number_from_suffix(suffix: str) -> int | None:
    """The doubleheader game number in an event suffix, or None.

    A non-None result means the team pair alone CANNOT identify the game —
    two Kalshi events share it. Resolvers that key on team names must
    disambiguate with ``unique_game_by_start`` (start time, fail closed) or
    drop the pair outright (``kalshi_rfi.discovery.drop_doubleheaders``).
    Answering such a lookup with ``LIMIT 1`` returns whichever row comes
    first — a wrong number, not a decline (#95 measured fairs off by
    0.05-0.11 from exactly this).
    """
    return split_game_number(suffix)[1]


# Kalshi's suffix minute and the Odds API's commence_time disagree by ~1
# minute on the same game (live 2026-09-01: suffix 26SEP012210STLLAD = 22:10
# ET vs Odds API 2026-09-02T02:11Z), so the match needs a tolerance. 30min is
# the leg surface's SURFACE_START_TOLERANCE_MIN — far wider than that skew and
# far narrower than any real doubleheader gap (2026-09-04 DET@CLE: 14:10 and
# 19:15, 5h apart).
SCHEDULE_START_TOLERANCE_MIN = 30.0


def as_naive_utc(dt: datetime | None) -> datetime | None:
    """Any datetime -> naive UTC, the convention every start time is compared in.

    ``parse_suffix_start_utc`` and ``mlb_target_lines.commence_time`` are
    already naive UTC; the Odds API's ``commence_time`` is tz-aware. Comparing
    the two raw raises TypeError, which a fail-safe ``except`` would swallow
    into a silent decline.
    """
    if dt is None:
        return None
    if dt.tzinfo is None:
        return dt
    return dt.astimezone(timezone.utc).replace(tzinfo=None)


def unique_game_by_start(kalshi_start: datetime | None,
                         candidates: list[tuple],
                         *,
                         tolerance_min: float = SCHEDULE_START_TOLERANCE_MIN,
                         ) -> tuple[object | None, str]:
    """Pick the ONE candidate game whose start matches a Kalshi event's.

    Inputs: the Kalshi event's first pitch (``parse_suffix_start_utc``) and
    ``candidates`` as ``(game_id, start_time)`` pairs already filtered to the
    event's team pair. Returns ``(game_id, "")`` on a unique match, else
    ``(None, reason)`` where reason is ``no_kalshi_start`` / ``unmatched`` /
    ``ambiguous``.

    AMBIGUITY FAILS CLOSED, the same rule as
    ``leg_surface.singles.match_book_game``. The team pair alone is not an
    identity: the Odds API ``events`` endpoint returns a multi-day window, so
    consecutive games of one series share it (10 of 25 events on 2026-09-01),
    and both games of a doubleheader share it on the SAME day. Answering
    either with an arbitrary row is the one failure a maker must never have.
    """
    kalshi_start = as_naive_utc(kalshi_start)
    if kalshi_start is None:
        return None, "no_kalshi_start"
    tolerance_sec = tolerance_min * 60.0
    matched: dict = {}
    for game_id, start in candidates:
        start = as_naive_utc(start)
        if start is None or game_id is None:
            continue
        if abs((start - kalshi_start).total_seconds()) > tolerance_sec:
            continue
        matched[game_id] = start
    if not matched:
        return None, "unmatched"
    if len(matched) > 1:
        return None, "ambiguous"
    return next(iter(matched)), ""


def _parse_event_suffix(suffix: str) -> tuple[str | None, str | None]:
    """Split a KXMLB* event suffix into (away_code, home_code).

    Format: YYMMMDDHHMM{AwayCode}{HomeCode}, optionally followed by a
    doubleheader marker G1/G2 (stripped here — see split_game_number). Date
    prefix is fixed at 11 chars. Each team code is 2 or 3 letters (KC/SF/SD/
    TB/AZ are 2-letter; the rest are 3-letter). Probes 3- then 2-letter home
    splits and returns the first where both codes are valid in
    _MLB_CODE_TO_TEAM. Returns (None, None) if no split matches — caller
    drops the event.
    """
    suffix, _game_number = split_game_number(suffix)
    if len(suffix) < 11 + 4:  # date prefix + at least 2+2 team chars
        return None, None
    team_block = suffix[11:]
    for home_len in (3, 2):
        if len(team_block) <= home_len:
            continue
        home = team_block[-home_len:]
        away = team_block[:-home_len]
        if home in _MLB_CODE_TO_TEAM and away in _MLB_CODE_TO_TEAM:
            return away, home
    return None, None


def _home_code_from_event_ticker(event_ticker: str) -> str | None:
    """Parse the home-team code from a Kalshi event ticker (2- or 3-letter)."""
    if "-" not in event_ticker:
        return None
    suffix = event_ticker.rsplit("-", 1)[-1]
    _, home = _parse_event_suffix(suffix)
    return home


def _typed_spread_from_leg(leg: dict):
    """SpreadLeg from any ``-{TEAM}{N}`` spread ticker (prefix-agnostic).

    Shared by the FG (KXMLBSPREAD) and F5 (KXMLBF5SPREAD) parse branches —
    the suffix grammar and sign semantics are identical (live-verified
    2026-08-11: KXMLBF5SPREAD-...-LAD3 is "LAD -2.5 first 5 innings").
    Callers gate on the series prefix; this helper never checks it.
    """
    mt = leg["market_ticker"]
    side = leg["side"]
    suffix = mt.rsplit("-", 1)[-1]
    n_chars = "".join(c for c in suffix if c.isdigit())
    team_chars = "".join(c for c in suffix if not c.isdigit())
    if not n_chars or not team_chars:
        return None
    n = int(n_chars)
    home_code = _home_code_from_event_ticker(leg.get("event_ticker", ""))
    team_is_home = (home_code is not None and team_chars == home_code)
    return fair_value.SpreadLeg(team_is_home=team_is_home, line_n=n, side=side)


def _typed_total_from_leg(leg: dict):
    """TotalLeg from any ``-{N}`` total ticker (prefix-agnostic).

    Shared by FG (KXMLBTOTAL) and F5 (KXMLBF5TOTAL): both use integer suffix
    N for line N - 0.5 (live-verified 2026-08-11: KXMLBF5TOTAL-...-7 has
    floor_strike 6.5). Callers gate on the series prefix.
    """
    side = leg["side"]
    try:
        n = int(leg["market_ticker"].rsplit("-", 1)[-1])
    except ValueError:
        return None
    return fair_value.TotalLeg(line_n=n, side=side)


def _leg_dict_to_typed(leg: dict, game_id: str):
    """Convert {market_ticker, event_ticker, side} to fair_value typed leg.

    Determines team_is_home by parsing the home code from the event_ticker
    (no DB lookup needed — the ticker self-encodes the home/away convention).

    Deliberately FG-only (issue #84): the taker calls this ungated on its
    candidate legs and has no period concept — teaching it F5 prefixes would
    let an F5 leg silently type as a full-game leg there. F5 parsing lives in
    legset.parse_leg, which carries the period.
    """
    mt = leg["market_ticker"]
    if mt.startswith("KXMLBSPREAD-"):
        return _typed_spread_from_leg(leg)
    if mt.startswith("KXMLBTOTAL-"):
        return _typed_total_from_leg(leg)
    return None


def _spread_line_from_legs(legs: list[dict]) -> float | None:
    """Signed home-perspective spread line of the first spread leg.

    Issue #70: sign follows the ticker's team — home margin -> -(n-0.5),
    away margin -> +(n-0.5). Returns 0.0 when no leg is a spread (legacy
    contract) and None when a spread leg's team cannot be resolved from its
    event ticker — a silently wrong sign is worse than a NULL in telemetry
    (this helper's only production caller is the taker's research events).
    """
    for l in legs:
        if l["market_ticker"].startswith("KXMLBSPREAD-"):
            try:
                typed = _leg_dict_to_typed(l, "")
            except (KeyError, TypeError, ValueError):
                return None
            if typed is None:
                return None
            home_code = _home_code_from_event_ticker(str(l.get("event_ticker", "")))
            if home_code is None:
                return None
            return (-(typed.line_n - 0.5) if typed.team_is_home
                    else (typed.line_n - 0.5))
    return 0.0


def _total_line_from_legs(legs: list[dict]) -> float:
    for l in legs:
        if l["market_ticker"].startswith("KXMLBTOTAL-"):
            try:
                n = int(l["market_ticker"].rsplit("-", 1)[-1])
                return n - 0.5
            except ValueError:
                continue
    return 0.0


def _moneyline_side(leg: dict) -> tuple[bool, str] | None:
    """(team_is_home, side) for a KXMLBGAME (moneyline) leg, or None.

    KXMLBGAME-{event_suffix}-{TEAMCODE}: YES = that team wins. team_is_home iff
    the ticker's team code equals the home code parsed from the event ticker.
    """
    mt = str(leg.get("market_ticker", ""))
    team = mt.rsplit("-", 1)[-1]
    home = _home_code_from_event_ticker(str(leg.get("event_ticker", "")))
    if not team or home is None:
        return None
    return (team == home, leg.get("side", "yes"))


def _event_codes_from_legs(legs: list[dict]) -> tuple[str | None, str | None]:
    """(away_code, home_code) parsed from any leg's event ticker. All legs of a
    game share the same event suffix, so the first resolvable one wins."""
    for l in legs:
        et = str(l.get("event_ticker", ""))
        if "-" in et:
            away, home = _parse_event_suffix(et.rsplit("-", 1)[-1])
            if away and home:
                return away, home
    return None, None


@dataclass(frozen=True)
class ComboDescriptor:
    """How to look a 2-leg Kalshi combo up in the mlb_sgp_odds grid.

    kind         : "spread_total" | "ml_total"
    spread_line  : home-perspective spread (None for ml_total — there's no spread)
    total_line   : the Over/Under line
    target_combo : the single grid cell matching the legs' STATED sides
                   (e.g. "Away ML + Under")
    combo_family : the 4 cells to n-way devig together
    away_code    : Kalshi away team code (for game resolution)
    home_code    : Kalshi home team code
    """
    kind: str
    spread_line: float | None
    total_line: float | None
    target_combo: str
    combo_family: tuple
    away_code: str
    home_code: str


def combo_descriptor(legs: list[dict]) -> "ComboDescriptor | None":
    """Classify a 2-leg combo into a grid-lookup descriptor, or None if it isn't
    a shape the SGP grid carries ({spread,total} or {moneyline,total}).

    The target cell is derived from the ACTUAL leg sides (not hardcoded), which
    is also the fix for the old maker's "Home Spread + Over"-always bug.
    """
    if not legs or len(legs) != 2:
        return None

    def _pfx(l):
        return str(l.get("market_ticker", "")).split("-")[0]

    by: dict[str, list[dict]] = {}
    for l in legs:
        by.setdefault(_pfx(l), []).append(l)

    away, home = _event_codes_from_legs(legs)
    if not away or not home:
        return None

    total_legs = by.get("KXMLBTOTAL")
    if not total_legs or len(total_legs) != 1:
        return None
    total = _leg_dict_to_typed(total_legs[0], "")        # TotalLeg
    if total is None:
        return None
    total_line = total.line_n - 0.5
    over_part = "Over" if total.side == "yes" else "Under"

    if len(by.get("KXMLBSPREAD", [])) == 1:
        spread = _leg_dict_to_typed(by["KXMLBSPREAD"][0], "")   # SpreadLeg
        if spread is None:
            return None
        # At a half-point line exactly one team covers, so spread 'no' == other.
        home_covers = ((spread.team_is_home and spread.side == "yes")
                       or (not spread.team_is_home and spread.side == "no"))
        part = "Home" if home_covers else "Away"
        # Issue #70: sign follows the ticker's team (home margin -> negative,
        # away margin -> positive), matching mlb_sgp_odds.spread_line.
        return ComboDescriptor(
            kind="spread_total",
            spread_line=(-(spread.line_n - 0.5) if spread.team_is_home
                         else (spread.line_n - 0.5)),
            total_line=total_line,
            target_combo=f"{part} Spread + {over_part}",
            combo_family=SPREAD_TOTAL_FAMILY,
            away_code=away, home_code=home)

    if len(by.get("KXMLBGAME", [])) == 1:
        ml = _moneyline_side(by["KXMLBGAME"][0])
        if ml is None:
            return None
        team_is_home, side = ml
        home_ml = ((team_is_home and side == "yes")
                   or (not team_is_home and side == "no"))
        part = "Home" if home_ml else "Away"
        return ComboDescriptor(
            kind="ml_total",
            spread_line=None,
            total_line=total_line,
            target_combo=f"{part} ML + {over_part}",
            combo_family=ML_TOTAL_FAMILY,
            away_code=away, home_code=home)

    return None
