#!/usr/bin/env python3
"""Rebuild novig_event_markets_query.json from the live novig.com app bundle.

Why: since 2026-09-22 Novig's GraphQL endpoint (api.novig.us/v1/graphql)
runs a Hasura allowlist. It executes only operations the novig.com app
ships and answers anything else with HTTP 200
{"errors": [{"message": "query is not allowed", ...}]}. Any Novig app
release that edits EventMarkets_Query, or a fragment it uses, breaks the
SGP structure fetch until this file is refreshed. The symptom in bot.log
is a BookTransportError at stage=structure reading
"GraphQL error: query is not allowed".

How: Hasura compares the parsed query with __typename stripped, so
whitespace and the __typename fields Apollo adds do not matter. The
selections and the ORDER of the definitions do: on 2026-10-02 the app's
order passed and reversed or sorted fragments were rejected. The app
bundle (graphql-codegen client preset) carries both pieces:
  * the source text of every operation and fragment, as the keys of its
    documents map:  "\\n  query EventMarkets_Query(...) {...}": _.EventMarkets_QueryDocument
  * the full document as a pre-parsed AST literal whose `definitions` list
    the operation, then each fragment it uses, in the order the app sends.
The rebuilt query is those sources joined in that order. The bundle is
only searched with regexes, never executed.

Usage:
    python3 mlb_sgp/refresh_novig_query.py

Side effects: reads novig.com and api.novig.us, POSTs the rebuilt query
for up to MAX_VERIFY_EVENTS posted games, then atomically REPLACES
mlb_sgp/novig_event_markets_query.json (keeping its `variables`) and
prints a diff of the query text. Running bots re-read that file on every
Novig structure fetch, so no restart is needed; commit the file so other
checkouts get it. Writes nothing if the query is unchanged or any check
fails. No DB access.
"""
from __future__ import annotations

import difflib
import json
import os
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

from curl_cffi import requests as cffi_requests

# Works as a script and as mlb_sgp.refresh_novig_query: the repo root makes
# the mlb_sgp package importable and mlb_sgp/ the top-level
# scraper_novig_sgp module the bots import (the verify_books.py pattern).
_REPO_ROOT = Path(__file__).resolve().parent.parent
for _path in (str(_REPO_ROOT), str(_REPO_ROOT / "mlb_sgp")):
    if _path not in sys.path:
        sys.path.insert(0, _path)

import scraper_novig_sgp                                          # noqa: E402
from mlb_sgp._shared import check_response, json_or_raise         # noqa: E402
from mlb_sgp.novig import build_line_structure                    # noqa: E402
from mlb_sgp.novig_client import (BOOK, NOVIG_LEAGUE_PAGE, Event,  # noqa: E402
                                  _parse_events_response)

OPERATION_NAME = "EventMarkets_Query"
APP_URL = "https://novig.com/"
# The allowlist check does not depend on the league, so any posted game can
# verify the query — these are tried in order until one has a game.
VERIFY_LEAGUES = ("MLB", "NFL", "NBA", "NHL")
VERIFY_WINDOW_HOURS = 14 * 24
# The soonest game may be at first pitch with its markets locked, so a game
# whose tree does not parse is skipped for the next one, up to this many.
MAX_VERIFY_EVENTS = 3

# Expo web build: <script src="/_expo/static/js/web/index-<hash>.js">
BUNDLE_SRC_RE = re.compile(r'src="(/_expo/static/js/web/index-[0-9a-f]+\.js)"')
# One documents-map entry: "<escaped GraphQL source>":<module>.<Name>Document
DOCUMENT_MAP_ENTRY_RE = re.compile(
    r'"((?:[^"\\]|\\.)*)":[A-Za-z_$][\w$]*\.\w+(?:Document|FragmentDoc)\b')
DEFINITION_NAME_RE = re.compile(
    r"^\s*(?:query|mutation|subscription|fragment)\s+(\w+)")
AST_DEFINITION_RE = re.compile(
    r'\{kind:"(?:OperationDefinition|FragmentDefinition)",'
    r'(?:operation:"\w+",)?name:\{kind:"Name",value:"(\w+)"\}')
DECLARED_VARIABLE_RE = re.compile(r"\$(\w+)\s*:")


def fetch_bundle(session: cffi_requests.Session) -> str:
    """The novig.com app's main JS bundle (~17 MB)."""
    page = session.get(APP_URL, timeout=30)
    check_response(BOOK, "bundle", page)
    match = BUNDLE_SRC_RE.search(page.text)
    if match is None:
        raise RuntimeError(f"no index-<hash>.js bundle referenced by {APP_URL} "
                           "— the app's build layout changed")
    bundle = session.get("https://novig.com" + match.group(1), timeout=60)
    check_response(BOOK, "bundle", bundle)
    return bundle.text


def _closing_bracket(text: str, open_index: int) -> int:
    """Index of the bracket closing the one at open_index (skips strings)."""
    depth = 0
    in_string = False
    i = open_index
    while i < len(text):
        ch = text[i]
        if in_string:
            if ch == "\\":
                i += 2
                continue
            if ch == '"':
                in_string = False
        elif ch == '"':
            in_string = True
        elif ch in "{[":
            depth += 1
        elif ch in "}]":
            depth -= 1
            if depth == 0:
                return i
        i += 1
    raise RuntimeError("unterminated AST literal in bundle")


def definition_order(bundle: str) -> list[str]:
    """Operation + fragment names, in the order the app's document lists them."""
    head = ('{kind:"Document",definitions:[{kind:"OperationDefinition",'
            f'operation:"query",name:{{kind:"Name",value:"{OPERATION_NAME}"}}')
    starts = [m.start() for m in re.finditer(re.escape(head), bundle)]
    if len(starts) != 1:
        raise RuntimeError(f"expected 1 {OPERATION_NAME} AST literal in the "
                           f"bundle, found {len(starts)}")
    literal = bundle[starts[0]:_closing_bracket(bundle, starts[0]) + 1]
    return AST_DEFINITION_RE.findall(literal)


def sources_by_name(bundle: str) -> dict[str, str]:
    """GraphQL source text per operation/fragment name, from the documents map.

    An entry that is not valid JSON once quoted is skipped: the minifier
    writes Latin-1 characters as JS-only \\xNN escapes, and one such
    character in a document we don't need must not block the refresh. If
    a needed document is the one skipped, build_query reports it missing.
    """
    sources: dict[str, str] = {}
    for escaped in DOCUMENT_MAP_ENTRY_RE.findall(bundle):
        try:
            text = json.loads(f'"{escaped}"').strip()
        except json.JSONDecodeError:
            continue
        name_match = DEFINITION_NAME_RE.match(text)
        if name_match is None:
            continue
        name = name_match.group(1)
        if name in sources and sources[name] != text:
            raise RuntimeError(f"two different sources for {name} in the bundle")
        sources[name] = text
    return sources


def build_query(bundle: str) -> str:
    order = definition_order(bundle)
    sources = sources_by_name(bundle)
    missing = [name for name in order if name not in sources]
    if not order or order[0] != OPERATION_NAME or missing:
        raise RuntimeError(
            f"cannot assemble {OPERATION_NAME}: order={order} missing "
            f"sources={missing} (a source using a JS-only escape such as "
            "\\xNN is skipped as undecodable)")
    return "\n\n".join(sources[name] for name in order)


def declared_variables(query: str) -> set[str]:
    """Variables the operation declares, e.g. {"eventId", "marketVisibleWhere"}."""
    operation_header = query.split("{", 1)[0]
    return set(DECLARED_VARIABLE_RE.findall(operation_header))


def has_two_sided_rungs(structure: dict) -> bool:
    """True when a parsed FG tree has a total rung AND a spread rung priced
    on both sides — an outcome id plus `available`, the fields
    build_line_structure and the parlay pricer actually read."""
    def priced(leg: dict | None) -> bool:
        return bool(leg and leg.get("id") and leg.get("available") is not None)

    total_rung = any(priced(structure["over"].get(line))
                     and priced(structure["under"].get(line))
                     for line in structure["over"])
    spread_rung = any(priced(structure["home_spread"].get(line))
                      and priced(structure["away_spread"].get(-line))
                      for line in structure["home_spread"])
    return total_rung and spread_rung


def find_verify_events(session: cffi_requests.Session,
                       now: datetime) -> list[Event]:
    """Posted pregame games from the first league in VERIFY_LEAGUES with any."""
    for league in VERIFY_LEAGUES:
        resp = session.get(NOVIG_LEAGUE_PAGE.format(league=league), timeout=20)
        check_response(BOOK, "events", resp)
        events = _parse_events_response(json_or_raise(BOOK, "events", resp),
                                        now=now,
                                        window_hours=VERIFY_WINDOW_HOURS)
        if events:
            return events
    raise RuntimeError(f"no posted pregame game in any of {VERIFY_LEAGUES} "
                       "to verify against — rerun once one is posted")


def verify_live(session: cffi_requests.Session, payload: dict,
                events: list[Event]) -> str:
    """Raise unless the rebuilt query yields a priceable market tree.

    Goes through scraper_novig_sgp._gql, the bots' own transport, so an
    allowlist rejection or HTTP error raises BookTransportError on the first
    game (the allowlist is not per-event). A game whose tree has no
    two-sided total and spread rung is skipped for the next one.
    """
    for event in events[:MAX_VERIFY_EVENTS]:
        variables = dict(payload["variables"])
        variables["eventId"] = event.event_id
        body = dict(payload)
        body["variables"] = variables
        data = scraper_novig_sgp._gql(session, json.dumps(body))
        event_rows = (data.get("data") or {}).get("event") or [{}]
        markets = event_rows[0].get("markets") or []
        structure = build_line_structure(markets, event.home_sym,
                                         event.away_sym)
        if has_two_sided_rungs(structure):
            return (f"{event.away_team} @ {event.home_team} "
                    f"({event.event_id}): {len(markets)} markets, two-sided "
                    "total and spread rungs parse")
    raise RuntimeError(
        f"Novig accepted the rebuilt query, but none of the first "
        f"{MAX_VERIFY_EVENTS} games returned a two-sided total and spread "
        "rung — the app may have dropped fields the parser reads (outcome "
        "id, available, competitor.symbol), or every game is locked; "
        "nothing written")


def write_payload(path: Path, payload: dict) -> None:
    """Replace ``path`` atomically: bots re-read it on every structure
    fetch, so they must never see a half-written file."""
    staging = path.with_name(path.name + ".tmp")
    staging.write_text(json.dumps(payload))
    os.replace(staging, path)


def main() -> int:
    path = scraper_novig_sgp.EVENT_MARKETS_PATH
    if not path.exists():
        raise RuntimeError(
            f"{path} is missing — restore it with `git checkout -- "
            f"mlb_sgp/{path.name}` first; its `variables` cannot be rebuilt "
            "from the bundle")
    current = json.loads(path.read_text())

    session = scraper_novig_sgp.init_session()
    query = build_query(fetch_bundle(session))
    declared = declared_variables(query)
    if declared != set(current["variables"]):
        raise RuntimeError(
            f"the app's {OPERATION_NAME} declares variables {sorted(declared)} "
            f"but {path.name} sends {sorted(current['variables'])} — edit "
            "`variables` to match, then rerun; nothing written")
    payload = {"operationName": OPERATION_NAME,
               "variables": current["variables"], "query": query}

    events = find_verify_events(session, datetime.now(timezone.utc))
    print(f"live check: {verify_live(session, payload, events)}")

    if query == current["query"]:
        print(f"{path.name} already matches the app — nothing written")
        return 0
    write_payload(path, payload)
    diff = difflib.unified_diff(current["query"].splitlines(), query.splitlines(),
                                "committed", "app", lineterm="", n=1)
    print("\n".join(diff))
    print(f"wrote {path} — running bots pick it up on their next Novig "
          "structure fetch; commit it")
    return 0


if __name__ == "__main__":
    sys.exit(main())
