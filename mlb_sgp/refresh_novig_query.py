#!/usr/bin/env python3
"""Rebuild novig_event_markets_query.json from the live novig.com app bundle.

Why: since 2026-09-22 Novig's GraphQL endpoint (api.novig.us/v1/graphql)
runs a Hasura allowlist. It executes only operations the novig.com app
ships and answers anything else with HTTP 200
{"errors": [{"message": "query is not allowed"}]}. Any Novig app release
that edits EventMarkets_Query, or a fragment it uses, breaks the SGP
structure fetch until this file is refreshed. The symptom is a
BookTransportError at stage=structure reading
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
The rebuilt query is those sources joined in that order.

Usage:
    python mlb_sgp/refresh_novig_query.py

Side effects: reads novig.com and api.novig.us, POSTs the rebuilt query
once against a live event, then REPLACES mlb_sgp/novig_event_markets_query.json
(keeping its `variables`) and prints a diff of the query text. Writes
nothing if the query is unchanged or the live check fails. No DB access.
"""
from __future__ import annotations

import difflib
import json
import re
import sys
from pathlib import Path

from curl_cffi import requests as cffi_requests

QUERY_PATH = Path(__file__).resolve().parent / "novig_event_markets_query.json"
OPERATION_NAME = "EventMarkets_Query"
APP_URL = "https://novig.com/"
NOVIG_GRAPHQL = "https://api.novig.us/v1/graphql"
TRADING_PAGE = "https://api.novig.us/nbx/v1/trading/{league}/page"
# The allowlist check does not depend on the league, so any posted game can
# verify the query — these are tried in order until one has a game.
VERIFY_LEAGUES = ("MLB", "NFL", "NBA", "NHL")

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


def fetch_bundle(session) -> str:
    """The novig.com app's main JS bundle (~17 MB)."""
    html = session.get(APP_URL, timeout=30).text
    match = BUNDLE_SRC_RE.search(html)
    if match is None:
        raise RuntimeError(f"no index-<hash>.js bundle referenced by {APP_URL} "
                           "— the app's build layout changed")
    return session.get("https://novig.com" + match.group(1), timeout=60).text


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
    """GraphQL source text per operation/fragment name, from the documents map."""
    sources: dict[str, str] = {}
    for escaped in DOCUMENT_MAP_ENTRY_RE.findall(bundle):
        text = json.loads(f'"{escaped}"').strip()
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
        raise RuntimeError(f"cannot assemble {OPERATION_NAME}: order={order} "
                           f"missing sources={missing}")
    return "\n\n".join(sources[name] for name in order)


def find_live_event_id(session) -> str:
    for league in VERIFY_LEAGUES:
        page = session.get(TRADING_PAGE.format(league=league), timeout=20).json()
        for section in page.get("sections") or []:
            for card in (section.get("content") or {}).get("components") or []:
                if card.get("type") == "game_event_card" and card.get("eventId"):
                    return card["eventId"]
    raise RuntimeError(f"no posted game in any of {VERIFY_LEAGUES} to verify against")


def count_live_markets(session, payload: dict, event_id: str) -> int:
    """POST the payload for one event; raise unless Novig returns a market tree."""
    body = dict(payload, variables=dict(payload["variables"], eventId=event_id))
    data = session.post(NOVIG_GRAPHQL, data=json.dumps(body),
                        headers={"Content-Type": "application/json"},
                        timeout=20).json()
    if data.get("errors"):
        raise RuntimeError(f"Novig rejected the rebuilt query: {data['errors'][0]}")
    events = (data.get("data") or {}).get("event") or []
    if not events or not events[0].get("markets"):
        raise RuntimeError(f"rebuilt query returned no markets for event {event_id}")
    return len(events[0]["markets"])


def main() -> int:
    session = cffi_requests.Session(impersonate="chrome")
    current = json.loads(QUERY_PATH.read_text())
    query = build_query(fetch_bundle(session))
    payload = {"operationName": OPERATION_NAME,
               "variables": current["variables"], "query": query}

    event_id = find_live_event_id(session)
    n_markets = count_live_markets(session, payload, event_id)
    print(f"live check: event {event_id} returned {n_markets} markets")

    if query == current["query"]:
        print(f"{QUERY_PATH.name} already matches the app — nothing written")
        return 0
    QUERY_PATH.write_text(json.dumps(payload))
    diff = difflib.unified_diff(current["query"].splitlines(), query.splitlines(),
                                "committed", "app", lineterm="", n=1)
    print("\n".join(diff))
    print(f"wrote {QUERY_PATH}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
