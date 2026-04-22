"""BetOnline NFL Draft markets adapter.

Consumes the envelope written by scrapers/recon_betonline.py
(one entry per bucket slug; each entry is the raw
api-offering.betonline.ag/get-contests-by-contest-type2 JSON).

Walks `ContestOfferings.DateGroup[].DescriptionGroup[].TimeGroup[].ContestExtended
.ContestGroupLine[].Contestants[]` per bucket; emits one OddsRow per
Contestant with `Name` + `Line.MoneyLine.Line` (signed int American odds).

Classifier dispatch is keyed off the bucket slug (not the Description
text) because BetOnline's `Description` phrasing is inconsistent
across buckets (e.g. 'Team to Draft Kenyon Sadiq' vs
'Arizona Cardinals 1st Drafted Player Position'), while the slug is
deterministic from the menu.
"""

from __future__ import annotations

import re
import sys
from datetime import datetime
from pathlib import Path
from typing import Iterator, List

from nfl_draft.scrapers._base import OddsRow


# ---------------------------------------------------------------------------
# Shared walker
# ---------------------------------------------------------------------------

def _iter_contests(bucket: dict) -> Iterator[tuple[dict, dict, dict]]:
    """Yield (description_group, contest_extended, contest_group_line) triples.

    Flattens the BetOnline offering envelope so classifiers only see a
    single Description + its runner lists. `ContestExtended` may ship as
    dict (single) or list (multi); normalize to list.
    """
    co = (bucket or {}).get("ContestOfferings") or {}
    for dg in (co.get("DateGroup") or []):
        for desc in (dg.get("DescriptionGroup") or []):
            for tg in (desc.get("TimeGroup") or []):
                ce = tg.get("ContestExtended")
                if not ce:
                    continue
                ces = ce if isinstance(ce, list) else [ce]
                for c in ces:
                    for cgl in (c.get("ContestGroupLine") or []):
                        yield desc, c, cgl


def _odds(c: dict) -> int | None:
    """Extract American odds (signed int) from a Contestant dict.

    BetOnline ships American, decimal, and fractional variants; we trust
    the American field. A MoneyLine.Line of 0 means 'no line' (market
    pulled or suspended).
    """
    line = (c.get("Line") or {}).get("MoneyLine") or {}
    val = line.get("Line")
    if val is None:
        return None
    try:
        iv = int(val)
    except (TypeError, ValueError):
        return None
    if iv == 0:
        return None
    return iv


# ---------------------------------------------------------------------------
# Bucket-specific classifiers (one per slug)
# ---------------------------------------------------------------------------

# "10th Overall Pick" / "2nd Overall Pick"
PICK_DESC_RE = re.compile(
    r"^(\d+)(?:st|nd|rd|th)\s+Overall\s+Pick\s*$", re.IGNORECASE,
)

# "First Wide Receiver Drafted" / "First Cornerback Drafted"
FIRST_POS_DESC_RE = re.compile(
    r"^First\s+(.+?)\s+Drafted\s*$", re.IGNORECASE,
)

# "Total Wide Receivers Drafted in 1st Round" / "Total ACC Players Drafted in 1st Round"
FIRST_ROUND_TOTAL_DESC_RE = re.compile(
    r"^Total\s+(.+?)\s+Drafted\s+in\s+1st\s+Round\s*$", re.IGNORECASE,
)

# "2nd Cornerback Selected" / "3rd Wide Receiver Selected"
NTH_POS_DESC_RE = re.compile(
    r"^(\d+)(?:st|nd|rd|th)\s+(.+?)\s+Selected\s*$", re.IGNORECASE,
)

# "Drafted Top 5" / "Drafted Top 10"
TOP_N_DESC_RE = re.compile(r"^Drafted\s+Top\s+(\d+)\s*$", re.IGNORECASE)
# "Drafted in Round 1" -> treat as top_32 (round 1 = first 32 picks)
ROUND_1_DESC_RE = re.compile(r"^Drafted\s+in\s+Round\s+1\s*$", re.IGNORECASE)

# "Team to Draft Kenyon Sadiq"
TEAM_TO_DRAFT_DESC_RE = re.compile(r"^Team\s+to\s+Draft\s+(.+?)\s*$", re.IGNORECASE)

# "Arizona Cardinals 1st Drafted Player Position"
TEAMS_1ST_POS_DESC_RE = re.compile(
    r"^(.+?)\s+1st\s+Drafted\s+Player\s+Position\s*$", re.IGNORECASE,
)

# BetOnline position words -> canonical abbreviations used by
# build_market_id('first_at_position', position=...).
POSITION_MAP = {
    "quarterback": "QB", "qb": "QB",
    "running back": "RB", "rb": "RB",
    "wide receiver": "WR", "wr": "WR",
    "tight end": "TE", "te": "TE",
    "cornerback": "CB", "cb": "CB",
    "safety": "S", "s": "S",
    "linebacker": "LB", "lb": "LB",
    "offensive lineman": "OL", "offensive line": "OL", "ol": "OL",
    "defensive back": "DB", "db": "DB",
    "defensive tackle": "DT", "dt": "DT",
    "defensive end": "DE", "de": "DE",
    "edge": "EDGE",
    "defensive line/edge": "EDGE",
    "defensive line / edge": "EDGE",
}


def _classify_1st_round(desc: dict, ce: dict, cgl: dict, now: datetime) -> Iterator[OddsRow]:
    label = (ce.get("Description") or "").strip()
    if not PICK_DESC_RE.match(label):
        # Unexpected label — yield a prop row so nothing is silently dropped.
        for c in (cgl.get("Contestants") or []):
            name = (c.get("Name") or "").strip()
            american = _odds(c)
            if name and american is not None:
                yield OddsRow(
                    book="betonline", book_label=label, book_subject=name,
                    american_odds=american, fetched_at=now,
                    market_group="prop_1st_round",
                )
        return
    for c in (cgl.get("Contestants") or []):
        name = (c.get("Name") or "").strip()
        american = _odds(c)
        if not name or american is None:
            continue
        yield OddsRow(
            book="betonline", book_label=label, book_subject=name,
            american_odds=american, fetched_at=now,
            market_group="pick_outright",
        )


def _classify_first_at_position(desc: dict, ce: dict, cgl: dict, now: datetime) -> Iterator[OddsRow]:
    """For `to-be-drafted-1st`: 'First Wide Receiver Drafted' with player runners."""
    label = (ce.get("Description") or "").strip()
    m = FIRST_POS_DESC_RE.match(label)
    if not m:
        for c in (cgl.get("Contestants") or []):
            name = (c.get("Name") or "").strip()
            american = _odds(c)
            if name and american is not None:
                yield OddsRow(
                    book="betonline", book_label=label, book_subject=name,
                    american_odds=american, fetched_at=now,
                    market_group="prop_first_at_position",
                )
        return
    pos_raw = m.group(1).strip().lower()
    pos = POSITION_MAP.get(pos_raw)
    if not pos:
        # Known shape but unknown position word — drop to prop so MARKET_MAP
        # surfaces the miss in quarantine rather than silently dropping.
        for c in (cgl.get("Contestants") or []):
            name = (c.get("Name") or "").strip()
            american = _odds(c)
            if name and american is not None:
                yield OddsRow(
                    book="betonline", book_label=label, book_subject=name,
                    american_odds=american, fetched_at=now,
                    market_group="prop_first_at_unknown_position",
                )
        return
    for c in (cgl.get("Contestants") or []):
        name = (c.get("Name") or "").strip()
        american = _odds(c)
        if not name or american is None:
            continue
        yield OddsRow(
            book="betonline", book_label=label, book_subject=name,
            american_odds=american, fetched_at=now,
            market_group="first_at_position",
        )


def _classify_nth_at_position(desc: dict, ce: dict, cgl: dict, now: datetime) -> Iterator[OddsRow]:
    """For `to-be-drafted-2nd`: '2nd Cornerback Selected' with player runners.

    Emits market_group `nth_at_position_N` (e.g. nth_at_position_2) so the
    MARKET_MAP builder can dispatch to build_market_id('nth_at_position',
    nth=N, position=..., player=...).
    """
    label = (ce.get("Description") or "").strip()
    m = NTH_POS_DESC_RE.match(label)
    if not m:
        for c in (cgl.get("Contestants") or []):
            name = (c.get("Name") or "").strip()
            american = _odds(c)
            if name and american is not None:
                yield OddsRow(
                    book="betonline", book_label=label, book_subject=name,
                    american_odds=american, fetched_at=now,
                    market_group="prop_nth_at_position",
                )
        return
    try:
        nth_val = int(m.group(1))
    except (ValueError, IndexError):
        return
    pos_raw = m.group(2).strip().lower()
    pos = POSITION_MAP.get(pos_raw)
    if not pos:
        for c in (cgl.get("Contestants") or []):
            name = (c.get("Name") or "").strip()
            american = _odds(c)
            if name and american is not None:
                yield OddsRow(
                    book="betonline", book_label=label, book_subject=name,
                    american_odds=american, fetched_at=now,
                    market_group="prop_nth_at_unknown_position",
                )
        return
    for c in (cgl.get("Contestants") or []):
        name = (c.get("Name") or "").strip()
        american = _odds(c)
        if not name or american is None:
            continue
        yield OddsRow(
            book="betonline", book_label=label, book_subject=name,
            american_odds=american, fetched_at=now,
            market_group=f"nth_at_position_{nth_val}",
        )


def _classify_to_be_selected(desc: dict, ce: dict, cgl: dict, now: datetime) -> Iterator[OddsRow]:
    """For `to-be-selected`: 'Drafted Top 5/10' or 'Drafted in Round 1'
    with player runners. Maps to canonical top_N_range (N in {5, 10, 32}).
    Round 1 = top 32 (all first-round picks).
    """
    label = (ce.get("Description") or "").strip()
    if ROUND_1_DESC_RE.match(label):
        n = 32
    else:
        m = TOP_N_DESC_RE.match(label)
        if not m:
            for c in (cgl.get("Contestants") or []):
                name = (c.get("Name") or "").strip()
                american = _odds(c)
                if name and american is not None:
                    yield OddsRow(
                        book="betonline", book_label=label, book_subject=name,
                        american_odds=american, fetched_at=now,
                        market_group="prop_to_be_selected",
                    )
            return
        try:
            n = int(m.group(1))
        except (ValueError, IndexError):
            return
    for c in (cgl.get("Contestants") or []):
        name = (c.get("Name") or "").strip()
        american = _odds(c)
        if not name or american is None:
            continue
        yield OddsRow(
            book="betonline", book_label=label, book_subject=name,
            american_odds=american, fetched_at=now,
            market_group=f"top_{n}_range",
        )


def _classify_mr_irrelevant(desc: dict, ce: dict, cgl: dict, now: datetime) -> Iterator[OddsRow]:
    """For `mr-irrelevant`: 'Draft Position of Mr. Irrelevant' where each
    contestant is a position word (e.g. 'Wide Receiver'). Maps to canonical
    mr_irrelevant_position (position-keyed, not player-keyed).
    """
    label = (ce.get("Description") or "").strip()
    for c in (cgl.get("Contestants") or []):
        name = (c.get("Name") or "").strip()  # position word
        american = _odds(c)
        if not name or american is None:
            continue
        yield OddsRow(
            book="betonline", book_label=label, book_subject=name,
            american_odds=american, fetched_at=now,
            market_group="mr_irrelevant_position",
        )


def _classify_team_to_draft(desc: dict, ce: dict, cgl: dict, now: datetime) -> Iterator[OddsRow]:
    """For `team-to-draft`: 'Team to Draft <Player>' where contestants are
    NFL team names. Maps to canonical team_first_pick(team, player) — the
    book_label carries the player; each contestant (team name) is the
    subject.
    """
    label = (ce.get("Description") or "").strip()
    if not TEAM_TO_DRAFT_DESC_RE.match(label):
        for c in (cgl.get("Contestants") or []):
            name = (c.get("Name") or "").strip()
            american = _odds(c)
            if name and american is not None:
                yield OddsRow(
                    book="betonline", book_label=label, book_subject=name,
                    american_odds=american, fetched_at=now,
                    market_group="prop_team_to_draft",
                )
        return
    for c in (cgl.get("Contestants") or []):
        name = (c.get("Name") or "").strip()  # team name
        american = _odds(c)
        if not name or american is None:
            continue
        yield OddsRow(
            book="betonline", book_label=label, book_subject=name,
            american_odds=american, fetched_at=now,
            market_group="team_drafts_player",
        )


def _classify_teams_1st_drafted_position(desc: dict, ce: dict, cgl: dict, now: datetime) -> Iterator[OddsRow]:
    """For `teams-1st-drafted-position`: '<Team> 1st Drafted Player Position'
    where contestants are position words. Team is in the book_label.
    Maps to team_first_pick_position(team, position).
    """
    label = (ce.get("Description") or "").strip()
    if not TEAMS_1ST_POS_DESC_RE.match(label):
        for c in (cgl.get("Contestants") or []):
            name = (c.get("Name") or "").strip()
            american = _odds(c)
            if name and american is not None:
                yield OddsRow(
                    book="betonline", book_label=label, book_subject=name,
                    american_odds=american, fetched_at=now,
                    market_group="prop_teams_1st_drafted_position",
                )
        return
    for c in (cgl.get("Contestants") or []):
        name = (c.get("Name") or "").strip()  # position word
        american = _odds(c)
        if not name or american is None:
            continue
        yield OddsRow(
            book="betonline", book_label=label, book_subject=name,
            american_odds=american, fetched_at=now,
            market_group="team_first_pick_position",
        )


def _classify_matchups(desc: dict, ce: dict, cgl: dict, now: datetime) -> Iterator[OddsRow]:
    """For `matchups`: 'To Be Drafted First' H2H pairs. Every description
    carries the same literal string, so we synthesize the label from the
    two contestants (sorted alphabetically) to disambiguate pairs. Sorting
    makes the label order-independent so cross-book joins hit regardless
    of which side each book listed first.
    """
    contestants = list(cgl.get("Contestants") or [])
    if len(contestants) != 2:
        # Odd shape — emit as prop so data isn't lost.
        label = (ce.get("Description") or "").strip() or "matchup"
        for c in contestants:
            name = (c.get("Name") or "").strip()
            american = _odds(c)
            if name and american is not None:
                yield OddsRow(
                    book="betonline", book_label=label, book_subject=name,
                    american_odds=american, fetched_at=now,
                    market_group="prop_matchup",
                )
        return
    names = sorted([(c.get("Name") or "").strip() for c in contestants])
    if not all(names):
        return
    label = f"{names[0]} vs {names[1]}"
    for c in contestants:
        name = (c.get("Name") or "").strip()
        american = _odds(c)
        if not name or american is None:
            continue
        yield OddsRow(
            book="betonline", book_label=label, book_subject=name,
            american_odds=american, fetched_at=now,
            market_group="matchup_before",
        )


def _classify_1st_round_props(desc: dict, ce: dict, cgl: dict, now: datetime) -> Iterator[OddsRow]:
    """For `1st-round-props`: 'Total X Drafted in 1st Round' with O/U + GroupLine.

    Group X is either a conference (ACC, Big Ten) or a position group
    (Quarterbacks, Wide Receivers). No canonical cross-book market today,
    so emit as a well-labeled prop with line encoded in subject
    ('Over 6.5' / 'Under 6.5') — matches the draft_position encoding so a
    future canonical can be added without data shape changes.
    """
    label = (ce.get("Description") or "").strip()
    line = cgl.get("GroupLine")
    if line is None:
        return
    try:
        line_val = float(line)
    except (TypeError, ValueError):
        return
    for c in (cgl.get("Contestants") or []):
        name = (c.get("Name") or "").strip()  # 'Over' or 'Under'
        american = _odds(c)
        if not name or american is None:
            continue
        subject = f"{name} {line_val:g}"
        yield OddsRow(
            book="betonline", book_label=label, book_subject=subject,
            american_odds=american, fetched_at=now,
            market_group="prop_first_round_total_ou",
        )


# Dispatch table: bucket slug -> classifier. Buckets not yet implemented
# fall through to a generic prop classifier so data is at least captured
# into draft_odds_unmapped.
CLASSIFIERS = {
    "1st-round": _classify_1st_round,
    "1st-round-props": _classify_1st_round_props,
    "to-be-drafted-1st": _classify_first_at_position,
    "to-be-drafted-2nd": _classify_nth_at_position,
    "to-be-selected": _classify_to_be_selected,
    "mr-irrelevant": _classify_mr_irrelevant,
    "team-to-draft": _classify_team_to_draft,
    "teams-1st-drafted-position": _classify_teams_1st_drafted_position,
    "matchups": _classify_matchups,
}


def _classify_generic_prop(desc: dict, ce: dict, cgl: dict, now: datetime,
                           slug: str) -> Iterator[OddsRow]:
    label = (ce.get("Description") or "").strip() or slug
    for c in (cgl.get("Contestants") or []):
        name = (c.get("Name") or "").strip()
        american = _odds(c)
        if not name or american is None:
            continue
        yield OddsRow(
            book="betonline", book_label=label, book_subject=name,
            american_odds=american, fetched_at=now,
            market_group=f"prop_{slug.replace('-', '_')}",
        )


# ---------------------------------------------------------------------------
# Public parser
# ---------------------------------------------------------------------------

def parse_response(raw: dict) -> List[OddsRow]:
    """Walk the envelope produced by recon_betonline.save_fixture and emit
    OddsRows for every Contestant in every bucket.

    Accepts either the full envelope ({"data": {slug: ...}}) or the inner
    bucket dict directly.
    """
    if not isinstance(raw, dict):
        return []
    buckets = raw.get("data") if "data" in raw and isinstance(raw["data"], dict) else raw

    now = datetime.now()
    rows: List[OddsRow] = []
    for slug, bucket in buckets.items():
        classifier = CLASSIFIERS.get(slug)
        for triple in _iter_contests(bucket):
            desc, ce, cgl = triple
            if classifier:
                rows.extend(classifier(desc, ce, cgl, now))
            else:
                rows.extend(_classify_generic_prop(desc, ce, cgl, now, slug))
    return rows


# ---------------------------------------------------------------------------
# Live fetch (delegates to recon_betonline.py)
# ---------------------------------------------------------------------------

def fetch_raw() -> dict:
    """Refresh BetOnline via the recon flow.

    Re-uses the module, not a subprocess: imports recon_betonline and calls
    its session + token + slug-loop helpers directly.
    """
    _THIS = Path(__file__).resolve()
    sys.path.insert(0, str(_THIS.parent.parent.parent))
    from nfl_draft.scrapers import recon_betonline as rb  # type: ignore[import-not-found]

    from curl_cffi import requests as cffi_requests
    session = cffi_requests.Session(impersonate="chrome")
    loaded = rb._load_cookies(session)
    if not loaded:
        raise RuntimeError(
            "No BetOnline cookies — run bet_logger/recon_betonline.py first."
        )
    session.get(rb.SITE_URL, timeout=20)
    token = session.get(rb.TOKEN_URL, timeout=20).json()["token"]
    headers = rb._build_auth_headers(token)
    slugs = rb._discover_slugs(session, headers)
    bundle: dict[str, dict] = {}
    for slug in slugs:
        headers = rb._build_auth_headers(token)
        data = rb._fetch_contest(session, headers, slug)
        if data is not None:
            bundle[slug] = data
    return {"book": "betonline", "data": bundle}


def fetch_draft_odds() -> List[OddsRow]:
    return parse_response(fetch_raw())
