# NFL Draft Portal — BetOnline Adapter Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add BetOnline as the 7th venue in the NFL Draft EV portal with parsing for all 11 NFL Draft market buckets (~870 mapped runners), plus 5 new canonical market_types that retroactively map existing Kalshi / Bookmaker / Wagerzon quarantine.

**Architecture:** New scraper `nfl_draft/scrapers/betonline.py` (parallel to `bookmaker.py`) emitting structured `OddsRow`s; 5 new market_types added to `nfl_draft/lib/market_map.py::build_market_id`; `_betonline_entries()` builds MARKET_MAP rows from the committed fixture (same pattern as `_bm_entries`); Kalshi's `_kalshi_market_id_for` extended for 4 previously-unmapped series.

**Tech Stack:** Python 3, DuckDB, `curl_cffi` (Chrome impersonation for Cloudflare), existing `kalshi_draft/venv`, pytest.

**Spec:** `docs/superpowers/specs/2026-04-21-nfl-draft-betonline-design.md`

**Working directory:** `.worktrees/nfl-draft-betonline` (branch `feature/nfl-draft-betonline`)

**Commands run from worktree:**
- Tests: `/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m pytest nfl_draft/tests/ -v`
- Single test: `/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m pytest nfl_draft/tests/unit/test_betonline.py::test_name -v`
- Recon: `/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python nfl_draft/scrapers/recon_betonline.py`

**NB on DuckDB lock isolation:** Smoke tests run against a temp DuckDB. Never point the worktree's `nfl_draft.run` at main's `nfl_draft.duckdb` — copy it or use a temp path. See Task 9.

---

## Task 1: Commit BetOnline recon script + fixture

Recon script and 3.4 MB fixture were produced during brainstorming. Verify they
load cleanly, add a regression test, and commit them as the first commit on
the branch.

**Files:**
- Already exists (untracked): `nfl_draft/scrapers/recon_betonline.py`
- Already exists (untracked): `nfl_draft/tests/fixtures/betonline/draft_markets.json`
- Create: `nfl_draft/tests/unit/test_betonline_fixture.py`
- Modify: `nfl_draft/scrapers/RECON_README.md` (add BetOnline section)

- [ ] **Step 1: Write the failing fixture sanity test**

Create `nfl_draft/tests/unit/test_betonline_fixture.py`:

```python
"""Fixture sanity checks — captured offline, parser subagent lives here."""

import json
from pathlib import Path

FIXTURES = Path(__file__).resolve().parent.parent / "fixtures"
BETONLINE_FIXTURE = FIXTURES / "betonline" / "draft_markets.json"


def test_betonline_fixture_present_and_complete():
    """Recon script should have captured all 11 NFL Draft buckets."""
    assert BETONLINE_FIXTURE.exists(), (
        f"Missing fixture at {BETONLINE_FIXTURE}. "
        "Run: python nfl_draft/scrapers/recon_betonline.py"
    )
    envelope = json.loads(BETONLINE_FIXTURE.read_text())
    assert envelope["book"] == "betonline"
    buckets = envelope["data"]
    assert isinstance(buckets, dict)

    # All 11 market buckets captured, each a full offering envelope.
    expected_slugs = {
        "1st-round", "1st-round-props", "draft-position", "matchups",
        "mr-irrelevant", "specials", "team-to-draft",
        "teams-1st-drafted-position", "to-be-drafted-1st",
        "to-be-drafted-2nd", "to-be-selected",
    }
    assert set(buckets.keys()) == expected_slugs, (
        f"fixture buckets != expected. got {sorted(buckets)}"
    )

    # Each bucket should have a ContestOfferings dict (not null/empty).
    for slug, payload in buckets.items():
        co = payload.get("ContestOfferings") or {}
        dg = co.get("DateGroup") or []
        assert dg, f"bucket {slug!r} has empty DateGroup"
```

- [ ] **Step 2: Run test to verify it passes now**

Run: `/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m pytest nfl_draft/tests/unit/test_betonline_fixture.py -v`

Expected: PASS (fixture already captured during brainstorming).

If FAIL with "Missing fixture", run the recon first:
`/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python nfl_draft/scrapers/recon_betonline.py`

- [ ] **Step 3: Add BetOnline section to RECON_README.md**

Append to `nfl_draft/scrapers/RECON_README.md`:

```markdown
## BetOnline

```
python nfl_draft/scrapers/recon_betonline.py
```

Uses `curl_cffi` with Chrome impersonation + the Cloudflare cookie jar
maintained by `bet_logger/recon_betonline.py`. Two-step auth:
`GET /get-token` (anonymous JWT) → authenticated POST to
`api-offering.betonline.ag/api/offering/Sports/get-contests-by-contest-type2`
per NFL Draft sub-market slug (discovered live via `get-menu`).

**When to re-run:**
- Player board rotates (new prospects listed).
- `draft_odds_unmapped` shows a spike in BetOnline rows (new market
  titles or subject spellings).
- First run after re-capturing cookies via `bet_logger/recon_betonline.py`.

**Cookie expiry:** On HTTP 403 (Cloudflare) or 401 (BetOnline app), refresh:
```
python bet_logger/recon_betonline.py
```
Then re-run recon_betonline.py. Refresh tokens rotate automatically on
each authenticated call, so cookies stay warm as long as the
`bet_logger/scraper_betonline.py` pipeline runs daily.
```

- [ ] **Step 4: Commit**

```bash
git add nfl_draft/scrapers/recon_betonline.py
git add nfl_draft/tests/fixtures/betonline/draft_markets.json
git add nfl_draft/tests/unit/test_betonline_fixture.py
git add nfl_draft/scrapers/RECON_README.md
git status --short   # verify only these 4 paths staged
git commit -m "$(cat <<'EOF'
feat(nfl_draft): add BetOnline recon script + fixture

Captures all 11 NFL Draft market buckets from api-offering.betonline.ag
via curl_cffi + Cloudflare cookie jar. 902 runners across 11 buckets in
the committed fixture. RECON_README entry documents the refresh path.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Extend `build_market_id` with 5 new canonical types

Every new canonical type gets one unit test and one implementation branch in
`build_market_id`. TDD one type at a time; one commit for the whole task.

**Files:**
- Modify: `nfl_draft/lib/market_map.py`
- Modify: `nfl_draft/tests/unit/test_market_map.py`

- [ ] **Step 1: Add failing tests for all 5 new types**

Append to `nfl_draft/tests/unit/test_market_map.py`:

```python
def test_build_market_id_nth_at_position():
    assert build_market_id(
        "nth_at_position", nth=2, position="WR", player="Jordyn Tyson"
    ) == "2_wr_jordyn-tyson"


def test_build_market_id_mr_irrelevant_position():
    assert build_market_id(
        "mr_irrelevant_position", position="Wide Receiver"
    ) == "mr_irrelevant_wide_receiver"


def test_build_market_id_team_first_pick_position():
    assert build_market_id(
        "team_first_pick_position", team="Arizona Cardinals", position="Wide Receiver"
    ) == "arizona_cardinals_first_pick_pos_wide_receiver"


def test_build_market_id_matchup_before_canonical_order():
    # Order-independent: whichever player is listed first, we get the same ID.
    a = build_market_id("matchup_before", player_a="Makai Lemon", player_b="Kenyon Sadiq")
    b = build_market_id("matchup_before", player_a="Kenyon Sadiq", player_b="Makai Lemon")
    assert a == b == "matchup_kenyon-sadiq_before_makai-lemon"


def test_build_market_id_draft_position_over_under_half_point():
    assert build_market_id(
        "draft_position_over_under", player="Caleb Downs", line=9.5, direction="under"
    ) == "draft_position_ou_caleb-downs_9p5_under"


def test_build_market_id_draft_position_over_under_whole_point():
    # Whole-point lines drop trailing .0.
    assert build_market_id(
        "draft_position_over_under", player="Abdul Carter", line=3.0, direction="over"
    ) == "draft_position_ou_abdul-carter_3_over"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m pytest nfl_draft/tests/unit/test_market_map.py -v -k "nth_at_position or mr_irrelevant or team_first_pick_position or matchup_before or draft_position_over_under"`

Expected: 6 FAIL with `ValueError: Unknown market_type`.

- [ ] **Step 3: Add helper + new cases to `build_market_id`**

Modify `nfl_draft/lib/market_map.py`. Above `def build_market_id`, add:

```python
def _slug_underscored(name: str) -> str:
    """Like slug() but uses underscores for spaces (team + position names).

    Used by canonical market_types where the subject is a multi-word team or
    position rather than a player. Keeps lowercase + strips punctuation like
    slug() does so 'Wide Receiver' -> 'wide_receiver' and 'Defensive Line/Edge'
    -> 'defensive_line_edge'.
    """
    out = name.lower().replace(".", "").replace("'", "")
    # Collapse any run of whitespace + slashes into a single underscore.
    import re
    out = re.sub(r"[\s/]+", "_", out).strip("_")
    return out


def _format_line(line: float) -> str:
    """Encode a pick line into the market_id safely.

    9.5 -> '9p5' (dot would be lexically noisy in the ID).
    3.0 -> '3' (drop trailing .0 so the ID doesn't gain a decimal needlessly).
    """
    s = f"{line:g}"  # 9.5 -> '9.5', 3.0 -> '3'
    return s.replace(".", "p")
```

Extend `build_market_id` by adding these cases BEFORE the final `raise`:

```python
    if market_type == "nth_at_position":
        return f"{kwargs['nth']}_{kwargs['position'].lower()}_{slug(kwargs['player'])}"
    if market_type == "mr_irrelevant_position":
        return f"mr_irrelevant_{_slug_underscored(kwargs['position'])}"
    if market_type == "team_first_pick_position":
        return (
            f"{_slug_underscored(kwargs['team'])}"
            f"_first_pick_pos_{_slug_underscored(kwargs['position'])}"
        )
    if market_type == "matchup_before":
        a, b = sorted([slug(kwargs["player_a"]), slug(kwargs["player_b"])])
        return f"matchup_{a}_before_{b}"
    if market_type == "draft_position_over_under":
        return (
            f"draft_position_ou_{slug(kwargs['player'])}"
            f"_{_format_line(kwargs['line'])}"
            f"_{kwargs['direction'].lower()}"
        )
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m pytest nfl_draft/tests/unit/test_market_map.py -v`

Expected: ALL tests PASS (existing 6 + new 6).

- [ ] **Step 5: Commit**

```bash
git add nfl_draft/lib/market_map.py nfl_draft/tests/unit/test_market_map.py
git commit -m "$(cat <<'EOF'
feat(nfl_draft): extend build_market_id with 5 new canonical types

Adds nth_at_position, mr_irrelevant_position, team_first_pick_position,
matchup_before, draft_position_over_under — unlocks canonical mapping of
Bookmaker nth-at-position, Kalshi team/matchup/OU series, and the 11
BetOnline draft buckets. matchup_before canonicalises player order; line
encoding uses 'p' for '.' to keep IDs safe.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Retroactive `nth_at_position` mapping for Bookmaker + Wagerzon

Both scrapers already emit `market_group="nth_at_position_N"` for N≥2 rows but
their `_<book>_market_id_for` helpers explicitly return None. Flip them on
now that the canonical type exists.

**Files:**
- Modify: `nfl_draft/config/markets.py` (extend `_bm_market_id_for` + `_wz_market_id_for`)
- Modify: `nfl_draft/tests/unit/test_scraper_parsing.py` (add regression test)

- [ ] **Step 1: Write the regression test**

Append to `nfl_draft/tests/unit/test_scraper_parsing.py`:

```python
def test_bookmaker_fixture_includes_nth_at_position_markets():
    """After v1 canonical-type expansion, BM's 2nd/3rd-WR-selected rows
    should round-trip through config._bm_entries() into MARKET_MAP instead
    of being silently dropped."""
    from nfl_draft.config.markets import _bm_entries
    entries = _bm_entries()
    # Any entry whose market_id starts with `2_` or `3_` followed by a
    # position slug is an nth_at_position canonical.
    nth = [e for e in entries if e[3].startswith(("2_", "3_", "4_", "5_"))]
    assert nth, "BM fixture should produce >=1 nth_at_position mapping"


def test_wagerzon_fixture_includes_nth_at_position_markets():
    from nfl_draft.config.markets import _wz_entries
    entries = _wz_entries()
    nth = [e for e in entries if e[3].startswith(("2_", "3_", "4_", "5_"))]
    # WZ may or may not post these consistently; only assert nonzero if the
    # fixture carries nth_at_position groups.
    from nfl_draft.scrapers.wagerzon import parse_response
    import json
    from pathlib import Path
    fix = Path(__file__).resolve().parent.parent / "fixtures" / "wagerzon" / "draft_markets.json"
    rows = parse_response(json.loads(fix.read_text()))
    has_nth = any(r.market_group.startswith("nth_at_position_") for r in rows)
    if has_nth:
        assert nth, "WZ fixture has nth rows but _wz_entries produced no mapping"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m pytest nfl_draft/tests/unit/test_scraper_parsing.py::test_bookmaker_fixture_includes_nth_at_position_markets -v`

Expected: FAIL (nth rows exist but are currently skipped in `_bm_market_id_for`).

- [ ] **Step 3: Extend `_bm_market_id_for` + `_wz_market_id_for`**

In `nfl_draft/config/markets.py`, find `_bm_market_id_for` and replace the
`# nth_at_position_N (N>=2) -> prop (no structured 2nd-at-position type)`
comment + return None with:

```python
    if group.startswith("nth_at_position_"):
        # e.g. group='nth_at_position_2'; the BM book_label carries the
        # position word (e.g. '2nd Wide Receiver Selected'). Re-parse to get
        # BOTH the nth and the position word.
        m = nth_pos_re.match(label)
        if not m:
            return None
        try:
            nth_val = int(m.group(1))
        except (ValueError, IndexError):
            return None
        pos = position_map.get(m.group(2).strip().lower())
        if not pos:
            return None
        return build_market_id(
            "nth_at_position", nth=nth_val, position=pos, player=subject,
        )
    return None
```

Apply the SAME block verbatim to `_wz_market_id_for` (it uses the same
`nth_pos_re`, `position_map` kwarg names).

- [ ] **Step 4: Run test to verify it passes**

Run: `/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m pytest nfl_draft/tests/unit/test_scraper_parsing.py -v`

Expected: ALL pass. Existing BM/WZ tests remain green (we only unlocked a
previously-None branch; nothing pre-existing should change shape).

- [ ] **Step 5: Commit**

```bash
git add nfl_draft/config/markets.py nfl_draft/tests/unit/test_scraper_parsing.py
git commit -m "$(cat <<'EOF'
feat(nfl_draft): rescue BM/WZ nth-at-position rows from quarantine

Now that build_market_id supports nth_at_position, the BM and WZ
market_id builders can emit canonical IDs for '2nd WR Selected' etc.
instead of returning None. Pure mapping change — no scraper output
change.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: BetOnline adapter — scaffold + 1st-round bucket (pick_outright)

TDD one bucket at a time. Task 4 covers the scaffolding + the first (and
largest) bucket; Tasks 5-12 each add one more bucket to `parse_response`.
Each bucket task has its own tiny test + commit so the diff stays
reviewable.

**Files:**
- Create: `nfl_draft/scrapers/betonline.py`
- Create: `nfl_draft/tests/unit/test_betonline.py`

- [ ] **Step 1: Write failing test for the `1st-round` bucket**

Create `nfl_draft/tests/unit/test_betonline.py`:

```python
"""Offline parser tests for scrapers/betonline.py.

Drives the parser off the committed fixture; never hits the network.
"""

import json
from pathlib import Path
import pytest

from nfl_draft.scrapers.betonline import parse_response

FIXTURE = Path(__file__).resolve().parent.parent / "fixtures" / "betonline" / "draft_markets.json"


@pytest.fixture(scope="module")
def rows():
    envelope = json.loads(FIXTURE.read_text())
    return parse_response(envelope)


def test_emits_pick_outright_rows(rows):
    pick = [r for r in rows if r.market_group == "pick_outright"]
    # At least the 9 Nth-Overall Pick markets × their contestants.
    assert len(pick) >= 100, f"expected >= 100 pick_outright rows, got {len(pick)}"
    # book_label is the market description (e.g. '10th Overall Pick').
    labels = {r.book_label for r in pick}
    assert "10th Overall Pick" in labels
    assert "2nd Overall Pick" in labels
    # Spot-check one: Caleb Downs at +350 for 10th Overall (from brainstorm capture).
    cd = [r for r in pick if r.book_label == "10th Overall Pick" and r.book_subject == "Caleb Downs"]
    assert cd and cd[0].american_odds == 350


def test_all_rows_have_betonline_book_and_valid_odds(rows):
    for r in rows:
        assert r.book == "betonline"
        assert isinstance(r.american_odds, int)
        assert r.american_odds != 0
        assert r.book_subject
        assert r.book_label
```

- [ ] **Step 2: Run test to verify it fails**

Run: `/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m pytest nfl_draft/tests/unit/test_betonline.py -v`

Expected: FAIL with `ModuleNotFoundError: No module named 'nfl_draft.scrapers.betonline'`.

- [ ] **Step 3: Create scraper scaffold + pick_outright classifier**

Create `nfl_draft/scrapers/betonline.py`:

```python
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
                    market_group=f"prop_1st_round",
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


# Dispatch table: bucket slug -> classifier.  Buckets not yet implemented
# fall through to a generic prop classifier so data is at least captured
# into draft_odds_unmapped.
CLASSIFIERS = {
    "1st-round": _classify_1st_round,
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
    its session + token + slug-loop helpers directly. Cleaner than shelling
    out and avoids re-serialising the fixture to disk when all we want is
    the live bundle.
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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m pytest nfl_draft/tests/unit/test_betonline.py -v`

Expected: `test_emits_pick_outright_rows` and `test_all_rows_have_betonline_book_and_valid_odds` PASS.

- [ ] **Step 5: Commit**

```bash
git add nfl_draft/scrapers/betonline.py nfl_draft/tests/unit/test_betonline.py
git commit -m "$(cat <<'EOF'
feat(nfl_draft): BetOnline adapter scaffold + 1st-round bucket

Walks get-contests-by-contest-type2 envelopes; dispatches per-bucket
classifiers keyed off slug. Task 4 implements 1st-round (pick_outright);
remaining 10 buckets follow in sequenced commits. Unimplemented buckets
fall through to a generic prop classifier so data is never dropped.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: BetOnline — `1st-round-props` + `to-be-drafted-1st` buckets

Both emit `first_at_position`. The descriptions look like "1st Quarterback
Selected" (1st-round-props) and "To Be Drafted First" with a position modifier.
Inspect the fixture first to get the exact description phrasings.

**Files:**
- Modify: `nfl_draft/scrapers/betonline.py`
- Modify: `nfl_draft/tests/unit/test_betonline.py`

- [ ] **Step 1: Inspect descriptions in the fixture**

Run:

```bash
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -c "
import json
f = json.load(open('nfl_draft/tests/fixtures/betonline/draft_markets.json'))['data']
for slug in ['1st-round-props', 'to-be-drafted-1st']:
    print('===', slug)
    descs = f[slug]['ContestOfferings']['DateGroup'][0]['DescriptionGroup']
    for d in descs[:12]:
        print(' ', repr(d['Description']))
"
```

Confirms the exact `Description` string format. Expect "1st Quarterback Selected"-shape in `1st-round-props`; `to-be-drafted-1st` likely uses the same phrasing.

- [ ] **Step 2: Write failing test**

Append to `nfl_draft/tests/unit/test_betonline.py`:

```python
def test_emits_first_at_position_from_1st_round_props(rows):
    fap = [r for r in rows if r.market_group == "first_at_position"]
    # Expect at least: QB, WR, RB, TE, OL, EDGE, CB, S, LB across both
    # 1st-round-props + to-be-drafted-1st. >= 30 contestants total.
    assert len(fap) >= 30, f"expected >= 30 first_at_position rows, got {len(fap)}"
    labels = {r.book_label for r in fap}
    # Spot check a well-known one.
    assert any("quarterback" in lbl.lower() for lbl in labels)
```

- [ ] **Step 3: Run test to verify it fails**

Run: `/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m pytest nfl_draft/tests/unit/test_betonline.py::test_emits_first_at_position_from_1st_round_props -v`

Expected: FAIL with `assert 0 >= 30` (no classifier registered yet).

- [ ] **Step 4: Add classifier + position map + regex**

Insert into `nfl_draft/scrapers/betonline.py`, above the `CLASSIFIERS` dict:

```python
# "1st Quarterback Selected" / "1st Offensive Lineman Selected" / etc.
FIRST_POS_DESC_RE = re.compile(
    r"^1st\s+(.+?)\s+Selected\s*$", re.IGNORECASE,
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
    "defensive line/edge": "EDGE",  # some BetOnline labels fold DL+EDGE
    "defensive line / edge": "EDGE",
}


def _classify_first_at_position(desc: dict, ce: dict, cgl: dict, now: datetime) -> Iterator[OddsRow]:
    label = (ce.get("Description") or "").strip()
    m = FIRST_POS_DESC_RE.match(label)
    if not m:
        # Not matchable — emit prop so data is captured, not lost.
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
        # Known market shape but unknown position word — drop to prop so the
        # MARKET_MAP builder surfaces the miss in quarantine, not silently.
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
```

Update the `CLASSIFIERS` dict:

```python
CLASSIFIERS = {
    "1st-round": _classify_1st_round,
    "1st-round-props": _classify_first_at_position,
    "to-be-drafted-1st": _classify_first_at_position,
}
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m pytest nfl_draft/tests/unit/test_betonline.py -v`

Expected: all pass.

- [ ] **Step 6: Commit**

```bash
git add nfl_draft/scrapers/betonline.py nfl_draft/tests/unit/test_betonline.py
git commit -m "$(cat <<'EOF'
feat(nfl_draft): BetOnline 1st-round-props + to-be-drafted-1st classifier

Both buckets emit first_at_position (same shape: '1st QB Selected' with
player contestants). Shares FIRST_POS_DESC_RE + POSITION_MAP which the
MARKET_MAP builder will re-use. Unmappable labels drop to prop fallback
so quarantine surfaces missing position words instead of silent dropouts.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: BetOnline — `to-be-drafted-2nd` bucket (nth_at_position_2)

Same shape as Task 5 but N=2. Emits `market_group="nth_at_position_2"`
so the MARKET_MAP builder can dispatch to `build_market_id("nth_at_position",
nth=2, ...)`.

**Files:**
- Modify: `nfl_draft/scrapers/betonline.py`
- Modify: `nfl_draft/tests/unit/test_betonline.py`

- [ ] **Step 1: Inspect fixture descriptions**

```bash
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -c "
import json
f = json.load(open('nfl_draft/tests/fixtures/betonline/draft_markets.json'))['data']
for d in f['to-be-drafted-2nd']['ContestOfferings']['DateGroup'][0]['DescriptionGroup']:
    print(' ', repr(d['Description']))
"
```

Confirms exact phrasing (likely "2nd Quarterback Selected" or similar).

- [ ] **Step 2: Write failing test**

Append to `test_betonline.py`:

```python
def test_emits_nth_at_position_2(rows):
    nth = [r for r in rows if r.market_group == "nth_at_position_2"]
    assert len(nth) >= 20, f"expected >= 20 nth_at_position_2 rows, got {len(nth)}"
    # book_label should contain '2nd' and a position.
    for r in nth:
        assert "2nd" in r.book_label.lower() or "second" in r.book_label.lower()
```

- [ ] **Step 3: Run test to verify it fails**

Run: `/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m pytest nfl_draft/tests/unit/test_betonline.py::test_emits_nth_at_position_2 -v`

Expected: FAIL.

- [ ] **Step 4: Add classifier**

Insert into `betonline.py`, above the `CLASSIFIERS` dict:

```python
# "2nd Quarterback Selected", "3rd Wide Receiver Selected", etc.
NTH_POS_DESC_RE = re.compile(
    r"^(\d+)(?:st|nd|rd|th)\s+(.+?)\s+Selected\s*$", re.IGNORECASE,
)


def _classify_nth_at_position(desc: dict, ce: dict, cgl: dict, now: datetime) -> Iterator[OddsRow]:
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
        # Unknown position word — emit prop so MARKET_MAP surfaces the miss.
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
```

Update `CLASSIFIERS`:

```python
CLASSIFIERS = {
    "1st-round": _classify_1st_round,
    "1st-round-props": _classify_first_at_position,
    "to-be-drafted-1st": _classify_first_at_position,
    "to-be-drafted-2nd": _classify_nth_at_position,
}
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m pytest nfl_draft/tests/unit/test_betonline.py -v`

Expected: all pass.

- [ ] **Step 6: Commit**

```bash
git add nfl_draft/scrapers/betonline.py nfl_draft/tests/unit/test_betonline.py
git commit -m "$(cat <<'EOF'
feat(nfl_draft): BetOnline to-be-drafted-2nd bucket (nth_at_position_2)

Same shape as Task 5 but parameterised N. NTH_POS_DESC_RE matches '2nd
Wide Receiver Selected' etc. market_group includes the N so the MARKET_MAP
builder can route it to build_market_id('nth_at_position', nth=N, ...).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 7: BetOnline — `to-be-selected` bucket (top_N_range)

Three descriptions: "Drafted in Round 1", "Drafted Top 10", "Drafted Top 5".
Each runner is a player; emit `top_N_range` with N in {5, 10, 32}. Round 1 = 32.

**Files:**
- Modify: `nfl_draft/scrapers/betonline.py`
- Modify: `nfl_draft/tests/unit/test_betonline.py`

- [ ] **Step 1: Write failing test**

Append to `test_betonline.py`:

```python
def test_emits_top_n_range_from_to_be_selected(rows):
    top_rows = [r for r in rows if r.market_group.startswith("top_") and r.market_group.endswith("_range")]
    assert len(top_rows) >= 30, f"expected >= 30 top_N_range rows, got {len(top_rows)}"
    groups = {r.market_group for r in top_rows}
    assert groups == {"top_5_range", "top_10_range", "top_32_range"}, (
        f"expected exactly top_{{5,10,32}}_range, got {groups}"
    )
```

- [ ] **Step 2: Run test to verify it fails**

Run: `/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m pytest nfl_draft/tests/unit/test_betonline.py::test_emits_top_n_range_from_to_be_selected -v`

Expected: FAIL.

- [ ] **Step 3: Add classifier**

Insert into `betonline.py`:

```python
# "Drafted Top 5" / "Drafted Top 10" / "Drafted in Round 1"
TOP_N_DESC_RE = re.compile(
    r"^Drafted\s+Top\s+(\d+)\s*$", re.IGNORECASE,
)
ROUND_1_RE = re.compile(r"^Drafted\s+in\s+Round\s+1\s*$", re.IGNORECASE)


def _classify_to_be_selected(desc: dict, ce: dict, cgl: dict, now: datetime) -> Iterator[OddsRow]:
    label = (ce.get("Description") or "").strip()
    if ROUND_1_RE.match(label):
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
        n = int(m.group(1))
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
```

Update `CLASSIFIERS`:

```python
    "to-be-selected": _classify_to_be_selected,
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m pytest nfl_draft/tests/unit/test_betonline.py -v`

Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add nfl_draft/scrapers/betonline.py nfl_draft/tests/unit/test_betonline.py
git commit -m "$(cat <<'EOF'
feat(nfl_draft): BetOnline to-be-selected bucket (top_N_range)

'Drafted Top 5/10' -> top_{5,10}_range; 'Drafted in Round 1' -> top_32_range.
Matches the canonical used by DK/FD/Kalshi so BetOnline prices cross-book
directly against those venues in the +EV detector.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 8: BetOnline — `mr-irrelevant` bucket (mr_irrelevant_position)

Single description ("Draft Position of Mr. Irrelevant") with position runners.

**Files:**
- Modify: `nfl_draft/scrapers/betonline.py`
- Modify: `nfl_draft/tests/unit/test_betonline.py`

- [ ] **Step 1: Write failing test**

```python
def test_emits_mr_irrelevant_position(rows):
    mi = [r for r in rows if r.market_group == "mr_irrelevant_position"]
    assert len(mi) >= 5, f"expected >= 5 mr_irrelevant rows, got {len(mi)}"
    # Subjects are position words, not player names.
    subjects = {r.book_subject.lower() for r in mi}
    assert any("receiver" in s or "lineman" in s or "back" in s for s in subjects)
```

- [ ] **Step 2: Run test to verify it fails**

Expected: FAIL.

- [ ] **Step 3: Add classifier**

```python
def _classify_mr_irrelevant(desc: dict, ce: dict, cgl: dict, now: datetime) -> Iterator[OddsRow]:
    label = (ce.get("Description") or "").strip()
    # Every contestant here is a position word, e.g. 'Wide Receiver'.
    for c in (cgl.get("Contestants") or []):
        name = (c.get("Name") or "").strip()
        american = _odds(c)
        if not name or american is None:
            continue
        yield OddsRow(
            book="betonline", book_label=label, book_subject=name,
            american_odds=american, fetched_at=now,
            market_group="mr_irrelevant_position",
        )
```

Update `CLASSIFIERS`:
```python
    "mr-irrelevant": _classify_mr_irrelevant,
```

- [ ] **Step 4: Run tests to verify they pass**

Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add nfl_draft/scrapers/betonline.py nfl_draft/tests/unit/test_betonline.py
git commit -m "feat(nfl_draft): BetOnline mr-irrelevant bucket (mr_irrelevant_position)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 9: BetOnline — `team-to-draft` + `teams-1st-drafted-position`

Two buckets, two classifiers:
- `team-to-draft`: Description "Team to Draft <Player>"; contestants are
  team names. Emit `team_drafts_player` (maps to `team_first_pick`
  canonical — uses `team=<team-abbr-or-full>`, `player=<player>`).
- `teams-1st-drafted-position`: Description "<Team> 1st Drafted Player
  Position"; contestants are position words. Emit
  `team_first_pick_position`.

**Files:**
- Modify: `nfl_draft/scrapers/betonline.py`
- Modify: `nfl_draft/tests/unit/test_betonline.py`

- [ ] **Step 1: Write failing tests**

Append to `test_betonline.py`:

```python
def test_emits_team_drafts_player(rows):
    tdp = [r for r in rows if r.market_group == "team_drafts_player"]
    assert len(tdp) >= 200, f"expected >= 200 team_drafts_player rows, got {len(tdp)}"
    # book_label contains 'Team to Draft'; subject is a team name.
    for r in tdp[:5]:
        assert "team to draft" in r.book_label.lower()


def test_emits_team_first_pick_position(rows):
    tfp = [r for r in rows if r.market_group == "team_first_pick_position"]
    assert len(tfp) >= 250, f"expected >= 250 team_first_pick_position rows, got {len(tfp)}"
    for r in tfp[:5]:
        assert "1st drafted" in r.book_label.lower() or "drafted player position" in r.book_label.lower()
```

- [ ] **Step 2: Run tests to verify they fail**

Expected: FAIL.

- [ ] **Step 3: Add classifiers**

```python
TEAM_TO_DRAFT_RE = re.compile(r"^Team\s+to\s+Draft\s+(.+?)\s*$", re.IGNORECASE)
TEAMS_1ST_POS_RE = re.compile(r"^(.+?)\s+1st\s+Drafted\s+Player\s+Position\s*$", re.IGNORECASE)


def _classify_team_to_draft(desc: dict, ce: dict, cgl: dict, now: datetime) -> Iterator[OddsRow]:
    label = (ce.get("Description") or "").strip()
    m = TEAM_TO_DRAFT_RE.match(label)
    if not m:
        # Label shape drift — emit prop but don't drop.
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
    label = (ce.get("Description") or "").strip()
    m = TEAMS_1ST_POS_RE.match(label)
    if not m:
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
```

Update `CLASSIFIERS`:
```python
    "team-to-draft": _classify_team_to_draft,
    "teams-1st-drafted-position": _classify_teams_1st_drafted_position,
```

- [ ] **Step 4: Run tests to verify they pass**

Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add nfl_draft/scrapers/betonline.py nfl_draft/tests/unit/test_betonline.py
git commit -m "feat(nfl_draft): BetOnline team-to-draft + teams-1st-drafted-position

team_drafts_player (256 runners) × team_first_pick_position (310 runners)
= 566 new rows. Both previously existed only at Kalshi and were unmapped;
Task 10 maps Kalshi's KXNFLDRAFTTEAM series so these cross-book.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 10: BetOnline — `matchups` bucket (matchup_before)

H2H "To Be Drafted First" between pairs of players. Each description has 2
contestants; emit two rows with `market_group="matchup_before"` and hint
fields that let the MARKET_MAP builder pair them via canonical ordering.

**Files:**
- Modify: `nfl_draft/scrapers/betonline.py`
- Modify: `nfl_draft/tests/unit/test_betonline.py`

- [ ] **Step 1: Write failing test**

```python
def test_emits_matchup_before_pairs(rows):
    mb = [r for r in rows if r.market_group == "matchup_before"]
    # 5 matchup descriptions × 2 contestants = 10 rows.
    assert len(mb) >= 8, f"expected >= 8 matchup_before rows, got {len(mb)}"
    # book_label must be consistent across the pair so the MARKET_MAP
    # builder can re-pair them; label carries both player names.
    labels = {r.book_label for r in mb}
    # Sanity: each label contains ' vs ' separator or both names.
    for r in mb[:4]:
        assert "vs" in r.book_label.lower() or " - " in r.book_label, (
            f"matchup label {r.book_label!r} needs both players; "
            f"MARKET_MAP builder uses label to pair opponents"
        )
```

- [ ] **Step 2: Run test to verify it fails**

Expected: FAIL.

- [ ] **Step 3: Add classifier**

The BetOnline fixture has `Description="To Be Drafted First"` for every
pair (same label for all 5 markets!), so we must rebuild the label from
the two contestant names to disambiguate. Insert:

```python
def _classify_matchups(desc: dict, ce: dict, cgl: dict, now: datetime) -> Iterator[OddsRow]:
    # BetOnline uses the same Description ('To Be Drafted First') for every
    # pair. Build a per-pair label from the contestants so the MARKET_MAP
    # builder can pair opponents. Sort so 'A vs B' and 'B vs A' collapse.
    contestants = list(cgl.get("Contestants") or [])
    if len(contestants) != 2:
        # Odd shape; fall through to prop so data isn't lost.
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
```

Update `CLASSIFIERS`:
```python
    "matchups": _classify_matchups,
```

- [ ] **Step 4: Run tests to verify they pass**

Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add nfl_draft/scrapers/betonline.py nfl_draft/tests/unit/test_betonline.py
git commit -m "feat(nfl_draft): BetOnline matchups bucket (matchup_before)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 11: BetOnline — `draft-position` bucket (draft_position_over_under)

Each description "Draft Position - <Player>" with GroupLine=N.5 and two
contestants ("Over", "Under"). Emit two rows per description with a
market_group carrying both the player (via book_label) and the line (via a
hint in book_subject). The MARKET_MAP builder re-parses the label+subject
to call `build_market_id("draft_position_over_under", player=..., line=..., direction=...)`.

**Files:**
- Modify: `nfl_draft/scrapers/betonline.py`
- Modify: `nfl_draft/tests/unit/test_betonline.py`

- [ ] **Step 1: Write failing test**

```python
def test_emits_draft_position_ou(rows):
    dpo = [r for r in rows if r.market_group == "draft_position_over_under"]
    # 23 player descriptions × 2 directions = 46 rows.
    assert len(dpo) >= 30, f"expected >= 30 draft_position_ou rows, got {len(dpo)}"
    # Subject encodes BOTH the direction and the line so the MARKET_MAP
    # builder can canonicalize; check format 'Over 9.5' or 'Under 9.5'.
    for r in dpo[:5]:
        assert r.book_subject.startswith(("Over ", "Under ")), (
            f"subject {r.book_subject!r} must start with 'Over '/'Under ' + line"
        )
```

- [ ] **Step 2: Run test to verify it fails**

Expected: FAIL.

- [ ] **Step 3: Add classifier**

```python
DRAFT_POSITION_DESC_RE = re.compile(
    r"^Draft\s+Position\s*[-:]\s*(.+?)\s*$", re.IGNORECASE,
)


def _classify_draft_position(desc: dict, ce: dict, cgl: dict, now: datetime) -> Iterator[OddsRow]:
    label = (ce.get("Description") or "").strip()
    m = DRAFT_POSITION_DESC_RE.match(label)
    if not m:
        # Shape drift -> prop.
        for c in (cgl.get("Contestants") or []):
            name = (c.get("Name") or "").strip()
            american = _odds(c)
            if name and american is not None:
                yield OddsRow(
                    book="betonline", book_label=label, book_subject=name,
                    american_odds=american, fetched_at=now,
                    market_group="prop_draft_position",
                )
        return
    # The line is on the ContestGroupLine; same line applies to both Over/Under.
    line = cgl.get("GroupLine")
    if line is None:
        return
    try:
        line_val = float(line)
    except (TypeError, ValueError):
        return
    # Encode line + direction into book_subject so the MARKET_MAP builder
    # can rebuild (player, line, direction) from (label, subject) alone.
    for c in (cgl.get("Contestants") or []):
        name = (c.get("Name") or "").strip()  # 'Over' or 'Under'
        american = _odds(c)
        if not name or american is None:
            continue
        subject = f"{name} {line_val:g}"  # 'Over 9.5'
        yield OddsRow(
            book="betonline", book_label=label, book_subject=subject,
            american_odds=american, fetched_at=now,
            market_group="draft_position_over_under",
        )
```

Update `CLASSIFIERS`:
```python
    "draft-position": _classify_draft_position,
```

- [ ] **Step 4: Run tests to verify they pass**

Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add nfl_draft/scrapers/betonline.py nfl_draft/tests/unit/test_betonline.py
git commit -m "feat(nfl_draft): BetOnline draft-position bucket (O/U pick line)

Encodes line + direction into book_subject ('Over 9.5') so the MARKET_MAP
builder can construct the canonical market_id with build_market_id(
'draft_position_over_under', player=..., line=..., direction=...).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 12: BetOnline — `specials` bucket (best-effort prop fallback)

Heterogeneous one-offs. For this bucket, we skip the classifier entry
entirely — the generic prop fallback in `parse_response` already handles it
(market_group = `prop_specials`). Just verify the test.

**Files:**
- Modify: `nfl_draft/tests/unit/test_betonline.py` (add an assertion that specials rows emerge as props)

- [ ] **Step 1: Write test**

```python
def test_specials_bucket_emits_prop_rows(rows):
    specials = [r for r in rows if r.market_group == "prop_specials"]
    # 26 contestants in the fixture; should round-trip as props.
    assert len(specials) >= 10, f"expected >= 10 prop_specials rows, got {len(specials)}"
    # Labels include things like 'How Many Trades Will Occur in Round 1 of Draft'.
    assert any("trades" in r.book_label.lower() for r in specials), (
        f"expected at least one trade-count special"
    )
```

- [ ] **Step 2: Run test to verify it passes (already handled by generic fallback)**

Run: `/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m pytest nfl_draft/tests/unit/test_betonline.py -v`

Expected: all pass. No scraper code change.

- [ ] **Step 3: Also add a global coverage assertion**

Append to `test_betonline.py`:

```python
def test_total_row_count_close_to_snapshot(rows):
    """Guard against regressions: fixture snapshot had 902 runners
    (brainstorm capture 2026-04-21). Allow ±15% drift for parser logic that
    skips zero-odds/empty-name edge cases."""
    assert 750 <= len(rows) <= 1050, f"unexpected total row count: {len(rows)}"


def test_v1_market_groups_are_allowlisted(rows):
    allowed = {
        "pick_outright", "first_at_position", "top_5_range", "top_10_range",
        "top_32_range", "nth_at_position_2", "mr_irrelevant_position",
        "team_drafts_player", "team_first_pick_position", "matchup_before",
        "draft_position_over_under",
    }
    actual = {r.market_group for r in rows}
    structured = actual & allowed
    props = {g for g in actual if g.startswith("prop_")}
    unexpected = actual - structured - props
    assert not unexpected, (
        f"unexpected market_groups not in v1 allowlist: {unexpected}"
    )
```

- [ ] **Step 4: Commit**

```bash
git add nfl_draft/tests/unit/test_betonline.py
git commit -m "test(nfl_draft): assert BetOnline specials quarantine + v1 coverage

Specials round-trip via the generic prop fallback (no classifier needed).
Coverage guard: total row count 750-1050; market_groups restricted to
v1 allowlist + prop_* fallbacks. Regression guard on future scraper
changes.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 13: Add `_betonline_entries()` to `config/markets.py`

Plugs BetOnline rows into the live MARKET_MAP. Mirrors the
`_bm_entries()` / `_h88_entries()` pattern.

**Files:**
- Modify: `nfl_draft/config/markets.py`
- Modify: `nfl_draft/tests/unit/test_scraper_parsing.py`

- [ ] **Step 1: Write failing test**

Append to `test_scraper_parsing.py`:

```python
def test_betonline_entries_cover_structured_groups():
    from nfl_draft.config.markets import _betonline_entries
    entries = _betonline_entries()
    assert len(entries) >= 500, f"expected >= 500 betonline entries, got {len(entries)}"
    # Every entry must be a 4-tuple (book, label, subject, market_id) with book='betonline'.
    for e in entries[:20]:
        assert len(e) == 4 and e[0] == "betonline"
    # Coverage: the canonical market_id prefixes we promised in the spec.
    ids = [e[3] for e in entries]
    assert any(m.startswith("pick_") for m in ids), "pick_outright missing"
    assert any(m.startswith("first_") for m in ids), "first_at_position missing"
    assert any(m.startswith("top_5_") or m.startswith("top_10_") or m.startswith("top_32_") for m in ids), "top_N missing"
    assert any(m.startswith("team_") and "_first_pick_" in m for m in ids), "team_first_pick missing"
    assert any(m.startswith("matchup_") for m in ids), "matchup_before missing"
    assert any(m.startswith("draft_position_ou_") for m in ids), "draft_position_ou missing"
    assert any(m.startswith("mr_irrelevant_") for m in ids), "mr_irrelevant_position missing"
    assert any("_first_pick_pos_" in m for m in ids), "team_first_pick_position missing"
```

- [ ] **Step 2: Run test to verify it fails**

Expected: FAIL with `ImportError: cannot import name '_betonline_entries'`.

- [ ] **Step 3: Add `_betonline_entries()` + `_betonline_market_id_for()`**

Add to `nfl_draft/config/markets.py` after `_h88_market_id_for`:

```python
# ---------------------------------------------------------------------------
# BetOnline
# ---------------------------------------------------------------------------

_BE_DRAFT_POSITION_OU_SUBJECT_RE = re.compile(
    r"^(Over|Under)\s+([\d.]+)\s*$", re.IGNORECASE,
)


def _betonline_entries() -> list[tuple[str, str, str, str]]:
    """Build BetOnline MARKET_MAP rows from the committed fixture.

    Each bucket's classifier emits a distinct market_group; the dispatch
    here maps group + label + subject to canonical market_ids using
    build_market_id (or returns None so the row quarantines).
    """
    raw = _load_fixture("betonline")
    if raw is None:
        return []
    from nfl_draft.scrapers.betonline import (
        parse_response, FIRST_POS_DESC_RE, NTH_POS_DESC_RE, PICK_DESC_RE,
        POSITION_MAP, TEAM_TO_DRAFT_RE, TEAMS_1ST_POS_RE, DRAFT_POSITION_DESC_RE,
    )
    rows = parse_response(raw)

    entries: list[tuple[str, str, str, str]] = []
    seen: set[tuple[str, str, str]] = set()
    for r in rows:
        key = (r.book, r.book_label, r.book_subject)
        if key in seen:
            continue
        mid = _betonline_market_id_for(
            r.market_group, r.book_label, r.book_subject,
            pick_re=PICK_DESC_RE, first_pos_re=FIRST_POS_DESC_RE,
            nth_pos_re=NTH_POS_DESC_RE, team_to_draft_re=TEAM_TO_DRAFT_RE,
            teams_1st_pos_re=TEAMS_1ST_POS_RE,
            draft_position_re=DRAFT_POSITION_DESC_RE,
            position_map=POSITION_MAP,
        )
        if mid is None:
            continue  # props + unmappable -> quarantine at runtime.
        entries.append((*key, mid))
        seen.add(key)
    return entries


def _betonline_market_id_for(
    group, label, subject, *,
    pick_re, first_pos_re, nth_pos_re,
    team_to_draft_re, teams_1st_pos_re, draft_position_re,
    position_map,
):
    """BetOnline-specific market_id builder.

    Dispatches on market_group (set by the scraper's classifier).
    """
    if group == "pick_outright":
        m = pick_re.match(label)
        if not m:
            return None
        return build_market_id(
            "pick_outright", pick_number=int(m.group(1)), player=subject,
        )
    if group == "first_at_position":
        m = first_pos_re.match(label)
        if not m:
            return None
        pos = position_map.get(m.group(1).strip().lower())
        if not pos:
            return None
        return build_market_id("first_at_position", position=pos, player=subject)
    if group.startswith("nth_at_position_"):
        m = nth_pos_re.match(label)
        if not m:
            return None
        try:
            nth_val = int(m.group(1))
        except (ValueError, IndexError):
            return None
        pos = position_map.get(m.group(2).strip().lower())
        if not pos:
            return None
        return build_market_id(
            "nth_at_position", nth=nth_val, position=pos, player=subject,
        )
    if group.endswith("_range") and group.startswith("top_"):
        try:
            n = int(group.split("_")[1])
        except (IndexError, ValueError):
            return None
        return build_market_id(
            "top_n_range", range_low=1, range_high=n, player=subject,
        )
    if group == "mr_irrelevant_position":
        return build_market_id("mr_irrelevant_position", position=subject)
    if group == "team_drafts_player":
        m = team_to_draft_re.match(label)
        if not m:
            return None
        player = m.group(1).strip()
        return build_market_id("team_first_pick", team=subject, player=player)
    if group == "team_first_pick_position":
        m = teams_1st_pos_re.match(label)
        if not m:
            return None
        team = m.group(1).strip()
        return build_market_id("team_first_pick_position", team=team, position=subject)
    if group == "matchup_before":
        # book_label format "A vs B" (built by _classify_matchups).
        if " vs " not in label:
            return None
        a, b = [p.strip() for p in label.split(" vs ", 1)]
        return build_market_id("matchup_before", player_a=a, player_b=b)
    if group == "draft_position_over_under":
        m_lbl = draft_position_re.match(label)
        m_sub = _BE_DRAFT_POSITION_OU_SUBJECT_RE.match(subject)
        if not m_lbl or not m_sub:
            return None
        player = m_lbl.group(1).strip()
        direction = m_sub.group(1).lower()
        try:
            line_val = float(m_sub.group(2))
        except (TypeError, ValueError):
            return None
        return build_market_id(
            "draft_position_over_under",
            player=player, line=line_val, direction=direction,
        )
    return None
```

Also append `_betonline_entries()` to the `MARKET_MAP` aggregator at the
bottom of the file:

```python
MARKET_MAP: list[tuple[str, str, str, str]] = list(
    dict.fromkeys(
        _dk_entries()
        + _fd_entries()
        + _bm_entries()
        + _wz_entries()
        + _h88_entries()
        + _kalshi_entries()
        + _betonline_entries()
    ).keys()
)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m pytest nfl_draft/tests/unit/test_scraper_parsing.py -v`

Expected: all pass. If `_betonline_entries` count < 500, debug by printing
per-group counts.

- [ ] **Step 5: Commit**

```bash
git add nfl_draft/config/markets.py nfl_draft/tests/unit/test_scraper_parsing.py
git commit -m "$(cat <<'EOF'
feat(nfl_draft): add _betonline_entries() to MARKET_MAP builder

~500+ new (book, label, subject, market_id) tuples feeding market_map at
seed time. Dispatches the 10 structured market_groups emitted by
scrapers.betonline to canonical market_ids via build_market_id.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 14: Map Kalshi `KXNFLDRAFTTEAM` / `KXNFLDRAFTMATCHUP` / `KXNFLDRAFTOU` series

Three new cases in `_kalshi_market_id_for`. With BetOnline now live, these
previously-unmappable series start cross-booking.

**Files:**
- Modify: `nfl_draft/config/markets.py`
- Modify: `nfl_draft/tests/unit/test_scraper_parsing.py`

- [ ] **Step 1: Write failing test**

Append to `test_scraper_parsing.py`:

```python
def test_kalshi_team_series_now_maps():
    """After Task 14, Kalshi's KXNFLDRAFTTEAM series rescues from quarantine.
    (KXNFLDRAFTMATCHUP is intentionally not mapped in v1 — tag->name
    resolver is post-v1 work.)"""
    from nfl_draft.config.markets import _kalshi_entries
    entries = _kalshi_entries()
    ids = [e[3] for e in entries]
    # team_first_pick canonical IDs have the form '{team}_first_pick_{player}'.
    # Exclude team_first_pick_position ('_first_pick_pos_') which doesn't exist in Kalshi fixture.
    team_first_pick = [m for m in ids if "_first_pick_" in m and "_first_pick_pos_" not in m]
    assert len(team_first_pick) >= 100, (
        f"Expected >= 100 Kalshi team_first_pick entries after mapping "
        f"KXNFLDRAFTTEAM (fixture has 512 markets); got {len(team_first_pick)}"
    )
```

- [ ] **Step 2: Run test to verify it fails**

Expected: FAIL.

- [ ] **Step 3: Extend `_kalshi_market_id_for`**

In `nfl_draft/config/markets.py`, find `_kalshi_market_id_for` and add the
following cases BEFORE the final `return None`:

```python
    # Team drafts player (e.g. KXNFLDRAFTTEAM-26-BAL-CWAR: Ravens drafts
    # Cam Ward). Ticker segment[2] is the 3-letter team code, segment[3]
    # is the player tag. subject is the player display name
    # ('Cam Ward') from yes_sub_title/custom_strike.Person.
    if series_ticker == "KXNFLDRAFTTEAM":
        parts = ticker.split("-")
        if len(parts) >= 4:
            team_code = parts[2]
            return build_market_id(
                "team_first_pick", team=team_code, player=subject,
            )
        return None

    # Matchups (e.g. KXNFLDRAFTMATCHUP-26-PLAYERA-PLAYERB).
    # Kalshi's market `yes_sub_title` is typically one of the two players;
    # the other is in the ticker's segments[3..4]. We need both names to
    # build the canonical ordered ID. If we only have one, drop.
    if series_ticker == "KXNFLDRAFTMATCHUP":
        parts = ticker.split("-")
        if len(parts) >= 5:
            # Both player tags live in parts[3] and parts[4]; they're
            # opaque short codes. The canonical ID uses the DISPLAY name,
            # not the code, so we'd need the other player's display name.
            # Fallback: use subject as player_a and parts[4] tag as
            # player_b — the adapter can still cross-book if it emits the
            # same (a, b) pair order under the same canonicalisation.
            # Prefer skipping when we can't resolve both names.
            return None  # Post-v1: add Kalshi player-tag -> name resolver.
        return None

    # Over/under pick lines (KXNFLDRAFTOU-26-CDOWN-9p5 or similar).
    # Segment[3] is the line encoded like '9p5'. subject is the player name
    # ('Caleb Downs'). market has a `yes_sub_title` that distinguishes
    # Over from Under — look for 'over' / 'under' in the ticker's tail
    # or in yes_sub_title. If neither is present, skip.
    if series_ticker == "KXNFLDRAFTOU":
        parts = ticker.split("-")
        if len(parts) < 4:
            return None
        line_tag = parts[3]  # 'T9p5' or '9p5' depending on Kalshi convention
        try:
            line_str = line_tag.lstrip("T").replace("p", ".")
            line_val = float(line_str)
        except (ValueError, IndexError):
            return None
        # Direction: Kalshi typically uses separate YES markets per side;
        # default to 'over' since YES='player drafted over pick X' is most
        # common. If the ticker encodes direction, refine here.
        direction = "over" if "OVER" in ticker.upper() else "under" if "UNDER" in ticker.upper() else "over"
        return build_market_id(
            "draft_position_over_under",
            player=subject, line=line_val, direction=direction,
        )
```

**Note on matchups:** The fixture's KXNFLDRAFTMATCHUP series has opaque
player tags that don't match the display names BetOnline uses. Mapping
them correctly requires a tag→name lookup that we don't have today. Skip
in v1; document as post-v1 work in the commit message.

- [ ] **Step 4: Run tests to verify they pass**

Expected: `test_kalshi_team_and_matchup_series_now_map` passes for the
`team_first_pick` half; the matchup half either passes if the fixture has
mappable markets or skips gracefully. Adjust the assertion to
`has_matchup or len([m for m in ids if '_first_pick_' in m]) >= 100` if
matchups prove unmappable.

- [ ] **Step 5: Commit**

```bash
git add nfl_draft/config/markets.py nfl_draft/tests/unit/test_scraper_parsing.py
git commit -m "$(cat <<'EOF'
feat(nfl_draft): map Kalshi KXNFLDRAFTTEAM + KXNFLDRAFTOU to new canonicals

Rescues ~500 Kalshi team-drafts-player rows and the O/U pick-line series
(0 live markets today; auto-joins BetOnline draft-position when reposted)
from quarantine. KXNFLDRAFTMATCHUP left unmapped in v1 — tag->name
resolver needed to match BetOnline's display-name-based matchup_before
canonicalisation.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 15: Register BetOnline in `run.py` orchestrator

One-line change.

**Files:**
- Modify: `nfl_draft/run.py`

- [ ] **Step 1: Add entry to `SCRAPERS` dict**

In `nfl_draft/run.py` find the `SCRAPERS = {...}` dict and add:

```python
    "betonline": "nfl_draft.scrapers.betonline",
```

So the dict now reads:

```python
SCRAPERS = {
    "kalshi": "nfl_draft.scrapers.kalshi",
    "draftkings": "nfl_draft.scrapers.draftkings",
    "fanduel": "nfl_draft.scrapers.fanduel",
    "bookmaker": "nfl_draft.scrapers.bookmaker",
    "wagerzon": "nfl_draft.scrapers.wagerzon",
    "hoop88": "nfl_draft.scrapers.hoop88",
    "betonline": "nfl_draft.scrapers.betonline",
}
```

- [ ] **Step 2: Verify module import works**

Run: `/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -c "from nfl_draft.scrapers import betonline; print(betonline.fetch_draft_odds.__doc__ or 'ok')"`

Expected: prints `ok` or a docstring — no ImportError.

- [ ] **Step 3: Commit**

```bash
git add nfl_draft/run.py
git commit -m "feat(nfl_draft): register BetOnline in orchestrator

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 16: Update `nfl_draft/README.md`

**Files:**
- Modify: `nfl_draft/README.md`

- [ ] **Step 1: Update current-venue-status block**

Open `nfl_draft/README.md`. Replace lines 4-8 (the current-venue-status block) with:

```markdown
Trader's-cockpit web portal for surfacing +EV NFL Draft bets across
Kalshi and 6 sportsbooks (DraftKings, FanDuel, Bookmaker, Wagerzon, Hoop88, BetOnline).

**Current venue status** (2026-04-21):
- Kalshi, DraftKings, Wagerzon, Hoop88, Bookmaker, BetOnline — posting draft markets, all scraping live.
- FanDuel — draft page temporarily offline as of 2026-04-20 (they pulled the "NFL Draft" tab from `customPageId=nfl`'s layout). Scraper is intact and will pick up automatically when FD reposts; the dashboard's staleness filter hides FD rows while they're gone.
```

- [ ] **Step 2: Add BetOnline to the `--mode scrape --book all` description**

Find line 15 ("pulls odds from all 6 venues...") and change to:

```markdown
- `--mode scrape --book all` — pulls odds from all 7 venues (kalshi, draftkings, fanduel, bookmaker, wagerzon, hoop88, betonline), devigs, writes to `draft_odds`, triggers legacy `kalshi_draft/edge_detector.py` + `consensus.py`
```

- [ ] **Step 3: Add a BetOnline credentials section**

Insert into the "Credentials" block after the Hoop88 section, before Item 3
("Run the one-time migration..."):

```markdown
   - **BetOnline** (`scrapers/betonline.py` -> `scrapers/recon_betonline.py`):
     no Keycloak login needed at runtime — the NFL Draft markets are
     anonymous-readable. The scraper uses a `curl_cffi` Chrome-impersonation
     session + the Cloudflare cookie jar at
     `bet_logger/recon_betonline_cookies.json` that `bet_logger/recon_betonline.py`
     maintains. Two-step auth: `GET /get-token` (anonymous JWT) →
     authenticated POST to `api-offering.betonline.ag/api/offering/Sports/get-contests-by-contest-type2`
     per NFL Draft bucket slug. On HTTP 403 (Cloudflare) or 401 (app), refresh
     cookies via `python bet_logger/recon_betonline.py`.
```

- [ ] **Step 4: Commit (combined docs commit)**

We'll batch Task 16 + Task 17 into one docs commit below.

---

## Task 17: Update `nfl_draft/scrapers/RECON_README.md` (if not already done in Task 1)

Task 1 already added the BetOnline section. Verify it's present; if not,
re-apply from Task 1 Step 3.

- [ ] **Step 1: Verify**

Run: `grep -q "^## BetOnline" nfl_draft/scrapers/RECON_README.md && echo ok || echo missing`

Expected: `ok`.

- [ ] **Step 2: Combined docs commit for Task 16 + 17**

```bash
git add nfl_draft/README.md nfl_draft/scrapers/RECON_README.md
git commit -m "$(cat <<'EOF'
docs(nfl_draft): document BetOnline venue, auth flow, refresh path

README venue count 6 -> 7; credentials section covers the two-step
anonymous-JWT + Cloudflare-cookie auth. RECON_README documents how to
re-run recon_betonline.py when cookies expire.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 18: Live smoke test (worktree-safe)

Confirm the end-to-end flow works before claiming the feature is done.
Runs from the worktree but points at a temp DuckDB so main's data is not
touched.

**Files:**
- None modified; this is a dry-run verification.

- [ ] **Step 1: Copy main's DuckDB into a temp path**

```bash
cp /Users/callancapitolo/NFLWork/nfl_draft/nfl_draft.duckdb /tmp/nfl_draft_betonline_smoke.duckdb
```

- [ ] **Step 2: Patch the DB_PATH via env shim**

Run the scrape with `NFL_DRAFT_DB_PATH` pointing at the temp DB IF the
codebase supports that env override. If not, temporarily edit
`nfl_draft/lib/db.py::DB_PATH` to `Path("/tmp/nfl_draft_betonline_smoke.duckdb")`
for the duration of this test (DO NOT COMMIT the edit).

```bash
# Verify db.py supports env override before editing:
grep -n "DB_PATH" nfl_draft/lib/db.py | head -5
```

If there's no env override, edit `db.py` line with `DB_PATH =` to:

```python
DB_PATH = Path("/tmp/nfl_draft_betonline_smoke.duckdb")
```

- [ ] **Step 3: Run the scrape from the worktree**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/nfl-draft-betonline
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m nfl_draft.run --mode scrape --book betonline 2>&1 | tail -30
```

Expected output includes a line like:
```
[scrape] betonline: mapped=850+ unmapped<=100
```

- [ ] **Step 4: Verify rows landed in the temp DB**

```bash
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -c "
import duckdb
con = duckdb.connect('/tmp/nfl_draft_betonline_smoke.duckdb', read_only=True)
total = con.execute(\"SELECT COUNT(*) FROM draft_odds WHERE book='betonline'\").fetchone()[0]
by_prefix = con.execute(\"SELECT split_part(market_id, '_', 1), COUNT(*) FROM draft_odds WHERE book='betonline' GROUP BY 1 ORDER BY 2 DESC\").fetchall()
print(f'betonline rows in temp db: {total}')
for p, n in by_prefix:
    print(f'  {p}: {n}')
"
```

Expected: `betonline rows in temp db: 850+` with at least these prefixes:
`pick`, `first`, `top`, `team`, `matchup`, `draft_position_ou`, `mr_irrelevant`.

- [ ] **Step 5: Cross-book spot check**

```bash
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -c "
import duckdb
con = duckdb.connect('/tmp/nfl_draft_betonline_smoke.duckdb', read_only=True)
# A market that should have rows from multiple books.
mid = con.execute(\"SELECT market_id FROM draft_odds WHERE book='betonline' AND market_id LIKE 'pick_1_overall_%' LIMIT 1\").fetchone()
if not mid:
    print('NO BETONLINE 1ST OVERALL ROW FOUND — investigate')
else:
    rows = con.execute('SELECT book, american_odds FROM draft_odds WHERE market_id = ?', [mid[0]]).fetchall()
    print(f'market_id={mid[0]}:')
    for b, a in rows: print(f'  {b}: {a}')
"
```

Expected: at least 2 books (betonline + one of kalshi/bm/dk/fd/hoop88) with
prices on the same market_id.

- [ ] **Step 6: Revert db.py change + clean up**

If you edited `db.py`, revert it:

```bash
git checkout -- nfl_draft/lib/db.py
```

Remove the temp DB:

```bash
rm /tmp/nfl_draft_betonline_smoke.duckdb
```

- [ ] **Step 7: No commit for this task (verification only).**

---

## Task 19: Pre-merge review + merge approval

**Files:** none — review-only task.

- [ ] **Step 1: Run the full test suite**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/nfl-draft-betonline
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m pytest nfl_draft/tests/ -v 2>&1 | tail -40
```

Expected: ALL GREEN. If any regression, STOP and fix before merging.

- [ ] **Step 2: Full diff audit**

```bash
git diff main..HEAD --stat
git diff main..HEAD | wc -l
```

Scan the full diff for the CLAUDE.md pre-merge review checklist:
- Data integrity: no duplicate writes, proper deduplication, incomplete records filtered
- Resource safety: all DB connections use context managers
- Edge cases: empty fixture, missing cookies, 0-line markets
- Dead code: no unused flags / imports
- Log/disk hygiene: no unbounded writes
- Security: no cookies/tokens / secrets in logs

Record findings in a review note (not committed) before proceeding.

- [ ] **Step 3: Request merge approval**

Ping the user:

> "All tests green; diff audit complete; no blockers found. Ready to merge
> `feature/nfl-draft-betonline` into `main`. Proceed?"

**DO NOT MERGE WITHOUT EXPLICIT USER APPROVAL.** Per CLAUDE.md branch
rules.

- [ ] **Step 4: On approval, merge**

```bash
cd /Users/callancapitolo/NFLWork
git checkout main
git merge --no-ff feature/nfl-draft-betonline -m "Merge feature/nfl-draft-betonline: add BetOnline venue to NFL Draft portal"
```

- [ ] **Step 5: Post-merge verification on main**

```bash
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m pytest nfl_draft/tests/ -v 2>&1 | tail -10
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m nfl_draft.lib.seed
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -c "
import duckdb
con = duckdb.connect('nfl_draft/nfl_draft.duckdb', read_only=True)
n = con.execute(\"SELECT COUNT(*) FROM market_map WHERE book='betonline'\").fetchone()[0]
print(f'BetOnline MARKET_MAP entries after seed: {n}')
"
```

Expected: MARKET_MAP has ≥500 BetOnline entries; pytest still green.

- [ ] **Step 6: Cleanup — worktree + branch**

```bash
cd /Users/callancapitolo/NFLWork
git worktree remove .worktrees/nfl-draft-betonline
git branch -d feature/nfl-draft-betonline
```

---

## Task 20 (conditional): Unmapped-markets audit

Only runs if there's spare time after Task 19. Time-box: 90 min.

- [ ] **Step 1: Query the unmapped table**

```bash
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -c "
import duckdb
con = duckdb.connect('nfl_draft/nfl_draft.duckdb', read_only=True)
rows = con.execute('''
  SELECT book, book_label, COUNT(*) AS n
  FROM draft_odds_unmapped
  GROUP BY 1, 2
  HAVING n >= 3
  ORDER BY 3 DESC
  LIMIT 30
''').fetchall()
for r in rows: print(f'  {r[0]:12s} | {r[2]:4d} | {r[1][:70]}')
"
```

- [ ] **Step 2: Decide per entry**

For each high-frequency (book, label) pair:
- Does it map to an existing canonical via `build_market_id`? If yes, add
  a MARKET_MAP entry (extend that book's `_<book>_market_id_for`).
- One-off or genuinely unique? Leave in quarantine.

- [ ] **Step 3: Commit (if any changes)**

If you added mappings, commit them:

```bash
git add nfl_draft/config/markets.py nfl_draft/tests/unit/test_scraper_parsing.py
git commit -m "feat(nfl_draft): fill cross-book mapping gaps from quarantine audit

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

If no changes needed, note in this task's checkbox: "Audit surfaced no
high-frequency gaps — skipped."

---

## Summary

- **Commits planned:** 14 (+1 conditional).
- **Files created:** `scrapers/betonline.py`, `tests/unit/test_betonline.py`,
  `tests/unit/test_betonline_fixture.py`.
- **Files modified:** `lib/market_map.py`, `config/markets.py`, `run.py`,
  `tests/unit/test_market_map.py`, `tests/unit/test_scraper_parsing.py`,
  `README.md`, `scrapers/RECON_README.md`.
- **Already in worktree (untracked):** `scrapers/recon_betonline.py`,
  `tests/fixtures/betonline/draft_markets.json`.
- **Pre-merge gate:** ALL tests green, full diff audit, user approval.
- **Post-merge cleanup:** remove worktree + branch.
