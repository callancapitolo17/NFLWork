# Wagerzon — dynamic league discovery + idgmtyp validation

**Branch:** `worktree-wz-dynamic-discovery` (worktree at
`.claude/worktrees/wz-dynamic-discovery/`)
**Status:** Drafted, awaiting user approval before implementation.

---

## Review Pack

**What we're building.** Replace Wagerzon's hardcoded `lg=` league-ID
list (25 magic numbers for MLB alone, most undocumented) with name-based
lookup against WZ's `ActiveLeaguesHelper.aspx` catalog at scrape time —
same pattern we just merged for BKM. Plus runtime validation of WZ's
`idgmtyp` enum so a future API change is caught loudly, not silently.

**Key decisions** (genuine choice points; not forced):

1. **All 5 sports in one PR.** Shared resolver plumbing means splitting
   by sport would duplicate infrastructure without separable value.
   Risk: bigger blast radius if something regresses. Mitigation: per-sport
   row-count parity checks before merge.
2. **`idgmtyp` validation runs warn-only, not fail-loud.** WZ's response
   shapes are slightly variable (NULL fields when no line posted yet),
   so strict assertions would drop legitimate data. Warnings catch drift
   without breaking the pipeline. Alternative considered: only fail when
   a critical field (e.g., team name) is missing; rejected as less
   informative.
3. **NFL/CBB get marked `# UNVERIFIED` (off-season at refactor time).**
   Catalog currently shows 0 active NFL leagues and 0 active CBB leagues
   so we cannot verify the Description patterns until next season. Same
   pattern as BKM CBB. The scraper will WARN-and-skip cleanly off-season
   and verify in-season.
4. **Match key is `(IndexName, Description)`** — `Description` alone
   would work for current MLB/NBA (no collisions) but the (IndexName,
   Description) pair is the real primary key and protects against future
   WZ data shape changes (e.g., if they add a NCAAB league with the same
   Description as an NBA one).
5. **Include `lg=2911` (1ST HALF TEAM TOTALS) but NOT `lg=4038`
   (PLAYER TO HIT 1ST HOME RUN).** Live spike revealed that the parser
   doesn't currently handle player-prop response shape: `lg=3908`,
   `lg=4717`, `lg=4718`, `lg=1986` already in config emit all-NULL
   rows because the parent dispatch doesn't have a player-prop branch.
   Adding `lg=4038` would just expand the junk surface. Instead, mark
   the 4 existing player-prop entries with a parser-limitation comment
   so the limitation is visible. Building player-prop parsing is a
   separate task.

**Risks / push back here:**

- **WZ is the primary book we trade against.** A regression here breaks
  parlay/single pricing, parlay placement, and the MLB dashboard. The
  blast radius is meaningfully larger than BKM was. Counter-argument:
  the parser logic itself is unchanged — only league resolution changes.
  Tests cover the resolver; the parser is exercised by the existing
  scrape against a real response. Want to do a 2-stage merge here
  (resolver only first, idgmtyp validation second)?
- **Off-season verification gap.** I cannot verify NFL/CBB patterns
  until the next season. We accept that the first in-season scrape
  could surface a "wrong Description pattern" WARN. The diagnostic
  output names what's available so the fix is trivial. Worth it?
- **Player-prop scraping is broken today and stays broken.** The
  spike confirmed the parser silently emits NULL rows for player
  props (lg=3908, 4717, 4718, 1986). The refactor preserves this
  existing behavior but adds a `# parser limitation` comment so
  future maintainers don't assume these markets are usable.
  Building real player-prop parsing is filed as a follow-up.

**Worth understanding** (opt-in learning, anchored to R where possible):

1. **Catalog endpoints as the "schema" of an external API.** Many
   sports books expose a metadata endpoint that lists their markets —
   `ActiveLeaguesHelper` for WZ, `GetRoutingInfo` for BKM. Treating
   that endpoint as the source of truth (instead of a hardcoded list)
   is like calling `colnames(df)` in R instead of hardcoding column
   indices: the code adapts when the schema evolves.
2. **`idgmtyp` as an unlabeled factor.** I taught this earlier — it's
   an integer code with no inline label in the response. Our parser
   maintains the translation. "Response shape validation" is the
   equivalent of `stopifnot(all(c("home", "away", "total") %in%
   names(row)))` — assert the columns you expect are present before
   trusting the parsing logic.

---

## Goal

Eliminate two failure modes in `wagerzon_odds/`:

1. **Hardcoded `lg=` IDs in `config.py`** — 25 MLB IDs, most marked `?`,
   with no validation that they map to the markets the parser expects.
   The bug we just fixed in BKM (league 503 mislabeled as F3) was
   exactly this pattern. WZ's surface is larger.
2. **`idgmtyp` enum drift** — parser routes each `GameChild` to a
   specific handler based on `idgmtyp`. If WZ ever reassigns an
   idgmtyp code, the parser silently parses wrong fields.

## Background

- WZ exposes `GET /wager/ActiveLeaguesHelper.aspx?WT=0` returning a
  514-row catalog of every active league with `IdLeague`, `Description`,
  `IdSport`, `IndexName`, `Active`. No existing scraper calls this.
- Schedule fetch is one POST to `NewScheduleHelper.aspx?lg=<csv>`. The
  response groups results by league inside `result.listLeagues[0]`, each
  with its own `Description` and `Games` array. The parser already does
  a partial dynamic period read at `scraper_v2.py:295-302` (substring
  matching on `Description`), which is exactly the shape this refactor
  formalises.
- `idgmtyp` is internal to each game's market tree (`GameChilds`); it
  is NOT in the catalog. Validation must happen at parse time.

## Pre-flight findings (already done, recon report above)

Cross-checked all current hardcoded `lg=` IDs against the live catalog:

| Sport | IDs in catalog | Notes |
|---|---|---|
| MLB | 25/25 active | All labels resolved; config comments were partly wrong (notably `417 ≠ Game Lines`, it's F5) |
| NBA | 5/5 active | All clean |
| NFL | 1/9 (only lg=124 Super Bowl coin toss visible) | Off-season — rest will reappear when NFL active leagues post |
| CBB | 0/3 | Off-season |
| college_baseball | 1/3 (lg=762 active) | Late-season — `1554, 4321` not currently in catalog |

Plus 2 active MLB leagues we're not yet scraping:
- `lg=2911` "MLB - 1ST HALF TEAM TOTALS" (parser handles via `idgmtyp=66`)
- `lg=4038` "MLB - PLAYER TO HIT 1ST HOME RUN" (parser handles via player-prop branch shared with lg=3908)

## Design

### Phase 1 — Dynamic `lg`-ID resolution

**New `config.py` shape (replaces `url_params: "lg=416,417,..."`):**

```python
SPORTS = {
    "mlb": {
        "sport_key": "baseball_mlb",
        "table_name": "mlb_odds",
        "wz_index_name": "MLB",
        "markets": [
            # Core game lines — main MLB markets the model consumes.
            {"description": "MLB - GAME LINES",                  "period": "fg",   "kind": "lines"},
            {"description": "MLB - 1ST 5 FULL INNINGS",          "period": "F5",   "kind": "lines"},
            {"description": "MLB - 3 INNINGS LINE",              "period": "F3",   "kind": "lines"},
            {"description": "MLB - 7 INNINGS LINE",              "period": "F7",   "kind": "lines"},
            {"description": "MLB - ALTERNATE RUNLINES & TOTALS", "period": "fg",   "kind": "alts"},
            # 3-way moneylines (special: routed via parse_3way_line).
            {"description": "MLB - 1ST INNING WINNER (3 WAY)",   "period": "f1",   "kind": "3way"},
            {"description": "MLB - 1ST 5 INN WINNER (3-WAY)",    "period": "F5",   "kind": "3way"},
            # Team totals.
            {"description": "MLB - TEAM TOTALS",                 "period": "fg",   "kind": "team_total"},
            {"description": "MLB - 1ST HALF TEAM TOTALS",        "period": "Half1","kind": "team_total"},  # NEW
            # Props + specials.
            {"description": "MLB - TEAM TO SCORE 1ST",                          "period": "fg", "kind": "ml_only"},
            {"description": "MLB - TEAM TO SCORE FIRST WINS THE GAME",          "period": "fg", "kind": "ml_only"},
            {"description": "MLB - SCORE 1ST INNING",                           "period": "f1", "kind": "ml_only"},
            {"description": "MLB - TEAM WITH THE HIGHEST SCORING INNING",       "period": "fg", "kind": "3way"},
            {"description": "MLB - DOUBLE RESULT",                              "period": "fg", "kind": "double_result"},
            {"description": "MLB - WINNING MARGIN",                             "period": "fg", "kind": "winning_margin"},
            {"description": "MLB - HITS +RUNS +ERRORS",                         "period": "fg", "kind": "hre"},
            {"description": "MLB - TOTAL HITS",                                 "period": "fg", "kind": "total_only"},
            {"description": "MLB - TOTAL BASES",                                "period": "fg", "kind": "total_only"},
            {"description": "MLB - TOTAL RUNS ODD/EVEN",                        "period": "fg", "kind": "ml_only"},
            # Pitcher props.
            {"description": "MLB - PITCHER STRIKEOUTS",          "period": "fg", "kind": "pitcher_prop"},
            {"description": "MLB - PITCHER HITS ALLOWED",        "period": "fg", "kind": "pitcher_prop"},
            {"description": "MLB - PITCHER WALKS ALLOWED",       "period": "fg", "kind": "pitcher_prop"},
            {"description": "MLB - PITCHER TOTAL OUTS",          "period": "fg", "kind": "pitcher_prop"},
            # Player + market props.
            {"description": "MLB - HOME RUN MARKET",                          "period": "fg", "kind": "player_prop"},
            {"description": "MLB - PLAYER TO HIT A HOME RUN",                 "period": "fg", "kind": "player_prop"},
            {"description": "MLB - PLAYER TO HIT 1ST HOME RUN",               "period": "fg", "kind": "player_prop"},  # NEW
            {"description": "MLB - 1ST PLATE APPEARANCE - EXACT RESULT",      "period": "fg", "kind": "exact_result"},
            {"description": "MLB - 1ST INNING EXACT HITS",                    "period": "fg", "kind": "exact_result"},
        ],
    },
    "nba": {
        "sport_key": "basketball_nba",
        "table_name": "nba_odds",
        "wz_index_name": "NBA",
        "markets": [
            {"description": "NBA - GAME LINES",       "period": "fg",    "kind": "lines"},
            {"description": "NBA - 1ST HALF LINES",   "period": "Half1", "kind": "lines"},
            {"description": "NBA - TEAM TOTALS",      "period": "fg",    "kind": "team_total"},
            {"description": "NBA - 1H TEAM TOTALS",   "period": "Half1", "kind": "team_total"},
            {"description": "NBA - ADJUSTED LINES",   "period": "fg",    "kind": "alts"},
        ],
    },
    "nfl": {
        "sport_key": "americanfootball_nfl",
        "table_name": "nfl_odds",
        "wz_index_name": "NFL",
        # UNVERIFIED (off-season at refactor time, 2026-05-20). First
        # in-season scrape may WARN-and-skip patterns that don't match
        # WZ's actual NFL Descriptions; the diagnostic output will name
        # what IS available. Update this list to match at that point.
        "markets": [
            {"description": "NFL - GAME LINES",                  "period": "fg",    "kind": "lines"},
            {"description": "NFL - 1ST HALF LINES",              "period": "Half1", "kind": "lines"},
            {"description": "NFL - 2ND HALF LINES",              "period": "Half2", "kind": "lines"},
            {"description": "NFL - QUARTER LINES",               "period": "Q",     "kind": "lines"},
            {"description": "NFL - TEAM TOTALS",                 "period": "fg",    "kind": "team_total"},
            {"description": "NFL - ALTERNATE LINES",             "period": "fg",    "kind": "alts"},
            # ... existing 9 lg IDs map to roughly these. Verify in-season.
        ],
    },
    "cbb": {
        "sport_key": "basketball_ncaab",
        "table_name": "cbb_odds",
        "wz_index_name": "CBB",
        # UNVERIFIED (off-season). See NFL note.
        "markets": [
            {"description": "CBB - GAME LINES",      "period": "fg",    "kind": "lines"},
            {"description": "CBB - 1ST HALF LINES",  "period": "Half1", "kind": "lines"},
            # Race-to-X props (race_to_10, race_to_20, race_to_40) — keep
            # separate URL params since they're fetched as standalone
            # leagues, not derivative children. See "Prop params" below.
        ],
        "prop_markets": [
            {"description": "CBB - RACE TO 10",  "key": "race_to_10",  "kind": "ml_only"},
            {"description": "CBB - RACE TO 20",  "key": "race_to_20",  "kind": "ml_only"},
            {"description": "CBB - RACE TO 40",  "key": "race_to_40",  "kind": "ml_only"},
        ],
    },
    "college_baseball": {
        "sport_key": "baseball_ncaa",
        "table_name": "college_baseball_odds",
        "wz_index_name": "NCAA BASEBALL",  # NOTE: differs from "MLB"
        "markets": [
            {"description": "COLLEGE BASEBALL",          "period": "fg", "kind": "lines"},
            {"description": "NCAA GAME",                 "period": "fg", "kind": "lines"},
            {"description": "1ST 5 INNINGS",             "period": "F5", "kind": "lines"},
            # Hard to verify late-season. Patterns are best-effort.
        ],
    },
}
```

**New helper functions in `scraper_v2.py`:**

```python
def fetch_active_leagues(session: requests.Session) -> Optional[list[dict]]:
    """GET /wager/ActiveLeaguesHelper.aspx?WT=0 → list of league dicts.
    Returns None on HTTP error / blocked request. Each row has:
        IdLeague, Description, IdSport, IndexName, Active, ...
    """

def resolve_leagues(catalog: list[dict] | None,
                    index_name: str,
                    markets: list[dict]) -> list[dict]:
    """Match each `markets` entry by (IndexName, Description) against the
    catalog. Returns a list of resolved configs with IdLeague filled in:
        [{"id": "416", "description": "MLB - GAME LINES",
          "period": "fg", "kind": "lines"}, ...]
    Markets that don't match are dropped with a clear WARN that lists
    what IS available under that IndexName (actionable diagnostic).
    """
```

**Wire into `fetch_odds_json`:**

```python
def fetch_odds_json(session, sport: str) -> dict:
    config = get_sport_config(sport)
    catalog = fetch_active_leagues(session)
    resolved = resolve_leagues(catalog, config["wz_index_name"], config["markets"])
    if not resolved:
        print(f"No {sport.upper()} markets resolved from WZ catalog "
              f"— off-season or wrong index_name. Saving empty.")
        return {"result": {"listLeagues": [[]]}}

    lg_csv = ",".join(lg["id"] for lg in resolved)
    print(f"Resolved {len(resolved)} WZ {sport.upper()} markets: "
          f"{', '.join(f'{lg[\"description\"]}→{lg[\"id\"]}' for lg in resolved[:5])}...")

    url = f"{WAGERZON_HELPER_URL}?WT=0&lg={lg_csv}"
    resp = session.get(url, timeout=30, headers={
        "Accept": "application/json, text/plain, */*",
        "X-Requested-With": "XMLHttpRequest",
    })
    resp.raise_for_status()
    return resp.json()
```

**Removed:** the substring-matching period detection at `parse_odds:295-302`.
The resolved markets already carry an explicit `period`, so we look it up
by the league's `IdLeague` in the response. One source of truth.

### Phase 2 — `idgmtyp` shape validation

**New constant in `scraper_v2.py`:**

```python
# Required fields per idgmtyp on each GameChild. If a child is missing
# any of these, log a warning and skip the row — never silently parse
# zeros where real odds should be. Add new idgmtyps here as they appear.
IDGMTYP_REQUIRED_FIELDS = {
    10: ["GameLines"],                  # FG parent (spread+total+ML inside)
    15: ["GameLines"],                  # 1st half
    19: ["GameLines"],                  # Hits / H+R+E totals
    25: ["GameLines"],                  # Alt lines / period totals
    29: ["GameLines", "vtm", "htm"],    # 3-way (parent has 3 outcomes)
    30: ["GameLines"],                  # Pitcher prop
    31: ["GameLines"],                  # Odd/even
    35: ["GameLines"],                  # FG team total
    44: ["GameLines"],                  # Score first
    47: ["GameLines"],                  # Score in 1st inning Y/N
    66: ["GameLines"],                  # 1H team total
}

def validate_idgmtyp_shape(child: dict, log_warn=True) -> bool:
    """Assert child's response shape matches what idgmtyp claims it has.
    Returns True if valid; warns + returns False otherwise.
    """
```

**Wire into `parse_odds`:** each branch that switches on `idgmtyp` first
calls `validate_idgmtyp_shape(child)` and skips on False. Aggregates a
per-scrape WARN summary so the operator sees "N child rows skipped due
to shape mismatch" at the end.

### New markets

- `lg=2911 = "MLB - 1ST HALF TEAM TOTALS"` — already supported by parser
  via `idgmtyp=66` (Team total 1H). Adding to MLB markets list is enough.
- `lg=4038 = "MLB - PLAYER TO HIT 1ST HOME RUN"` — **NOT included.**
  Spike (2026-05-20) confirmed parser doesn't handle player-prop shape;
  adding 4038 would expand the existing NULL-row junk surface. Existing
  player-prop lg entries (1986, 3908, 4717, 4718) get a clear
  parser-limitation comment in the new config. Real player-prop parsing
  is a separate follow-up task.

## Version control

**Branch:** `worktree-wz-dynamic-discovery` (already created).

**Commits** (three, for reviewability):

1. **C1:** `feat(wagerzon): add ActiveLeaguesHelper-driven league resolver`
   - New: `fetch_active_leagues`, `resolve_leagues`, tests.
   - No behavior change to scraper output — pure additive infrastructure.
2. **C2:** `feat(wagerzon): migrate all 5 sport configs to dynamic discovery`
   - Rewrite `config.py` `SPORTS` shape.
   - Wire `fetch_odds_json` to use resolver.
   - Add the 2 new MLB markets.
   - Update parser to use resolved period (drop substring matching).
   - End-to-end verification: row-count parity per sport before/after.
3. **C3:** `feat(wagerzon): add idgmtyp shape validation to parser`
   - Add `IDGMTYP_REQUIRED_FIELDS` + `validate_idgmtyp_shape`.
   - Wire into `parse_odds`, warn-only.
   - Tests.

**Pre-merge:** executive review of `git diff main..HEAD` per project rule.

## Worktree

Already created at `.claude/worktrees/wz-dynamic-discovery/`. Synced to
local main. After merge: `git worktree remove
.claude/worktrees/wz-dynamic-discovery && git branch -d
worktree-wz-dynamic-discovery`.

## Documentation

- **`wagerzon_odds/README.md`** — rewrite to describe ActiveLeaguesHelper
  flow, new SPORT_CONFIG shape, how to add a market. Mirror BKM
  README structure.
- **`wagerzon_odds/CLAUDE.md`** — note the dynamic-discovery contract
  in the "Pitfalls" / "Quick map" section. Update the `idgmtyp` doc
  pointer to the new validation constant.
- **MEMORY.md** — update `scraper_id_validation_gaps.md`: move WZ
  from "HIGH risk" to "safe pattern". Update top-line summary.

## Testing strategy

**Unit tests (`wagerzon_odds/tests/test_resolve_leagues.py`):**

Mirror BKM's `test_resolve_leagues.py` structure. Cover:
- Happy-path: resolves multiple markets in order
- Partial-match: 2/3 patterns resolve, 1 warns
- Missing market: WARN with diagnostic listing what's available
- Wrong index_name: WARN with "no leagues under this IndexName"
- Empty / None catalog: returns empty, no crash
- Malformed catalog entry (missing IdLeague): skip with warning
- Off-season simulation: catalog has rows but none match → empty resolved

**Unit tests (`wagerzon_odds/tests/test_idgmtyp_validation.py`):**

- Validates each idgmtyp shape against fixture children.
- Asserts WARN-not-crash on missing fields.
- Asserts unknown idgmtyp is tolerated (forward-compat) with a different
  WARN.

**Integration verification (manual, before merge):**

1. Run `scraper_v2.py mlb` against live WZ.
2. Compare `wagerzon.duckdb mlb_odds` row count + market distribution to
   a pre-refactor snapshot. Expected: same or more (we added 2 markets).
3. Repeat for NBA (in-season).
4. Run `scraper_v2.py nfl/cbb` — expected to WARN cleanly and write 0 rows.
5. Run `scraper_v2.py college_baseball` — partial resolution; verify
   parser handles whatever resolves.

## Files modified

| File | Change |
|---|---|
| `wagerzon_odds/config.py` | Rewrite `SPORTS` to declare markets by Description pattern |
| `wagerzon_odds/scraper_v2.py` | Add `fetch_active_leagues`, `resolve_leagues`, `validate_idgmtyp_shape`, `IDGMTYP_REQUIRED_FIELDS`. Wire resolver into `fetch_odds_json`. Drop substring period detection in `parse_odds`. |
| `wagerzon_odds/README.md` | Rewrite for new flow |
| `wagerzon_odds/CLAUDE.md` | Update contract notes |
| `wagerzon_odds/tests/test_resolve_leagues.py` | NEW |
| `wagerzon_odds/tests/test_idgmtyp_validation.py` | NEW |
| `memory/scraper_id_validation_gaps.md` | WZ → safe pattern |
| `memory/MEMORY.md` | Update top-line summary |

## Sequence of operations during implementation

1. Spike: hit `lg=4038` live, confirm parser shape match. ~5 min.
2. C1 — resolver infra. ~1 hr. Run tests.
3. C2 — config rewrite + wire. ~2 hr. Run scraper for MLB, diff DB.
4. C3 — idgmtyp validation. ~1 hr. Run scraper, verify no spurious warns.
5. Docs (README/CLAUDE.md/memory). ~30 min.
6. Pre-merge review per CLAUDE.md checklist.
7. Ask user for merge approval.
8. Merge to main, clean up worktree.

Total estimated: 4-5 hours of focused work.

## Open questions for user approval

1. **Phase 2 in this PR vs follow-up?** User already said in this PR.
   Plan reflects that.
2. **Add unscraped markets in this PR vs follow-up?** User already said
   in this PR. Plan reflects that.
3. **Should I do the lg=4038 spike before writing the code (Step 1
   above), or pre-flight as part of the plan approval?** Recommend
   doing it now during plan review since it's 5 minutes and could
   change a small piece of the design.
4. **Two-stage merge** (resolver only, then idgmtyp validation in a
   second merge)? Lower risk if you're concerned about WZ blast radius.
   I lean against — the changes are independent in code but reviewing
   one big diff is easier than two. Your call.

