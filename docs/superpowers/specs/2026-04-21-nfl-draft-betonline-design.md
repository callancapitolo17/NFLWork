# NFL Draft Portal — BetOnline Adapter

**Status:** Draft (design only; awaiting user review)
**Author:** Callan + Claude
**Date:** 2026-04-21
**Target merge:** before NFL Draft kickoff 2026-04-23 20:00 ET

## 1. Goal

Add BetOnline as the 7th venue in the NFL Draft EV portal so its prices flow
into `draft_odds` alongside Kalshi, DraftKings, FanDuel, Bookmaker, Wagerzon,
and Hoop88. Doing this also unlocks cross-book +EV detection against Kalshi's
`top_N` / `team-drafts-player` / `matchup` series that were previously stored
but unmapped, since BetOnline mirrors those markets with sportsbook odds.

## 2. Non-goals

- No dashboard UI changes. Once rows land in `draft_odds`, the existing
  Cross-Book Grid and +EV Candidates tabs pick them up via the 2-hour
  staleness filter.
- No change to the trade-tape, edge detector, or consensus writer.
- No bet-placement integration (existing `bet_placer/navigator_betonline.py`
  is unrelated to this feature).
- No conversion of BetOnline's O/U `draft-position` lines to Kalshi's
  discrete `top_N` probabilities. That math (interpolating a pick-number
  distribution) is post-draft scope.

## 3. Context and constraints

**BetOnline API discovery** (confirmed 2026-04-21 via live XHR capture + curl probes):

- NFL Draft markets are keyed under `sportsbook/futures-and-props/nfl-draft/*`,
  **not** `football/nfl-draft`. The site URL shows `football` but the offering
  API uses `futures-and-props` → `ContestType: "nfl-draft"`.
- Two-step auth: `GET https://www.betonline.ag/get-token` returns an anonymous
  JWT. Include as `Authorization: Bearer <jwt>` + seven custom headers
  (`gsetting: bolnasite`, `contests: na`, `gmt-offset`, `utc-offset`,
  `actual-time`, `iso-time`, `utc-time`) — identical to the header set already
  used by `bet_logger/scraper_betonline.py`.
- Cloudflare sits in front of `api-offering.betonline.ag`. We reuse the cookie
  jar maintained by `bet_logger/recon_betonline.py`
  (`bet_logger/recon_betonline_cookies.json`), which carries `cf_clearance`.
- Market odds endpoint:
  `POST https://api-offering.betonline.ag/api/offering/Sports/get-contests-by-contest-type2`
  with `{"ContestType":"nfl-draft","ContestType2":"<slug>","filterTime":0}`.
  One call per bucket slug. Menu endpoint `get-menu` enumerates the 11 live slugs.
- Response shape:
  `ContestOfferings.DateGroup[].DescriptionGroup[].TimeGroup[].ContestExtended.ContestGroupLine[].Contestants[]`.
  Each `Contestants[]` entry has `Name` + `Line.MoneyLine.Line` (American odds
  as signed int).
- Total captured: 902 runners across 11 buckets (fixture: 3.4 MB).

**Portal conventions**

- Scrapers emit `OddsRow(book, book_label, book_subject, american_odds,
  fetched_at, market_group)`.
- `nfl_draft/config/markets.py` builds `MARKET_MAP` by re-parsing each book's
  committed fixture and emitting `(book, book_label, book_subject, market_id)`
  tuples. Only structured `market_group`s map; `prop_*` rows intentionally
  fall into `draft_odds_unmapped`.
- `market_id` is constructed deterministically by
  `nfl_draft.lib.market_map.build_market_id(market_type, **kwargs)`. Current
  types: `pick_outright`, `first_at_position`, `top_n_range`, `team_first_pick`,
  `prop`. Adding new types requires a matching case in that function — this
  is the portal's single source of truth for cross-book joins.

## 4. Scope

### 4.1 BetOnline adapter — 10 of 11 buckets

| # | Bucket (slug)                | Runners | Market_group emitted            | Canonical market_type used      | Cross-book overlap |
|---|------------------------------|--------:|---------------------------------|---------------------------------|--------------------|
| 1 | `1st-round`                  |     116 | `pick_outright`                 | `pick_outright`                 | Kalshi, BM, DK, FD, Hoop88 |
| 2 | `1st-round-props`            |      24 | `first_at_position`             | `first_at_position`             | Kalshi, BM, Hoop88 |
| 3 | `to-be-drafted-1st`          |      18 | `first_at_position`             | `first_at_position`             | Kalshi, BM, Hoop88 |
| 4 | `to-be-drafted-2nd`          |      26 | `nth_at_position_2`             | **`nth_at_position` (new)**     | BM                 |
| 5 | `to-be-selected`             |      60 | `top_N_range` (N∈{5,10,32})     | `top_n_range`                   | Kalshi, DK, FD     |
| 6 | `mr-irrelevant`              |      10 | `mr_irrelevant_position`        | **`mr_irrelevant_position` (new)** | BM (currently prop) |
| 7 | `team-to-draft`              |     256 | `team_drafts_player`            | `team_first_pick`               | Kalshi (currently unmapped) |
| 8 | `teams-1st-drafted-position` |     310 | `team_first_pick_position`      | **`team_first_pick_position` (new)** | — (BetOnline-unique) |
| 9 | `matchups`                   |      10 | `matchup_before`                | **`matchup_before` (new)**      | Kalshi (currently unmapped) |
| 10 | `specials`                  |    0-26 | varies (parse what we can)      | `prop` (fallback)               | mostly BetOnline-only |

Total mapped runners targeted: **~830** (most of 902; the rest stay as `prop_*`
in quarantine by design — same pattern as Wagerzon/Bookmaker props).

### 4.2 Deliberately out of scope (v1)

**`draft-position` (46 runners)** — Over/Under lines with a numeric pick
threshold (e.g. "Caleb Downs O/U 9.5"). The existing `draft_odds` schema
stores binary outcomes indexed by `market_id`; encoding the line in the ID
(e.g. `draft_position_ou:caleb-downs:9.5`) is trivial but the cross-book
join against Kalshi's discrete `top_N_*` requires a pick-distribution
interpolation (post-draft work). The recon script still captures this
bucket into the fixture so nothing is lost — we just don't parse it to
`OddsRow` in v1.

**`specials`** — Heterogeneous one-offs (e.g. "How many trades in Round 1
O/U 4.5"). Parse any that match an existing canonical market_type; emit
everything else as `prop_<label>` so it quarantines cleanly. No engineering
time sunk into special-casing single-market novelties 24 hours before
kickoff.

### 4.3 Cross-book gap-fill — Kalshi

BetOnline's arrival makes three currently-unmapped Kalshi series suddenly
cross-bookable:

- `KXNFLDRAFTTEAM-*` (team drafts player) → map to `team_first_pick` ← joins
  BetOnline `team-to-draft`.
- `KXNFLDRAFTMATCHUP-*` (player X drafted before player Y) → map to new
  `matchup_before` ← joins BetOnline `matchups`.
- `KXNFLDRAFTTOP-26-R1-*` (already mapped to `top_32_range`) — verify it
  joins BetOnline `to-be-selected/Drafted in Round 1`. If the BetOnline
  label comes through as "Drafted in Round 1" while Kalshi emits
  `top_32_...`, we align by using `build_market_id("top_n_range",
  range_high=32, range_low=1, player=...)` from both sides.

Update `nfl_draft/config/markets.py::_kalshi_market_id_for` to add cases
for these series. No Kalshi scraper changes required — we already ingest
the tickers.

**Bonus retroactive mapping:** adding the `nth_at_position` canonical type
in §5.2 also lets us rescue Bookmaker's and Wagerzon's currently-quarantined
"2nd/3rd WR Selected" rows by extending `_bm_market_id_for` and
`_wz_market_id_for` to emit that canonical. Same branch, small delta —
folded into commit 2 of §9.

### 4.4 Cross-book gap-fill — other books (audit-driven, not speculative)

Before writing code we query `draft_odds_unmapped` for high-frequency
unmapped `(book, book_label, book_subject)` keys and decide:

1. If ≥3 rows of the same shape exist with a clear canonical analog, add a
   `MARKET_MAP` entry.
2. If one-off, leave quarantined.

This audit becomes one commit. Time-box: 90 minutes. If no high-frequency
gaps exist (likely, since the portal shipped clean on 2026-04-17), skip
the commit and document that in the merge PR.

## 5. Architecture

### 5.1 Files created

- `nfl_draft/scrapers/recon_betonline.py` **(already written in this worktree)**
  — curl_cffi + Cloudflare cookie jar + JWT flow; discovers slugs live via
  `get-menu`; POSTs `get-contests-by-contest-type2` per slug; saves one
  bundle fixture keyed by slug to
  `nfl_draft/tests/fixtures/betonline/draft_markets.json`.
- `nfl_draft/scrapers/betonline.py` — adapter exporting:
  - `parse_response(raw: dict) -> List[OddsRow]` — walks the bundle,
    dispatches per-slug classifiers, yields one OddsRow per Contestant.
  - `fetch_raw() -> dict` — runs the recon flow and returns the envelope.
    Called by `fetch_draft_odds()`. Re-uses `recon_betonline` helpers —
    no duplicated HTTP code.
  - `fetch_draft_odds() -> List[OddsRow]` — composition.
  - Module-level regexes / position maps that `_betonline_entries()` in
    `config/markets.py` re-uses (mirrors the BM/WZ pattern).
- `nfl_draft/tests/unit/test_betonline.py` — parses the committed fixture;
  asserts:
  - Total row count ≥ 800 (allows ±10% drift from the 902-runner snapshot).
  - Per-bucket row counts are within ±10% of the captured snapshot.
  - Every emitted `market_group` is in the v1 allowlist.
  - Spot checks: specific (market_group, book_subject, american_odds)
    triples for 3 well-known markets (1st Overall Pick, 1st QB Selected,
    Drafted Top 10).

### 5.2 Files modified

- `nfl_draft/lib/market_map.py::build_market_id`
  - Add new market_types:
    - `nth_at_position(nth, position, player)` →
      `{nth}_{position_lower}_{slug(player)}` (e.g. `2_wr_jordyn-tyson`).
    - `mr_irrelevant_position(position)` →
      `mr_irrelevant_{position_lower_underscored}` (e.g.
      `mr_irrelevant_wide_receiver`). Position-keyed, not player-keyed —
      matches Bookmaker's "Last Pick Position".
    - `team_first_pick_position(team, position)` →
      `{team_lower_underscored}_first_pick_pos_{position_lower_underscored}`
      (e.g. `arizona_cardinals_first_pick_pos_wide_receiver`). Team-keyed,
      no player.
    - `matchup_before(player_a, player_b)` →
      `matchup_{min(slug_a, slug_b)}_before_{max(slug_a, slug_b)}` — player
      order canonicalised so the two books emit the same ID regardless of
      which side was listed first.
- `nfl_draft/config/markets.py`
  - Add `_betonline_entries()` mirroring the `_bm_entries()` / `_h88_entries()`
    pattern: load fixture, re-parse via `scrapers.betonline.parse_response`,
    emit MARKET_MAP tuples for structured groups.
  - Add `_betonline_market_id_for(group, label, subject, hints)` dispatching
    on the 10 v1 market_groups.
  - Extend `_kalshi_market_id_for` to handle `KXNFLDRAFTTEAM` and
    `KXNFLDRAFTMATCHUP` prefixes (section 4.3).
  - Append `_betonline_entries()` output to the `MARKET_MAP` build list.
- `nfl_draft/run.py`
  - Add `"betonline": "nfl_draft.scrapers.betonline"` to the `SCRAPERS` dict.
- `nfl_draft/README.md`
  - Add BetOnline to the venue list and current-venue-status block.
  - Document auth model: reuses `bet_logger/recon_betonline_cookies.json`.
  - Document refresh path: `python bet_logger/recon_betonline.py`.
- `nfl_draft/scrapers/RECON_README.md`
  - Add BetOnline entry.
- `docs/superpowers/specs/2026-04-21-nfl-draft-betonline-design.md`
  - This file.

### 5.3 Data flow (unchanged at system level)

```
cron → run.py --mode scrape --book all
     → scrapers.betonline.fetch_draft_odds()
        → recon_betonline.fetch_flow()
          → GET /get-token (JWT)
          → POST get-menu                    (discover 11 slugs)
          → POST get-contests-by-contest-type2  × per slug
          ← bundle dict keyed by slug
        → parse_response(bundle)
          → per-slug classifier → OddsRow[]
     → write_or_quarantine(rows)
        → market_map lookup → draft_odds  (mapped)
                             → draft_odds_unmapped  (props / specials one-offs)
```

## 6. Unit boundaries

Each of the following units has ONE clear purpose, a stable interface, and
is independently testable:

- **`recon_betonline`** — authenticated HTTP only. No OddsRow construction.
  Input: nothing. Output: raw bundle dict (or exit code on auth failure).
- **`betonline.parse_response`** — raw bundle → OddsRow[]. Pure function.
  No I/O. Fully testable off the committed fixture.
- **Per-slug classifier functions** (`_classify_pick_outright`,
  `_classify_first_position`, `_classify_team_drafts_player`, etc.) — each
  takes a `DescriptionGroup` + `Contestant` and returns a `market_group`
  + any hints needed to build the canonical ID. Small, local, testable.
- **`_betonline_market_id_for`** — `(group, label, subject, hints) → str |
  None`. Never touches I/O; never imports scrapers. Pure dispatch.
- **MARKET_MAP builder** (`config/markets.py::_betonline_entries`) — loads
  fixture, runs `parse_response`, calls `_betonline_market_id_for` per row.

Files stay small: adapter target ~300 LOC, config extension ~80 LOC,
market_map.py additions ~30 LOC.

## 7. Error handling

- Auth failure (403 / 401 on `get-contests-by-contest-type2`):
  `recon_betonline` exits non-zero with a pointer to
  `bet_logger/recon_betonline.py` for cookie refresh.
- Per-book error isolation already guaranteed by `run.py`'s try/except —
  a BetOnline blowup doesn't stop the other 6 venues.
- Unmapped markets flow into `draft_odds_unmapped` — not an error,
  just "didn't cross-book".

## 8. Testing

- **Unit**: `tests/unit/test_betonline.py` as described in 5.1. Run via
  `/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m pytest
  nfl_draft/tests/ -v`. Target: green.
- **MARKET_MAP regression**: after re-seeding, confirm ≥ 250 new
  `(betonline, *)` rows in the `market_map` table and that a spot-checked
  player (e.g. Caleb Downs) has matching `market_id`s from BetOnline and
  Kalshi for the `first_s` / `top_10` / `pick_10_overall` markets.
- **Live smoke**: `python -m nfl_draft.run --mode scrape --book betonline`
  from the worktree (reading `bet_logger/recon_betonline_cookies.json` from
  main via `_main_repo_root()`). Expect ≥ 800 mapped rows on first run.
  No writes to main DB unless explicitly run from main; worktree tests
  should point at a temp DuckDB.
- **Dashboard smoke** (post-merge, from main): visit Cross-Book Grid,
  confirm BetOnline column populates for overlapping markets within 15
  minutes of the next cron tick.

## 9. Version control

- **Worktree:** `.worktrees/nfl-draft-betonline` **(already created)**
- **Branch:** `feature/nfl-draft-betonline` **(already active)**
- **Commit sequence** (each is independently runnable):
  1. `feat(nfl_draft): add BetOnline recon script + fixture` — recon script
     (already written) + captured 3.4 MB fixture. `RECON_README.md` entry.
  2. `feat(nfl_draft): extend build_market_id with 4 new types` — adds
     `nth_at_position`, `mr_irrelevant_position`, `team_first_pick_position`,
     `matchup_before`. Includes `test_market_map.py` additions covering
     each new type's determinism.
  3. `feat(nfl_draft): add BetOnline adapter (10 of 11 buckets)` — scraper,
     per-slug classifiers, `_betonline_entries()`, unit tests.
  4. `feat(nfl_draft): map Kalshi team/matchup series to new canonicals` —
     `_kalshi_market_id_for` extensions. Re-uses already-existing Kalshi
     fixture. Unit test for the new mappings.
  5. `feat(nfl_draft): register BetOnline in orchestrator` — one-line
     SCRAPERS dict addition in run.py.
  6. `docs(nfl_draft): add BetOnline to README + RECON_README` — same
     commit as the feature per the CLAUDE.md documentation rule.
  7. *(conditional)* `feat(nfl_draft): fill cross-book mapping gaps from
     audit` — only if the `draft_odds_unmapped` audit turns up ≥ 3 unmapped
     markets with clear canonical analogs.

- **Pre-merge review:** full `git diff main..HEAD` audit per CLAUDE.md
  checklist (data integrity, resource safety, edge cases, dead code, log
  hygiene, security). Document findings before merging.
- **Approval gate:** merge to main requires explicit user approval.
  No exceptions.
- **Worktree cleanup:** after merge, `git worktree remove
  .worktrees/nfl-draft-betonline` + `git branch -d
  feature/nfl-draft-betonline`.

## 10. Documentation updates

- `nfl_draft/README.md` — BetOnline venue listed, auth model documented.
- `nfl_draft/scrapers/RECON_README.md` — how to re-run recon_betonline.py.
- This spec (`docs/superpowers/specs/2026-04-21-nfl-draft-betonline-design.md`).
- No CLAUDE.md changes needed; the existing NFL Draft Portal reference in
  the root CLAUDE.md already covers the portal abstractly.

## 11. Risks and mitigations

| Risk | Likelihood | Mitigation |
|---|---|---|
| BetOnline cookies expire mid-draft-weekend | Medium | `bet_logger/recon_betonline.py` is already well-tested; refresh is a 2-min manual run. Document this in README. |
| New canonical types in `build_market_id` introduce collisions with existing IDs | Low | All new IDs are namespaced (`matchup_*`, `mr_irrelevant_*`, `{team}_first_pick_pos_*`); collision-proof by construction. Unit test each. |
| Kalshi team/matchup mapping introduces false cross-books (different player-name spellings) | Medium | The `slug()` canonicalizer already handles case/punctuation; but if Kalshi uses "Shedeur Sanders" and BetOnline uses "Shedeur-Sanders-Jr", mapping diverges silently. Mitigation: the MARKET_MAP unit test prints a diff of subject-slugs per book for the top 30 prospect names, and we eyeball it pre-merge. Anything that looks like the same human with different slugs gets an alias entry in `config/players.py`. Kept manual because a name-matching CI rule costs more to maintain than the error it prevents. |
| `team-to-draft` × 256 rows inflates the `market_map` table size significantly | Low | 256 rows is trivial; total map size still < 5k entries. No concern. |
| Worktree runs hit main's DuckDB | High if careless | Tests use temp DuckDB via the standard pattern. Live smoke from worktree writes to the worktree's copy only after a `cp nfl_draft.duckdb`. |

## 12. Success criteria

- `python -m nfl_draft.run --mode scrape --book betonline` prints
  `mapped=~800+ unmapped<=150` on first run.
- After one cron tick on `main`, the Cross-Book Grid shows BetOnline values
  for at least:
  - 1st-10th Overall Pick outrights
  - 1st WR/QB/RB/CB Selected
  - Drafted Top 5 / 10 players
  - Team-drafts-player for the top 10 players
- `pytest nfl_draft/tests/ -v` is green.
- No stale `draft_odds_unmapped` growth > 200 rows per BetOnline scrape
  cycle.

## 13. Questions for the user

None remaining. All design-space decisions resolved in-conversation 2026-04-21:
- v1 scope = 10 of 11 buckets (skip `draft-position`).
- 4 new canonical types approved conceptually; naming per §5.2.
- Cross-book Kalshi gap-fill included in this branch, not a separate ticket.
- Generic unmapped-markets audit time-boxed to 90 min (one commit if
  anything found, skipped otherwise).
