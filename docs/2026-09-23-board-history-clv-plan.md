# Board history + settled-bet grading (CLV) — design

2026-09-23. Status: design for review, no code. Follows the BFA / Wagerzon sources merged the same day (9ee2622).

## 1. Goal

Record the Unabated board so every settled bet the bets service holds (Kalshi, Novig, BetOnline, BFA, Wagerzon) can be matched to its game and scored against the market's close: CLV per bet, rolled up by book, league, market and period. Open-bet flagging is untouched; this reuses its matcher.

## 2. Facts this rests on (verified 2026-09-23)

- A league snapshot is `content.unabated.com/markets/v2/league/{id}/odds.json?t=<30s bucket>`, gzipped JSON: `odds` keyed `lg{league}:pt{period}:pregame` → rows `{eventId, eventStart, betTypeId, eventTeams[{rotationNumber, …}], sides: {"si0:tid6": {"ms4": {points, americanPrice, bacr, ge, liquidity, sequenceNumber, alternateLines[]}}}}`; `teams` `{id: {name, abbreviation}}`; `marketSources` `[{id, name, isActive, isEnabledForGameOdds, hasLiquidity}]`. Sizes today: NFL 3.3 MB, CFB 8.8 MB, MLB 1.3 MB compressed.
- The board carries current events only; nothing historical is served, and the changes stream is incomplete anonymously.
- `unabated_edge/feed.py::parse_v2` reads the same file but full-game rows only, no periods, no rotations, and keeps raw line dicts; `unabated_edge` records closes for its own two sports (soccer totals, MLB run totals) in `line_snapshots`.
- The matcher (`extension/bets.js` `gameMatches`) is pure: same league, then team pair or rotation with every known team in the game, then time within 30 min of the start. `isMatchable` gates on `status === "open"`. `teams.js` resolves spellings from the registered team lists, `ALIASES` and the learned crosswalk.
- `bets.duckdb::bets` keeps every record ever seen, settled included, and the service is its only writer. DuckDB allows one writer per file, so no second process may open a file another process holds open for writing.

## 3. Components

### 3.1 Recorder — `unabated_ticket/board_history/recorder.py`

- One process, run by launchd every 5 min. Covers the seven US team-sport leagues in `feed.LEAGUES` (NFL 1, CFB 2, NBA 3, CBB 4, MLB 5, NHL 6, WNBA 7); soccer later.
- Per league, a rule decides whether this tick fetches: hourly when no event starts within 24 h, every 15 min otherwise, every 5 min when an event starts within 60 min. A rule, not a setting.
- Its own small parser mirrors `feed.js ingestSnapshotRow`: events with both rotations, every period key, main lines per book with `bacr`. Alternate lines are not recorded in v1. Fixture-tested on a trimmed real snapshot. Fails loudly on a file without an `odds` map, the `feed.js` rule.
- Writes `unabated_ticket/board_history/board_history.duckdb`, opening, writing and closing per poll so the file is free between polls:
  - `events(league, event_id PK, start_time TIMESTAMPTZ, away_team_id, home_team_id, away_name, home_name, away_rotation, home_rotation, first_seen_at, last_seen_at, closed_at)` — upsert.
  - `teams(league, team_id PK, name, abbreviation, last_seen_at)` — upsert; the list `teams.js` registers at runtime.
  - `lines(league, event_id, period, bet_type, book_id, side_index, points, price, fair_prob, ge, sequence_number, changed_at)` — PK on (event_id, period, bet_type, book_id, side_index); rewritten only when points, price or fair changed (the store's #125 compare).
  - `closes` — same columns plus `captured_at`, `close_age_sec`; append-once: at the first poll after `start_time` the event's `lines` rows are copied and `events.closed_at` set. If the recorder was down at kickoff the last lines seen before the start are copied and `close_age_sec` says how old they are.
  - `recorder_runs(league, started_at, ok, error, n_events, n_lines_changed)` — append, pruned like `source_runs`.
- Growth: `lines` stays the size of the live board; `closes` adds roughly 30k rows on a full Saturday, of the order of 10 M rows a year, hundreds of MB. Fine for DuckDB.

### 3.2 Grader — `board_history/grade.py` + `board_history/grade.js`

- Python orchestrates, node matches, so the panel's matcher stays the only matcher.
- Input: `GET 127.0.0.1:8094/bets.json?days=N` for the settled records and the crosswalk (and pins, once the manual-attach feature exists); `board_history.duckdb` read-only for the events of the bets' dates, the teams and the closes, with a retry on the recorder's brief write lock.
- `grade.js` requires `../extension/bets.js` and `teams.js`: registers the stored teams, applies the crosswalk, runs the game rule for settled records through a new exported `matchSettled(records, events)` that lifts the open-only gate, and prints `{betId → eventId, matchedBy: pair | rotation | pin}` or the reason. A settled BFA college row carries no league: it is tried against the CFB and the CBB events of its date; exactly one game whose two names resolve wins, two fit is "ambiguous league", none is unmatched. Never a guess.
- CLV, in Python, from `closes`: the consensus row's `fair_prob` (Unabated's devigged `bacr`) for the bet's side at the main number. v1 scores a bet only when its number equals the close's main number; otherwise it records `close_points` and `points_delta` (signed by side) and leaves `clv_pts` null with reason `off the main number`. Moneylines always score. `clv_pts = fair_close_prob − implied_prob(bet_price)`, positive when the bet beat the market. An alt-number ladder from recorded `alternateLines` is v2.
- Output goes through the service, `POST /grades.json {rows: [...]}`, so the service stays its DB's only writer: `bets.duckdb::bet_grades(bet_id PK, event_id, matched_by, league, bet_points, bet_price, close_points, close_fair_prob, points_delta, clv_pts, reason, graded_at)`; a regrade upserts.
- Crosswalk lessons: a settled bet matched by both team names teaches (venue, league, spelling → Unabated team id) through the existing `POST /crosswalk.json`, same conflict rules. That is where the obscure college spellings will come from, and it improves open-bet matching.
- Cadence: launchd hourly; about a second of work.

### 3.3 Service — `bets_service`

- `POST /grades.json` with the crosswalk routes' guards (Host allowlist, JSON content type, row validation) and `DELETE /grades.json?betId=` for a regrade.
- `/bets.json` embeds `grade` on every record that has one, and a `gradeSummary` rollup by venue, league, bet type and period — count, mean CLV, hit rate, P&L — over 30, 90 and 365 days. The panel does no math.
- `bet_grades` is never pruned; it is the history.

### 3.4 Panel — the Bets tab's "Settled" section

- A new section under the open bets, collapsed by default. Header strip: `Settled · 30 d · CLV +0.8 pts · 61 bets · 54% hit`, with venue and league filters.
- One row per settled bet: date, venue, description, result, close and fair, CLV coloured green or red, the matched game in the tooltip. A row the grader could not match shows the reason, and carries the Attach control when that feature lands.
- Pure builders `settledRows` and `settledSummary` in `betsview.js`, node-tested; rendering in `panel.js`; CSS in `panel.css`.
- Later, not v1: the Edges tab shows your CLV at that book and market beside a line.

### 3.5 Report — `board_history/report.py`

A CLI table of the same rollups for the terminal, useful before the panel section lands and for the numbers the panel does not show.

## 4. What v1 does not do

- No alt-line ladder: a bet off the main number gets `points_delta`, not CLV.
- No intraday line history: current lines plus closes only.
- No backfill: recording starts the day the recorder lands. Every week without it is a week of settled bets with no close.
- No pins UI: manual attach is its own feature; the grader reads pins once they exist.
- No bet placement, no changes-stream use.

## 5. Decisions taken as defaults

1. Baseline: Unabated's consensus fair at close, not one book's price.
2. Cadence: the kickoff-aware rule in 3.1, no setting.
3. Scope: every venue the service holds.
4. Leagues: the seven US team-sport leagues; soccer later.

## 6. Version control, worktree, documentation

- Branch `feature/board-history-clv`, worktree `.claude/worktrees/board-history-clv`; this document is its first commit.
- Commits, in order: (1) this document; (2) recorder, schema, parser fixture, tests; (3) service `bet_grades`, the two routes, the `/bets.json` embedding, tests; (4) grader (Python and node), CLV tests, report; (5) panel Settled section, node tests; (6) launchd plists, docs.
- Docs in the same commits: `unabated_ticket/README.md` (a "Board history and CLV" section, the Bets tab section), root `CLAUDE.md` bullet, `board_history/README.md` (launchd install, paths, troubleshooting).
- Tests: `./unabated_ticket/check.sh` with `board_history/tests` added to its pytest path. Pre-merge review per `CLAUDE.md`; Cal merges.

## 7. Rollout

1. Land the recorder first and alone, and let it run. A week of closes is worth more than the grader.
2. Grader and service routes second; grade once by hand and compare a handful of bets against the closes on the Unabated site.
3. Panel section last.

## 8. Risks

- DuckDB locks: the recorder writes briefly, the grader reads with retry, the service never opens `board_history.duckdb`.
- Snapshot schema drift: the parser fails loudly, and `feed.js` breaking first would show in the panel anyway.
- Bandwidth: 1 to 1.5 GB a day in football season at the 15-min cadence (CFB is 8.8 MB per pull); the rule skips leagues with nothing in the next 24 h.
- Clocks: `eventStart` is ISO from the feed, closes stamp UTC; the BFA plus-7-hours clock affects only the bet side and is already handled.

Effort: recorder about 350 lines plus tests, service 200, grader 300 plus 100 of JavaScript, panel 200, docs. Two to three days.
