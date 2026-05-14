# FanDuel recon — 2026-05-13

Probed via `/Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python` against the
live FanDuel event-page endpoint using `FanDuelClient` from
`mlb_sgp/fd_client.py`. Cross-checked against `fd_odds/fd.duckdb::mlb_odds`.

Reference events probed (all returned identical F7 market set):
- `35604278` — Kansas City Royals @ Chicago White Sox
- `35604277` — St. Louis Cardinals @ Athletics
- `35604276` — Chicago Cubs @ Atlanta Braves

## F7 markets

**FanDuel does NOT post the F7 line markets that exist for F5.** Across every
event probed, the *only* "First 7 Innings" market FD returns is:

- `'First 7 Innings Result'` — a 3-way moneyline (Team A / Tie / Team B).
  Example runners: `Kansas City Royals -118`, `Tie +630`, `Chicago White Sox +160`.

For comparison, FD posts these F5 markets on the same event (none of which
have F7 analogues):

- `'First 5 Innings Run Line'`
- `'First 5 Innings Total Runs'`
- `'First 5 Innings Money Line'`
- `'First 5 Innings Alternate Run Lines'`
- `'First 5 Innings Alternate Total Runs'`

So `mlb_bets_book_prices` rows that originate from a model F7 spread / total /
ML line have **no FD price to match** — there is no F7 market on FD beyond the
3-way result. The dashboard will continue to render "no FD pill" for F7 rows
unless we (a) start parsing `First 7 Innings Result` as a 2-way ML by dropping
the Tie leg, or (b) accept that FD F7 coverage is intentionally absent.

The prior probe showing zero F7 rows in `fd_odds/fd.duckdb` today is **not** a
whitelist bug — it's faithful reflection of FD's market offering. The
whitelist could be extended with `'First 7 Innings Result'` if we want a 2-way
ML pill (treating the Tie as devig); no other F7 strings exist to add.

## Alt-run-line capture

Verified on event `35600618` (`Kansas City Royals` @ `Chicago White Sox`):

- **FD API `Alternate Run Lines` (FG) returns 12 runners** spanning 6 distinct
  paired lines:
  - `Kansas City Royals +3.5` (-1100) / `Chicago White Sox -3.5` (+600)
  - `Kansas City Royals +2.5` (-620)  / `Chicago White Sox -2.5` (+400)
  - `Kansas City Royals +1.5` (-330)  / `Chicago White Sox -1.5` (+240)
  - `Kansas City Royals -1.5` (+120)  / `Chicago White Sox +1.5` (-154)
  - `Kansas City Royals -2.5` (+190)  / `Chicago White Sox +2.5` (-250)
  - `Kansas City Royals -3.5` (+285)  / `Chicago White Sox +3.5` (-400)

  i.e. each team is offered at both favored and underdog handicaps, yielding
  6 distinct *bets-pairs* (12 runners total).

- **`fd_odds/fd.duckdb::mlb_odds` persisted only 3 rows** for this game's FG
  alt_spreads (`home_spread` ∈ {1.5, 2.5, 3.5}; the symmetric -1.5/-2.5/-3.5
  rows for the other team are missing). Same pattern across all 13 events
  probed today.

- F5 alt_spreads shows the same shape — 4 rows where FD posts 8 runners.

(Skipped Step 5 R trace since the gap is already evident at the DuckDB layer
upstream of `get_fd_odds()`.)

## Bug location

`mlb_sgp/scraper_fanduel_singles.py::parse_runners_to_wide_rows`, line ~154:

```python
elif market_type == "alternate_spreads" and effective_line is not None:
    bucket_line = abs(effective_line)
```

The `abs()` on the bucket key collapses `KC -2.5 / CHW +2.5` (one bet pair)
and `KC +2.5 / CHW -2.5` (the OPPOSITE bet pair) into the same bucket
`(FG, alt_spreads, 2.5)`. Whichever pair the loop sees second silently
overwrites the first via the `row["home_spread"] = effective_line` / etc.
assignments.

The dashboard correspondingly sees only the favorite-favored direction of
each alt magnitude, and no FD pill is rendered when the model's bet is on
the dog-favored side (e.g. dashboard wants a CHW -2.5 line; DB only has
KC -2.5).

This is **not** a Tools.R collapse and **not** an FD payload gap — FD returns
the data; the scraper drops half of it during wide-row bucketing.

## Fix path

In `parse_runners_to_wide_rows`, change the alt_spreads bucket key from
`abs(effective_line)` to a *signed* representation that preserves which team
is favored. Simplest: bucket by `(period, "alternate_spreads", home_team_signed_line)`
where `home_team_signed_line` is `effective_line` when the runner is the
home team and `-effective_line` when the runner is the away team. The two
runners of a paired line (KC -2.5 + CHW +2.5) will share that key while the
opposite-direction pair (KC +2.5 + CHW -2.5) gets a distinct one.

Concretely: compute a `signed_home_line` before bucketing for alt_spreads
rows, e.g.

- if `home_team in r.name`: `signed_home_line = effective_line`
- elif `away_team in r.name`: `signed_home_line = -effective_line`

then `bucket_line = signed_home_line` (replacing `abs(...)`). The downstream
home/away assignment logic on lines 177–185 already places each side
correctly, so no further changes needed there.

This is the entire Task 4 scope; F7 (Task 3) reduces to a one-line whitelist
add of `'First 7 Innings Result'` ONLY IF we want a 2-way ML pill — otherwise
Task 3 is closed as "FD doesn't post F7 line markets, by design".
