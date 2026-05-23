# TZ Audit Findings — 2026-05-22 (UTC 2026-05-23 02:04)

Run: `python tests/timezone_audit.py`

## Sign convention

Δ = `(odds_api_commence_time_UTC)` − `(scraper_recorded_dt_treated_as_UTC)`

So **positive Δ = the scraper's clock lags UTC by Δ hours**, which is the
expected sign for any TZ west of Greenwich.

| Δ      | Implied scraper TZ                |
|--------|-----------------------------------|
| 0h     | UTC                               |
| +4h    | America/New_York (EDT)            |
| +5h    | America/New_York (EST)            |
| +7h    | America/Los_Angeles (PDT)         |
| +8h    | America/Los_Angeles (PST)         |

NOTE: the source plan text predicted "+7h" for BKM with the opposite sign
convention (scraper − api). Our script's sign is **flipped from that note**
but the magnitudes are identical. Treat the magnitude as load-bearing.

## Modal Δ per scraper

| Scraper   | Modal Δ   | Inferred TZ                    | Matched games | Notes                                                   |
|-----------|-----------|--------------------------------|---------------|---------------------------------------------------------|
| wagerzon  | +4.02h    | America/New_York (EDT)         | 6 / 6         | Matches CLAIMED. (12 unmatched are pitcher-name props)  |
| hoop88    | +7.02h    | America/Los_Angeles (PDT)      | 6 / 6         | **BUG**: Tools.R CLAIMS ET — actual is PT.              |
| bookmaker | +7.02h    | America/Los_Angeles (PDT)      | 22 / 25       | **BUG confirmed**: CLAIMED ET, actually PT (live bug).  |
| bet105    | +0.02h    | UTC                            | 6 / 6         | Matches CLAIMED.                                        |

The +0.02h drift in modal Δ (and +7.02h for BKM/Hoop88, +4.02h for WZ) is
explained by the Odds API publishing first-pitch times with ~1min precision
(`23:06`, `23:11`, `23:41`) while scrapers round to 5-min grid (`19:40`,
`16:40`). Every modal Δ is exactly the TZ offset; the `.02h` ≈ 1min slop is
expected.

## Surprises

1. **Hoop88 is on PT, not ET**. The plan brief CLAIMED Hoop88 → America/New_York;
   the audit shows +7.02h (PDT) across all 6 matched games. This means **two**
   scrapers (BKM + Hoop88) are misinterpreted by `.parse_wz_game_dt` in
   Tools.R, not just BKM. Phase 3f (Hoop88) needs the same PT→UTC conversion
   logic that Phase 3g (BKM) needs.

2. **Wagerzon "unmatched" rows are pitcher-name props**, not real games — the
   `home_team` / `away_team` columns contain values like
   `"J DEGROM @ Texas Rangers"`. These shouldn't be in the same `mlb_odds`
   table as game-level rows; flag for Phase 3e investigation (likely a
   separate market_type the WZ parser is dumping into the wrong column or
   a side-effect of the alt-line scraper).

3. **Bookmaker 3 outliers** in the +7h modal:
   - `New York Mets @ Miami Marlins  05/22 16:10 → +28.02h`: scraper has a
     row dated 05/22 but the matching API row is 05/23 evening. Likely a
     stale row from a postponement/reschedule, OR BKM lists doubleheaders
     differently than Odds API.
   - `Cleveland Guardians @ Philadelphia Phillies 05/22 15:40 → +28.43h`:
     same pattern — stale 05/22 row matched to 05/23 API.
   - `St. Louis Cardinals @ Cincinnati Reds 05/23 16:15 → +0.93h`: BKM lists
     this as a 4:15pm PT game; Odds API says 5:11pm UTC. **A 56-minute
     real-clock discrepancy**, NOT a TZ issue. Probably a game-time update
     that BKM hasn't refreshed yet, or a doubleheader-game-2 ambiguity.

   None of these change the inferred TZ; they're consistent with stale data
   or doubleheader ambiguity in the matching heuristic, not a TZ signal.

## Implications for Phase 2 (`tools/migrate_scraper_schemas.py`)

Migration map for `commence_time_utc` derivation per scraper:

```python
SCRAPER_TZ = {
    "wagerzon":  "America/New_York",  # confirmed EDT
    "hoop88":    "America/Los_Angeles",  # confirmed PDT (was misclassified as ET)
    "bookmaker": "America/Los_Angeles",  # confirmed PDT (active bug)
    "bet105":    "UTC",                  # confirmed UTC
}
```

The two Group-A bugs to escalate:
- Hoop88 is the SECOND scraper (after BKM) that needs PT→UTC conversion.
- WZ "unmatched" rows containing pitcher names are a data-shape concern
  for Phase 3e but don't affect TZ inference.

## Raw output

```
Odds API returned 26 MLB rows across 15 distinct matchups
Now UTC: 2026-05-23T02:04:22.531195+00:00
Sample API entries:
  Tampa Bay Rays @ New York Yankees: ['2026-05-22T23:06:00+00:00', '2026-05-23T17:36:00+00:00']
  Pittsburgh Pirates @ Toronto Blue Jays: ['2026-05-22T23:08:00+00:00', '2026-05-23T19:08:00+00:00']
  Minnesota Twins @ Boston Red Sox: ['2026-05-22T23:11:00+00:00']

=== wagerzon (6 matched / 18 distinct rows in DB) ===
  12 DB rows had no Odds API match
    unmatched: J DEGROM @ Texas Rangers @ Los Angeles Angels  (05/22 21:38)
    unmatched: T SUGANO @ Colorado Rockies @ Arizona Diamondbacks  (05/22 21:40)
    unmatched: N CAMERON @ Seattle Mariners @ Kansas City Royals  (05/22 19:40)
    unmatched: J SPRINGS @ Athletics @ San Diego Padres  (05/22 21:40)
    unmatched: L HENDERSON @ Los Angeles Dodgers @ Milwaukee Brewers  (05/22 19:40)
  Colorado Rockies @ Arizona Diamondbacks  scraper=05/22 21:40  api=2026-05-23T01:41:00+00:00  Δ = +4.02h
  Texas Rangers @ Los Angeles Angels  scraper=05/22 21:38  api=2026-05-23T01:39:00+00:00  Δ = +4.02h
  Los Angeles Dodgers @ Milwaukee Brewers  scraper=05/22 19:40  api=2026-05-22T23:41:00+00:00  Δ = +4.02h
  Seattle Mariners @ Kansas City Royals  scraper=05/22 19:40  api=2026-05-22T23:41:00+00:00  Δ = +4.02h
  Athletics @ San Diego Padres  scraper=05/22 21:40  api=2026-05-23T01:41:00+00:00  Δ = +4.02h
  Chicago White Sox @ San Francisco Giants  scraper=05/22 22:15  api=2026-05-23T02:16:00+00:00  Δ = +4.02h
  → MODAL Δ: +4.0h  (6/6 games)  ⇒  America/New_York (EDT)

=== hoop88 (6 matched / 6 distinct rows in DB) ===
  Chicago White Sox @ San Francisco Giants  scraper=05/22 19:15  api=2026-05-23T02:16:00+00:00  Δ = +7.02h
  Los Angeles Dodgers @ Milwaukee Brewers  scraper=05/22 16:40  api=2026-05-22T23:41:00+00:00  Δ = +7.02h
  Seattle Mariners @ Kansas City Royals  scraper=05/22 16:40  api=2026-05-22T23:41:00+00:00  Δ = +7.02h
  Colorado Rockies @ Arizona Diamondbacks  scraper=05/22 18:40  api=2026-05-23T01:41:00+00:00  Δ = +7.02h
  Texas Rangers @ Los Angeles Angels  scraper=05/22 18:38  api=2026-05-23T01:39:00+00:00  Δ = +7.02h
  Athletics @ San Diego Padres  scraper=05/22 18:40  api=2026-05-23T01:41:00+00:00  Δ = +7.02h
  → MODAL Δ: +7.0h  (6/6 games)  ⇒  America/Los_Angeles (PDT)

=== bookmaker (25 matched / 25 distinct rows in DB) ===
  Washington Nationals @ Atlanta Braves  scraper=05/22 16:20  api=2026-05-22T23:15:41+00:00  Δ = +6.93h
  Colorado Rockies @ Arizona Diamondbacks  scraper=05/22 18:40  api=2026-05-23T01:41:00+00:00  Δ = +7.02h
  Cleveland Guardians @ Philadelphia Phillies  scraper=05/23 13:05  api=2026-05-23T20:06:00+00:00  Δ = +7.02h
  Athletics @ San Diego Padres  scraper=05/23 18:40  api=2026-05-24T01:41:00+00:00  Δ = +7.02h
  Chicago White Sox @ San Francisco Giants  scraper=05/22 19:15  api=2026-05-23T02:16:00+00:00  Δ = +7.02h
  Washington Nationals @ Atlanta Braves  scraper=05/23 13:10  api=2026-05-23T20:11:00+00:00  Δ = +7.02h
  Houston Astros @ Chicago Cubs  scraper=05/23 11:20  api=2026-05-23T18:21:00+00:00  Δ = +7.02h
  Pittsburgh Pirates @ Toronto Blue Jays  scraper=05/23 12:07  api=2026-05-23T19:08:00+00:00  Δ = +7.02h
  Chicago White Sox @ San Francisco Giants  scraper=05/23 13:05  api=2026-05-23T20:06:00+00:00  Δ = +7.02h
  Cleveland Guardians @ Philadelphia Phillies  scraper=05/22 15:40  api=2026-05-23T20:06:00+00:00  Δ = +28.43h
  Pittsburgh Pirates @ Toronto Blue Jays  scraper=05/22 16:07  api=2026-05-22T23:08:00+00:00  Δ = +7.02h
  Los Angeles Dodgers @ Milwaukee Brewers  scraper=05/22 16:40  api=2026-05-22T23:41:00+00:00  Δ = +7.02h
  Los Angeles Dodgers @ Milwaukee Brewers  scraper=05/23 16:15  api=2026-05-23T23:16:00+00:00  Δ = +7.02h
  Colorado Rockies @ Arizona Diamondbacks  scraper=05/23 19:10  api=2026-05-24T02:11:00+00:00  Δ = +7.02h
  Texas Rangers @ Los Angeles Angels  scraper=05/23 19:05  api=2026-05-24T02:06:00+00:00  Δ = +7.02h
  Detroit Tigers @ Baltimore Orioles  scraper=05/22 16:20  api=2026-05-22T23:16:00+00:00  Δ = +6.93h
  Seattle Mariners @ Kansas City Royals  scraper=05/22 16:40  api=2026-05-22T23:41:00+00:00  Δ = +7.02h
  St. Louis Cardinals @ Cincinnati Reds  scraper=05/23 16:15  api=2026-05-23T17:11:00+00:00  Δ = +0.93h
  New York Mets @ Miami Marlins  scraper=05/22 16:10  api=2026-05-23T20:11:00+00:00  Δ = +28.02h
  Tampa Bay Rays @ New York Yankees  scraper=05/22 16:05  api=2026-05-22T23:06:00+00:00  Δ = +7.02h
  Athletics @ San Diego Padres  scraper=05/22 18:40  api=2026-05-23T01:41:00+00:00  Δ = +7.02h
  Minnesota Twins @ Boston Red Sox  scraper=05/22 16:10  api=2026-05-22T23:11:00+00:00  Δ = +7.02h
  Texas Rangers @ Los Angeles Angels  scraper=05/22 18:38  api=2026-05-23T01:39:00+00:00  Δ = +7.02h
  Detroit Tigers @ Baltimore Orioles  scraper=05/23 13:05  api=2026-05-23T20:06:00+00:00  Δ = +7.02h
  Seattle Mariners @ Kansas City Royals  scraper=05/23 13:10  api=2026-05-23T20:11:00+00:00  Δ = +7.02h
  → MODAL Δ: +7.0h  (22/25 games)  ⇒  America/Los_Angeles (PDT)
    other rounded-to-0.5h Δ values: {7.0: 22, 28.5: 1, 1.0: 1, 28.0: 1}

=== bet105 (6 matched / 6 distinct rows in DB) ===
  Colorado Rockies @ Arizona Diamondbacks  scraper=05/23 01:40  api=2026-05-23T01:41:00+00:00  Δ = +0.02h
  Athletics @ San Diego Padres  scraper=05/23 01:40  api=2026-05-23T01:41:00+00:00  Δ = +0.02h
  Seattle Mariners @ Kansas City Royals  scraper=05/22 23:40  api=2026-05-22T23:41:00+00:00  Δ = +0.02h
  Chicago White Sox @ San Francisco Giants  scraper=05/23 02:15  api=2026-05-23T02:16:00+00:00  Δ = +0.02h
  Los Angeles Dodgers @ Milwaukee Brewers  scraper=05/22 23:40  api=2026-05-22T23:41:00+00:00  Δ = +0.02h
  Texas Rangers @ Los Angeles Angels  scraper=05/23 01:38  api=2026-05-23T01:39:00+00:00  Δ = +0.02h
  → MODAL Δ: +0.0h  (6/6 games)  ⇒  UTC
```
