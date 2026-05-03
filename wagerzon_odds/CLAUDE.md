# wagerzon_odds module

Wagerzon-specific scrapers, pricer, placer, and shared auth/account/balance
helpers used by the MLB correlated parlay dashboard.

## Quick map

- **Scraping odds:** `scraper_v2.py`, `scraper_specials.py` write to
  `wagerzon_odds/wagerzon.duckdb`.
- **Pricing parlays:** `parlay_pricer.py::get_parlay_price()` calls
  `ConfirmWagerHelper` with `RiskWin=2` (skips balance validation).
- **Placing parlays:** `parlay_placer.py::place_parlays(specs, account)`
  calls `PostWagerMultipleHelper`. **Requires** an explicit account.
- **Balance:** `wagerzon_balance.py::fetch_available_balance(account)`
  returns a `BalanceSnapshot`. The `available` field is
  `RealAvailBalance` (cash + credit) — the amount WZ will actually let
  you wager, not raw cash. `cash` (`CurrentBalance`) is exposed for
  tooltip / debug.
- **Account registry:** `wagerzon_accounts.py` discovers accounts by
  env-var suffix (`WAGERZON{X}_USERNAME` / `WAGERZON{X}_PASSWORD`).
  See `README.md` for the pattern.
- **Auth:** `wagerzon_auth.py::get_session(account)` caches one logged-in
  `requests.Session` per account label, in-memory only. Per-label locks
  so concurrent logins for different accounts run in parallel.
- **3-way F5 ML scraping:** `scraper_v2.py::parse_odds()` routes `idgmtyp=29` parent
  games (league `lg=1280` — "MLB - 1ST 5 INN WINNER (3-WAY)") through
  `parse_3way_line()`. Records have `market = "h2h_3way_1st_5_innings"`,
  `period = "f5"`, and the new `draw_ml INTEGER` column populated with
  the third (draw) outcome's American odds. `*_odds` tables include
  `draw_ml` via idempotent `ALTER TABLE ADD COLUMN IF NOT EXISTS`.

## Adding a new Wagerzon account

1. Pick a single uppercase suffix letter (e.g. `D`).
2. Add `WAGERZOND_USERNAME` and `WAGERZOND_PASSWORD` to `bet_logger/.env`.
3. Add the same pair to `bet_logger/.env.example` (placeholder values).
4. Restart any process that uses the registry (MLB dashboard server).
5. Optional: add `--account d` support to `bet_logger/scraper_wagerzon.py`
   if you want the bet-logging side to also pick up tickets from the
   account (separate concern; see `bet_logger/CLAUDE.md`).

## Pitfalls

- The auth flow scrapes `__VIEWSTATE` etc. from the Default.aspx form.
  If WZ changes the login page, both `wagerzon_auth.py` and any
  dependent test fixtures need updating. The login flow now raises a
  loud `RuntimeError` with the account label if `__VIEWSTATE` isn't
  found, so a markup change fails fast instead of silently posting
  empty fields.
- The balance endpoint URL is hardcoded in `wagerzon_balance.py` — if
  WZ moves it, update the constant + retest.
- The balance response body contains `result.Player` and
  `result.Password` (plaintext password) — never log the raw response.
  The parser extracts only the two numeric balance fields. The
  `WagerzonAccount` dataclass has `password` marked `repr=False` so
  it's excluded from any default representation.
- Pricing always uses the primary account regardless of which account
  will eventually place the bet. This is intentional (`RiskWin=2`) but
  means a primary-account login failure breaks pricing for all
  accounts.
