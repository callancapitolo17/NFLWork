"""DuckDB persistence for the bets service.

Owns unabated_ticket/bets_service/bets.duckdb (this service is its only
writer). Tables:
  bets         one row per normalised record, UPSERT on `id` (the venue-native
               bet id) — every record ever seen, never pruned: the future CLV
               work needs the full history. `record` holds the JSON the panel
               receives; the other columns exist for the retention window and
               ad-hoc queries. A row whose content is unchanged is not
               rewritten (see upsert_bets; `content_hash` is the compare key),
               so `last_seen_at` and the record's own `sourceFetchedAt` are
               the last poll that CHANGED the record, not the last poll that
               saw it — no column records the latter.
  source_runs  APPEND, one row per source poll whether it succeeded or not —
               the per-source freshness the panel shows. Pruned to
               `source_runs_retention_days` (see _prune_source_runs_locked),
               always keeping each source's latest run and latest OK run.
  team_crosswalk  one row per (venue, league, venue team) -> Unabated team id,
               learned by the panel from a bet the board joined by its venue
               id (#118 step 4; the panel POSTs them, see service.py). INSERT
               only: a held key is never rewritten — the same id is a no-op,
               a different id is a conflict the caller is told about. Served
               with every /bets.json; DELETE /crosswalk.json empties it.
               The one exception is a manual attach (pin_bet): Cal is the
               authority, so its rows REPLACE a held key and carry
               `pinned_bet_id`; unpin_bet deletes exactly those rows.
  bet_pins     one row per bet Cal attached to a board event by hand (the
               panel's Attach control), UPSERT on `bet_id`. Never pruned: a
               few a week, and the history of manual matches.
  bet_fill_fairs  one row per bet: the fair (Unabated's `bacr`, a whole
               American price) the bet's own line had when it was placed, read
               by the panel off its in-memory line history
               (extension/fillfair.js) and POSTed to /fill_fairs.json. INSERT
               only, first capture wins (ON CONFLICT DO NOTHING): a saved fair
               is never rewritten, and nothing is ever backfilled. Never pruned
               (~14 bets a day). Served with /bets.json for the bets in its
               window (see load_fill_fairs).
  bet_fill_books  one row per bet: what the other books offered at its number
               the moment it was placed — how many books posted it, how many
               beat the price taken, and the best of them (server/fillbooks.js,
               read off the runner's line history; Cal 2026-10-13). POSTed to
               /fill_books.json by the server runner. INSERT only, first
               capture wins (ON CONFLICT DO NOTHING); nothing is backfilled.
               Stored only: nothing serves it back yet. Never pruned.
  bet_closing_fairs  one row per bet: Unabated's fair for the bet's own line on
               the newest snapshot before its game started — the close CLV is
               measured against. POSTed to /closing_fairs.json by the server
               runner (server/closefair.js) while the game is still to start;
               UPSERT keeps the row with the LATEST fair_observed_at, so a late
               or repeated POST never moves it backwards, and every row is
               observed before its event_start (service.py validates). Never
               pruned (one row per straight bet). Served with /bets.json for
               the bets in its window (see load_closing_fairs).
  bet_exclusions  one row per bet Cal removed from the Bet Tracker (its Bets
               page: bets at BFA or Wagerzon that are not his). The bet stays in
               `bets`; the tracker leaves an excluded bet, and a parlay any of
               whose legs is excluded, out of every number. Restore deletes the
               row. Never pruned.
  bet_dismissals  one row per open bet Cal dismissed from the unmatched list
               (the panel's or the phone's Dismiss: "stop flagging this bet").
               INSERT on Dismiss (a dismissed bet keeps its first time), DELETE
               on Restore. Served only while the bet is open — the flag ends
               when it settles — and never pruned (a few a week).
  teaser_blocks  one row per market (`<eventId>:bt<2|3>`, teaser.js
               marketKeyOf) Cal marked "Can't tease" (Buckeye keeps some CFB
               games off its teaser menu), UPSERT on the key, DELETE on
               Restore. Served only until the game starts; never pruned.
  edge_skips   one row per Edges line Cal skipped (the Skip chip on a row a
               book beats; pricecheck.js skipKeyOf: the line at every book,
               '<rung>'; rows written before 2026-10-16 end '@<bookId>'),
               UPSERT on the line, DELETE of the line (old per-book rows
               too) on Restore. Served only until the game starts; never
               pruned.
  edge_settings  at most ONE row (settings_id = 1): the Edges settings the
               server runner (unabated_ticket/server/runner.js, phone page plan
               step 1) reads every cycle, written by PUT /settings.json and
               UPSERTed whole on every write. Each column is an override; NULL
               means "the panel's default", and the defaults live only in
               extension/edgerows.js, never here. Never pruned (one row).
  manual_orders  one row per Kalshi or Novig order Cal placed from the Ticket's
               Place button (manual_orders.py), on Place's own key: INSERTED
               before it is sent (status 'sending'), UPDATED with what the
               venue answered and what filled. The client_order_id ('uman-' or
               'nman-' + the click's id) is the key, so one click can never send
               twice. Never pruned (a few a day); the bet itself reaches `bets`
               through the venue's poll.
A failed poll writes a source_runs row and touches nothing in `bets`, so a
dark source keeps serving its previous records.

The DuckDB connection is shared between the poll thread and the HTTP handler
threads behind one lock (DuckDB connections are not thread-safe).
"""
import hashlib
import json
import logging
import threading
from datetime import datetime, timedelta, timezone
from pathlib import Path

import duckdb

from unabated_ticket.bets_service.sources.kalshi import PRIMARY_SUBACCOUNT, subaccount_of_record_id

log = logging.getLogger(__name__)

_SCHEMA = """
CREATE TABLE IF NOT EXISTS bets (
    id             VARCHAR PRIMARY KEY,
    source         VARCHAR,
    venue          VARCHAR,
    league         VARCHAR,
    status         VARCHAR,
    placed_at      TIMESTAMPTZ,
    closed_at      TIMESTAMPTZ,
    event_date     VARCHAR,
    record         VARCHAR,      -- JSON, the normalised record as served
    first_seen_at  TIMESTAMPTZ,
    last_seen_at   TIMESTAMPTZ
);
-- #125: sha256 of the record content the UPSERT compares (see _content_hash).
-- Rows from before the column carry NULL and are rewritten once.
ALTER TABLE bets ADD COLUMN IF NOT EXISTS content_hash VARCHAR;
CREATE TABLE IF NOT EXISTS source_runs (
    source       VARCHAR,
    started_at   TIMESTAMPTZ,
    finished_at  TIMESTAMPTZ,
    ok           BOOLEAN,
    error        VARCHAR,
    n_records    INTEGER
);
CREATE TABLE IF NOT EXISTS team_crosswalk (
    venue               VARCHAR NOT NULL,
    league              VARCHAR NOT NULL,
    venue_team_key      VARCHAR NOT NULL,   -- Novig team id; Kalshi event-title name
    venue_team_name     VARCHAR,
    unabated_team_id    VARCHAR NOT NULL,   -- the id the board's lines carry; key = "<league>:<id>"
    unabated_team_name  VARCHAR,
    learned_from        VARCHAR,            -- "<bet id> on board event <event id>"
    learned_at          TIMESTAMPTZ NOT NULL,
    PRIMARY KEY (venue, league, venue_team_key)
);
-- The bet whose manual attach wrote the row; NULL on rows learned from an id join.
ALTER TABLE team_crosswalk ADD COLUMN IF NOT EXISTS pinned_bet_id VARCHAR;
CREATE TABLE IF NOT EXISTS bet_pins (
    bet_id          VARCHAR PRIMARY KEY,   -- bets.id, the venue-native bet id
    venue           VARCHAR NOT NULL,
    league          VARCHAR NOT NULL,
    event_id        VARCHAR NOT NULL,      -- Unabated's eventId, as a string like the team ids
    event_start     TIMESTAMPTZ,
    away_team_id    VARCHAR,
    home_team_id    VARCHAR,
    away_team_name  VARCHAR,
    home_team_name  VARCHAR,
    pinned_at       TIMESTAMPTZ NOT NULL
);
CREATE TABLE IF NOT EXISTS bet_fill_fairs (
    bet_id            VARCHAR PRIMARY KEY,   -- bets.id, the venue-native bet id
    line_key          VARCHAR NOT NULL,      -- the Unabated feed line the fair was read from
    points            DOUBLE,                -- that line's number; NULL on a moneyline
    fair_american     INTEGER NOT NULL,      -- Unabated's fair (bacr) for the bet's side at the fill
    fair_observed_at  TIMESTAMPTZ NOT NULL,  -- when the panel saw that fair: at or before placed_at
    placed_at         TIMESTAMPTZ NOT NULL,  -- the bet's placedAt
    captured_at       TIMESTAMPTZ NOT NULL   -- when this service stored the row
);
CREATE TABLE IF NOT EXISTS bet_fill_books (
    bet_id        VARCHAR PRIMARY KEY,   -- bets.id, the venue-native bet id
    placed_at     TIMESTAMPTZ NOT NULL,  -- the bet's placedAt
    fill_price    DOUBLE NOT NULL,       -- the bet's American price on Unabated's basis (Kalshi with its fee)
    book_count    INTEGER NOT NULL,      -- other books with a price at its number at the fill
    books_better  INTEGER NOT NULL,      -- of those, how many beat the price taken
    best_book     VARCHAR,               -- the best other book; NULL when book_count = 0
    best_price    DOUBLE,                -- its American price at the fill; NULL when book_count = 0
    captured_at   TIMESTAMPTZ NOT NULL   -- when this service stored the row
);
CREATE TABLE IF NOT EXISTS bet_closing_fairs (
    bet_id            VARCHAR PRIMARY KEY,   -- bets.id, the venue-native bet id
    line_key          VARCHAR NOT NULL,      -- the Unabated feed line the fair was read from
    points            DOUBLE,                -- that line's number; NULL on a moneyline
    fair_american     INTEGER NOT NULL,      -- Unabated's fair (bacr) for the bet's side
    fair_observed_at  TIMESTAMPTZ NOT NULL,  -- the league snapshot it was read off; before event_start
    event_start       TIMESTAMPTZ NOT NULL,  -- the board game's start
    updated_at        TIMESTAMPTZ NOT NULL   -- when this service last wrote the row
);
CREATE TABLE IF NOT EXISTS bet_exclusions (
    bet_id       VARCHAR PRIMARY KEY,   -- bets.id, the venue-native bet id
    venue        VARCHAR NOT NULL,
    excluded_at  TIMESTAMPTZ NOT NULL
);
CREATE TABLE IF NOT EXISTS bet_dismissals (
    bet_id        VARCHAR PRIMARY KEY,   -- bets.id, the venue-native bet id
    dismissed_at  TIMESTAMPTZ NOT NULL
);
CREATE TABLE IF NOT EXISTS teaser_blocks (
    market_key   VARCHAR PRIMARY KEY,    -- teaser.js marketKeyOf: '<eventId>:bt2' spread, ':bt3' total
    event_start  TIMESTAMPTZ NOT NULL,   -- the game's start: the block is served until then
    blocked_at   TIMESTAMPTZ NOT NULL
);
CREATE TABLE IF NOT EXISTS edge_skips (
    skip_key     VARCHAR PRIMARY KEY,    -- pricecheck.js skipKeyOf: '<eventId>:<period>:<betType>:<side>|<points|ml>' (older rows: '...@<bookId>')
    event_start  TIMESTAMPTZ NOT NULL,   -- the game's start: the skip is served until then
    skipped_at   TIMESTAMPTZ NOT NULL
);
CREATE TABLE IF NOT EXISTS edge_settings (
    settings_id           INTEGER PRIMARY KEY CHECK (settings_id = 1),  -- one row
    bankroll              DOUBLE,      -- dollars; NULL = the default (edgerows.js)
    kelly_multiplier      DOUBLE,      -- 0 < m <= 1
    league_ids            INTEGER[],   -- Unabated league ids (feed.LEAGUES)
    period_type_ids       INTEGER[],   -- feed.PERIODS ids, 1 = full game
    bet_type_ids          INTEGER[],   -- 1 moneyline, 2 spread, 3 total
    book_mode             VARCHAR,     -- 'default' (the default books) | 'all' (every live book) | 'custom'
    book_ids              INTEGER[],   -- the books ticked; set exactly when book_mode = 'custom'
    min_edge_pct          DOUBLE,
    min_stake             DOUBLE,      -- Min suggested bet, dollars; 0 = off
    max_line_age_hours    DOUBLE,
    min_liquidity_to_win  DOUBLE,      -- dollars; 0 = off
    include_alts          BOOLEAN,
    sort_by               VARCHAR,     -- 'edge' | 'stake' | 'start' | 'exposure'
    group_by_market       BOOLEAN,
    updated_at            TIMESTAMPTZ NOT NULL
);
CREATE TABLE IF NOT EXISTS manual_orders (
    client_order_id  VARCHAR PRIMARY KEY,   -- 'uman-' + the click's id
    placed_at        TIMESTAMPTZ NOT NULL,  -- when the row was inserted, just before the send
    event_id         BIGINT NOT NULL,       -- Unabated's eventId
    ticker           VARCHAR NOT NULL,
    contract_side    VARCHAR NOT NULL,      -- 'yes' | 'no'
    label            VARCHAR,               -- the line in words, for the log
    ticket_all_in    DOUBLE NOT NULL,       -- the Ticket's price, all-in, as a probability
    ask_price        DOUBLE NOT NULL,       -- the limit sent (the live ask)
    all_in_price     DOUBLE NOT NULL,       -- that ask with Kalshi's taker fee
    stake            DOUBLE NOT NULL,       -- dollars Cal asked to bet
    count_requested  INTEGER NOT NULL,
    status           VARCHAR NOT NULL,      -- sending | filled | partial | unfilled | rejected | unknown
    kalshi_order_id  VARCHAR,
    count_filled     DOUBLE,
    cost_filled      DOUBLE,                -- dollars paid, fees included
    error            VARCHAR,
    updated_at       TIMESTAMPTZ NOT NULL
);
-- Novig Place: rows from before it are all Kalshi. On a Novig row `ticker` is
-- the Novig outcome id, `contract_side` 'buy' (every Novig order buys its
-- outcome), `ask_price` = `all_in_price` (a pregame take pays no fee) and a
-- contract pays 1c, not $1.
ALTER TABLE manual_orders ADD COLUMN IF NOT EXISTS venue VARCHAR DEFAULT 'kalshi';
ALTER TABLE manual_orders ADD COLUMN IF NOT EXISTS novig_order_id VARCHAR;
-- The board line a Novig Place order bought (taker.js novigLineOf, as the bot
-- desk keeps it in taker.duckdb::novig_taker_orders): Novig's catalog names
-- neither the board's teams nor its periods, so Place's sub-account read
-- (sources/novig_subaccount.py) names its bets with these. NULL on Kalshi rows.
ALTER TABLE manual_orders ADD COLUMN IF NOT EXISTS market_id VARCHAR;
ALTER TABLE manual_orders ADD COLUMN IF NOT EXISTS event_start TIMESTAMPTZ;
ALTER TABLE manual_orders ADD COLUMN IF NOT EXISTS league VARCHAR;
ALTER TABLE manual_orders ADD COLUMN IF NOT EXISTS bet_type VARCHAR;      -- spread | total
ALTER TABLE manual_orders ADD COLUMN IF NOT EXISTS period VARCHAR;
ALTER TABLE manual_orders ADD COLUMN IF NOT EXISTS side VARCHAR;          -- away | home | over | under
ALTER TABLE manual_orders ADD COLUMN IF NOT EXISTS points DOUBLE;
ALTER TABLE manual_orders ADD COLUMN IF NOT EXISTS away_team VARCHAR;
ALTER TABLE manual_orders ADD COLUMN IF NOT EXISTS home_team VARCHAR;
ALTER TABLE manual_orders ADD COLUMN IF NOT EXISTS away_team_id VARCHAR;
ALTER TABLE manual_orders ADD COLUMN IF NOT EXISTS home_team_id VARCHAR;
-- When Place's sub-account read first saw the order's market SETTLED: the
-- bet's fixed closedAt, kept so a restart never moves it.
ALTER TABLE manual_orders ADD COLUMN IF NOT EXISTS settled_at TIMESTAMPTZ;
"""

MANUAL_ORDER_INSERT_COLUMNS = (
    "client_order_id", "placed_at", "venue", "event_id", "ticker", "contract_side", "label", "ticket_all_in",
    "ask_price", "all_in_price", "stake", "count_requested", "status", "updated_at",
)
# A Novig row's board line; a Kalshi row leaves them NULL.
MANUAL_ORDER_NOVIG_LINE_COLUMNS = (
    "market_id", "event_start", "league", "bet_type", "period", "side", "points", "away_team", "home_team",
    "away_team_id", "home_team_id",
)
MANUAL_ORDER_UPDATE_COLUMNS = ("status", "kalshi_order_id", "novig_order_id", "count_filled", "cost_filled", "error")
_INSERT_MANUAL_ORDER_COLUMNS = MANUAL_ORDER_INSERT_COLUMNS + MANUAL_ORDER_NOVIG_LINE_COLUMNS
_INSERT_MANUAL_ORDER = (f"INSERT INTO manual_orders ({', '.join(_INSERT_MANUAL_ORDER_COLUMNS)}) "
                        f"VALUES ({', '.join(['?'] * len(_INSERT_MANUAL_ORDER_COLUMNS))})")
# Every Novig Place order Novig gave an id, as the sub-account read's
# contexts (sources/novig_subaccount.py, taker_store.NOVIG_CONTEXT_COLUMNS).
_SELECT_NOVIG_PLACE_CONTEXTS = """
SELECT novig_order_id, league, bet_type, period, side, points, away_team, home_team, away_team_id, home_team_id,
       epoch_ms(event_start), event_id, market_id, ticker, ask_price, epoch_ms(settled_at)
FROM manual_orders
WHERE venue = 'novig' AND novig_order_id IS NOT NULL
"""
# A Place order on the game that may hold contracts the stake's read missed:
# sent and filled, partly filled, or no definite answer, since the read. A
# Novig 'unfilled' counts too (its fills are read once, right after the order
# goes final, and Novig's fills feed can lag that read). A Novig order Novig
# never answered and never named (no order id) counts whatever the read time:
# its fill can reach no read until it is named (`lost_novig_place_orders`).
_SELECT_MANUAL_MAYBE_FILLED = """
SELECT client_order_id, ticker, status, epoch(placed_at)
FROM manual_orders
WHERE event_id = ?
  AND (? IS NULL OR coalesce(venue, 'kalshi') = ?)
  AND ((placed_at >= ? AND (status IN ('sending', 'unknown', 'filled', 'partial')
                            OR (coalesce(venue, 'kalshi') = 'novig' AND status = 'unfilled')))
       OR (coalesce(venue, 'kalshi') = 'novig' AND status = 'unknown' AND novig_order_id IS NULL))
ORDER BY placed_at
"""
# Novig Place orders Novig never answered (no order id): looked up by clientId.
_SELECT_LOST_NOVIG_PLACE = """
SELECT client_order_id, ticker
FROM manual_orders
WHERE venue = 'novig' AND status = 'unknown' AND novig_order_id IS NULL AND placed_at <= ?
ORDER BY placed_at
"""
# Novig Place orders Novig named but never showed final to the desk: settled by
# Place's sub-account read once Novig shows them final.
_SELECT_UNCONFIRMED_NOVIG_PLACE = """
SELECT client_order_id, novig_order_id, count_requested
FROM manual_orders
WHERE venue = 'novig' AND status = 'unknown' AND novig_order_id IS NOT NULL AND placed_at <= ? AND placed_at >= ?
ORDER BY placed_at
"""
# An unconfirmed order is re-read for this long, then left 'unknown' (a row
# Novig never shows final must not slow every novig_place poll for good).
UNCONFIRMED_NOVIG_PLACE_WINDOW = timedelta(hours=24)

# An open straight Novig bet (app or sub-account) on one outcome.
_SELECT_OPEN_NOVIG_ON_OUTCOME = """
SELECT 1 FROM bets
WHERE venue = 'novig' AND status = 'open'
  AND json_extract_string(record, '$.venueIds.outcomeId') = ?
  AND coalesce(json_extract(record, '$.isParlayLeg')::BOOLEAN, false) = false
LIMIT 1
"""

_INSERT_CROSSWALK = """
INSERT INTO team_crosswalk (venue, league, venue_team_key, venue_team_name, unabated_team_id,
                            unabated_team_name, learned_from, learned_at)
VALUES (?, ?, ?, ?, ?, ?, ?, ?)
"""

_SELECT_CROSSWALK_KEYS = "SELECT venue, league, venue_team_key, unabated_team_id FROM team_crosswalk"

_SELECT_CROSSWALK = """
SELECT venue, league, venue_team_key, venue_team_name, unabated_team_id, unabated_team_name,
       learned_from, epoch(learned_at), pinned_bet_id
FROM team_crosswalk
ORDER BY learned_at DESC, venue, league, venue_team_key
"""

# A manual attach replaces whatever is held for the key (see the module docstring).
_REPLACE_CROSSWALK = """
INSERT INTO team_crosswalk (venue, league, venue_team_key, venue_team_name, unabated_team_id,
                            unabated_team_name, learned_from, learned_at, pinned_bet_id)
VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
ON CONFLICT (venue, league, venue_team_key) DO UPDATE SET
    venue_team_name = excluded.venue_team_name, unabated_team_id = excluded.unabated_team_id,
    unabated_team_name = excluded.unabated_team_name, learned_from = excluded.learned_from,
    learned_at = excluded.learned_at, pinned_bet_id = excluded.pinned_bet_id
"""

_UPSERT_PIN = """
INSERT INTO bet_pins (bet_id, venue, league, event_id, event_start, away_team_id, home_team_id,
                      away_team_name, home_team_name, pinned_at)
VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
ON CONFLICT (bet_id) DO UPDATE SET
    venue = excluded.venue, league = excluded.league, event_id = excluded.event_id,
    event_start = excluded.event_start, away_team_id = excluded.away_team_id,
    home_team_id = excluded.home_team_id, away_team_name = excluded.away_team_name,
    home_team_name = excluded.home_team_name, pinned_at = excluded.pinned_at
"""

_SELECT_PINS = """
SELECT bet_id, venue, league, event_id, epoch(event_start), away_team_id, home_team_id,
       away_team_name, home_team_name, epoch(pinned_at)
FROM bet_pins
ORDER BY pinned_at DESC, bet_id
"""

# One statement per UPSERT_ROWS_PER_STATEMENT rows ({values} is that many
# 12-placeholder tuples): DuckDB runs executemany one statement per row, which
# measured 440 ms for a 528-row poll against 20-28 ms for multi-row VALUES.
_UPSERT_BETS = """
INSERT INTO bets (id, source, venue, league, status, placed_at, closed_at, event_date,
                  record, content_hash, first_seen_at, last_seen_at)
VALUES {values}
ON CONFLICT (id) DO UPDATE SET
    source = excluded.source, venue = excluded.venue, league = excluded.league,
    status = excluded.status, placed_at = excluded.placed_at, closed_at = excluded.closed_at,
    event_date = excluded.event_date, record = excluded.record, content_hash = excluded.content_hash,
    last_seen_at = excluded.last_seen_at
WHERE bets.content_hash IS DISTINCT FROM excluded.content_hash
"""

_UPSERT_ROW_PLACEHOLDERS = "(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"
UPSERT_ROWS_PER_STATEMENT = 200

# Insert-only, first capture wins: a bet that has a saved fair keeps it. The
# same multi-row VALUES shape as _UPSERT_BETS ({values} is that many
# 7-placeholder tuples, at most UPSERT_ROWS_PER_STATEMENT).
_INSERT_FILL_FAIRS = """
INSERT INTO bet_fill_fairs (bet_id, line_key, points, fair_american, fair_observed_at, placed_at, captured_at)
VALUES {values}
ON CONFLICT (bet_id) DO NOTHING
"""

_FILL_FAIR_ROW_PLACEHOLDERS = "(?, ?, ?, ?, ?, ?, ?)"

# Insert-only, first capture wins, as bet_fill_fairs.
_INSERT_FILL_BOOKS = """
INSERT INTO bet_fill_books (bet_id, placed_at, fill_price, book_count, books_better, best_book, best_price, captured_at)
VALUES (?, ?, ?, ?, ?, ?, ?, ?)
ON CONFLICT (bet_id) DO NOTHING
"""

# The saved fairs of the bets /bets.json serves — the same window as
# _SELECT_WINDOW. The table is never pruned, so serving all of it would grow
# every poll for good; the panel only reads the fairs of the bets it holds.
_SELECT_FILL_FAIRS = """
SELECT f.bet_id, f.line_key, f.points, f.fair_american,
       epoch(f.fair_observed_at), epoch(f.placed_at), epoch(f.captured_at)
FROM bet_fill_fairs f
JOIN bets b ON b.id = f.bet_id
WHERE b.status = 'open' OR b.closed_at IS NULL OR b.closed_at >= ?
ORDER BY f.placed_at, f.bet_id
"""

# Newest observation wins: a row only moves forward in time, so a POST that
# lands late (or twice) never replaces a later reading of the same bet.
_UPSERT_CLOSING_FAIR = """
INSERT INTO bet_closing_fairs (bet_id, line_key, points, fair_american, fair_observed_at, event_start, updated_at)
VALUES (?, ?, ?, ?, ?, ?, ?)
ON CONFLICT (bet_id) DO UPDATE SET
    line_key = excluded.line_key, points = excluded.points, fair_american = excluded.fair_american,
    fair_observed_at = excluded.fair_observed_at, event_start = excluded.event_start,
    updated_at = excluded.updated_at
WHERE excluded.fair_observed_at > bet_closing_fairs.fair_observed_at
"""

# The same window as _SELECT_FILL_FAIRS.
_SELECT_CLOSING_FAIRS = """
SELECT c.bet_id, c.line_key, c.points, c.fair_american, epoch(c.fair_observed_at), epoch(c.event_start)
FROM bet_closing_fairs c
JOIN bets b ON b.id = c.bet_id
WHERE b.status = 'open' OR b.closed_at IS NULL OR b.closed_at >= ?
ORDER BY c.event_start, c.bet_id
"""

_SELECT_EXCLUSIONS = "SELECT bet_id, venue, epoch(excluded_at) FROM bet_exclusions ORDER BY excluded_at DESC, bet_id"

# A dismissal ends when its bet settles, so only open bets' rows are served.
_SELECT_DISMISSALS = """
SELECT d.bet_id, epoch(d.dismissed_at)
FROM bet_dismissals d
JOIN bets b ON b.id = d.bet_id
WHERE b.status = 'open'
ORDER BY d.dismissed_at DESC, d.bet_id
"""

_UPSERT_TEASER_BLOCK = """
INSERT INTO teaser_blocks (market_key, event_start, blocked_at) VALUES (?, ?, ?)
ON CONFLICT (market_key) DO UPDATE SET event_start = excluded.event_start, blocked_at = excluded.blocked_at
"""

# A block is moot once its game has started (the panel's liveTeaserBlocks rule).
_SELECT_TEASER_BLOCKS = """
SELECT market_key, epoch_ms(event_start), epoch(blocked_at)
FROM teaser_blocks
WHERE event_start > ?
ORDER BY event_start, market_key
"""

_UPSERT_EDGE_SKIP = """
INSERT INTO edge_skips (skip_key, event_start, skipped_at) VALUES (?, ?, ?)
ON CONFLICT (skip_key) DO UPDATE SET event_start = excluded.event_start, skipped_at = excluded.skipped_at
"""

def edge_skip_line(skip_key: str) -> str:
    """The line a skip key names (pricecheck.js skipLineOf): a key from a
    panel not yet reloaded carries a '@<bookId>' suffix, dropped here."""
    return skip_key.split("@", 1)[0]


# A skip is moot once its game has started, like a teaser block.
_SELECT_EDGE_SKIPS = """
SELECT skip_key, epoch_ms(event_start), epoch(skipped_at)
FROM edge_skips
WHERE event_start > ?
ORDER BY event_start, skip_key
"""

# The settings fields in the API's camelCase, in column order: the shape
# GET/PUT /settings.json speaks and the order of every statement below.
EDGE_SETTINGS_FIELDS = (
    "bankroll", "multiplier", "leagues", "periods", "betTypes", "bookMode", "bookIds", "minEdgePct",
    "minStake", "maxLineAgeHours", "minLiquidityToWin", "includeAlts", "sortBy", "groupByMarket",
)

_UPSERT_EDGE_SETTINGS = """
INSERT INTO edge_settings (settings_id, bankroll, kelly_multiplier, league_ids, period_type_ids, bet_type_ids,
                           book_mode, book_ids, min_edge_pct, min_stake, max_line_age_hours,
                           min_liquidity_to_win, include_alts, sort_by, group_by_market, updated_at)
VALUES (1, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
ON CONFLICT (settings_id) DO UPDATE SET
    bankroll = excluded.bankroll, kelly_multiplier = excluded.kelly_multiplier,
    league_ids = excluded.league_ids, period_type_ids = excluded.period_type_ids,
    bet_type_ids = excluded.bet_type_ids, book_mode = excluded.book_mode, book_ids = excluded.book_ids,
    min_edge_pct = excluded.min_edge_pct, min_stake = excluded.min_stake,
    max_line_age_hours = excluded.max_line_age_hours, min_liquidity_to_win = excluded.min_liquidity_to_win,
    include_alts = excluded.include_alts, sort_by = excluded.sort_by,
    group_by_market = excluded.group_by_market, updated_at = excluded.updated_at
"""

_SELECT_EDGE_SETTINGS = """
SELECT bankroll, kelly_multiplier, league_ids, period_type_ids, bet_type_ids, book_mode, book_ids,
       min_edge_pct, min_stake, max_line_age_hours, min_liquidity_to_win, include_alts, sort_by,
       group_by_market, epoch(updated_at)
FROM edge_settings
WHERE settings_id = 1
"""

# Keeps every source's latest run AND latest successful run whatever their
# age: source_status() reads exactly those two, and a source failing for
# longer than the window must keep showing its last good read, not blank out.
_PRUNE_SOURCE_RUNS = """
DELETE FROM source_runs
WHERE started_at < ?
  AND rowid NOT IN (SELECT rowid FROM source_runs
                    QUALIFY row_number() OVER (PARTITION BY source, ok ORDER BY started_at DESC) = 1)
"""

# Open bets always; settled/closed ones within the window; a non-open bet with
# no closedAt is kept (dropping it would be silent) — same rule as bets.js.
_SELECT_KALSHI_IDS = """
SELECT id FROM bets WHERE venue = 'kalshi'
"""

# Cal's Novig app feed (sources/novig.py, ids "novig:<order id>") also lists
# the orders of his API sub-accounts, which their own reads hold as
# "novig:<order id>:place" / ":bot" (sources/novig_subaccount.py): the app's
# copy of such an order is left out, so every order is served (sized,
# tracked) once, under the sub-account read that knows its board line. A
# sub-account copy with no board line (unmatchable: an order no desk names,
# or a lost Place send not yet adopted) sizes nothing, so then the app's copy
# is served instead of it. Either way one copy.
_SELECT_WINDOW = """
SELECT b.record FROM bets AS b
WHERE (b.status = 'open' OR b.closed_at IS NULL OR b.closed_at >= ?)
  AND NOT (b.venue = 'novig' AND EXISTS (
      SELECT 1 FROM bets AS sub
      WHERE sub.id IN (b.id || ':place', b.id || ':bot')
        AND json_extract_string(sub.record, '$.unmatchable') IS NULL))
  AND NOT (b.venue = 'novig' AND regexp_matches(b.id, ':(place|bot)$')
           AND json_extract_string(b.record, '$.unmatchable') IS NOT NULL
           AND EXISTS (SELECT 1 FROM bets AS app WHERE app.id = regexp_replace(b.id, ':(place|bot)$', '')))
ORDER BY b.placed_at DESC, b.id
"""

# Timestamps come back as epoch seconds: reading a TIMESTAMPTZ into Python
# needs pytz, which the service's venv does not carry.
_SELECT_LATEST_RUNS = """
SELECT source, epoch(started_at), ok, error, n_records
FROM source_runs
QUALIFY row_number() OVER (PARTITION BY source ORDER BY started_at DESC) = 1
"""

_SELECT_LATEST_OK_RUNS = """
SELECT source, epoch(started_at), n_records
FROM source_runs
WHERE ok
QUALIFY row_number() OVER (PARTITION BY source ORDER BY started_at DESC) = 1
"""


# `sourceFetchedAt` is the poll's own clock — every source restamps every
# record with it on every poll (sources/kalshi.py:182) — so an unchanged-row
# check that compared it would never skip anything. Excluding it means the
# STORED stamp is the last poll that changed the record. That is the right
# value to serve: the panel breaks merge ties on it (bets.js dedupeByNativeId,
# newest wins), and "content current as of its last change" is the honest
# claim — restamping it with the latest poll would also mark records the
# source no longer returns (Novig/BetOnline pull a rolling window; `bets`
# serves open rows forever) as fresh, and a stale copy would then beat a
# newer content-script read for good.
VOLATILE_RECORD_FIELDS = ("sourceFetchedAt",)

# source_runs is pruned at most this often: source_status() only ever reads the
# latest run per source, so a DELETE on every poll would scan the table under
# the store lock for nothing.
SOURCE_RUNS_PRUNE_INTERVAL_SEC = 3600


def _content_hash(record: dict) -> str:
    """sha256 of the record with the volatile fields removed, keys sorted —
    the UPSERT's compare key: an exact text compare (no Python `==`, where
    1 == 1.0 and True == 1). The hash covers the record only — a change to
    how the derived columns (status, closed_at, ...) are extracted does not
    rewrite existing rows, and a build from before the column leaves the
    hash stale while rewriting `record`; after either, force a one-time
    rewrite with `UPDATE bets SET content_hash = NULL`."""
    content = {key: value for key, value in record.items() if key not in VOLATILE_RECORD_FIELDS}
    return hashlib.sha256(json.dumps(content, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _last_record_per_id(records: list[dict]) -> list[dict]:
    """The records with any duplicate id collapsed to its LAST occurrence.
    No source emits one today; if one did, a multi-row INSERT ... ON CONFLICT
    silently keeps the FIRST row for a key repeated within the statement
    (verified on DuckDB 1.4.4) — the old per-row executemany applied them in
    order, so last-wins is the behaviour kept, and the duplicate is logged."""
    by_id: dict[str, dict] = {}
    duplicates: set[str] = set()
    for record in records:
        if record["id"] in by_id:
            duplicates.add(record["id"])
        by_id[record["id"]] = record
    if duplicates:
        log.warning("upsert_bets: %d record id(s) repeated in one poll, last occurrence kept: %s",
                    len(duplicates), sorted(duplicates))
    return list(by_id.values())


def _first_row_per_bet(rows: list[dict]) -> list[dict]:
    """The fill-fair rows with any repeated betId collapsed to its FIRST
    occurrence: the first capture wins inside one request as it does against
    the table. A multi-row INSERT would silently keep the first too (DuckDB
    1.4.4); collapsing here makes that the stated rule and logs the repeat."""
    by_bet: dict[str, dict] = {}
    duplicates: set[str] = set()
    for row in rows:
        if row["betId"] in by_bet:
            duplicates.add(row["betId"])
            continue
        by_bet[row["betId"]] = row
    if duplicates:
        log.warning("save_fill_fairs: %d bet id(s) repeated in one request, first row kept: %s",
                    len(duplicates), sorted(duplicates))
    return list(by_bet.values())


def _parse_iso(value: object) -> datetime | None:
    if not isinstance(value, str):
        return None
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        return None
    return parsed if parsed.tzinfo else parsed.replace(tzinfo=timezone.utc)


def _epoch_to_iso(epoch_seconds: float | None) -> str | None:
    if epoch_seconds is None:
        return None
    return datetime.fromtimestamp(epoch_seconds, timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _utc_of_epoch_ms(epoch_ms: int | None) -> datetime | None:
    return None if epoch_ms is None else datetime.fromtimestamp(epoch_ms / 1000, timezone.utc)


class BetsStore:
    def __init__(self, db_path: str | Path, source_runs_retention_days: int):
        """`source_runs_retention_days` 0 disables the prune (the repo's
        convention for switching a guard off); negative is a config error."""
        if source_runs_retention_days < 0:
            raise ValueError(
                f"source_runs_retention_days must be 0 (no prune) or more, got {source_runs_retention_days}")
        self._lock = threading.Lock()
        self._source_runs_retention_days = source_runs_retention_days
        self._last_source_runs_prune_at: datetime | None = None
        self._con = duckdb.connect(str(db_path))
        self._con.execute(_SCHEMA)

    def close(self) -> None:
        with self._lock:
            self._con.close()

    def upsert_bets(self, records: list[dict], seen_at: datetime) -> int:
        """UPSERT the records; a row whose content_hash is unchanged is left
        alone by the UPSERT's WHERE (measured: zero WAL bytes over 20
        identical polls of 500 rows, against 1.3 MB unconditionally). Returns
        how many rows were inserted or updated. A source re-sends every record
        it has ever seen on every poll (Kalshi: ~111 records every 60 s, of
        which ~109 are unchanged), and a DuckDB update is copy-on-write —
        rewriting an identical row appends a new version and turns the old
        one into garbage, which is what grew the WAL to 9.8 MB against a
        4.7 MB file holding 500 rows (#125). `last_seen_at` is written only
        on those rows — it is the last poll that CHANGED the record."""
        rows = [(
            record["id"], record.get("source"), record.get("venue"), record.get("league"),
            record.get("status"), _parse_iso(record.get("placedAt")), _parse_iso(record.get("closedAt")),
            record.get("eventDate"), json.dumps(record, separators=(",", ":")), _content_hash(record),
            seen_at, seen_at,
        ) for record in _last_record_per_id(records)]
        written = 0
        with self._lock:
            for start in range(0, len(rows), UPSERT_ROWS_PER_STATEMENT):
                chunk = rows[start:start + UPSERT_ROWS_PER_STATEMENT]
                statement = _UPSERT_BETS.format(values=", ".join([_UPSERT_ROW_PLACEHOLDERS] * len(chunk)))
                # DuckDB's INSERT result is the rows it inserted or updated —
                # the ones the WHERE let through (verified on 1.4.4).
                [written_in_chunk] = self._con.execute(statement, [value for row in chunk for value in row]).fetchone()
                written += written_in_chunk
        return written

    def log_source_run(self, source: str, started_at: datetime, finished_at: datetime,
                       ok: bool, error: str | None, n_records: int) -> None:
        with self._lock:
            self._con.execute(
                "INSERT INTO source_runs (source, started_at, finished_at, ok, error, n_records) "
                "VALUES (?, ?, ?, ?, ?, ?)",
                [source, started_at, finished_at, ok, error, n_records])
            self._prune_source_runs_locked(finished_at)

    def _prune_source_runs_locked(self, now: datetime) -> int:
        """DELETE source_runs rows older than the retention window — never a
        source's latest run of either outcome, so its latest run and latest
        OK run always survive — at most once per
        SOURCE_RUNS_PRUNE_INTERVAL_SEC; returns how many went. The store lock
        is already taken. Nothing reads a run but source_status(), which
        full-scans this table twice under the lock on every /bets.json and
        /health, so an unbounded table (~2,000 rows/day measured) slows every
        poll and every panel refresh. Never raises: the run that triggered it
        is already logged, and an error here must not re-enter
        run_source_once's failure path and log that poll a second time as
        failed. The throttle advances first, so a persistent error is
        reported hourly, not every poll."""
        if self._source_runs_retention_days == 0:
            return 0
        if (self._last_source_runs_prune_at is not None
                and (now - self._last_source_runs_prune_at).total_seconds() < SOURCE_RUNS_PRUNE_INTERVAL_SEC):
            return 0
        self._last_source_runs_prune_at = now
        cutoff = now - timedelta(days=self._source_runs_retention_days)
        try:
            [deleted] = self._con.execute(_PRUNE_SOURCE_RUNS, [cutoff]).fetchone()
        except Exception:  # noqa: BLE001 — a failed prune is logged, never a failed poll
            log.exception("source_runs: prune of rows older than %s failed; retrying in %ds",
                          cutoff.isoformat(), SOURCE_RUNS_PRUNE_INTERVAL_SEC)
            return 0
        if deleted:
            log.info("source_runs: pruned %d row(s) older than %s", deleted, cutoff.isoformat())
        return deleted

    def load_kalshi_tickers(self, subaccount: int = PRIMARY_SUBACCOUNT) -> set[str]:
        """Every Kalshi market ticker a stored record of `subaccount` names (ids
        are kalshi:<ticker>:<side>, plus :s<N> in numbered subaccount N), so the
        Kalshi source reads that subaccount's pre-cutoff fills of them."""
        with self._lock:
            rows = self._con.execute(_SELECT_KALSHI_IDS).fetchall()
        tickers: set[str] = set()
        for (record_id,) in rows:
            try:
                if subaccount_of_record_id(record_id) == subaccount:
                    tickers.add(record_id.split(":")[1])
            except ValueError as error:
                log.warning("load_kalshi_tickers: skipping stored id %r: %s", record_id, error)
        return tickers

    def load_bets(self, days: int, now: datetime) -> list[dict]:
        """Records in the retention window: open, or closed within `days` of
        `now`; an order the Novig app feed and a sub-account read both hold,
        once (_SELECT_WINDOW)."""
        cutoff = now - timedelta(days=days)
        with self._lock:
            rows = self._con.execute(_SELECT_WINDOW, [cutoff]).fetchall()
        return [json.loads(row[0]) for row in rows]

    def learn_crosswalk(self, rows: list[dict], learned_at: datetime) -> dict:
        """INSERT the crosswalk rows not yet held (validated shape: see
        service.validate_crosswalk_rows). A held key is never rewritten: the
        same Unabated id is a no-op, a different one is returned in
        `conflicts` with what is held. Returns {learned, conflicts}."""
        learned = 0
        conflicts: list[dict] = []
        with self._lock:
            held = {(venue, league, key): team_id
                    for venue, league, key, team_id in self._con.execute(_SELECT_CROSSWALK_KEYS).fetchall()}
            for row in rows:
                key = (row["venue"], row["league"], row["venueTeamKey"])
                if key in held:
                    if held[key] != row["unabatedTeamId"]:
                        conflicts.append({"venue": row["venue"], "league": row["league"],
                                          "venueTeamKey": row["venueTeamKey"],
                                          "held": held[key], "proposed": row["unabatedTeamId"]})
                    continue
                self._con.execute(_INSERT_CROSSWALK, [
                    row["venue"], row["league"], row["venueTeamKey"], row.get("venueTeamName"),
                    row["unabatedTeamId"], row.get("unabatedTeamName"), row.get("learnedFrom"), learned_at])
                held[key] = row["unabatedTeamId"]
                learned += 1
        if conflicts:
            log.warning("crosswalk: %d row(s) refused, a different Unabated id is held: %s", len(conflicts), conflicts)
        return {"learned": learned, "conflicts": conflicts}

    def load_crosswalk(self) -> list[dict]:
        """Every crosswalk row, newest first, in the panel's camelCase shape."""
        with self._lock:
            rows = self._con.execute(_SELECT_CROSSWALK).fetchall()
        return [{
            "venue": venue, "league": league, "venueTeamKey": venue_team_key, "venueTeamName": venue_team_name,
            "unabatedTeamId": unabated_team_id, "unabatedTeamName": unabated_team_name,
            "learnedFrom": learned_from, "learnedAt": _epoch_to_iso(learned_at), "pinnedBetId": pinned_bet_id,
        } for venue, league, venue_team_key, venue_team_name, unabated_team_id, unabated_team_name,
              learned_from, learned_at, pinned_bet_id in rows]

    def clear_crosswalk(self) -> int:
        """DELETE every crosswalk row; returns how many were held."""
        with self._lock:
            count = self._con.execute("SELECT count(*) FROM team_crosswalk").fetchone()[0]
            self._con.execute("DELETE FROM team_crosswalk")
        log.info("crosswalk: cleared %d row(s)", count)
        return count

    def bet_venue(self, bet_id: str) -> str | None:
        """The venue of the stored bet with this id, or None when `bets` holds
        no such record (a pin must name a real bet, of its own venue)."""
        with self._lock:
            row = self._con.execute("SELECT venue FROM bets WHERE id = ?", [bet_id]).fetchone()
        return row[0] if row else None

    def pin_bet(self, pin: dict, crosswalk_rows: list[dict], pinned_at: datetime) -> None:
        """Cal's manual attach, in one transaction: UPSERT the pin (validated
        shape: see service.validate_pin_request), drop the crosswalk rows an
        earlier attach of the same bet taught, then write this attach's rows,
        REPLACING a held key (a manual lesson overrides an automatic one)."""
        with self._lock:
            self._con.execute("BEGIN TRANSACTION")
            try:
                self._con.execute(_UPSERT_PIN, [
                    pin["betId"], pin["venue"], pin["league"], pin["eventId"], _parse_iso(pin.get("eventStart")),
                    pin.get("awayTeamId"), pin.get("homeTeamId"), pin.get("awayTeamName"), pin.get("homeTeamName"),
                    pinned_at])
                self._con.execute("DELETE FROM team_crosswalk WHERE pinned_bet_id = ?", [pin["betId"]])
                for row in crosswalk_rows:
                    self._con.execute(_REPLACE_CROSSWALK, [
                        row["venue"], row["league"], row["venueTeamKey"], row.get("venueTeamName"),
                        row["unabatedTeamId"], row.get("unabatedTeamName"), row.get("learnedFrom"), pinned_at,
                        pin["betId"]])
                self._con.execute("COMMIT")
            except Exception:
                self._con.execute("ROLLBACK")
                raise
        log.info("pin: bet %s -> %s event %s, %d crosswalk row(s) taught",
                 pin["betId"], pin["league"], pin["eventId"], len(crosswalk_rows))

    def unpin_bet(self, bet_id: str) -> dict:
        """Undo an attach: DELETE the pin and the crosswalk rows it taught.
        A row it replaced is not restored. Returns {removedPin, removedRows}."""
        with self._lock:
            self._con.execute("BEGIN TRANSACTION")
            try:
                [removed_pin] = self._con.execute("DELETE FROM bet_pins WHERE bet_id = ?", [bet_id]).fetchone()
                [removed_rows] = self._con.execute(
                    "DELETE FROM team_crosswalk WHERE pinned_bet_id = ?", [bet_id]).fetchone()
                self._con.execute("COMMIT")
            except Exception:
                self._con.execute("ROLLBACK")
                raise
        log.info("pin: bet %s unpinned, %d crosswalk row(s) removed", bet_id, removed_rows)
        return {"removedPin": removed_pin > 0, "removedRows": removed_rows}

    def load_pins(self) -> list[dict]:
        """Every pin, newest first, in the panel's camelCase shape."""
        with self._lock:
            rows = self._con.execute(_SELECT_PINS).fetchall()
        return [{
            "betId": bet_id, "venue": venue, "league": league, "eventId": event_id,
            "eventStart": _epoch_to_iso(event_start), "awayTeamId": away_team_id, "homeTeamId": home_team_id,
            "awayTeamName": away_team_name, "homeTeamName": home_team_name, "pinnedAt": _epoch_to_iso(pinned_at),
        } for bet_id, venue, league, event_id, event_start, away_team_id, home_team_id,
              away_team_name, home_team_name, pinned_at in rows]

    def save_fill_fairs(self, rows: list[dict], captured_at: datetime) -> int:
        """INSERT the fill fairs of bets that have none (validated shape: see
        service.validate_fill_fair_rows). A bet that already has one keeps it
        — ON CONFLICT DO NOTHING, the first capture wins and a saved fair is
        never rewritten. Returns how many rows were inserted."""
        values = [(
            row["betId"], row["lineKey"], row["points"], row["fairAmerican"],
            _parse_iso(row["fairObservedAt"]), _parse_iso(row["placedAt"]), captured_at,
        ) for row in _first_row_per_bet(rows)]
        saved = 0
        with self._lock:
            for start in range(0, len(values), UPSERT_ROWS_PER_STATEMENT):
                chunk = values[start:start + UPSERT_ROWS_PER_STATEMENT]
                statement = _INSERT_FILL_FAIRS.format(values=", ".join([_FILL_FAIR_ROW_PLACEHOLDERS] * len(chunk)))
                # The INSERT's result row is the rows it inserted: the ones DO NOTHING let through.
                [inserted] = self._con.execute(statement, [value for row in chunk for value in row]).fetchone()
                saved += inserted
        if saved:
            log.info("fill fairs: saved %d of %d row(s)", saved, len(values))
        return saved

    def save_fill_books(self, rows: list[dict], captured_at: datetime) -> int:
        """INSERT the other-books rows of bets that have none (validated shape:
        see service.validate_fill_book_rows). A bet that already has one keeps
        it — ON CONFLICT DO NOTHING. Returns how many rows were inserted."""
        saved = 0
        with self._lock:
            for row in rows:
                [inserted] = self._con.execute(_INSERT_FILL_BOOKS, [
                    row["betId"], _parse_iso(row["placedAt"]), row["fillPrice"], row["bookCount"],
                    row["booksBetter"], row["bestBook"], row["bestPrice"], captured_at,
                ]).fetchone()
                saved += inserted
        if saved:
            log.info("fill books: saved %d of %d row(s)", saved, len(rows))
        return saved

    def load_fill_fairs(self, days: int, now: datetime) -> list[dict]:
        """The saved fill fairs of the bets load_bets(days, now) serves, oldest
        placement first, in the panel's camelCase shape."""
        cutoff = now - timedelta(days=days)
        with self._lock:
            rows = self._con.execute(_SELECT_FILL_FAIRS, [cutoff]).fetchall()
        return [{
            "betId": bet_id, "lineKey": line_key, "points": points, "fairAmerican": fair_american,
            "fairObservedAt": _epoch_to_iso(fair_observed_at), "placedAt": _epoch_to_iso(placed_at),
            "capturedAt": _epoch_to_iso(captured_at),
        } for bet_id, line_key, points, fair_american, fair_observed_at, placed_at, captured_at in rows]

    def save_closing_fairs(self, rows: list[dict], updated_at: datetime) -> int:
        """UPSERT the closing fairs (validated shape: see
        service.validate_closing_fair_rows), each kept only when it is newer
        than the stored one. Returns how many rows were written."""
        written = 0
        with self._lock:
            for row in rows:
                [count] = self._con.execute(_UPSERT_CLOSING_FAIR, [
                    row["betId"], row["lineKey"], row["points"], row["fairAmerican"],
                    _parse_iso(row["fairObservedAt"]), _parse_iso(row["eventStart"]), updated_at,
                ]).fetchone()
                written += count
        return written

    def load_closing_fairs(self, days: int, now: datetime) -> list[dict]:
        """The closing fairs of the bets load_bets(days, now) serves, earliest game first."""
        cutoff = now - timedelta(days=days)
        with self._lock:
            rows = self._con.execute(_SELECT_CLOSING_FAIRS, [cutoff]).fetchall()
        return [{
            "betId": bet_id, "lineKey": line_key, "points": points, "fairAmerican": fair_american,
            "fairObservedAt": _epoch_to_iso(fair_observed_at), "eventStart": _epoch_to_iso(event_start),
        } for bet_id, line_key, points, fair_american, fair_observed_at, event_start in rows]

    def set_exclusions(self, bet_ids: list[str], venue_by_bet: dict[str, str], excluded: bool,
                       at: datetime) -> int:
        """Exclude (INSERT, an already-excluded bet keeps its time) or restore
        (DELETE) the bets. Returns how many rows changed."""
        changed = 0
        with self._lock:
            for bet_id in bet_ids:
                if excluded:
                    statement = ("INSERT INTO bet_exclusions (bet_id, venue, excluded_at) VALUES (?, ?, ?) "
                                 "ON CONFLICT (bet_id) DO NOTHING")
                    [count] = self._con.execute(statement, [bet_id, venue_by_bet[bet_id], at]).fetchone()
                else:
                    [count] = self._con.execute("DELETE FROM bet_exclusions WHERE bet_id = ?", [bet_id]).fetchone()
                changed += count
        log.info("exclusions: %s %d of %d bet(s)", "excluded" if excluded else "restored", changed, len(bet_ids))
        return changed

    def load_exclusions(self) -> list[dict]:
        with self._lock:
            rows = self._con.execute(_SELECT_EXCLUSIONS).fetchall()
        return [{"betId": bet_id, "venue": venue, "excludedAt": _epoch_to_iso(excluded_at)}
                for bet_id, venue, excluded_at in rows]

    def set_dismissals(self, bet_ids: list[str], dismissed: bool, at: datetime) -> int:
        """Dismiss (INSERT, an already-dismissed bet keeps its time) or
        restore (DELETE) the bets. Returns how many rows changed."""
        changed = 0
        with self._lock:
            for bet_id in bet_ids:
                if dismissed:
                    statement = ("INSERT INTO bet_dismissals (bet_id, dismissed_at) VALUES (?, ?) "
                                 "ON CONFLICT (bet_id) DO NOTHING")
                    [count] = self._con.execute(statement, [bet_id, at]).fetchone()
                else:
                    [count] = self._con.execute("DELETE FROM bet_dismissals WHERE bet_id = ?", [bet_id]).fetchone()
                changed += count
        log.info("dismissals: %s %d of %d bet(s)", "dismissed" if dismissed else "restored", changed, len(bet_ids))
        return changed

    def load_dismissals(self) -> list[dict]:
        """[{betId, dismissedAt}] of the dismissed bets still open, newest first."""
        with self._lock:
            rows = self._con.execute(_SELECT_DISMISSALS).fetchall()
        return [{"betId": bet_id, "dismissedAt": _epoch_to_iso(dismissed_at)} for bet_id, dismissed_at in rows]

    def set_teaser_block(self, market_key: str, event_start: datetime, blocked: bool, at: datetime) -> int:
        """Mark (UPSERT) or clear (DELETE) one market's Can't tease. Returns
        how many rows changed (an UPSERT always counts 1)."""
        with self._lock:
            if blocked:
                [count] = self._con.execute(_UPSERT_TEASER_BLOCK, [market_key, event_start, at]).fetchone()
            else:
                [count] = self._con.execute("DELETE FROM teaser_blocks WHERE market_key = ?", [market_key]).fetchone()
        log.info("teaser blocks: %s %s", "blocked" if blocked else "restored", market_key)
        return count

    def load_teaser_blocks(self, now: datetime) -> list[dict]:
        """[{marketKey, eventStartMs, blockedAt}] of the blocks whose game is still to start."""
        with self._lock:
            rows = self._con.execute(_SELECT_TEASER_BLOCKS, [now]).fetchall()
        return [{"marketKey": market_key, "eventStartMs": int(event_start_ms), "blockedAt": _epoch_to_iso(blocked_at)}
                for market_key, event_start_ms, blocked_at in rows]

    def set_edge_skip(self, skip_key: str, event_start: datetime, skipped: bool, at: datetime) -> int:
        """Skip (UPSERT) or restore (DELETE) one Edges line at every book.
        Returns how many rows changed (an UPSERT always counts 1). Restore
        also deletes the line's old per-book rows ('<line>@<bookId>')."""
        line = edge_skip_line(skip_key)
        with self._lock:
            if skipped:
                [count] = self._con.execute(_UPSERT_EDGE_SKIP, [line, event_start, at]).fetchone()
            else:
                [count] = self._con.execute(
                    "DELETE FROM edge_skips WHERE skip_key = ? OR starts_with(skip_key, ?)", [line, f"{line}@"]
                ).fetchone()
        log.info("edge skips: %s %s", "skipped" if skipped else "restored", line)
        return count

    def load_edge_skips(self, now: datetime) -> list[dict]:
        """[{skipKey, eventStartMs, skippedAt}] of the skips whose game is still to start."""
        with self._lock:
            rows = self._con.execute(_SELECT_EDGE_SKIPS, [now]).fetchall()
        return [{"skipKey": skip_key, "eventStartMs": int(event_start_ms), "skippedAt": _epoch_to_iso(skipped_at)}
                for skip_key, event_start_ms, skipped_at in rows]

    def load_edge_settings(self) -> dict:
        """{settings: {field: value or None}, updatedAt: ISO or None} — every
        field of EDGE_SETTINGS_FIELDS, None where nothing overrides the
        default; all None and updatedAt None before the first write."""
        with self._lock:
            row = self._con.execute(_SELECT_EDGE_SETTINGS).fetchone()
        if row is None:
            return {"settings": dict.fromkeys(EDGE_SETTINGS_FIELDS), "updatedAt": None}
        *values, updated_at = row
        return {"settings": dict(zip(EDGE_SETTINGS_FIELDS, values)), "updatedAt": _epoch_to_iso(updated_at)}

    def save_edge_settings(self, settings: dict, updated_at: datetime) -> None:
        """UPSERT the one settings row whole (validated shape: see
        service.validate_settings_update). `settings` must carry every field
        of EDGE_SETTINGS_FIELDS; None stores NULL, i.e. the default."""
        missing = [field for field in EDGE_SETTINGS_FIELDS if field not in settings]
        if missing:
            raise ValueError(f"save_edge_settings: expected every settings field, missing {missing}")
        with self._lock:
            self._con.execute(_UPSERT_EDGE_SETTINGS, [settings[field] for field in EDGE_SETTINGS_FIELDS] + [updated_at])
        log.info("edge settings: saved %s", {field: value for field, value in settings.items() if value is not None})

    def bet_is_open(self, bet_id: str) -> bool:
        """Whether `bets` holds this id as an open record."""
        with self._lock:
            row = self._con.execute("SELECT status FROM bets WHERE id = ?", [bet_id]).fetchone()
        return bool(row) and row[0] == "open"

    def novig_outcome_is_held(self, outcome_id: str) -> bool:
        """Whether an open straight Novig bet is on this outcome id."""
        with self._lock:
            row = self._con.execute(_SELECT_OPEN_NOVIG_ON_OUTCOME, [outcome_id]).fetchone()
        return row is not None

    def insert_manual_order(self, order: dict) -> bool:
        """INSERT a Place-button order before it is sent (status 'sending').
        False when its client_order_id is already held: that click was sent."""
        with self._lock:
            held = self._con.execute("SELECT 1 FROM manual_orders WHERE client_order_id = ?",
                                     [order["client_order_id"]]).fetchone()
            if held:
                return False
            values = ([order[column] for column in MANUAL_ORDER_INSERT_COLUMNS]
                      + [order.get(column) for column in MANUAL_ORDER_NOVIG_LINE_COLUMNS])
            self._con.execute(_INSERT_MANUAL_ORDER, values)
        return True

    def novig_place_order_contexts(self) -> dict[str, dict]:
        """Novig order id -> the board line a Novig Place order bought, in the
        bot desk's context shape (taker_store.NOVIG_CONTEXT_COLUMNS), for every
        order Novig gave an id. Place's sub-account read names its bets with
        these; `settled_at` is the bet's fixed closedAt once set."""
        with self._lock:
            rows = self._con.execute(_SELECT_NOVIG_PLACE_CONTEXTS).fetchall()
        contexts = {}
        for (order_id, league, bet_type, period, side, points, away_team, home_team, away_team_id, home_team_id,
             event_start_ms, event_id, market_id, outcome_id, price, settled_at_ms) in rows:
            contexts[order_id] = {
                "league": league, "bet_type": bet_type, "period": period, "side": side, "points": points,
                "away_team": away_team, "home_team": home_team, "away_team_id": away_team_id,
                "home_team_id": home_team_id, "event_start": _utc_of_epoch_ms(event_start_ms), "event_id": event_id,
                "market_id": market_id, "outcome_id": outcome_id, "price": price,
                "settled_at": _utc_of_epoch_ms(settled_at_ms),
            }
        return contexts

    def adopt_novig_place_order(self, client_order_id: str, novig_order_id: str | None, at: datetime) -> None:
        """A lost Novig Place order (no answer, no id) found on Novig's order
        lists gets its id, so its fills are named from now on; one listed
        nowhere is marked rejected: it never reached Novig."""
        if novig_order_id is None:
            fields = {"status": "rejected", "error": "never reached Novig: in none of its order lists"}
        else:
            fields = {"novig_order_id": novig_order_id, "error": "answer lost; found on Novig by clientId"}
        self.update_manual_order(client_order_id, fields, at)

    def unconfirmed_novig_place_orders(self, placed_before: datetime) -> list[dict]:
        """Novig Place orders with a Novig id still 'unknown' (never read final,
        or their fills unread), placed at or before `placed_before` and within
        UNCONFIRMED_NOVIG_PLACE_WINDOW of it: [{client_order_id, novig_order_id, count_requested}]."""
        with self._lock:
            rows = self._con.execute(_SELECT_UNCONFIRMED_NOVIG_PLACE,
                                     [placed_before, placed_before - UNCONFIRMED_NOVIG_PLACE_WINDOW]).fetchall()
        return [{"client_order_id": client_order_id, "novig_order_id": novig_order_id, "count_requested": count_requested}
                for client_order_id, novig_order_id, count_requested in rows]

    def confirm_novig_place_order(self, client_order_id: str, status: str, filled: int, paid: float,
                                  at: datetime) -> None:
        """An 'unknown' Novig Place order Novig now shows final: what filled."""
        self.update_manual_order(client_order_id, {"status": status, "count_filled": filled, "cost_filled": paid,
                                                   "error": "confirmed later by the Place sub-account read"}, at)

    def note_novig_place_settled(self, novig_order_id: str, at: datetime) -> datetime:
        """Keep when Place's sub-account read first saw this order's market
        SETTLED (only the first sighting is kept); returns the kept time."""
        with self._lock:
            self._con.execute("UPDATE manual_orders SET settled_at = ? "
                              "WHERE novig_order_id = ? AND venue = 'novig' AND settled_at IS NULL", [at, novig_order_id])
            row = self._con.execute("SELECT epoch_ms(settled_at) FROM manual_orders WHERE novig_order_id = ?",
                                    [novig_order_id]).fetchone()
        return _utc_of_epoch_ms(row[0]) if row and row[0] is not None else at

    def update_manual_order(self, client_order_id: str, fields: dict, at: datetime) -> None:
        """UPDATE the named columns of one manual order (and updated_at)."""
        unknown = set(fields) - set(MANUAL_ORDER_UPDATE_COLUMNS)
        if unknown:
            raise ValueError(f"update_manual_order: unknown columns {sorted(unknown)}")
        assignments = ", ".join(f"{column} = ?" for column in fields)
        with self._lock:
            self._con.execute(f"UPDATE manual_orders SET {assignments}, updated_at = ? WHERE client_order_id = ?",
                              [*fields.values(), at, client_order_id])

    def manual_orders_maybe_filled_since(self, event_id: int, since: datetime, venue: str | None = None) -> list[dict]:
        """Manual orders on this game that may hold contracts the read at
        `since` missed (_SELECT_MANUAL_MAYBE_FILLED), at one venue ('kalshi' |
        'novig') or, with venue None, at either."""
        with self._lock:
            rows = self._con.execute(_SELECT_MANUAL_MAYBE_FILLED, [event_id, venue, venue, since]).fetchall()
        return [{"client_order_id": client_order_id, "ticker": ticker, "status": status,
                 "placed_at": _epoch_to_iso(placed_at)} for client_order_id, ticker, status, placed_at in rows]

    def lost_novig_place_orders(self, placed_before: datetime) -> list[dict]:
        """Novig Place orders with no definite answer and no Novig order id,
        placed at or before `placed_before`: [{client_order_id, outcome_id}]."""
        with self._lock:
            rows = self._con.execute(_SELECT_LOST_NOVIG_PLACE, [placed_before]).fetchall()
        return [{"client_order_id": client_order_id, "outcome_id": ticker} for client_order_id, ticker in rows]

    def source_status(self) -> dict[str, dict]:
        """{source: {fetchedAt, ok, error, count}} — fetchedAt/count from the
        latest SUCCESSFUL run (freshness), ok/error from the latest run of any kind."""
        with self._lock:
            latest = self._con.execute(_SELECT_LATEST_RUNS).fetchall()
            latest_ok = {row[0]: row for row in self._con.execute(_SELECT_LATEST_OK_RUNS).fetchall()}
        status = {}
        for source, started_at, ok, error, _n_records in latest:
            ok_row = latest_ok.get(source)
            status[source] = {
                "fetchedAt": _epoch_to_iso(ok_row[1]) if ok_row else None,
                "ok": bool(ok),
                "error": error,
                "count": ok_row[2] if ok_row else 0,
            }
        return status
