# Answer Keys/MLB Answer Key/market_edge.R
# Market-consensus edge detection for the MLB bets tab.
#
# find_market_edges() reads the canonical per-book odds list MLB.R already
# builds (book_odds_by_book) and flags any book whose offered price beats the
# LEAVE-ONE-OUT median devigged fair of the OTHER books at the same wager.
# Pure: no DuckDB, no model dependency. Output rows merge onto mlb_bets_combined.
#
# Depends on (source first): Tools.R (devig_american, compute_ev, kelly_stake),
# odds_screen.R (compute_bet_row_id). Tested by tests/test_market_edge.R.

library(dplyr)
library(tibble)

# American odds -> raw implied probability (no devig).
.implied_prob <- function(o) ifelse(o > 0, 100 / (o + 100), -o / (-o + 100))

# Devig a single 2-outcome pair, returning the fair prob for each input side.
.pair_fairs <- function(o1, o2) {
  d <- devig_american(o1, o2)   # data.frame(p1, p2)
  c(d$p1, d$p2)
}

.empty_market_edges <- function() {
  tibble(
    bet_row_id = character(), id = character(), game_id = character(),
    market = character(), market_type = character(), period = character(),
    line = numeric(), bet_on = character(), bookmaker_key = character(),
    odds = integer(), prob = numeric(), ev = numeric(), n_books = integer(),
    bet_size = numeric(), edge_source = character(),
    model_ev = numeric(), market_ev = numeric(),
    home_team = character(), away_team = character(),
    pt_start_time = as.POSIXct(character(), tz = "UTC")
  )
}

#' @param book_odds_by_book named list of canonical frames
#'        (game_id, market, period, side, line, american_odds, fetch_time)
#' @param game_info optional tibble(game_id, home_team, away_team, pt_start_time)
#' @param threshold minimum market EV to flag (default 0.02)
#' @param min_others minimum OTHER books to form the leave-one-out yardstick
#' @param staleness_min drop quotes older than this many minutes before `now`
#' @param now reference time for staleness (default Sys.time())
#' @param bankroll,kelly_mult sizing inputs for kelly_stake()
find_market_edges <- function(book_odds_by_book,
                              game_info     = NULL,
                              threshold     = 0.02,
                              min_others    = 1,
                              staleness_min = 30,
                              now           = Sys.time(),
                              bankroll      = 100,
                              kelly_mult    = 0.25) {

  if (length(book_odds_by_book) == 0) return(.empty_market_edges())

  # 1. Stack books into one long frame tagged with bookmaker_key.
  long <- bind_rows(lapply(names(book_odds_by_book), function(bk) {
    df <- book_odds_by_book[[bk]]
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df$bookmaker_key <- bk
    df
  }))
  if (is.null(long) || nrow(long) == 0) return(.empty_market_edges())

  # 2. Drop stale quotes.
  long <- long %>%
    filter(!is.na(fetch_time),
           as.numeric(difftime(now, fetch_time, units = "mins")) <= staleness_min)
  if (nrow(long) == 0) return(.empty_market_edges())

  # base market type collapses alt + main into one bucket (like
  # .related_market_types); line_key normalizes 8.5 vs 8.50.
  long <- long %>%
    mutate(base_mt  = gsub("^alternate_", "", market),
           line_key = ifelse(is.na(line), NA_real_, round(as.numeric(line), 1)),
           absline  = ifelse(is.na(line_key), -1, abs(line_key)))

  # 3. Validate each book quotes BOTH sides (a proper 2-outcome pair):
  #    Over+Under for totals, the two teams for spreads / h2h, etc.
  #    Books that quote only one side are dropped (no pair → excluded from
  #    both the consensus AND the candidate set).
  #    After pair validation, store the DEVIGGED fair probability per row —
  #    this is the vig-removed true probability, which is the correct
  #    yardstick for measuring edge against consensus.
  paired <- long %>%
    group_by(game_id, base_mt, period, bookmaker_key, absline) %>%
    filter(n() == 2) %>%
    mutate(fair = .pair_fairs(american_odds[1], american_odds[2])) %>%
    ungroup() %>%
    filter(!is.na(fair))
  if (nrow(paired) == 0) return(.empty_market_edges())

  # 4. Leave-one-out consensus: each book judged vs the median DEVIGGED fair
  #    probability of the OTHER books at the same (game, base_mt, period, side, line).
  #    Using devigged fairs removes bookmaker margin from the yardstick so that
  #    edge is measured against the true market probability, not the vig-inflated
  #    implied price. A soft book offering +110 vs a consensus of -110/-110 reads
  #    as +5.0% EV (fair=50%, book implied=47.6%, EV=50%-47.6%=+5.0%).
  judged <- paired %>%
    group_by(game_id, base_mt, period, side, line_key) %>%
    mutate(
      n_books  = n(),
      loo_fair = vapply(seq_len(n()),
                        function(i) stats::median(fair[-i]), numeric(1))
    ) %>%
    ungroup() %>%
    filter(n_books - 1 >= min_others, !is.na(loo_fair))
  if (nrow(judged) == 0) return(.empty_market_edges())

  # 5. EV of each book's price vs its leave-one-out consensus; keep flags.
  flagged <- judged %>%
    mutate(market_ev = compute_ev(loo_fair, .implied_prob(american_odds))) %>%
    filter(market_ev >= threshold)
  if (nrow(flagged) == 0) return(.empty_market_edges())

  # Best book per (side, line), then best line per side (mirror model dedup).
  best <- flagged %>%
    group_by(game_id, base_mt, period, side, line_key) %>%
    slice_max(market_ev, n = 1, with_ties = FALSE) %>%
    group_by(game_id, base_mt, period, side) %>%
    slice_max(market_ev, n = 1, with_ties = FALSE) %>%
    ungroup()

  out <- best %>%
    transmute(
      id            = game_id,
      game_id       = game_id,
      market        = base_mt,
      market_type   = base_mt,
      period        = period,
      line          = line_key,
      bet_on        = side,
      bookmaker_key = bookmaker_key,
      odds          = as.integer(american_odds),
      prob          = loo_fair,
      ev            = market_ev,
      n_books       = as.integer(n_books),
      edge_source   = "market",
      model_ev      = NA_real_,
      market_ev     = market_ev,
      bet_size      = kelly_stake(market_ev, .implied_prob(american_odds),
                                  bankroll, kelly_mult)
    ) %>%
    # NOTE: this bet_row_id only joins against mlb_bets_combined AFTER MLB.R is
    # migrated to compute_bet_row_id() (Task 3). The legacy inline hash in MLB.R
    # uses a different recipe; until Task 3, a join would silently match zero rows.
    mutate(bet_row_id = compute_bet_row_id(id, market_type, period, line, bet_on))

  # Attach game header info if provided.
  if (!is.null(game_info) && nrow(game_info) > 0) {
    out <- out %>% left_join(game_info, by = "game_id")
  } else {
    out <- out %>% mutate(home_team = NA_character_, away_team = NA_character_,
                          pt_start_time = as.POSIXct(NA, tz = "UTC"))
  }
  out
}
