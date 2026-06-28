# =============================================================================
# CBB Extreme-Samples Shrinkage A/B Backtest
# =============================================================================
# PURPOSE: Measure whether Beta-Binomial shrinkage of the Elihu/Feustel
# matched-sample probability toward the (devigged) market price improves
# out-of-sample probability accuracy and betting outcomes WITHOUT overfitting.
#
# Estimator (Section A of the research):
#   p_hat = (x + n0 * p_mkt) / (n + n0)  =  w*(x/n) + (1-w)*p_mkt,  w = n/(n+n0)
#   - x/n   : raw matched-sample empirical frequency (current production)
#   - p_mkt : devigged market price at that exact line (the shrink target)
#   - n0    : prior strength, FIT ON A TEMPORAL TRAIN SPLIT (never hand-tuned)
#
# Kelly-under-uncertainty (Section C, Baker & McHale 2013):
#   alpha = edge^2 / (edge^2 + Var(p_hat)),   Var(p_hat) = p(1-p)/(n+n0+1)
#   stake = 0.25 * alpha * full_kelly
#
# Validation (Section D): the metric we report comes from games the n0 fit
# never saw. Primary evidence = held-out log-loss / Brier / calibration with a
# Diebold-Mariano significance test. ROI is corroboration only (high variance).
#
# This script does NOT modify Tools.R or production. It is a measurement
# instrument. Reads cbb.duckdb read-only from the MAIN checkout.
# =============================================================================

# --- paths: source Tools.R from THIS worktree, read DB from main checkout ----
WORKTREE_AK <- "/Users/callancapitolo/NFLWork/.claude/worktrees/extreme-samples-shrinkage/Answer Keys"
MAIN_DB     <- "/Users/callancapitolo/NFLWork/Answer Keys/cbb.duckdb"

setwd(WORKTREE_AK)
suppressMessages({
  library(data.table); library(duckdb); library(dplyr); library(tidyverse)
  library(DBI); library(parallel)
})
source("Tools.R")

# =============================================================================
# CONFIG
# =============================================================================
N_TEST_GAMES <- as.integer(Sys.getenv("N_TEST_GAMES", unset = NA))  # NA = all
SAMPLE_PCT   <- 0.02     # production: CBBPrepare.R
KELLY_MULT   <- 0.25
BANKROLL     <- 1000
EV_THRESHOLD <- 0.05     # production: CBBCombine.R
N_CORES      <- max(1, detectCores() - 1)
N0_GRID      <- c(0, 5, 10, 15, 20, 30, 40, 50, 75, 100, 150, 200, 300, 500, 1000, Inf)
MIN_NONPUSH  <- 10       # same floor as production backtest

cat("Cores:", N_CORES, "| N_TEST_GAMES:", ifelse(is.na(N_TEST_GAMES), "ALL", N_TEST_GAMES), "\n")

# =============================================================================
# HELPERS
# =============================================================================
american_to_prob <- function(odds) ifelse(odds > 0, 100/(odds+100), abs(odds)/(abs(odds)+100))
devig_pair <- function(o1, o2) { p1<-american_to_prob(o1); p2<-american_to_prob(o2); t<-p1+p2; list(prob1=p1/t, prob2=p2/t) }
calc_logloss <- function(p, a) -(a*log(pmax(p,1e-10)) + (1-a)*log(pmax(1-p,1e-10)))
calc_brier   <- function(p, a) (p - a)^2
calc_ev      <- function(p, book_odds) { d <- ifelse(book_odds>0, 1+book_odds/100, 1+100/abs(book_odds)); p*d - 1 }

# Beta-Binomial shrinkage toward per-row market prob p_mkt with strength n0.
shrink_prob <- function(x, n, p_mkt, n0) {
  if (is.infinite(n0)) return(p_mkt)
  (x + n0 * p_mkt) / (n + n0)
}
# Posterior variance of the shrunk proportion.
shrink_var <- function(p_hat, n, n0) {
  if (is.infinite(n0)) return(rep(0, length(p_hat)))
  p_hat * (1 - p_hat) / (n + n0 + 1)
}
# Baker & McHale (2013) Kelly fraction multiplier from estimation variance.
bm_alpha <- function(edge, var_p) {
  a <- edge^2 / (edge^2 + var_p)
  ifelse(is.finite(a), a, 1)
}
# Kelly stake from a (already-final) probability + book odds.
kelly_stake_p <- function(p, book_odds, bankroll, kelly_mult, alpha = 1) {
  d <- ifelse(book_odds>0, 1+book_odds/100, 1+100/abs(book_odds))
  b <- d - 1
  q <- 1 - p
  f <- (b*p - q)/b
  f <- pmax(0, f)
  bankroll * kelly_mult * alpha * f
}

# =============================================================================
# LOAD + PREP (mirror CBB_Backtest_Parallel.R, add game_date)
# =============================================================================
cat("Loading data (read-only) from main checkout...\n")
con <- dbConnect(duckdb(), dbdir = MAIN_DB, read_only = TRUE)
pbp <- dbGetQuery(con, "SELECT * FROM cbb_betting_pbp")
closing_odds <- dbGetQuery(con, "
  SELECT id as odds_id, home_spread, total_line,
         spread_home_odds, spread_away_odds, tot_over_odds, tot_under_odds
  FROM cbb_closing_odds WHERE spread_home_odds IS NOT NULL")
deriv_odds <- dbGetQuery(con, "SELECT * FROM cbb_derivative_closing_odds")
dbDisconnect(con)
cat("PBP rows:", nrow(pbp), "| deriv rows:", nrow(deriv_odds), "\n")

closing_odds <- closing_odds %>% mutate(
  implied_home=american_to_prob(spread_home_odds), implied_away=american_to_prob(spread_away_odds),
  implied_over=american_to_prob(tot_over_odds),   implied_under=american_to_prob(tot_under_odds),
  devig_home=implied_home/(implied_home+implied_away), devig_over=implied_over/(implied_over+implied_under))

consensus <- closing_odds %>% group_by(odds_id) %>% summarize(
  home_spread=median(home_spread,na.rm=TRUE), total_line=median(total_line,na.rm=TRUE),
  consensus_home_cover_prob=mean(devig_home,na.rm=TRUE), consensus_over_prob=mean(devig_over,na.rm=TRUE),
  .groups="drop")

pbp <- pbp %>% select(-any_of(c("home_spread","total_line")))  # avoid join collision (pbp has its own)
DT <- pbp %>% inner_join(consensus, by="odds_id") %>%
  filter(!is.na(home_spread), !is.na(total_line), !is.na(game_home_margin_h1), !is.na(game_home_margin_h2)) %>%
  mutate(actual_cover=ifelse(game_home_margin_fg > -home_spread,1,0),
         actual_over =ifelse(game_total_fg > total_line,1,0)) %>%
  as.data.table()
cat("Games with full data:", nrow(DT), "\n")

# Derivative odds -> wide (H1 spreads / totals / ml)
spreads_h1 <- deriv_odds %>% filter(market=="spreads_h1") %>%
  mutate(side=ifelse(outcome_name==home_team,"home","away"),
         line=ifelse(side=="home",outcome_point,-outcome_point)) %>%
  group_by(event_id,line) %>% summarize(
    best_home_odds=max(outcome_price[side=="home"],na.rm=TRUE),
    best_away_odds=max(outcome_price[side=="away"],na.rm=TRUE), .groups="drop") %>%
  filter(is.finite(best_home_odds), is.finite(best_away_odds))
totals_h1 <- deriv_odds %>% filter(market=="totals_h1") %>%
  mutate(side=tolower(outcome_name)) %>% group_by(event_id,outcome_point) %>% summarize(
    best_over_odds=max(outcome_price[side=="over"],na.rm=TRUE),
    best_under_odds=max(outcome_price[side=="under"],na.rm=TRUE), .groups="drop") %>%
  rename(line=outcome_point) %>% filter(is.finite(best_over_odds), is.finite(best_under_odds))
ml_h1 <- deriv_odds %>% filter(market=="h2h_h1") %>%
  mutate(side=ifelse(outcome_name==home_team,"home","away")) %>% group_by(event_id) %>% summarize(
    best_home_odds=max(outcome_price[side=="home"],na.rm=TRUE),
    best_away_odds=max(outcome_price[side=="away"],na.rm=TRUE), .groups="drop") %>%
  filter(is.finite(best_home_odds), is.finite(best_away_odds))

games_with_derivs <- unique(c(spreads_h1$event_id, totals_h1$event_id, ml_h1$event_id))
DT_test <- DT %>% filter(odds_id %in% games_with_derivs)
cat("Games with derivative odds:", nrow(DT_test), "\n")

if (!is.na(N_TEST_GAMES) && N_TEST_GAMES < nrow(DT_test)) {
  set.seed(42); test_ids <- sample(unique(DT_test$odds_id), N_TEST_GAMES)
} else test_ids <- unique(DT_test$odds_id)
cat("Pricing", length(test_ids), "games\n")

disp <- compute_dispersion(DT, moneyline=FALSE, spread_col="home_spread", total_col="total_line")
ss <- disp$ss; st <- disp$st
N  <- round(nrow(DT) * SAMPLE_PCT, 0)
cat("Sample size N:", N, "| ss:", round(ss,3), "st:", round(st,3), "\n\n")

# Support guard reference: 99th pct of historical lines
P99_SPREAD <- quantile(abs(DT$home_spread), 0.99, na.rm=TRUE)
P99_TOTAL  <- quantile(DT$total_line, 0.99, na.rm=TRUE)

# =============================================================================
# WORKER: price one game, emit per-(line,side) rows with x, n, book_prob, etc.
# =============================================================================
process_single_game <- function(test_id) {
  out <- list()
  tg <- DT_test[DT_test$odds_id==test_id, ][1, ]
  parent_spread <- tg$home_spread; parent_total <- tg$total_line
  target_cover  <- tg$consensus_home_cover_prob; target_over <- tg$consensus_over_prob
  if (is.na(parent_spread)||is.na(parent_total)||is.na(target_cover)||is.na(target_over)) return(NULL)

  DT_train <- DT[DT$odds_id != test_id, ]
  sr <- tryCatch(run_answer_key_sample(id=test_id, parent_spread=parent_spread, parent_total=parent_total,
            target_cover=target_cover, target_over=target_over, DT=DT_train, ss=ss, st=st, N=N,
            use_spread_line=TRUE), error=function(e) NULL)
  if (is.null(sr)) return(NULL)
  smp <- sr$sample; final_N <- sr$final_N
  gdate <- as.character(tg$game_date)
  is_extreme <- (abs(parent_spread) > P99_SPREAD) || (parent_total > P99_TOTAL) || (final_N < 0.5*N)

  add_row <- function(market, line, side, x, n, p_mkt, book_odds, actual) {
    out[[length(out)+1]] <<- data.frame(game_id=test_id, game_date=gdate, market=market, line=line,
      side=side, x=x, n=n, algo_prob=x/n, book_prob=p_mkt, book_odds=book_odds, actual=actual,
      final_N=final_N, is_extreme=is_extreme,
      # production guard signal: matched-sample COLLAPSE only (final_N<0.5*N),
      # NOT the p99-line terms in is_extreme. Mirrors Tools.R low_confidence.
      low_conf=(final_N < 0.5*N), stringsAsFactors=FALSE)
  }

  # H1 SPREADS
  gs <- spreads_h1[spreads_h1$event_id==test_id, ]
  if (nrow(gs)>0) {
    m <- smp$game_home_margin_h1; am <- tg$game_home_margin_h1
    for (j in 1:nrow(gs)) {
      line<-gs$line[j]; np <- m != -line; nn<-sum(np)
      if (nn<MIN_NONPUSH || am==-line) next
      x <- sum(m[np] > -line); bk <- devig_pair(gs$best_home_odds[j], gs$best_away_odds[j])
      ah <- ifelse(am > -line,1,0)
      add_row("spreads_h1", line, "home", x, nn, bk$prob1, gs$best_home_odds[j], ah)
      add_row("spreads_h1", line, "away", nn-x, nn, bk$prob2, gs$best_away_odds[j], 1-ah)
    }
  }
  # H1 TOTALS
  gt <- totals_h1[totals_h1$event_id==test_id, ]
  if (nrow(gt)>0) {
    tt <- smp$game_total_h1; at <- tg$game_total_h1
    for (j in 1:nrow(gt)) {
      line<-gt$line[j]; np <- tt != line; nn<-sum(np)
      if (nn<MIN_NONPUSH || at==line) next
      x <- sum(tt[np] > line); bk <- devig_pair(gt$best_over_odds[j], gt$best_under_odds[j])
      ao <- ifelse(at > line,1,0)
      add_row("totals_h1", line, "over",  x, nn, bk$prob1, gt$best_over_odds[j], ao)
      add_row("totals_h1", line, "under", nn-x, nn, bk$prob2, gt$best_under_odds[j], 1-ao)
    }
  }
  # H1 MONEYLINE
  gm <- ml_h1[ml_h1$event_id==test_id, ]
  if (nrow(gm)>0) {
    m <- smp$game_home_margin_h1; am <- tg$game_home_margin_h1
    np <- m != 0; nn <- sum(np)
    if (nn>=MIN_NONPUSH && am!=0) {
      x <- sum(m[np] > 0); bk <- devig_pair(gm$best_home_odds[1], gm$best_away_odds[1])
      ah <- ifelse(am>0,1,0)
      add_row("h2h_h1", NA_real_, "home", x, nn, bk$prob1, gm$best_home_odds[1], ah)
      add_row("h2h_h1", NA_real_, "away", nn-x, nn, bk$prob2, gm$best_away_odds[1], 1-ah)
    }
  }
  if (length(out)==0) return(NULL)
  do.call(rbind, out)
}

# =============================================================================
# RUN PARALLEL
# =============================================================================
cat("Running parallel pricing...\n")
t0 <- Sys.time()
cl <- makeCluster(N_CORES)
clusterExport(cl, c("DT","DT_test","spreads_h1","totals_h1","ml_h1","ss","st","N",
                    "devig_pair","american_to_prob","process_single_game",
                    "MIN_NONPUSH","P99_SPREAD","P99_TOTAL"))
clusterEvalQ(cl, { library(data.table); library(dplyr)
  setwd("/Users/callancapitolo/NFLWork/.claude/worktrees/extreme-samples-shrinkage/Answer Keys")
  source("Tools.R") })
res_list <- parLapplyLB(cl, test_ids, function(id)
  tryCatch(process_single_game(id), error=function(e) list(error=TRUE, id=id, msg=conditionMessage(e))))
stopCluster(cl)
elapsed <- as.numeric(difftime(Sys.time(), t0, units="mins"))
cat("Priced in", round(elapsed,1), "min\n")

errs <- sapply(res_list, function(x) is.list(x) && !is.null(x$error))
cat("Errors:", sum(errs), "\n")
if (sum(errs)>0) for (e in res_list[errs][1:min(5,sum(errs))]) cat("  ", e$id, ":", e$msg, "\n")
res <- bind_rows(res_list[!errs & !sapply(res_list, is.null)])
cat("Total predictions:", nrow(res), "| games:", n_distinct(res$game_id), "\n\n")
if (nrow(res)==0) { cat("No results.\n"); quit(save="no", status=1) }

# Persist raw predictions so analysis can be re-run without re-pricing.
scratch <- dbConnect(duckdb(), dbdir="cbb_backtest_shrinkage.duckdb")
dbWriteTable(scratch, "preds", res, overwrite=TRUE)
dbDisconnect(scratch, shutdown=TRUE)
cat("Saved preds to cbb_backtest_shrinkage.duckdb::preds\n\n")

# =============================================================================
# ANALYSIS  (sourced separately so it can be iterated without re-pricing)
# =============================================================================
source("CBB Answer Key/CBB_Backtest_Shrinkage_Analysis.R")
