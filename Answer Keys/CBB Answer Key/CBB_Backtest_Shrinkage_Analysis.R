# =============================================================================
# CBB Shrinkage A/B — ANALYSIS (reads preds, fits n0 on TRAIN, evaluates on TEST)
# Sourced by CBB_Backtest_Shrinkage_AB.R, or run standalone after pricing:
#   setwd(".../Answer Keys"); source("CBB Answer Key/CBB_Backtest_Shrinkage_Analysis.R")
# =============================================================================
if (!exists("res")) {
  suppressMessages({ library(duckdb); library(dplyr); library(DBI) })
  con <- dbConnect(duckdb(), dbdir="cbb_backtest_shrinkage.duckdb", read_only=TRUE)
  res <- dbGetQuery(con, "SELECT * FROM preds"); dbDisconnect(con, shutdown=TRUE)
  KELLY_MULT<-0.25; BANKROLL<-1000; EV_THRESHOLD<-0.05
  N0_GRID <- c(0,5,10,15,20,30,40,50,75,100,150,200,300,500,1000,Inf)
  american_to_prob <- function(o) ifelse(o>0,100/(o+100),abs(o)/(abs(o)+100))
  calc_logloss <- function(p,a) -(a*log(pmax(p,1e-10))+(1-a)*log(pmax(1-p,1e-10)))
  calc_brier   <- function(p,a) (p-a)^2
  calc_ev      <- function(p,bo){d<-ifelse(bo>0,1+bo/100,1+100/abs(bo)); p*d-1}
  shrink_prob  <- function(x,n,pm,n0) if(is.infinite(n0)) pm else (x+n0*pm)/(n+n0)
  shrink_var   <- function(p,n,n0) if(is.infinite(n0)) rep(0,length(p)) else p*(1-p)/(n+n0+1)
  bm_alpha     <- function(e,v){a<-e^2/(e^2+v); ifelse(is.finite(a),a,1)}
  kelly_stake_p<- function(p,bo,br,km,al=1){d<-ifelse(bo>0,1+bo/100,1+100/abs(bo));b<-d-1;f<-pmax(0,(b*p-(1-p))/b);br*km*al*f}
}

res <- res %>% mutate(game_date = as.Date(game_date))
res <- res[!is.na(res$game_date), ]

# ---- Temporal split: earlier games fit n0, later games are the holdout -------
cut_date <- as.Date(quantile(unclass(res$game_date), 0.5, type=1), origin="1970-01-01")
train <- res %>% filter(game_date <= cut_date)
test  <- res %>% filter(game_date >  cut_date)
cat("=====================================================================\n")
cat("TEMPORAL SPLIT  cut =", as.character(cut_date), "\n")
cat("  TRAIN games:", n_distinct(train$game_id), " preds:", nrow(train), "\n")
cat("  TEST  games:", n_distinct(test$game_id),  " preds:", nrow(test),  "\n")
cat("  extreme preds — train:", sum(train$is_extreme), " test:", sum(test$is_extreme), "\n\n")

mean_ll_for_n0 <- function(df, n0) {
  p <- shrink_prob(df$x, df$n, df$book_prob, n0); mean(calc_logloss(p, df$actual))
}
# ---- (1) Grid-search n0 on TRAIN (minimize shrunk held-out log-loss) ---------
grid_train <- sapply(N0_GRID, function(n0) mean_ll_for_n0(train, n0))
n0_cv <- N0_GRID[which.min(grid_train)]

# ---- (2) Empirical-Bayes n0 (1-D MLE of Beta-Binomial, per-row prior=book) ---
nll_bb <- function(n0) {
  if (n0 <= 0) return(Inf)
  a <- n0*train$book_prob; b <- n0*(1-train$book_prob)
  ll <- lchoose(train$n, train$x) + lbeta(train$x+a, train$n-train$x+b) - lbeta(a,b)
  -sum(ll[is.finite(ll)])
}
n0_eb <- tryCatch(optimize(nll_bb, c(0.5, 2000))$minimum, error=function(e) NA)

cat("---------------------------------------------------------------------\n")
cat("n0 FIT ON TRAIN\n")
cat("  CV grid-search (min log-loss): n0 =", n0_cv, "\n")
cat("  Empirical-Bayes MLE          : n0 =", round(n0_eb,1), "\n")
cat("  TRAIN log-loss vs n0:\n")
print(data.frame(n0=N0_GRID, train_logloss=round(grid_train,5)))
cat("\n")

# Principled choice: empirical Bayes (marginal likelihood) is the PRIMARY n0 —
# it measures how much the model's matched-sample counts actually deviate from
# the market beyond binomial noise, is fit once on TRAIN and frozen, and does
# NOT degenerate to "defer to book / place no bets" the way pure log-loss CV can
# (the closing line is near-unbeatable on AVERAGE proper score by construction).
# CV grid is reported as a diagnostic / sanity cross-check only.
N0 <- if (is.finite(n0_eb)) round(n0_eb) else n0_cv
cat(">>> FROZEN n0 for evaluation:", N0, "(empirical Bayes; CV grid is diagnostic)\n\n")

# ---- Build the three probability arms on TEST --------------------------------
test <- test %>% mutate(
  p_raw    = algo_prob,
  p_shrunk = shrink_prob(x, n, book_prob, N0),
  p_book   = book_prob,
  var_shrunk = shrink_var(p_shrunk, n, N0)
)

# ---- (3) PRIMARY EVIDENCE: held-out proper scores ----------------------------
score_tbl <- test %>% summarize(
  raw_logloss   = mean(calc_logloss(p_raw,    actual)),
  shrunk_logloss= mean(calc_logloss(p_shrunk, actual)),
  book_logloss  = mean(calc_logloss(p_book,   actual)),
  raw_brier     = mean(calc_brier(p_raw,    actual)),
  shrunk_brier  = mean(calc_brier(p_shrunk, actual)),
  book_brier    = mean(calc_brier(p_book,   actual)))
cat("=====================================================================\n")
cat("PRIMARY: HELD-OUT PROPER SCORES (TEST, lower=better)\n")
cat("=====================================================================\n")
cat(sprintf("  log-loss:  raw=%.5f  shrunk=%.5f  book=%.5f   (shrunk-raw=%+.5f)\n",
    score_tbl$raw_logloss, score_tbl$shrunk_logloss, score_tbl$book_logloss,
    score_tbl$shrunk_logloss-score_tbl$raw_logloss))
cat(sprintf("  Brier:     raw=%.5f  shrunk=%.5f  book=%.5f   (shrunk-raw=%+.5f)\n\n",
    score_tbl$raw_brier, score_tbl$shrunk_brier, score_tbl$book_brier,
    score_tbl$shrunk_brier-score_tbl$raw_brier))

# ---- (4) Diebold-Mariano, clustered by game (within-game preds correlated) ---
dm_clustered <- function(df, p_a, p_b) {
  d_game <- df %>% mutate(la=calc_logloss(.data[[p_a]],actual), lb=calc_logloss(.data[[p_b]],actual)) %>%
    group_by(game_id) %>% summarize(d = mean(la - lb), .groups="drop")
  dbar <- mean(d_game$d); se <- sd(d_game$d)/sqrt(nrow(d_game))
  list(dbar=dbar, se=se, t=dbar/se, n_games=nrow(d_game))
}
dm <- dm_clustered(test, "p_raw", "p_shrunk")  # positive t => shrunk better (raw loss higher)
cat("Diebold-Mariano (game-clustered) raw vs shrunk log-loss:\n")
cat(sprintf("  mean per-game loss diff (raw-shrunk)=%+.5f  SE=%.5f  t=%+.2f  (|t|>1.96 sig)  games=%d\n\n",
    dm$dbar, dm$se, dm$t, dm$n_games))

# ---- (5) Calibration / reliability (TEST) ------------------------------------
cal <- function(df, pcol) df %>% mutate(b=cut(.data[[pcol]], seq(0,1,0.1), include.lowest=TRUE)) %>%
  group_by(b) %>% summarize(n=n(), pred=mean(.data[[pcol]]), obs=mean(actual), .groups="drop") %>%
  mutate(cal_err=pred-obs)
cat("CALIBRATION (TEST) — RAW model:\n");   print(cal(test,"p_raw") %>% mutate(across(c(pred,obs,cal_err),~round(.,3))))
cat("\nCALIBRATION (TEST) — SHRUNK model:\n"); print(cal(test,"p_shrunk") %>% mutate(across(c(pred,obs,cal_err),~round(.,3))))
ece <- function(df,pcol){ c<-cal(df,pcol); sum(c$n*abs(c$cal_err))/sum(c$n) }
cat(sprintf("\n  ECE raw=%.4f  shrunk=%.4f  book=%.4f  (lower=better)\n\n", ece(test,"p_raw"), ece(test,"p_shrunk"), ece(test,"p_book")))

# ---- (6) ROI corroboration: raw vs shrunk vs shrunk+BakerMcHale (TEST) -------
grade <- function(df, pcol, use_bm) {
  d <- df %>% mutate(ev=calc_ev(.data[[pcol]], book_odds)) %>% filter(ev > EV_THRESHOLD)
  if (nrow(d)==0) return(data.frame(arm=pcol, n_bets=0, staked=0, pnl=0, roi=NA, win=NA, n_extreme=0))
  al <- if (use_bm) bm_alpha(d[[pcol]] - 1/ifelse(d$book_odds>0,1+d$book_odds/100,1+100/abs(d$book_odds)), shrink_var(d[[pcol]], d$n, N0)) else 1  # prob-space edge (dim-consistent Baker-McHale)
  d <- d %>% mutate(alpha=al,
    stake=kelly_stake_p(.data[[pcol]], book_odds, BANKROLL, KELLY_MULT, alpha),
    dec=ifelse(book_odds>0,1+book_odds/100,1+100/abs(book_odds)),
    pnl=ifelse(actual==1, stake*(dec-1), -stake))
  data.frame(arm=pcol, n_bets=nrow(d), staked=sum(d$stake), pnl=sum(d$pnl),
             roi=sum(d$pnl)/sum(d$stake)*100, win=mean(d$actual)*100, n_extreme=sum(d$is_extreme))
}
roi <- bind_rows(
  cbind(label="RAW (production)",        grade(test,"p_raw",   FALSE)),
  cbind(label="SHRUNK (Kelly plug-in)",  grade(test,"p_shrunk",FALSE)),
  cbind(label="SHRUNK + Baker-McHale",   grade(test,"p_shrunk",TRUE)))
cat("=====================================================================\n")
cat("ROI CORROBORATION (TEST, +EV>5% bets) — high variance, NOT primary\n")
cat("=====================================================================\n")
print(roi %>% mutate(staked=round(staked), pnl=round(pnl), roi=round(roi,2), win=round(win,1)))
cat("\n")

# ---- (7) EXTREME subset: does shrinkage specifically tame phantom edges? -----
ext <- test %>% filter(is_extreme)
if (nrow(ext)>0) {
  cat("EXTREME-SAMPLE SUBSET (is_extreme: |spread|>p99 | total>p99 | final_N<0.5N)\n")
  cat(sprintf("  preds=%d  games=%d\n", nrow(ext), n_distinct(ext$game_id)))
  cat(sprintf("  log-loss: raw=%.4f shrunk=%.4f book=%.4f\n",
      mean(calc_logloss(ext$p_raw,ext$actual)), mean(calc_logloss(ext$p_shrunk,ext$actual)),
      mean(calc_logloss(ext$p_book,ext$actual))))
  re<-grade(ext,"p_raw",FALSE); se<-grade(ext,"p_shrunk",TRUE)
  cat(sprintf("  +EV bets: raw n=%d roi=%.1f%%  | shrunk+BM n=%d roi=%.1f%%\n\n",
      re$n_bets, re$roi, se$n_bets, se$roi))
}

# ---- (7b) ROI vs n0 sweep on TEST (betting objective — complements log-loss) -
roi_for_n0 <- function(n0, use_bm=TRUE) {
  p <- shrink_prob(test$x, test$n, test$book_prob, n0)
  ev <- calc_ev(p, test$book_odds); sel <- ev > EV_THRESHOLD
  if (sum(sel)==0) return(data.frame(n0=n0, n_bets=0, roi=NA, pnl=0, win=NA))
  d <- test[sel,]; pp <- p[sel]
  dec <- ifelse(d$book_odds>0,1+d$book_odds/100,1+100/abs(d$book_odds))
  al <- if (use_bm) bm_alpha(pp - 1/dec, shrink_var(pp,d$n,n0)) else 1  # prob-space edge
  stake <- kelly_stake_p(pp, d$book_odds, BANKROLL, KELLY_MULT, al)
  pnl <- ifelse(d$actual==1, stake*(dec-1), -stake)
  data.frame(n0=n0, n_bets=sum(sel), roi=sum(pnl)/sum(stake)*100, pnl=round(sum(pnl)), win=round(mean(d$actual)*100,1))
}
cat("ROI vs n0 sweep (TEST, shrunk + Baker-McHale):\n")
print(bind_rows(lapply(N0_GRID, roi_for_n0)))
cat("\n")

# ---- (8) Overfit check: TEST log-loss vs n0 (is TRAIN-chosen n0 near TEST opt?)
grid_test <- sapply(N0_GRID, function(n0) mean(calc_logloss(shrink_prob(test$x,test$n,test$book_prob,n0), test$actual)))
cat("OVERFIT CHECK — log-loss vs n0 on TRAIN vs TEST (optimum should be stable):\n")
print(data.frame(n0=N0_GRID, train_ll=round(grid_train,5), test_ll=round(grid_test,5)))
cat(sprintf("  TRAIN-opt n0=%s | TEST-opt n0=%s  (close => not overfit)\n",
    N0_GRID[which.min(grid_train)], N0_GRID[which.min(grid_test)]))
cat("=====================================================================\n")
cat("ANALYSIS COMPLETE\n")
