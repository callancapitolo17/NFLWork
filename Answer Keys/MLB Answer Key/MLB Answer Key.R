library(data.table)
library(oddsapiR)
Sys.setenv(ODDS_API_KEY = "4c490c81e9b826798390440845c89d50")
mlb_hist <- toa_sports_odds_history(
  sport_key   = "baseball_mlb", 
  regions     = "us",
  markets     = "h2h,totals",
  odds_format = "american",
  date_format = "iso",
  date_from   = "2015-01-01",
  date_to     = "2024-12-31"
)
# Helper: de-vig two American-style odds, return data.frame
devig_american <- function(odd1, odd2) {
  p1_raw <- ifelse(odd1 > 0,
                   100 / (odd1 + 100),
                   -odd1 / (-odd1 + 100))
  p2_raw <- ifelse(odd2 > 0,
                   100 / (odd2 + 100),
                   -odd2 / (-odd2 + 100))
  total_raw <- p1_raw + p2_raw
  data.frame(p1 = p1_raw/total_raw, p2 = p2_raw/total_raw)
}

# 2. Fetch historical odds and compute de-vig probabilities
DT <- load_schedules(1999:2024) %>%
  select(game_id, season, week,
         home_score, away_score,
         spread_line, total_line,
         home_spread_odds, away_spread_odds,
         over_odds, under_odds) %>%
  filter(!is.na(spread_line), !is.na(total_line),
         !is.na(home_spread_odds), !is.na(away_spread_odds),
         !is.na(over_odds), !is.na(under_odds)) %>%
  as.data.table()
qs <- dt[, quantile(spread_line, probs = c(0.05, 0.95))]
ss <- qs[2] - qs[1]
qt <- dt[, quantile(total_line,  probs = c(0.05, 0.95))]
st <- qt[2] - qt[1]
# Compute de-vigged probabilities for spread and over
spr <- devig_american(DT$home_spread_odds, DT$away_spread_odds)
DT[, `:=`(home_spread_prob = spr$p1,
          away_spread_prob = spr$p2)]
ou <- devig_american(DT$over_odds, DT$under_odds)
DT[, `:=`(over_prob = ou$p1,
          under_prob = ou$p2)]

distance_index <- function(dt, ps, pt, ss, st) {
  dt[, index := ((spread_line - ps)/ss)^2 + ((total_line - pt)/st)^2]
}

# 1) Mean‐matching loop to select initial N games
library(data.table)

# Function 1: mean‐matching to select initial N games
mean_match <- function(dt, N, parent_spread, parent_total,
                       ss, st, max_iter_mean, tol_mean) {
  dt <- copy(dt)
  adj_spread <- parent_spread
  adj_total  <- parent_total
  for (iter in seq_len(max_iter_mean)) {
    distance_index(dt, adj_spread, adj_total, ss, st)
    setorder(dt, index)
    dt[, included := FALSE]
    dt[1:N, included := TRUE]
    
    # compute current means and errors
    mean_s <- dt[included == TRUE, mean(spread_line)]
    mean_t <- dt[included == TRUE, mean(total_line)]
    err_s  <- mean_s   - parent_spread
    err_t  <- mean_t   - parent_total
    
    # stop if within tolerance
    if (abs(err_s) < tol_mean && abs(err_t) < tol_mean) break
    
    # adjust the parent targets
    adj_spread <- adj_spread - err_s
    adj_total  <- adj_total  - err_t
  }
  # return updated dt and targets
  list(dt = dt,
       parent_spread = parent_spread,
       parent_total  = parent_total)
}

# Function 2: balance sample by greedy remove/add + shrink‐on‐stall
balance_sample <- function(dt, N, target_cover, target_over, tol_error) {
  dt <- copy(dt)
  # initialize errors
  cover_error <- dt[included == TRUE, sum(actual_cover)] -
    round(target_cover * N)
  over_error  <- dt[included == TRUE, sum(actual_over )] -
    round(target_over  * N)
  M <- nrow(dt)
  
  repeat {
    # stop if within tolerance
    if (abs(cover_error) <= tol_error && abs(over_error) <= tol_error) break
    
    # mark failures
    removal_failed  <- TRUE
    addition_failed <- TRUE
    
    # ---- remove one game from 1..N ----
    in1toN    <- dt[1:N,   included]
    cov_help1 <- (cover_error > 0 & dt[1:N, actual_cover] == 1) |
      (cover_error < 0 & dt[1:N, actual_cover] == 0)
    ov_help1  <- (over_error  > 0 & dt[1:N, actual_over ] == 1) |
      (over_error  < 0 & dt[1:N, actual_over ] == 0)
    
    both1   <- which(in1toN & cov_help1 & ov_help1)
    either1 <- which(in1toN & (cov_help1 | ov_help1))
    
    if (length(both1) > 0) {
      i <- max(both1)
      dt[i, included := FALSE]
      removal_failed <- FALSE
    } else if (length(either1) > 0) {
      i <- max(either1)
      dt[i, included := FALSE]
      removal_failed <- FALSE
    }
    
    if (!removal_failed) {
      cover_error <- dt[included == TRUE, sum(actual_cover)] -
        round(target_cover * N)
      over_error  <- dt[included == TRUE, sum(actual_over )] -
        round(target_over  * N)
    }
    
    # ---- add one game from (N+1)..M ----
    inNplus   <- dt[(N+1):M, included]
    cov_help2 <- (cover_error > 0 & dt[(N+1):M, actual_cover] == 0) |
      (cover_error < 0 & dt[(N+1):M, actual_cover] == 1)
    ov_help2  <- (over_error  > 0 & dt[(N+1):M, actual_over ] == 0) |
      (over_error  < 0 & dt[(N+1):M, actual_over ] == 1)
    
    both2   <- which(!inNplus & cov_help2 & ov_help2)
    either2 <- which(!inNplus & (cov_help2 | ov_help2))
    
    if (length(both2) > 0) {
      j <- max(both2) + N
      dt[j, included := TRUE]
      addition_failed <- FALSE
    } else if (length(either2) > 0) {
      j <- max(either2) + N
      dt[j, included := TRUE]
      addition_failed <- FALSE
    }
    
    if (!addition_failed) {
      cover_error <- dt[included == TRUE, sum(actual_cover)] -
        round(target_cover * N)
      over_error  <- dt[included == TRUE, sum(actual_over )] -
        round(target_over  * N)
    }
    
    # ---- if neither helped, shrink sample size ----
    if (removal_failed && addition_failed) {
      worst_i <- dt[1:N, which.max(index)]
      dt[worst_i, included := FALSE]
      N <- N - 1
      cover_error <- dt[included == TRUE, sum(actual_cover)] -
        round(target_cover * N)
      over_error  <- dt[included == TRUE, sum(actual_over )] -
        round(target_over  * N)
      next
    }
  }
  
  # return final dt and stats
  list(dt          = dt,
       final_N     = N,
       cover_error = cover_error,
       over_error  = over_error)
}


parent_spread <- -8.5
parent_total  <- 49.0
home_sp_odds <- +110; away_sp_odds <- -130
over_odds    <- -110; under_odds    <- +100
N <- round(nrow(dt)*0.025,0)

matched_mean<- mean_match(dt,N,parent_spread,parent_total,ss,st,max_iter_mean = 500, tol_mean = 0.005)
mean(matched_mean$dt[included == TRUE,total_line])

desired_over_prob<- devig_american(over_odds,under_odds)[[1]]
desired_cover_prob<- devig_american(home_sp_odds,away_sp_odds)[[1]]
final_sample <- balance_sample(matched_mean$dt,N,desired_cover_prob,desired_over_prob,tol_error = 1)
mean(matched_mean$dt[included == TRUE,total_line])

dt_final   <- final_sample$dt
final_N    <- final_sample$final_N
cover_err  <- final_sample$cover_error
over_err   <- final_sample$over_error

# 2) Compute final metrics
final_mean_s  <- dt_final[included == TRUE, mean(spread_line)]
final_mean_t  <- dt_final[included == TRUE, mean(total_line)]
final_cover_r <- dt_final[included == TRUE, sum(actual_cover)] / final_N
final_over_r  <- dt_final[included == TRUE, sum(actual_over )] / final_N

# 3) Print everything alongside your targets
cat("===== Sanity Check =====\n")
cat(sprintf("Sample size:      %d  (target was %d)\n", final_N, N))
cat(sprintf("Spread mean:      %.3f  (target %.3f; err %.3f)\n",
            final_mean_s, parent_spread, final_mean_s - parent_spread))
cat(sprintf("Total mean:       %.3f  (target %.3f; err %.3f)\n",
            final_mean_t, parent_total, final_mean_t - parent_total))
cat(sprintf("Cover rate:       %.3f  (target %.3f; err %d games)\n",
            final_cover_r, desired_cover_prob, cover_err))
cat(sprintf("Over rate:        %.3f  (target %.3f; err %d games)\n",
            final_over_r, desired_over_prob, over_err))


