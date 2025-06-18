simulate_series_distribution <- function(
    home_advantage,              # e.g. "Grizzlies"
    challenger,                  # e.g. "Thunder"
    p_home_win_at_home = NULL,          # market P(home_advantage wins G1 at home)
    p_challenger_win_at_home = NULL,  # optional override
    HCA = 0.088,                 # league home‑court edge (NBA ≈0.088)
    pattern = c(2,2,1,1,1),      # 2‑2‑1‑1‑1 for best‑of‑7
    max_wins = 4,                # wins needed to clinch
    n_iter = 1e5,
    current_home_adv_wins = 0,
    current_challenger_wins = 0
    # Monte Carlo iterations
) {
  # 1) back out challenger’s home‑win % if not provided
  if (is.null(p_challenger_win_at_home)) {
    neutral_home_p      <- p_home_win_at_home - HCA
    neutral_challenger_p <- 1 - neutral_home_p
    p_challenger_win_at_home <- neutral_challenger_p + HCA
  }
  if (is.null(p_home_win_at_home)) {
    neutral_challenger_p      <- p_challenger_win_at_home - HCA
    neutral_home_p <- 1 - neutral_challenger_p
    p_home_win_at_home <- neutral_home_p + HCA
  }
  
  # 2) helpers
  prob_to_american <- function(p) {
    if (p >= 0.5) -round((p/(1-p)) * 100) else round(((1-p)/p) * 100)
  }
  
  # 3) build who’s at home each game
  hosts      <- rep(c(home_advantage, challenger), length.out = length(pattern))
  home_order <- rep(hosts, times = pattern)[(current_home_adv_wins + current_challenger_wins+1):sum(pattern)]
  
  # 4) simulate one series
  single_run <- function() {
    hw <- current_home_adv_wins; cw <- current_challenger_wins
    for (loc in home_order) {
      if (hw == max_wins || cw == max_wins) break
      if (loc == home_advantage) {
        if (runif(1) < p_home_win_at_home) hw <- hw + 1L else cw <- cw + 1L
      } else {
        if (runif(1) < p_challenger_win_at_home) cw <- cw + 1L else hw <- hw + 1L
      }
    }
    c(hw, cw)
  }
  
  # 5) run Monte Carlo
  res_mat <- replicate(n_iter, single_run())
  hws <- res_mat[1,]; cws <- res_mat[2,]
  
  # 6a) series‑win table
  series_tbl <- data.frame(
    team  = c(home_advantage, challenger),
    count = c(sum(hws == max_wins), sum(cws == max_wins)),
    stringsAsFactors = FALSE
  )
  series_tbl$prob          <- series_tbl$count / n_iter
  series_tbl$american_odds <- vapply(series_tbl$prob, prob_to_american, numeric(1))
  
  # 6b) exact‑outcome table
  outcome_str <- ifelse(
    hws == max_wins,
    paste0(home_advantage, " ", max_wins, "–", challenger, " ", cws),
    paste0(challenger,     " ", max_wins, "–", home_advantage, " ", hws)
  )
  exact_tbl <- as.data.frame(table(outcome_str), stringsAsFactors = FALSE)
  names(exact_tbl) <- c("outcome", "count")
  exact_tbl$prob          <- exact_tbl$count / n_iter
  exact_tbl$american_odds <- vapply(exact_tbl$prob, prob_to_american, numeric(1))
  
  # 6c) series‑length distribution
  lengths   <- hws + cws
  length_tbl <- as.data.frame(table(lengths), stringsAsFactors = FALSE)
  names(length_tbl) <- c("length", "count")
  length_tbl$length       <- as.integer(length_tbl$length)
  length_tbl$prob         <- length_tbl$count / n_iter
  length_tbl$american_odds<- vapply(length_tbl$prob, prob_to_american, numeric(1))
  
  # 6d) spread (margin) distribution
  diffs    <- hws - cws
  spread0  <- as.data.frame(table(diffs), stringsAsFactors = FALSE)
  names(spread0) <- c("diff", "count")
  spread0$diff   <- as.integer(spread0$diff)
  spread0$margin <- abs(spread0$diff)
  spread0$team   <- ifelse(spread0$diff > 0, home_advantage, challenger)
  spread0$spread <- paste0(spread0$team, " by ", spread0$margin)
  spread0$prob   <- spread0$count / n_iter
  spread0$american_odds <- vapply(spread0$prob, prob_to_american, numeric(1))
  spread_tbl <- spread0[order(spread0$team, -spread0$margin),
                        c("spread", "team", "margin", "prob", "american_odds")]
  
  # 6e) auto‑generate every margin prop (1.5, 2.5, …)
  margin_thresholds <- seq(1.5, max_wins - 0.5, by = 1)
  margin_props <- do.call(rbind, lapply(margin_thresholds, function(th) {
    min_marg <- ceiling(th)
    df_home <- data.frame(
      prop = paste0(home_advantage, " -", th),
      prob = sum(spread0$prob[spread0$team == home_advantage &
                                spread0$margin >= min_marg]),
      stringsAsFactors = FALSE
    )
    df_away <- data.frame(
      prop = paste0(challenger, " -", th),
      prob = sum(spread0$prob[spread0$team == challenger &
                                spread0$margin >= min_marg]),
      stringsAsFactors = FALSE
    )
    rbind(df_home, df_away)
  }))
  margin_props$american_odds <- vapply(margin_props$prob, prob_to_american, numeric(1))
  
  # 6f) auto‑generate total‑games props (4.5, 5.5, …)
  total_thresholds <- seq(max_wins + 0.5, sum(pattern) - 0.5, by = 1)
  total_props <- do.call(rbind, lapply(total_thresholds, function(th) {
    under_max <- floor(th)
    over_min  <- ceiling(th)
    df_under <- data.frame(
      prop = paste0("Under ", th, " games"),
      prob = sum(length_tbl$prob[length_tbl$length <= under_max]),
      stringsAsFactors = FALSE
    )
    df_over <- data.frame(
      prop = paste0("Over ", th, " games"),
      prob = sum(length_tbl$prob[length_tbl$length >= over_min]),
      stringsAsFactors = FALSE
    )
    rbind(df_under, df_over)
  }))
  total_props$american_odds <- vapply(total_props$prob, prob_to_american, numeric(1))
  
  # 7) return all tables
  list(
    series       = series_tbl,
    exact        = exact_tbl[order(-exact_tbl$prob), ],
    length       = length_tbl[order(length_tbl$length), ],
    spread       = spread_tbl,
    margin_props = margin_props,
    total_props  = total_props
  )
}

simulate_series_distribution(
  home_advantage         = "Lakers",
  challenger             = "Wolves",
  p_challenger_win_at_home  =  0.3,
  n_iter                 = 1e5,
  current_home_adv_wins = 1,
  current_challenger_wins = 3,
  #HCA = 0.03
  
)
