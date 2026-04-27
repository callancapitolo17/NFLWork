# conditional_kelly.R — Kelly residuals for the singles when a combined parlay
# is already placed. Used by the MLB Dashboard's Parlay tab to update the
# source rows' Kelly column after a combined ticket is locked in.
#
# Approach: for two independent legs A and B and a fixed combo stake s_combo,
# numerically maximize E[log(1 + total return)] over (s_a, s_b) using the
# 4-outcome joint distribution. Reuses the same Kelly-portfolio framing as
# parlay_multivariate_kelly() but with one position held fixed.

conditional_kelly_residuals <- function(p_a, d_a,
                                         p_b, d_b,
                                         s_combo, d_combo,
                                         bankroll, kelly_mult = 1.0) {
  # Input guards — degenerate cases return zero residuals rather than crash
  if (!is.finite(bankroll) || bankroll <= 0) return(list(s_a = 0, s_b = 0))
  if (s_combo < 0 || s_combo >= bankroll)    return(list(s_a = 0, s_b = 0))

  # Per-bet returns: profit per $1 wagered if win; -$1 if lose
  b_a <- d_a - 1
  b_b <- d_b - 1
  b_c <- d_combo - 1

  # Combo as bankroll fraction (held fixed during optimization)
  f_c <- s_combo / bankroll

  # Joint probabilities under cross-game independence
  p_ab      <- p_a * p_b
  p_a_only  <- p_a * (1 - p_b)
  p_b_only  <- (1 - p_a) * p_b
  p_neither <- (1 - p_a) * (1 - p_b)

  # Clamp returns at -0.999999 so log1p never sees a value <= -1
  # (the optimizer probes corners where total exposure can exceed 100%)
  safe_log1p <- function(r) log1p(pmax(r, -0.999999))

  # Negative expected log return — minimize this = maximize log growth
  neg_e_log <- function(f) {
    f_a <- f[1]; f_b <- f[2]
    r_both    <- f_a * b_a + f_b * b_b + f_c * b_c
    r_a_only  <- f_a * b_a - f_b           - f_c
    r_b_only  <- -f_a       + f_b * b_b    - f_c
    r_neither <- -f_a       - f_b           - f_c
    -(p_ab      * safe_log1p(r_both)    +
      p_a_only  * safe_log1p(r_a_only)  +
      p_b_only  * safe_log1p(r_b_only)  +
      p_neither * safe_log1p(r_neither))
  }

  result <- optim(
    par    = c(0.01, 0.01),
    fn     = neg_e_log,
    method = "L-BFGS-B",
    lower  = c(0, 0),
    upper  = c(0.5, 0.5)  # safe upper bound — Kelly never recommends > 50% bankroll
  )

  # IMPORTANT: convergence guard — optim returns 0 on success, non-zero on failure
  if (result$convergence != 0) {
    warning("conditional_kelly_residuals: optimizer did not converge (code ",
            result$convergence, "). Returning zero residuals.")
    return(list(s_a = 0, s_b = 0))
  }

  s_a <- result$par[1] * kelly_mult * bankroll
  s_b <- result$par[2] * kelly_mult * bankroll

  list(s_a = round(max(s_a, 0), 2),
       s_b = round(max(s_b, 0), 2))
}
