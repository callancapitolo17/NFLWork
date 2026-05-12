#!/usr/bin/env Rscript
# Extract sample-based correlation factor for FG fav+over combos,
# for the games we placed bets on, then compare to realized outcomes.

suppressPackageStartupMessages({
  library(duckdb)
  library(dplyr)
  library(jsonlite)
})

ARGV <- commandArgs(trailingOnly = TRUE)
OUT <- if (length(ARGV) >= 1) ARGV[1] else "/tmp/fg_fav_leak/sample_corr.json"

con <- dbConnect(duckdb(), dbdir = "/Users/callancapitolo/NFLWork/Answer Keys/mlb_mm.duckdb", read_only = TRUE)
samples <- tbl(con, "mlb_game_samples") %>% collect()
dbDisconnect(con)

if (!"home_margin" %in% names(samples)) {
  stop("mlb_game_samples missing home_margin column")
}
if (!"total_final_score" %in% names(samples)) {
  stop("mlb_game_samples missing total_final_score column")
}

# Sample-based per-game correlation factor for fav+over.
# Conditioning lines: spread -1.5 (home covers when margin > 1.5; away covers when margin < -1.5)
# Test multiple total lines to mirror the H3 stratification: 6.5, 7.5, 8.5, 9.5
per_game_lines <- expand.grid(
  game_id = unique(samples$game_id),
  total_line = c(6.5, 7.5, 8.5, 9.5)
)

compute_corr <- function(s, total_line) {
  # Returns a list of metrics for one (game_id, total_line) combination
  p_home_cov <- mean(s$home_margin > 1.5, na.rm = TRUE)
  p_away_cov <- mean(s$home_margin < -1.5, na.rm = TRUE)
  p_over     <- mean(s$total_final_score > total_line, na.rm = TRUE)
  p_home_AND_over <- mean(s$home_margin > 1.5 & s$total_final_score > total_line, na.rm = TRUE)
  p_away_AND_over <- mean(s$home_margin < -1.5 & s$total_final_score > total_line, na.rm = TRUE)
  list(
    p_home_cov = p_home_cov,
    p_away_cov = p_away_cov,
    p_over = p_over,
    joint_home = p_home_AND_over,
    joint_away = p_away_AND_over,
    corr_home = if (p_home_cov > 0 && p_over > 0) p_home_AND_over / (p_home_cov * p_over) else NA,
    corr_away = if (p_away_cov > 0 && p_over > 0) p_away_AND_over / (p_away_cov * p_over) else NA
  )
}

results <- list()
for (gid in unique(samples$game_id)) {
  s <- samples %>% filter(game_id == gid)
  for (tl in c(6.5, 7.5, 8.5, 9.5)) {
    m <- compute_corr(s, tl)
    results[[length(results) + 1]] <- c(list(game_id = gid, total_line = tl), m)
  }
}
df_corr <- bind_rows(lapply(results, function(x) as.data.frame(x, stringsAsFactors = FALSE)))

# Summary stats per total_line
summary_by_total <- df_corr %>%
  group_by(total_line) %>%
  summarise(
    n_games = n(),
    mean_corr_home = mean(corr_home, na.rm = TRUE),
    median_corr_home = median(corr_home, na.rm = TRUE),
    mean_corr_away = mean(corr_away, na.rm = TRUE),
    median_corr_away = median(corr_away, na.rm = TRUE),
    mean_joint_home = mean(joint_home, na.rm = TRUE),
    mean_joint_away = mean(joint_away, na.rm = TRUE)
  ) %>%
  as.data.frame()

out <- list(
  n_games = length(unique(samples$game_id)),
  n_samples_per_game = nrow(samples) / length(unique(samples$game_id)),
  by_total_line = lapply(seq_len(nrow(summary_by_total)), function(i) as.list(summary_by_total[i, ]))
)

writeLines(toJSON(out, auto_unbox = TRUE, pretty = TRUE), OUT)
cat("Wrote", OUT, "\n")
