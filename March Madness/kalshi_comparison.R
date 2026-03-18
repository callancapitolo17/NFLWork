# =============================================================================
# Kalshi March Madness Comparison
# =============================================================================
# Compares simulator probabilities to Kalshi market prices.
# Identifies +EV opportunities across championship, advancement, and seed props.
#
# Run after props_analysis.R or standalone.

source("shared.R")
source("espn_bracket.R")
library(httr)
library(jsonlite)

# =============================================================================
# 1. Run simulation
# =============================================================================
cat("Loading data and running simulation...\n")
br <- fetch_espn_bracket()
ts <- get_teams_std()
fb <- br$bracket %>% select(team, seed, region, play_in)
bwr <- fetch_bracket_with_ratings(fb, ts)
bracket_64 <- as.data.frame(resolve_first_four(bwr, br$games))
region_order <- get_region_order(bracket_64)

n_sims <- 10000
t0 <- Sys.time()
raw_sims <- map_dfr(1:n_sims, function(i) {
  result <- simulate_tournament_fast(bracket_64, region_order)
  result$sim_id <- i
  result
})
cat(sprintf("Sim done in %.1f sec.\n\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# Team-level results
team_probs <- raw_sims %>%
  group_by(team, seed) %>%
  summarise(
    R32 = mean(Round_32), S16 = mean(Sweet_16), E8 = mean(Elite_8),
    F4 = mean(Final_4), TG = mean(Title_Game), Champ = mean(Champion),
    .groups = "drop"
  )

# Seed-level upset props
seed_upset_props <- raw_sims %>%
  filter(seed >= 11) %>%
  group_by(sim_id, seed) %>%
  summarise(n_r32_wins = sum(Round_32), .groups = "drop") %>%
  group_by(seed) %>%
  summarise(p_at_least_one = mean(n_r32_wins >= 1), .groups = "drop")

# =============================================================================
# 2. Fetch Kalshi markets
# =============================================================================
cat("Fetching Kalshi markets...\n")

fetch_kalshi_event <- function(event_ticker) {
  resp <- GET(sprintf("https://api.elections.kalshi.com/trade-api/v2/markets?event_ticker=%s&limit=200", event_ticker),
              add_headers("Accept" = "application/json"))
  data <- fromJSON(content(resp, "text", encoding = "UTF-8"))
  if (is.data.frame(data$markets) && nrow(data$markets) > 0) {
    data$markets %>%
      mutate(
        yes_ask = as.numeric(yes_ask_dollars),
        yes_bid = as.numeric(yes_bid_dollars),
        last_price = as.numeric(last_price_dollars)
      ) %>%
      select(ticker, title, subtitle, yes_ask, yes_bid, last_price)
  } else {
    tibble()
  }
}

# Fetch all round-specific markets
kalshi <- bind_rows(
  fetch_kalshi_event("KXMARMAD-26") %>% mutate(round = "Champ"),
  fetch_kalshi_event("KXMARMADROUND-26T2") %>% mutate(round = "TG"),
  fetch_kalshi_event("KXMARMADROUND-26F4") %>% mutate(round = "F4"),
  fetch_kalshi_event("KXMARMADROUND-26E8") %>% mutate(round = "E8"),
  fetch_kalshi_event("KXMARMADROUND-26S16") %>% mutate(round = "S16"),
  fetch_kalshi_event("KXMARMADROUND-26RO32") %>% mutate(round = "R32")
)

# Also get seed props
seed_events <- GET("https://api.elections.kalshi.com/trade-api/v2/events?series_ticker=KXMARMADSEEDWIN&status=open&limit=50",
                   add_headers("Accept" = "application/json"))
seed_evts <- fromJSON(content(seed_events, "text", encoding = "UTF-8"))$events
if (is.data.frame(seed_evts) && nrow(seed_evts) > 0) {
  seed_markets <- map_dfr(seed_evts$event_ticker, function(et) {
    fetch_kalshi_event(et) %>% mutate(round = "seed_prop")
  })
  kalshi <- bind_rows(kalshi, seed_markets)
}

cat(sprintf("Kalshi markets: %d\n\n", nrow(kalshi)))

# =============================================================================
# 3. Match Kalshi markets to sim probabilities
# =============================================================================

# Extract team name from Kalshi title
extract_team <- function(title) {
  # "Will Duke win the..." -> "Duke"
  # "Will Michigan St. qualify for..." -> "Michigan St."
  team <- str_replace(title, "^Will\\s+", "")
  team <- str_replace(team, "\\s+(win|qualify|make|advance|beat).*$", "")
  trimws(team)
}

kalshi <- kalshi %>%
  filter(!is.na(yes_ask), yes_ask > 0, yes_ask < 1) %>%
  mutate(kalshi_team = extract_team(title))

# Map Kalshi team names to our bracket names
kalshi <- kalshi %>%
  mutate(standard_team = map_chr(kalshi_team, ~ get_standard_team(.x, teams_std = ts)))

# Match with sim probabilities
comparison <- kalshi %>%
  left_join(team_probs, by = c("standard_team" = "team")) %>%
  mutate(
    sim_prob = case_when(
      round == "Champ" ~ Champ,
      round == "TG" ~ TG,
      round == "F4" ~ F4,
      round == "E8" ~ E8,
      round == "S16" ~ S16,
      round == "R32" ~ R32,
      TRUE ~ NA_real_
    ),
    kalshi_prob = yes_ask,
    edge = sim_prob - kalshi_prob,
    ev_per_dollar = ifelse(!is.na(sim_prob) & sim_prob > 0,
                           (sim_prob / kalshi_prob - 1) * 100, NA)
  ) %>%
  filter(!is.na(sim_prob), !is.na(kalshi_prob))

# =============================================================================
# 4. Output: Best edges
# =============================================================================

cat("=== BEST YES EDGES (sim > Kalshi price) ===\n\n")
yes_edges <- comparison %>%
  filter(edge > 0) %>%
  arrange(desc(ev_per_dollar)) %>%
  head(25)

cat(sprintf("%-45s %5s %6s %6s %6s %6s\n", "Market", "Round", "Sim", "Kalshi", "Edge", "EV%"))
cat(paste(rep("-", 90), collapse = ""), "\n")
for (i in 1:nrow(yes_edges)) {
  r <- yes_edges[i, ]
  cat(sprintf("%-45s %5s %5.1f%% %5.1f%% %+5.1f%% %+5.0f%%\n",
              substr(r$title, 1, 45), r$round,
              r$sim_prob * 100, r$kalshi_prob * 100,
              r$edge * 100, r$ev_per_dollar))
}

cat("\n=== BEST NO EDGES (sim < Kalshi price) ===\n\n")
no_edges <- comparison %>%
  filter(edge < -0.03) %>%
  mutate(
    no_sim = 1 - sim_prob,
    no_kalshi = 1 - kalshi_prob,
    no_edge = no_sim - no_kalshi,
    no_ev = (no_sim / no_kalshi - 1) * 100
  ) %>%
  arrange(desc(no_ev)) %>%
  head(25)

cat(sprintf("%-45s %5s %6s %6s %6s %6s\n", "Market (NO side)", "Round", "NoSim", "NoKal", "Edge", "EV%"))
cat(paste(rep("-", 90), collapse = ""), "\n")
for (i in 1:nrow(no_edges)) {
  r <- no_edges[i, ]
  cat(sprintf("%-45s %5s %5.1f%% %5.1f%% %+5.1f%% %+5.0f%%\n",
              substr(r$title, 1, 45), r$round,
              r$no_sim * 100, r$no_kalshi * 100,
              r$no_edge * 100, r$no_ev))
}

# Seed props comparison
cat("\n=== SEED UPSET PROPS ===\n\n")
seed_props <- kalshi %>% filter(round == "seed_prop")
if (nrow(seed_props) > 0) {
  for (i in 1:nrow(seed_props)) {
    sp <- seed_props[i, ]
    cat(sprintf("  %-60s kalshi=%.0f¢\n", substr(sp$title, 1, 60), sp$yes_ask * 100))
    # Try to match with our upset props
    if (grepl("#16", sp$title)) {
      our_prob <- seed_upset_props$p_at_least_one[seed_upset_props$seed == 16]
      if (length(our_prob) > 0) cat(sprintf("    Our sim: %.1f%% | Edge: %+.1f%%\n", our_prob * 100, (our_prob - sp$yes_ask) * 100))
    }
  }
}

cat("\n=== SUMMARY ===\n")
cat(sprintf("Total markets compared: %d\n", nrow(comparison)))
cat(sprintf("YES edges (sim > kalshi): %d\n", sum(comparison$edge > 0.02)))
cat(sprintf("NO edges (kalshi > sim): %d\n", sum(comparison$edge < -0.03)))
cat(sprintf("Within 2%%: %d\n", sum(abs(comparison$edge) <= 0.02)))
