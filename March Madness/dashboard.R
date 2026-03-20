# =============================================================================
# March Madness Simulator Dashboard
# =============================================================================
# Interactive query tool for simulation results.
# Run: Rscript dashboard.R (opens in browser)
#
# Answers questions like:
#   - How many wins does the SEC have?
#   - What's P(1-seed wins championship)?
#   - What's P(exactly 2 #1-seeds in Final Four)?
#   - What's P(a 12-seed beats a 5-seed)?
#   - Conference win totals, seed performance, team advancement

library(shiny)
library(dplyr)
library(DT)

# =============================================================================
# 1. Run simulation and prepare data
# =============================================================================

cat("Loading simulation data...\n")
source("shared.R")
source("espn_bracket.R")

br <- fetch_espn_bracket()
ts <- get_teams_std()
fb <- br$bracket %>% select(team, seed, region, play_in)
bwr <- fetch_bracket_with_ratings(fb, ts)
bracket_64 <- as.data.frame(resolve_first_four(bwr, br$games))
region_order <- get_region_order(bracket_64)

# Filter out eliminated teams (respect completed game results)
games_played <- br$games
eliminated <- games_played %>%
  filter(status == "final", round != "First Four", !is.na(winner)) %>%
  mutate(loser = ifelse(winner == team1, team2, team1)) %>%
  pull(loser)

current_bracket <- bracket_64 %>%
  filter(!team %in% eliminated) %>%
  mutate(status = "pending")

n_remaining <- nrow(current_bracket)
tourney_state <- br$tournament_state
cat(sprintf("Tournament: %s | Round: %s | %d teams remaining\n",
            tourney_state$state, tourney_state$current_round %||% "N/A", n_remaining))

# Add conference
conf_map <- tryCatch({
  cbd_bpi_ratings() %>%
    transmute(team_short = team, conference = conf) %>%
    mutate(team = map_chr(team_short, ~ get_standard_team(.x, teams_std = ts))) %>%
    select(team, conference) %>% distinct(team, .keep_all = TRUE)
}, error = function(e) NULL)

# Use Dynamic Simulator approach: handles partial rounds, auto-advances byes
simulate_round_dynamic <- function(teams, game_number = 1, region_order_auto = NULL) {
  teams_ordered <- get_bracket_matchups(teams, region_order_auto)
  team_pairs <- split(teams_ordered, rep(1:(nrow(teams_ordered)/2), each = 2))
  winners <- map(team_pairs, function(pair) {
    if (all(pair$status == "advanced")) { pair[1, ] }
    else if (any(pair$status == "advanced") && any(pair$status == "pending")) {
      pair %>% filter(status == "advanced") %>% slice(1)
    } else { simulate_game(pair[1, ], pair[2, ], game_number)$winner }
  })
  bind_rows(winners)
}

simulate_remaining <- function(bracket, region_order_auto = NULL) {
  if (is.null(region_order_auto)) region_order_auto <- get_region_order(bracket)
  round_names <- get_remaining_rounds(nrow(bracket))
  team_progress <- bracket %>% select(team, seed)
  for (r in round_names) team_progress[[r]] <- 0
  round_num <- 1; teams_round <- bracket
  while (nrow(teams_round) > 1) {
    teams_round <- simulate_round_dynamic(teams_round, game_number = round_num, region_order_auto = region_order_auto)
    if (round_num <= length(round_names)) {
      rn <- round_names[round_num]
      team_progress <- team_progress %>%
        mutate("{rn}" := ifelse(team %in% teams_round$team, 1, .data[[rn]]))
    }
    round_num <- round_num + 1
  }
  team_progress
}

# Pad missing round columns so all sims have consistent columns
pad_round_cols <- function(df) {
  all_rounds <- c("Round_32", "Sweet_16", "Elite_8", "Final_4", "Title_Game", "Champion")
  # Map get_remaining_rounds names to standard column names
  name_map <- c("Round 32"="Round_32", "Sweet 16"="Sweet_16", "Elite 8"="Elite_8",
                "Final 4"="Final_4", "Title Game"="Title_Game", "Champion"="Champion")
  # Rename any columns that use space format
  for (old_name in names(name_map)) {
    new_name <- name_map[old_name]
    if (old_name %in% names(df) && !new_name %in% names(df)) {
      names(df)[names(df) == old_name] <- new_name
    }
  }
  # Add missing columns as 1 (team already passed those rounds)
  for (r in all_rounds) {
    if (!r %in% names(df)) df[[r]] <- 1L
  }
  df
}

n_sims <- 50000
cat(sprintf("Running %dk simulations (%d teams)...\n", n_sims/1000, n_remaining))
t0 <- Sys.time()
raw <- map_dfr(1:n_sims, function(i) {
  r <- simulate_remaining(current_bracket, region_order)
  r <- pad_round_cols(r)
  r$region <- current_bracket$region[match(r$team, current_bracket$team)]
  r$sim_id <- i
  r
})
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("Done in %.0f seconds.\n", elapsed))

# Add conference
if (!is.null(conf_map)) {
  raw <- raw %>% left_join(conf_map, by = "team")
}

# Pre-compute totals per sim
raw <- raw %>% mutate(
  total_wins = Round_32 + Sweet_16 + Elite_8 + Final_4 + Title_Game + Champion
)

cat("Dashboard ready.\n")

# =============================================================================
# 2. Kalshi Edge Finder
# =============================================================================

KALSHI_BASE <- "https://api.elections.kalshi.com/trade-api/v2"
TAKER_FEE_RATE <- 0.07

kalshi_request <- function(path) {
  tryCatch(fromJSON(paste0(KALSHI_BASE, path), simplifyDataFrame = FALSE),
           error = function(e) { warning(sprintf("Kalshi API: %s", e$message)); NULL })
}

fetch_kalshi_markets <- function(series_ticker = NULL, event_ticker = NULL) {
  markets <- list(); cursor <- NULL
  repeat {
    path <- "/markets?"
    if (!is.null(series_ticker)) path <- paste0(path, "series_ticker=", series_ticker, "&")
    if (!is.null(event_ticker)) path <- paste0(path, "event_ticker=", event_ticker, "&")
    path <- paste0(path, "status=open&limit=200")
    if (!is.null(cursor) && nchar(cursor) > 0) path <- paste0(path, "&cursor=", cursor)
    data <- kalshi_request(path)
    if (is.null(data)) break
    batch <- data$markets
    if (length(batch) == 0) break
    markets <- c(markets, batch)
    cursor <- data$cursor
    if (is.null(cursor) || nchar(cursor) == 0) break
    Sys.sleep(0.1)
  }
  markets
}

kalshi_fee_cents <- function(price_cents) {
  p <- price_cents / 100
  TAKER_FEE_RATE * p * (1 - p) * 100
}

kalshi_compute_ev <- function(sim_prob, yes_ask_cents, yes_bid_cents) {
  fair <- sim_prob * 100; no_ask <- 100 - yes_bid_cents
  yes_fee <- kalshi_fee_cents(yes_ask_cents)
  yes_edge <- fair - yes_ask_cents - yes_fee
  yes_ev <- if (yes_ask_cents > 0) yes_edge / yes_ask_cents else -Inf
  no_fee <- kalshi_fee_cents(no_ask)
  no_edge <- (100 - fair) - no_ask - no_fee
  no_ev <- if (no_ask > 0) no_edge / no_ask else -Inf
  list(yes_ev = yes_ev, no_ev = no_ev, yes_fee = yes_fee, no_fee = no_fee,
       best_side = if (yes_ev >= no_ev) "YES" else "NO",
       best_ev = max(yes_ev, no_ev))
}

# Team name matching: Kalshi title → sim team name
match_kalshi_team <- function(title, team_probs) {
  m <- str_match(title, "^Will (.+?) qualify for")
  if (is.na(m[1,2])) return(NA_character_)
  kn <- tolower(str_replace_all(m[1,2], "['\\.\\-]", ""))
  kn <- str_replace_all(kn, "\\bst\\.?\\b", "st")
  kn <- str_squish(kn)
  cn <- tolower(str_replace_all(team_probs$team, "['\\.\\-]", ""))
  cn <- str_replace_all(cn, "\\bst\\.?\\b", "st")
  cn <- str_squish(cn)
  idx <- which(cn == kn)
  if (length(idx) == 1) return(team_probs$team[idx])
  idx <- which(str_detect(cn, fixed(kn)) | str_detect(kn, cn))
  if (length(idx) == 1) return(team_probs$team[idx])
  NA_character_
}

fetch_all_kalshi_edges <- function(raw, team_probs) {
  cat("Fetching Kalshi tournament props...\n")
  library(stringr)
  all_edges <- tibble()

  # Detect rows per sim
  first_team <- raw$team[1]
  n_per_sim <- which(raw$team == first_team)[2] - 1L
  n_sims_local <- nrow(raw) / n_per_sim
  raw$sim_id <- rep(1:n_sims_local, each = n_per_sim)

  # --- Round Advancement ---
  round_events <- list(
    list(event="KXMARMADROUND-26RO32", col="Round_32", label="Round 32"),
    list(event="KXMARMADROUND-26S16", col="Sweet_16", label="Sweet 16"),
    list(event="KXMARMADROUND-26E8", col="Elite_8", label="Elite 8"),
    list(event="KXMARMADROUND-26F4", col="Final_4", label="Final 4"),
    list(event="KXMARMADROUND-26T2", col="Title_Game", label="Title Game"))

  for (re in round_events) {
    markets <- fetch_kalshi_markets(event_ticker = re$event)
    cat(sprintf("  %s: %d markets\n", re$label, length(markets)))
    for (m in markets) {
      sim_team <- match_kalshi_team(m$title %||% "", team_probs)
      if (is.na(sim_team)) next
      sim_prob <- team_probs[[re$col]][team_probs$team == sim_team]
      if (length(sim_prob) == 0) next
      ya <- round(as.numeric(m$yes_ask_dollars %||% 0) * 100)
      yb <- round(as.numeric(m$yes_bid_dollars %||% 0) * 100)
      ev <- kalshi_compute_ev(sim_prob, ya, yb)
      all_edges <- bind_rows(all_edges, tibble(
        ticker=m$ticker %||% "", category=paste0("Round: ", re$label),
        title=m$title %||% "", team=sim_team, sim_prob=sim_prob,
        yes_ask=ya, yes_bid=yb, yes_fee=ev$yes_fee, no_fee=ev$no_fee,
        yes_ev=ev$yes_ev, no_ev=ev$no_ev,
        best_side=ev$best_side, best_ev=ev$best_ev))
    }
  }

  # --- Seed-level props helper ---
  add_seed_props <- function(series, category, compute_fn) {
    markets <- fetch_kalshi_markets(series_ticker = series)
    cat(sprintf("  %s: %d markets\n", series, length(markets)))
    for (m in markets) {
      fs <- as.numeric(m$floor_strike %||% NA)
      cs <- as.numeric(m$cap_strike %||% NA)
      st <- m$strike_type %||% "greater_or_equal"
      et <- m$event_ticker %||% ""
      ttl <- m$title %||% ""
      sp <- compute_fn(raw, fs, cs, st, et, ttl, n_sims_local)
      if (is.na(sp)) next
      ya <- round(as.numeric(m$yes_ask_dollars %||% 0) * 100)
      yb <- round(as.numeric(m$yes_bid_dollars %||% 0) * 100)
      ev <- kalshi_compute_ev(sp, ya, yb)
      all_edges <- bind_rows(all_edges, tibble(
        ticker=m$ticker %||% "", category=category, title=ttl, team=NA_character_,
        sim_prob=sp, yes_ask=ya, yes_bid=yb, yes_fee=ev$yes_fee, no_fee=ev$no_fee,
        yes_ev=ev$yes_ev, no_ev=ev$no_ev,
        best_side=ev$best_side, best_ev=ev$best_ev))
    }
    all_edges
  }

  # Seed Sum
  all_edges <- add_seed_props("KXMARMADSEEDSUM", "Seed Sum", function(raw, fs, cs, st, et, ttl, ns) {
    rc <- if (grepl("F4", et)) "Final_4" else if (grepl("T2", et)) "Title_Game" else return(NA_real_)
    adv <- raw %>% filter(.data[[rc]]==1) %>% group_by(sim_id) %>% summarise(ss=sum(seed), .groups="drop")
    adv <- tibble(sim_id=1:ns) %>% left_join(adv, by="sim_id") %>% mutate(ss=coalesce(ss, 0L))
    if (st=="between" && !is.na(cs)) mean(adv$ss>=fs & adv$ss<=cs) else mean(adv$ss>=fs)
  })

  # Seed Count Per Round
  all_edges <- add_seed_props("KXMARMADSEEDROUND", "Seed Count", function(raw, fs, cs, st, et, ttl, ns) {
    sm <- str_match(et, "S(\\d+)"); if (is.na(sm[1,2])) return(NA_real_)
    ts <- as.integer(sm[1,2])
    rc <- case_when(grepl("R8",et)~"Elite_8", grepl("R16",et)~"Sweet_16",
                    grepl("R32",et)~"Round_32", grepl("F4",et)~"Final_4", TRUE~NA_character_)
    if (is.na(rc)) return(NA_real_)
    sc <- raw %>% filter(seed==ts, .data[[rc]]==1) %>% group_by(sim_id) %>% summarise(n=n(), .groups="drop")
    sc <- tibble(sim_id=1:ns) %>% left_join(sc, by="sim_id") %>% mutate(n=coalesce(n, 0L))
    if (grepl("exactly",ttl,ignore.case=TRUE) || st=="between") mean(sc$n==fs) else mean(sc$n>=fs)
  })

  # Highest Seed
  all_edges <- add_seed_props("KXMARMADSEED", "Highest Seed", function(raw, fs, cs, st, et, ttl, ns) {
    rc <- case_when(grepl("-26R32$",et)~"Round_32", grepl("-26R16$",et)~"Sweet_16",
                    grepl("-26R8$",et)~"Elite_8", grepl("-26F4$",et)~"Final_4",
                    grepl("-26$",et)~"Champion", TRUE~NA_character_)
    if (is.na(rc)) return(NA_real_)
    ms <- raw %>% filter(.data[[rc]]==1) %>% group_by(sim_id) %>% summarise(mx=max(seed), .groups="drop")
    ms <- tibble(sim_id=1:ns) %>% left_join(ms, by="sim_id") %>% mutate(mx=coalesce(mx, 0L))
    if (st=="between" && !is.na(cs)) mean(ms$mx>=fs & ms$mx<=cs)
    else if (grepl("exactly",ttl,ignore.case=TRUE)) mean(ms$mx==fs) else mean(ms$mx>=fs)
  })

  # Upsets
  all_edges <- add_seed_props("KXMARMADUPSET", "Upsets", function(raw, fs, cs, st, et, ttl, ns) {
    rc <- case_when(grepl("R64",et)~"Round_32", grepl("R32",et)~"Sweet_16",
                    grepl("R16",et)~"Elite_8", grepl("R8",et)~"Final_4", TRUE~NA_character_)
    if (is.na(rc)) return(NA_real_)
    us <- if (rc=="Round_32") 9:16 else if (rc=="Sweet_16") 5:16 else if (rc=="Elite_8") 3:16 else 2:16
    uc <- raw %>% filter(.data[[rc]]==1, seed %in% us) %>% group_by(sim_id) %>% summarise(n=n(), .groups="drop")
    uc <- tibble(sim_id=1:ns) %>% left_join(uc, by="sim_id") %>% mutate(n=coalesce(n, 0L))
    if (grepl("exactly",ttl,ignore.case=TRUE)) mean(uc$n==fs) else mean(uc$n>=fs)
  })

  # Seed Wins
  all_edges <- add_seed_props("KXMARMADSEEDWIN", "Seed Wins", function(raw, fs, cs, st, et, ttl, ns) {
    sm <- str_match(et, "S(\\d+)$")
    if (!is.na(sm[1,2])) { ts <- as.integer(sm[1,2]) }
    else { nums <- as.integer(str_extract_all(et, "\\d+")[[1]]); nums <- nums[nums>=10 & nums<=16]
           if (length(nums)==0) return(NA_real_); ts <- nums }
    w <- raw %>% filter(seed %in% ts, Round_32==1) %>% group_by(sim_id) %>% summarise(n=n(), .groups="drop")
    w <- tibble(sim_id=1:ns) %>% left_join(w, by="sim_id") %>% mutate(n=coalesce(n, 0L))
    if (grepl("exactly",ttl,ignore.case=TRUE)) mean(w$n==fs) else mean(w$n>=fs)
  })

  all_edges %>% arrange(desc(best_ev))
}

# Pre-compute team probs for Kalshi matching
team_probs <- raw %>%
  group_by(team, seed) %>%
  summarise(Round_32=mean(Round_32), Sweet_16=mean(Sweet_16), Elite_8=mean(Elite_8),
            Final_4=mean(Final_4), Title_Game=mean(Title_Game), Champion=mean(Champion),
            .groups="drop")

# Fetch Kalshi edges at startup
kalshi_edges <- tryCatch(fetch_all_kalshi_edges(raw, team_probs), error = function(e) {
  cat(sprintf("Kalshi edge error: %s\n", e$message)); tibble()
})
cat(sprintf("Kalshi: %d markets, %d +EV\n", nrow(kalshi_edges), sum(kalshi_edges$best_ev > 0, na.rm=TRUE)))

cat("Starting server...\n")

# =============================================================================
# 3. Shiny App
# =============================================================================

prob_to_american <- function(prob) {
  ifelse(prob >= 0.5, round(-100 * prob / (1 - prob)), round(100 / prob - 100))
}

ui <- fluidPage(
  titlePanel(sprintf("March Madness Simulator (%dk sims | %d teams remaining)",
                     n_sims/1000, n_remaining)),

  tabsetPanel(
    # --- Tab 1: Team Advancement ---
    tabPanel("Team Advancement",
      br(),
      DT::dataTableOutput("team_table")
    ),

    # --- Tab 2: Seed Props ---
    tabPanel("Seed Props",
      br(),
      fluidRow(
        column(4, selectInput("seed_num", "Seed:", choices = 1:16, selected = 1)),
        column(4, selectInput("seed_round", "Round:",
          choices = c("Round_32", "Sweet_16", "Elite_8", "Final_4", "Title_Game", "Champion"),
          selected = "Final_4"))
      ),
      h4(textOutput("seed_title")),
      tableOutput("seed_counts"),
      h4("Individual Teams at This Seed:"),
      tableOutput("seed_teams")
    ),

    # --- Tab 3: Conference Props ---
    tabPanel("Conference Props",
      br(),
      h4("Conference Win Totals & Championship Probability"),
      DT::dataTableOutput("conf_table"),
      br(),
      fluidRow(
        column(4, selectInput("conf_sel", "Conference:",
          choices = if(!is.null(conf_map)) sort(unique(conf_map$conference)) else "N/A")),
        column(4, numericInput("conf_wins_threshold", "Wins >=", value = 10, min = 0, max = 30))
      ),
      h4(textOutput("conf_over_title")),
      textOutput("conf_over_result")
    ),

    # --- Tab 4: Custom Query ---
    tabPanel("Custom Query",
      br(),
      h4("Build a custom probability query"),
      fluidRow(
        column(3, selectInput("q_filter_type", "Filter by:",
          choices = c("Seed", "Team", "Conference", "Region"))),
        column(3, uiOutput("q_filter_value_ui")),
        column(3, selectInput("q_round", "Round:",
          choices = c("Round_32", "Sweet_16", "Elite_8", "Final_4", "Title_Game", "Champion", "total_wins"),
          selected = "Final_4")),
        column(3, selectInput("q_condition", "Condition:",
          choices = c("at least", "exactly", "at most")))
      ),
      fluidRow(
        column(3, numericInput("q_value", "Value:", value = 1, min = 0, max = 10)),
        column(3, actionButton("q_run", "Calculate", class = "btn-primary"))
      ),
      br(),
      h3(textOutput("q_result")),
      h4(textOutput("q_odds")),
      br(),
      h4("Distribution:"),
      tableOutput("q_distribution")
    ),

    # --- Tab 5: Final Four Seed Sum ---
    tabPanel("F4 Seed Sum",
      br(),
      h4("Final Four Seed Sum Distribution"),
      tableOutput("f4_seed_sum"),
      br(),
      fluidRow(
        column(4, numericInput("f4_line", "Over/Under Line:", value = 10, min = 4, max = 40))
      ),
      h4(textOutput("f4_over_result"))
    ),

    # --- Tab 6: Upset Props ---
    tabPanel("Upset Props",
      br(),
      h4("Round of 64 Upset Probabilities"),
      tableOutput("upset_table"),
      br(),
      h4("Upset Count Distribution (total upsets in Round of 64)"),
      tableOutput("upset_count_table")
    ),

    # --- Tab 7: Highest Seed by Round ---
    tabPanel("Highest Seed",
      br(),
      fluidRow(
        column(4, selectInput("hs_round", "Round:",
          choices = c("Round_32", "Sweet_16", "Elite_8", "Final_4", "Title_Game", "Champion"),
          selected = "Round_32"))
      ),
      h4(textOutput("hs_title")),
      tableOutput("hs_table")
    ),

    # --- Tab 8: Kalshi Edges ---
    tabPanel("Kalshi Edges",
      br(),
      fluidRow(
        column(4, selectInput("ke_category", "Category:",
          choices = c("All", sort(unique(kalshi_edges$category))), selected = "All")),
        column(4, selectInput("ke_min_ev", "Min EV%:",
          choices = c("All", "1%", "3%", "5%", "10%", "15%"), selected = "All")),
        column(4, actionButton("ke_refresh", "Refresh Kalshi Prices", class = "btn-primary"))
      ),
      h4(textOutput("ke_summary")),
      DT::dataTableOutput("ke_table")
    ),

    # --- Tab 9: Seed Sum by Round ---
    tabPanel("Seed Sum",
      br(),
      fluidRow(
        column(4, selectInput("ss_round", "Round:",
          choices = c("Final_4", "Elite_8", "Sweet_16", "Title_Game"),
          selected = "Final_4"))
      ),
      h4(textOutput("ss_title")),
      tableOutput("ss_buckets"),
      br(),
      h4("Exact Distribution:"),
      tableOutput("ss_exact"),
      br(),
      fluidRow(
        column(4, numericInput("ss_line", "Custom Over/Under:", value = 10, min = 2, max = 60))
      ),
      h4(textOutput("ss_custom_result"))
    )
  )
)

server <- function(input, output, session) {

  # --- Team Advancement ---
  output$team_table <- DT::renderDataTable({
    raw %>%
      group_by(team, seed, region) %>%
      summarise(
        R32 = sprintf("%.1f%%", mean(Round_32)*100),
        S16 = sprintf("%.1f%%", mean(Sweet_16)*100),
        E8 = sprintf("%.1f%%", mean(Elite_8)*100),
        F4 = sprintf("%.1f%%", mean(Final_4)*100),
        Title = sprintf("%.1f%%", mean(Title_Game)*100),
        Champ = sprintf("%.1f%%", mean(Champion)*100),
        `Champ Odds` = prob_to_american(mean(Champion)),
        .groups = "drop"
      ) %>%
      arrange(desc(as.numeric(gsub("%", "", Champ))))
  }, options = list(pageLength = 25))

  # --- Seed Props ---
  seed_data <- reactive({
    s <- as.integer(input$seed_num)
    r <- input$seed_round
    raw %>% filter(seed == s) %>%
      group_by(sim_id) %>%
      summarise(count = sum(.data[[r]]), .groups = "drop")
  })

  output$seed_title <- renderText({
    sprintf("#%s seeds reaching %s:", input$seed_num, gsub("_", " ", input$seed_round))
  })

  output$seed_counts <- renderTable({
    d <- seed_data()
    max_count <- max(d$count)
    data.frame(
      Count = 0:max_count,
      Probability = sapply(0:max_count, function(k) sprintf("%.1f%%", mean(d$count == k)*100)),
      `American Odds` = sapply(0:max_count, function(k) {
        p <- mean(d$count == k)
        if (p > 0 && p < 1) prob_to_american(p) else "N/A"
      }),
      check.names = FALSE
    )
  })

  output$seed_teams <- renderTable({
    s <- as.integer(input$seed_num)
    r <- input$seed_round
    raw %>% filter(seed == s) %>%
      group_by(team) %>%
      summarise(
        Probability = sprintf("%.1f%%", mean(.data[[r]])*100),
        .groups = "drop"
      ) %>%
      arrange(desc(Probability))
  })

  # --- Conference Props ---
  output$conf_table <- DT::renderDataTable({
    if (!"conference" %in% names(raw)) return(data.frame(Message = "No conference data"))
    raw %>%
      filter(!is.na(conference)) %>%
      group_by(sim_id, conference) %>%
      summarise(wins = sum(total_wins), champ = max(Champion), .groups = "drop") %>%
      group_by(conference) %>%
      summarise(
        `Avg Wins` = round(mean(wins), 1),
        `P(Champ)` = sprintf("%.1f%%", mean(champ >= 1)*100),
        .groups = "drop"
      ) %>%
      arrange(desc(`Avg Wins`))
  }, options = list(pageLength = 20))

  output$conf_over_title <- renderText({
    sprintf("%s wins >= %d:", input$conf_sel, input$conf_wins_threshold)
  })

  output$conf_over_result <- renderText({
    if (!"conference" %in% names(raw)) return("No conference data")
    conf_wins <- raw %>%
      filter(conference == input$conf_sel) %>%
      group_by(sim_id) %>%
      summarise(wins = sum(total_wins), .groups = "drop")
    p <- mean(conf_wins$wins >= input$conf_wins_threshold)
    sprintf("%.1f%% (%s)", p*100, prob_to_american(p))
  })

  # --- Custom Query ---
  output$q_filter_value_ui <- renderUI({
    choices <- switch(input$q_filter_type,
      "Seed" = 1:16,
      "Team" = sort(unique(raw$team)),
      "Conference" = if("conference" %in% names(raw)) sort(unique(na.omit(raw$conference))) else "N/A",
      "Region" = sort(unique(raw$region))
    )
    selectInput("q_filter_value", input$q_filter_type, choices = choices)
  })

  q_result <- eventReactive(input$q_run, {
    req(input$q_filter_value)
    filtered <- switch(input$q_filter_type,
      "Seed" = raw %>% filter(seed == as.integer(input$q_filter_value)),
      "Team" = raw %>% filter(team == input$q_filter_value),
      "Conference" = raw %>% filter(conference == input$q_filter_value),
      "Region" = raw %>% filter(region == input$q_filter_value)
    )

    by_sim <- filtered %>%
      group_by(sim_id) %>%
      summarise(val = sum(.data[[input$q_round]]), .groups = "drop")

    v <- input$q_value
    p <- switch(input$q_condition,
      "at least" = mean(by_sim$val >= v),
      "exactly" = mean(by_sim$val == v),
      "at most" = mean(by_sim$val <= v)
    )

    # Distribution
    max_val <- max(by_sim$val)
    dist <- data.frame(
      Count = 0:max_val,
      Probability = sapply(0:max_val, function(k) sprintf("%.2f%%", mean(by_sim$val == k)*100)),
      Odds = sapply(0:max_val, function(k) {
        pk <- mean(by_sim$val == k)
        if (pk > 0 & pk < 1) prob_to_american(pk) else "N/A"
      })
    )

    list(p = p, dist = dist)
  })

  output$q_result <- renderText({
    r <- q_result()
    sprintf("Probability: %.2f%%", r$p * 100)
  })

  output$q_odds <- renderText({
    r <- q_result()
    if (r$p > 0 & r$p < 1) sprintf("American Odds: %s", prob_to_american(r$p)) else ""
  })

  output$q_distribution <- renderTable({ q_result()$dist })

  # --- F4 Seed Sum ---
  f4_sums <- reactive({
    raw %>% filter(Final_4 == 1) %>%
      group_by(sim_id) %>%
      summarise(seed_sum = sum(seed), .groups = "drop")
  })

  output$f4_seed_sum <- renderTable({
    sums <- f4_sums()
    breaks <- c(4, 6, 8, 10, 12, 14, 16, 18, 20, 25, 30, 40)
    data.frame(
      `Seed Sum` = paste0("Over ", breaks),
      Probability = sapply(breaks, function(b) sprintf("%.1f%%", mean(sums$seed_sum > b)*100)),
      Odds = sapply(breaks, function(b) {
        p <- mean(sums$seed_sum > b)
        if (p > 0 & p < 1) prob_to_american(p) else "N/A"
      }),
      check.names = FALSE
    )
  })

  output$f4_over_result <- renderText({
    sums <- f4_sums()
    p_over <- mean(sums$seed_sum > input$f4_line)
    p_under <- mean(sums$seed_sum < input$f4_line)
    sprintf("Over %.1f: %.1f%% (%s) | Under %.1f: %.1f%% (%s)",
      input$f4_line, p_over*100, prob_to_american(p_over),
      input$f4_line, p_under*100, prob_to_american(p_under))
  })

  # --- Upset Props ---
  output$upset_table <- renderTable({
    matchups <- data.frame(fav=c(1,2,3,4,5,6,7,8), dog=c(16,15,14,13,12,11,10,9))
    result <- lapply(1:nrow(matchups), function(i) {
      fav_seed <- matchups$fav[i]; dog_seed <- matchups$dog[i]
      by_sim <- raw %>% filter(seed == dog_seed) %>%
        group_by(sim_id) %>% summarise(wins = sum(Round_32), .groups = "drop")
      p_indiv <- mean(raw$Round_32[raw$seed == dog_seed])
      p_at_least_1 <- mean(by_sim$wins >= 1)
      data.frame(
        Matchup = sprintf("#%d vs #%d", fav_seed, dog_seed),
        `Indiv Upset` = sprintf("%.1f%%", p_indiv*100),
        `P(1+ Upset)` = sprintf("%.1f%%", p_at_least_1*100),
        `1+ Odds` = prob_to_american(p_at_least_1),
        check.names = FALSE
      )
    })
    bind_rows(result)
  })

  output$upset_count_table <- renderTable({
    # Total upsets in R64: count how many underdogs (seed > opponent seed) win
    # Underdogs are seeds 9-16
    upsets_per_sim <- raw %>% filter(seed >= 9) %>%
      group_by(sim_id) %>% summarise(n_upsets = sum(Round_32), .groups = "drop")
    max_u <- min(max(upsets_per_sim$n_upsets), 20)
    data.frame(
      Upsets = 0:max_u,
      Probability = sapply(0:max_u, function(k) sprintf("%.2f%%", mean(upsets_per_sim$n_upsets == k)*100)),
      `Cumulative (>=)` = sapply(0:max_u, function(k) sprintf("%.1f%%", mean(upsets_per_sim$n_upsets >= k)*100)),
      check.names = FALSE
    )
  })

  # --- Highest Seed by Round ---
  output$hs_title <- renderText({
    sprintf("Highest numerical seed to qualify for %s:", gsub("_", " ", input$hs_round))
  })

  output$hs_table <- renderTable({
    r <- input$hs_round
    highest <- raw %>% filter(.data[[r]] == 1) %>%
      group_by(sim_id) %>% summarise(highest = max(seed), .groups = "drop")

    results <- data.frame(
      `Highest Seed` = paste0("#", 1:16),
      `P(Exactly)` = sapply(1:16, function(s) sprintf("%.1f%%", mean(highest$highest == s)*100)),
      `P(>= this seed)` = sapply(1:16, function(s) sprintf("%.1f%%", mean(highest$highest >= s)*100)),
      Odds = sapply(1:16, function(s) {
        p <- mean(highest$highest == s)
        if (p > 0.001 & p < 0.999) prob_to_american(p) else "N/A"
      }),
      check.names = FALSE
    )
    results %>% filter(as.numeric(gsub("[^0-9]", "", `P(Exactly)`)) > 0)
  })

  # --- Seed Sum by Round ---
  output$ss_title <- renderText({
    sprintf("Sum of seeds in %s:", gsub("_", " ", input$ss_round))
  })

  output$ss_buckets <- renderTable({
    r <- input$ss_round
    sums <- raw %>% filter(.data[[r]] == 1) %>%
      group_by(sim_id) %>% summarise(seed_sum = sum(seed), .groups = "drop")

    # Kalshi-style buckets (adjust based on round)
    n_teams <- switch(r, "Final_4" = 4, "Elite_8" = 8, "Sweet_16" = 16, "Title_Game" = 2, 4)
    min_sum <- n_teams  # all 1-seeds

    if (n_teams == 4) {
      buckets <- list("4"=c(4,4), "5-6"=c(5,6), "7-8"=c(7,8), "9-10"=c(9,10),
                       "11-14"=c(11,14), "15-20"=c(15,20), "21-30"=c(21,30), "31+"=c(31,100))
    } else if (n_teams == 2) {
      buckets <- list("2"=c(2,2), "3-4"=c(3,4), "5-6"=c(5,6), "7-10"=c(7,10), "11+"=c(11,40))
    } else {
      # Generic
      max_s <- max(sums$seed_sum)
      step <- max(1, round((max_s - min_sum) / 8))
      buckets <- list()
      lo <- min_sum
      while (lo <= max_s) {
        hi <- min(lo + step - 1, max_s + 10)
        buckets[[sprintf("%d-%d", lo, hi)]] <- c(lo, hi)
        lo <- hi + 1
      }
    }

    data.frame(
      Bucket = names(buckets),
      Probability = sapply(buckets, function(b) sprintf("%.1f%%", mean(sums$seed_sum >= b[1] & sums$seed_sum <= b[2])*100)),
      Odds = sapply(buckets, function(b) {
        p <- mean(sums$seed_sum >= b[1] & sums$seed_sum <= b[2])
        if (p > 0.001 & p < 0.999) prob_to_american(p) else "N/A"
      }),
      check.names = FALSE
    )
  })

  output$ss_exact <- renderTable({
    r <- input$ss_round
    sums <- raw %>% filter(.data[[r]] == 1) %>%
      group_by(sim_id) %>% summarise(seed_sum = sum(seed), .groups = "drop")

    vals <- sort(unique(sums$seed_sum))
    data.frame(
      `Seed Sum` = vals,
      Probability = sapply(vals, function(v) sprintf("%.1f%%", mean(sums$seed_sum == v)*100)),
      `Cumulative (<=)` = sapply(vals, function(v) sprintf("%.1f%%", mean(sums$seed_sum <= v)*100)),
      check.names = FALSE
    ) %>% filter(as.numeric(gsub("[^0-9.]", "", Probability)) >= 0.1)
  })

  output$ss_custom_result <- renderText({
    r <- input$ss_round
    sums <- raw %>% filter(.data[[r]] == 1) %>%
      group_by(sim_id) %>% summarise(seed_sum = sum(seed), .groups = "drop")
    p_over <- mean(sums$seed_sum > input$ss_line)
    p_under <- mean(sums$seed_sum < input$ss_line)
    p_exact <- mean(sums$seed_sum == input$ss_line)
    sprintf("Over %d: %.1f%% (%s) | Under %d: %.1f%% (%s) | Exactly %d: %.1f%%",
      input$ss_line, p_over*100, prob_to_american(p_over),
      input$ss_line, p_under*100, prob_to_american(p_under),
      input$ss_line, p_exact*100)
  })

  # --- Kalshi Edges ---
  ke_data <- reactiveVal(kalshi_edges)

  observeEvent(input$ke_refresh, {
    showNotification("Refreshing Kalshi prices...", type = "message")
    new_edges <- tryCatch(fetch_all_kalshi_edges(raw, team_probs), error = function(e) {
      showNotification(sprintf("Error: %s", e$message), type = "error"); ke_data()
    })
    ke_data(new_edges)
    showNotification(sprintf("Refreshed: %d markets", nrow(new_edges)), type = "message")
  })

  output$ke_summary <- renderText({
    d <- ke_data()
    if (nrow(d) == 0) return("No Kalshi markets found.")
    pos <- sum(d$best_ev > 0, na.rm = TRUE)
    best <- if (pos > 0) max(d$best_ev[d$best_ev > 0], na.rm = TRUE) * 100 else 0
    sprintf("%d markets | %d +EV edges | Best: +%.1f%%", nrow(d), pos, best)
  })

  output$ke_table <- DT::renderDataTable({
    d <- ke_data()
    if (nrow(d) == 0) return(data.frame(Message = "No markets"))

    # Filter by category
    if (input$ke_category != "All") d <- d %>% filter(category == input$ke_category)

    # Filter by min EV
    min_ev <- switch(input$ke_min_ev,
      "1%" = 0.01, "3%" = 0.03, "5%" = 0.05, "10%" = 0.10, "15%" = 0.15, -Inf)
    d <- d %>% filter(best_ev >= min_ev)

    d %>%
      transmute(
        Category = category,
        Market = title,
        `Fair` = sprintf("%.1f%%", sim_prob * 100),
        `YES Ask` = paste0(yes_ask, "\u00A2"),
        `NO Ask` = paste0(100 - yes_bid, "\u00A2"),
        `YES Fee` = sprintf("%.1f\u00A2", yes_fee),
        `NO Fee` = sprintf("%.1f\u00A2", no_fee),
        `YES EV` = sprintf("%+.1f%%", yes_ev * 100),
        `NO EV` = sprintf("%+.1f%%", no_ev * 100),
        Side = best_side,
        `Best EV` = sprintf("%+.1f%%", best_ev * 100),
        .sort_ev = best_ev
      ) %>%
      arrange(desc(.sort_ev)) %>%
      select(-.sort_ev)
  }, options = list(pageLength = 50))
}

shinyApp(ui, server, options = list(port = 8085, host = "0.0.0.0", launch.browser = TRUE))
