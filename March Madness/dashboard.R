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

# Add conference
conf_map <- tryCatch({
  cbd_bpi_ratings() %>%
    transmute(team_short = team, conference = conf) %>%
    mutate(team = map_chr(team_short, ~ get_standard_team(.x, teams_std = ts))) %>%
    select(team, conference) %>% distinct(team, .keep_all = TRUE)
}, error = function(e) NULL)

# Sim with fitted params
sim_one <- function() {
  df <- bracket_64
  all_teams <- df$team; all_seeds <- df$seed; n_all <- length(all_teams)
  progress <- matrix(0L, nrow = n_all, ncol = 6)
  seed_order <- c(1,16,8,9,5,12,4,13,6,11,3,14,7,10,2,15)
  teams <- df$team; seeds <- df$seed; regions <- df$region
  ratings <- ifelse(is.na(df$composite_rating), 0, df$composite_rating)
  n <- length(teams); round_num <- 1
  while (n > 1) {
    if (n == 4) { reg_idx <- match(regions, region_order); ord <- order(reg_idx); ord <- ord[c(1,3,2,4)]
    } else if (n == 2) { ord <- 1:2
    } else { seed_pos <- match(seeds, seed_order); ord <- order(match(regions, region_order), seed_pos) }
    teams <- teams[ord]; seeds <- seeds[ord]; regions <- regions[ord]; ratings <- ratings[ord]
    n_games <- n %/% 2; i1 <- seq(1,n,by=2); i2 <- seq(2,n,by=2)
    expected_diff <- ratings[i1] - ratings[i2]
    actual_margin <- rnorm(n_games, mean = expected_diff, sd = 11.2)
    w <- ifelse(actual_margin > 0, i1, i2)
    teams <- teams[w]; seeds <- seeds[w]; regions <- regions[w]
    ratings <- ratings[w] + rnorm(n_games, 0, 0.17)
    n <- length(teams); col <- min(round_num, 6)
    progress[match(teams, all_teams), col] <- 1L; round_num <- round_num + 1
  }
  data.frame(team=all_teams, seed=all_seeds, region=df$region,
    Round_32=progress[,1], Sweet_16=progress[,2], Elite_8=progress[,3],
    Final_4=progress[,4], Title_Game=progress[,5], Champion=progress[,6],
    stringsAsFactors=FALSE)
}

n_sims <- 50000
cat(sprintf("Running %dk simulations...\n", n_sims/1000))
t0 <- Sys.time()
raw <- map_dfr(1:n_sims, function(i) {
  r <- sim_one()
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

cat("Dashboard ready. Starting server...\n")

# =============================================================================
# 2. Shiny App
# =============================================================================

prob_to_american <- function(prob) {
  ifelse(prob >= 0.5, round(-100 * prob / (1 - prob)), round(100 / prob - 100))
}

ui <- fluidPage(
  titlePanel(sprintf("March Madness Simulator (%dk sims)", n_sims/1000)),

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

    # --- Tab 8: Seed Sum by Round ---
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
}

shinyApp(ui, server, options = list(port = 8085, host = "0.0.0.0", launch.browser = TRUE))
