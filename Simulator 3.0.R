# Load required libraries
library(rvest)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(hoopR)
library(googlesheets4)
library(jsonlite)
library(readr)
library(gt)
library(httr)

# Uncomment if needed:
# gs4_auth()

######################################
### 1. Acquire Bracket Data (Scraping) ###
######################################
bracket_url <- "https://www.ncaa.com/march-madness-live/bracket"  # update if needed
bracket_page <- read_html(bracket_url)
region_nodes <- bracket_page %>% html_nodes("div.bracket-container div.region")

scrape_region <- function(region_node) {
  region_name <- region_node %>% 
    html_node("span.subtitle") %>% 
    html_text(trim = TRUE)
  
  round_nodes <- region_node %>% html_nodes("div.region-round")
  
  if (length(round_nodes) > 0) {
    map_df(round_nodes, function(round_node) {
      round_class <- round_node %>% html_attr("class")
      round_label <- str_extract(round_class, "round-\\d+")
      game_nodes <- round_node %>% html_nodes("div.game-pod")
      
      if (length(game_nodes) > 0) {
        map_df(game_nodes, function(game_node) {
          top_seed <- game_node %>% 
            html_node("div.team.team-top span.overline") %>% 
            html_text(trim = TRUE)
          top_team <- game_node %>% 
            html_node("div.team.team-top p.body_2") %>% 
            html_text(trim = TRUE)
          bottom_seed <- game_node %>% 
            html_node("div.team.team-bottom span.overline") %>% 
            html_text(trim = TRUE)
          bottom_team <- game_node %>% 
            html_node("div.team.team-bottom p.body_2") %>% 
            html_text(trim = TRUE)
          
          tibble(
            region   = region_name,
            round    = round_label,  # "round-1" means Round of 64, "round-2" = Round of 32, etc.
            team     = c(top_team, bottom_team),
            seed     = c(top_seed, bottom_seed),
            play_in  = 0
          )
        })
      } else {
        team_nodes <- round_node %>% html_nodes("div.team")
        teams <- team_nodes %>% map_chr(~ .x %>% html_node("p.body_2") %>% html_text(trim = TRUE))
        seeds <- team_nodes %>% map_chr(~ .x %>% html_node("span.overline") %>% html_text(trim = TRUE))
        tibble(
          region   = region_name,
          round    = round_label,
          team     = teams,
          seed     = seeds,
          play_in  = 0
        )
      }
    })
  } else {
    team_nodes <- region_node %>% html_nodes("div.team")
    teams <- team_nodes %>% map_chr(~ .x %>% html_node("p.body_2") %>% html_text(trim = TRUE))
    seeds <- team_nodes %>% map_chr(~ .x %>% html_node("span.overline") %>% html_text(trim = TRUE))
    tibble(
      region   = region_name,
      round    = "Round of 64",
      team     = teams,
      seed     = seeds,
      play_in  = 0
    )
  }
}

bracket_data <- map_df(region_nodes, scrape_region)

# --- Functions to scrape Final Four and Final ---
scrape_final_four <- function(page) {
  node <- page %>% html_node(xpath = "//*[@id='602']")
  if(is.null(node)) return(tibble())
  top_seed <- node %>% html_node("div.team.team-top span.overline") %>% html_text(trim = TRUE)
  top_team <- node %>% html_node("div.team.team-top p.body_2") %>% html_text(trim = TRUE)
  bottom_seed <- node %>% html_node("div.team.team-bottom span.overline") %>% html_text(trim = TRUE)
  bottom_team <- node %>% html_node("div.team.team-bottom p.body_2") %>% html_text(trim = TRUE)
  tibble(
    region = "Final Four",
    round  = "Final Four",
    team   = c(top_team, bottom_team),
    seed   = c(top_seed, bottom_seed),
    play_in = 0
  )
}

scrape_final <- function(page) {
  node <- page %>% html_node(xpath = "//*[@id='701']")
  if(is.null(node)) return(tibble())
  top_seed <- node %>% html_node("div.team.team-top span.overline") %>% html_text(trim = TRUE)
  top_team <- node %>% html_node("div.team.team-top p.body_2") %>% html_text(trim = TRUE)
  bottom_seed <- node %>% html_node("div.team.team-bottom span.overline") %>% html_text(trim = TRUE)
  bottom_team <- node %>% html_node("div.team.team-bottom p.body_2") %>% html_text(trim = TRUE)
  tibble(
    region = "Final",
    round  = "Final",
    team   = c(top_team, bottom_team),
    seed   = c(top_seed, bottom_seed),
    play_in = 0
  )
}

final_bracket <- bind_rows(
  bracket_data,
  scrape_final_four(bracket_page),
  scrape_final(bracket_page)
)

######################################
### 2. Acquire Power Ratings         ###
######################################
teams_std <- espn_mbb_teams(2025) %>% 
  mutate(team = ifelse(team == "McNeese", "McNeese State", team))

clean_text <- function(x) {
  sub("St\\.$", "State", x, ignore.case = TRUE)
}

get_standard_team <- function(team, teams_std) {
  if(is.na(team)) return(NA_character_)
  team_clean <- clean_text(team)
  for(i in 1:nrow(teams_std)) {
    variants <- teams_std[i, c("abbreviation", "display_name", "short_name", "mascot", "nickname", "team")]
    for(variant in variants) {
      if(!is.na(variant) && team_clean == clean_text(variant)) {
        return(teams_std$display_name[i])
      }
    }
  }
  team
}

# ESPN BPI Ratings
base_url <- "https://site.web.api.espn.com/apis/fitt/v3/sports/basketball/mens-college-basketball/powerindex?region=us&lang=en&groups=50&limit=50&page="
all_teams <- list()
for (page in 1:8) {
  url <- paste0(base_url, page)
  response <- GET(url)
  if(status_code(response) != 200) {
    cat("Error fetching page", page, "\n")
    next
  }
  content_text <- content(response, "text", encoding = "UTF-8")
  data <- fromJSON(content_text)
  team_names <- data$teams$team$displayName
  bpi_values <- map_dbl(data$teams$categories, function(cat) {
    if(is.data.frame(cat) && "name" %in% names(cat)) {
      bpi_row <- cat %>% filter(name == "bpi")
      if(nrow(bpi_row) > 0 && !is.null(bpi_row$values[[1]][1])) {
        return(bpi_row$values[[1]][1])
      }
    }
    NA_real_
  })
  team_data <- tibble(
    team = team_names,
    bpi = bpi_values
  )
  all_teams[[page]] <- team_data
}
final_bpi_data <- bind_rows(all_teams)
clean_bpi_data <- final_bpi_data %>% 
  mutate(standard_team = map_chr(team, ~ get_standard_team(.x, teams_std = teams_std)))

sheet_url <- "https://docs.google.com/spreadsheets/d/10o9dwZeyREliM8iOIevpHIhNeCHYeVGE-y50N3KHujQ/edit?gid=0#gid=0"
kenpom_data <- read_sheet(sheet_url) %>% 
  mutate(Team = str_replace(Team, "\\s\\d+$", "")) %>% 
  as.data.frame() %>% 
  mutate(Team = ifelse(Team == "Connecticut", "UConn", 
                       ifelse(Team == "Mississippi", "Ole Miss",
                              ifelse(Team == "Nebraska Omaha", "Omaha", Team))))
clean_kenpom_data <- kenpom_data %>%
  mutate(standard_team = map_chr(Team, ~ get_standard_team(.x, teams_std = teams_std))) %>% 
  select(standard_team, NetRtg) %>% 
  mutate(NetRtg = map_dbl(NetRtg, ~ ifelse(is.null(.x), NA_real_, as.numeric(.x))))

torvik_url <- "https://barttorvik.com/trank.php?year=2025&t=0&json=1"
response <- GET(torvik_url, add_headers("User-Agent" = "Mozilla/5.0"))
if(status_code(response) == 200) {
  torvik_data <- as.data.frame(fromJSON(content(response, "text", encoding = "UTF-8")))
  colnames(torvik_data)[c(1,2,3,4)] <- c("Team", "OffEff", "DefEff", "TorvikPower")
  torvik_data <- torvik_data %>%
    mutate(
      OffEff = as.numeric(OffEff),
      DefEff = as.numeric(DefEff),
      TorvikPower = as.numeric(TorvikPower),
      TorvikMargin = ((OffEff - DefEff) / 100) * 70
    ) %>% 
    select(Team, TorvikMargin)
} else {
  stop("Failed to retrieve Torvik data. HTTP status:", status_code(response))
}
clean_torvik_data <- torvik_data %>% 
  mutate(Team = ifelse(Team == "Connecticut", "UConn", 
                       ifelse(Team == "Mississippi", "Ole Miss",
                              ifelse(Team == "Nebraska Omaha", "Omaha", Team)))) %>% 
  mutate(standard_team = map_chr(Team, ~ get_standard_team(.x, teams_std = teams_std))) %>% 
  select(standard_team, TorvikMargin)

evan_miya <- read_sheet("https://docs.google.com/spreadsheets/d/1zaDMGo6bRMis_VD7yH4uZY9aFXbp0QUPp3HRoE8x3EY/edit?usp=drive_web&pli=1&authuser=0")
group_size <- 21
evan_miya$group <- rep(1:(nrow(evan_miya)/group_size), each = group_size)[1:nrow(evan_miya)]
evan_miya <- evan_miya %>%
  group_by(group) %>%
  mutate(row_number = row_number()) %>%
  ungroup() %>% as.data.frame()
evan_miya_wide <- evan_miya %>%
  pivot_wider(names_from = row_number, values_from = 1) %>%  
  select(-group) %>% 
  mutate(`2` = `2` %>% str_replace_all("[^A-Za-z0-9 ]", "") %>% str_squish())
colnames(evan_miya_wide) <- c(
  "Relative Ranking", "Team", "O-Rate", "D-Rate", "Relative Rating",
  "Opponent Adjust", "Pace Adjust", "Off Rank", "Def Rank", "True Tempo",
  "Tempo Rank", "Injury Rank", "Home Rank", "Roster Rank",
  "Kill Shots Per Game", "Kill Shots Conceded Per Game", "Kill Shots Margin Per Game",
  "Total Kill Shots", "Total Kill Shots Conceded", "D1 Wins", "D1 Losses"
)
clean_evan_miya_wide <- evan_miya_wide %>% 
  mutate(Team = ifelse(Team == "Connecticut", "UConn", 
                       ifelse(Team == "Mississippi", "Ole Miss", Team))) %>% 
  mutate(Team = sub("AM", "A&M", Team)) %>% 
  mutate(standard_team = map_chr(Team, ~ get_standard_team(.x, teams_std = teams_std %>% 
                                                             mutate(nickname = sub("'", "", nickname),
                                                                    nickname = sub("St.", "Saint ", nickname)
                                                             )))) %>% 
  select(standard_team, `Relative Rating`) %>% 
  mutate(`Relative Rating` = map_dbl(`Relative Rating`, ~ ifelse(is.null(.x), NA_real_, as.numeric(.x))))

power_ratings <- clean_bpi_data %>%
  left_join(clean_kenpom_data %>% rename(KenPomRating = NetRtg), by = c("standard_team")) %>%
  left_join(clean_torvik_data, by = c("standard_team")) %>%
  left_join(clean_evan_miya_wide %>% rename(EvanMiyaRating = `Relative Rating`), by = c("standard_team")) %>%  
  mutate(
    EvanMiyaRating = as.numeric(unlist(EvanMiyaRating)),
    KenPomMargin = (as.numeric(KenPomRating) / 100) * 70,
    EvanMiyaMargin = (EvanMiyaRating / 100) * 70,
    BPIMargin = bpi
  ) %>%
  select(team, KenPomMargin, BPIMargin, EvanMiyaMargin, TorvikMargin)

######################################
### 3. Build Current Bracket         ###
######################################
# Merge bracket data with power ratings.
final_bracket <- final_bracket %>%
  mutate(standard_team = map_chr(team, ~ get_standard_team(.x, teams_std = teams_std))) %>%
  mutate(seed = as.numeric(seed))

bracket_with_ratings <- left_join(final_bracket, power_ratings, 
                                  by = c("standard_team" = "team")) %>% 
  rowwise() %>%
  mutate(composite_rating = mean(c_across(KenPomMargin:TorvikMargin), na.rm = TRUE)) %>%
  ungroup()

# --- Mark round levels ---
# "round-1" = Round of 64, "round-2" = Round of 32, "round-3" = Sweet 16, "round-4" = Elite 8,
# "Final Four" = 5, "Final" = 6.
bracket_with_ratings <- bracket_with_ratings %>%
  mutate(round_num = case_when(
    grepl("round-", round) ~ as.numeric(str_extract(round, "\\d+")),
    round == "Final Four" ~ 5,
    round == "Final" ~ 6,
    TRUE ~ 1
  ))

# For each team, keep only the deepest row.
largest_round_per_team <- bracket_with_ratings %>%
  group_by(team) %>%
  filter(round_num == max(round_num, na.rm = TRUE)) %>%
  ungroup()

# Build the current bracket from these rows.
current_bracket <- largest_round_per_team %>%
  filter(!is.na(team) & team != "") %>%
  arrange(seed)  # Order by seed globally (not by region)

# --- New Logic: Determine Matchup Outcomes Based on Stage ---
# In a mid-round scenario, teams in a pairing may be at different stages.
# For each matchup (based on seeding order), if one team’s round_num is higher (i.e. advanced),
# that team is automatically the winner, and its opponent is eliminated.
# Otherwise, if both teams are at the same stage, the game is simulated.

######################################
### 4. Simulation Functions          ###
######################################
# Helper: Order teams by seeding (ignoring region)
get_bracket_matchups <- function(teams) {
  # Define a matchup order based on standard seeding.
  match_order <- tibble(
    seed = c(1,16,8,9,5,12,4,13,6,11,3,14,7,10,2,15),
    matchup_order = 1:16
  )
  teams %>%
    left_join(match_order, by = "seed") %>%
    arrange(matchup_order) %>%
    select(-matchup_order)
}

simulate_game <- function(team1, team2, game_number = 1, beta1 = 0.1, beta2 = 0.05, sd_margin = 11.2) {
  # If both teams are at the same stage, simulate the game.
  if(team1$round_num == team2$round_num) {
    expected_diff <- team1$composite_rating - team2$composite_rating
    actual_margin <- rnorm(1, mean = expected_diff, sd = sd_margin)
    if(actual_margin > 0) return(team1) else return(team2)
  } else {
    # Otherwise, the team with the higher round_num (i.e. further advanced) wins automatically.
    if(team1$round_num > team2$round_num) return(team1) else return(team2)
  }
}

simulate_round <- function(teams, game_number = 1) {
  teams <- get_bracket_matchups(teams)
  # If odd number of teams, give the last one a bye.
  if(nrow(teams) %% 2 == 1) {
    bye_team <- teams[nrow(teams), ]
    main_teams <- teams[-nrow(teams), ]
  } else {
    bye_team <- NULL
    main_teams <- teams
  }
  team_pairs <- split(main_teams, rep(1:(nrow(main_teams)/2), each = 2))
  winners <- map_df(team_pairs, function(pair) {
    simulate_game(pair[1, ], pair[2, ], game_number)
  })
  if(!is.null(bye_team)) winners <- bind_rows(winners, bye_team)
  winners
}

# Determine remaining rounds dynamically.
get_remaining_rounds <- function(current_N) {
  rounds_needed <- ceiling(log2(current_N))
  round_names <- c("Round of 64", "Round of 32", "Sweet 16", "Elite 8", "Final Four", "Title Game", "Champion")
  list(round_names = round_names[(length(round_names) - rounds_needed + 1):length(round_names)])
}

simulate_remaining_tournament <- function(bracket) {
  current_N <- nrow(bracket)
  rounds_info <- get_remaining_rounds(current_N)
  round_names <- rounds_info$round_names
  
  team_progress <- bracket %>% select(team, seed)
  for(r in round_names) {
    team_progress[[r]] <- 0
  }
  
  round_num <- 1
  teams_round <- bracket
  while(nrow(teams_round) > 1) {
    teams_round <- simulate_round(teams_round, game_number = round_num)
    if(round_num <= length(round_names)) {
      current_round_name <- round_names[round_num]
      team_progress <- team_progress %>%
        mutate("{current_round_name}" := ifelse(team %in% teams_round$team, 1, .data[[current_round_name]]))
    }
    round_num <- round_num + 1
  }
  team_progress
}

######################################
### 5. Monte Carlo Simulation        ###
######################################
n_simulations <- 10
sim_results <- map_dfr(1:n_simulations, ~ simulate_remaining_tournament(current_bracket))

# Convert probabilities to American odds.
prob_to_american <- function(prob, eps = 1e-6) {
  prob <- pmax(pmin(prob, 1 - eps), eps)
  ifelse(prob >= 0.5, -100 * prob / (1 - prob), 100 / prob - 100)
}

team_results <- sim_results %>%
  group_by(team, seed) %>%
  summarise(across(everything(), mean), .groups = "drop") %>%
  arrange(desc(Champion)) %>%
  filter(team %in% current_bracket$team) %>%  # Only include teams still in contention.
  mutate(across(-c(team, seed), prob_to_american))

######################################
### 6. Display Results               ###
######################################
bets <- team_results %>%
  gt() %>%
  tab_header(title = "Monte Carlo Simulation: In-Progress Tournament") %>%
  fmt_number(columns = vars(-team, -seed), decimals = 0)

bets
