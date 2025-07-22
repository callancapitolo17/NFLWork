library(understatr)
#https://abhiamishra.github.io/ggshakeR/articles/Guide_to_Pitch_Plots.html
#devtools::install_github("statsbomb/StatsBombR")
# remotes::install_github('ewenme/understatr')
#Wages----
man_utd_url <- "https://fbref.com/en/squads/19538871/Manchester-United-Stats"
df <- fb_squad_wages(team_urls = man_utd_url)


df %>% 
  mutate(leaving = ifelse(Player %in% c("Christian Eriksen","Victor Lindelöf", "Jonny Evans"), "Confirmed Departure","Part of Squad")) %>% 
  group_by(Age,leaving) %>% 
  summarize(total_wages = sum(AnnualWageUSD)) %>% 
  ungroup() %>% 
  mutate(pct_of_total = total_wages/sum(total_wages)) %>% 
  ggplot(aes(x = Age, y = pct_of_total, fill = leaving))+
  geom_bar(stat = "identity")+
  labs(y = "% of Total Wages", title = "What Ages are Manchester United's Wages Allocated From 2024-2025 Squad?")+
  scale_fill_manual(values = c("Confirmed Departure" = "white", 
                               "Part of Squad" = "red"))+
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color="white"),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "bold"),
        plot.title = element_text(hjust = .5, colour = "white", face = "bold", size = 32),
        plot.subtitle = element_text(hjust = .5, colour = "white", size = 24),
        plot.caption = element_text(colour = "white", size = 12),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "black", color="black"),
        panel.background = element_rect(fill = "black", color="black"),
        axis.ticks = element_line(color = "white"),
        axis.text = element_text(face = "bold", colour = "white",size = 24),
        axis.title = element_text(color = "white", size = 14),
        panel.border = element_rect(colour = "white", fill = NA, size = 1))
ggsave("UnitedWages.png", width = 14, height = 10)

#Shot Radar
shot_data <- get_player_shots(player_id = 8995) # Salah's `player_id` on Understat is `1250`
plot_shot(shot_data %>% filter(year == 2024),highlight_goals = TRUE, average_location = FALSE)
ggsave("ShotRadar.png", width = 14, height = 10)
man_city_shots <- understat_team_season_shots(team_url = "https://understat.com/team/Manchester_City/2020")
plot_shot(man_city_shots,highlight_goals = TRUE, average_location = TRUE)

#passing
library(StatsBombR)
Comp <- FreeCompetitions() %>%
  filter(competition_id == 11 & season_name == "2014/2015")
Matches <- FreeMatches(Comp)
StatsBombData <- free_allevents(MatchesDF = Matches, Parallel = TRUE)
plotting_data  <- allclean(StatsBombData)
plotting_data <- plotting_data %>%
  rename("x" = "location.x",
         "y" = "location.y",
         "finalX" = "pass.end_location.x",
         "finalY" = "pass.end_location.y")
heat_data <- plotting_data %>%
  select(x, y)

heatPlot <- plot_heatmap(data = heat_data, type = "bin")

heatPlot
plotting_data_alba  <- plotting_data %>%
  filter(type.name == "Pass" & team.name == "Barcelona" & player.name == "Jordi Alba Ramos")
plot_sonar(data = plotting_data_alba)
plotting_data_alba_single <- plotting_data_alba 

passPlot <- plot_pass(data = plotting_data_alba_single, 
                      progressive_pass = TRUE, type = "all")
passPlot

passflowPlot <- plot_passflow(data = plotting_data_alba)

passflowPlot

unique(plotting_data$match_id) # Find all match ID's from the data set

convexPlot <- plotting_data %>%
  filter(match_id == 266631,
         team.name == "Barcelona") %>%
  plot_convexhull()

convexPlot


finalData <- plotting_data %>%
  filter(team.name == "Barcelona") %>%
  group_by(player.name) %>%
  summarise_at(vars(x, y, minute), list(name = mean), na.rm = TRUE) %>%
  na.omit() %>%
  rename("x" = "x_name") %>%
  rename("y" = "y_name")
plotVoronoi <- plot_voronoi(finalData)
plotVoronoi

threatData <- calculate_epv(plotting_data)
finalData <- threatData %>%
  filter(team.name == "Barcelona") %>%
  group_by(player.name) %>%
  filter(EPVStart != -1) %>% #Getting rid of -1 as that as how we identify null values
  #We get the mean and sum of these three columns
  summarise_at(vars(x, y, EPVStart), list(name = mean, sum), na.rm = TRUE) %>%
  na.omit() %>%
  rename("x" = "x_name") %>%
  rename("y" = "y_name") %>%
  rename("EPV_Start" =  "EPVStart_fn1") #Renaming the relevant columns

#Plot voronoi plot

plotVoronoi <- plot_voronoi(data = finalData, 
                            fill = "EPV_Start", 
                            title = "Barcelona Voronoi Diagram")
plotVoronoi

unique(plotting_data$match_id) # Find all match ID's from the data set

passnetPlot <- plotting_data %>%
  filter(match_id == 266631) %>%
  plot_passnet(team_name = "Barcelona")

passnetPlot

optadata <- ggshakeR::SampleEventData

optadata <- optadata %>%
  mutate(teamId = case_when(teamId == 2151 ~ "Team 1",
                            teamId == 3070 ~ "Team 2")) %>%
  mutate(playerId = case_when(playerId == "1" ~ "Cameron",
                              playerId == "2" ~ "Sahaj",
                              playerId == "3" ~ "Trevor",
                              playerId == "4" ~ "Jolene",
                              playerId == "5" ~ "David",
                              playerId == "6" ~ "Marcus",
                              playerId == "7" ~ "Ryo",
                              playerId == "8" ~ "Anthony",
                              playerId == "9" ~ "Aabid",
                              playerId == "10" ~ "Kylian",
                              playerId == "11" ~ "Harry",
                              playerId == "12" ~ "Victor",
                              playerId == "13" ~ "Jesse",
                              playerId == "14" ~ "Borges",
                              playerId == "15" ~ "Miguel",
                              playerId == "16" ~ "Alex",
                              playerId == "17" ~ "Abdul",
                              playerId == "18" ~ "Eric",
                              playerId == "19" ~ "Marten",
                              playerId == "20" ~ "Robert",
                              playerId == "21" ~ "Lionel",
                              playerId == "22" ~ "Dean",
                              playerId == "23" ~ "Aaron",
                              playerId == "24" ~ "Benjamin",
                              playerId == "25" ~ "Diego",
                              playerId == "26" ~ "Abhishek",
                              playerId == "27" ~ "Samuel",
                              playerId == "28" ~ "Edward",
                              playerId == "29" ~ "Malik",
                              playerId == "30" ~ "Albert",
                              playerId == "31" ~ "Paul",
                              playerId == "32" ~ "Smith",
                              playerId == "33" ~ "Seth",
                              playerId == "34" ~ "Mohammed",
                              playerId == "35" ~ "Eder",
                              playerId == "36" ~ "Adam",
                              playerId == "37" ~ "Harsh",
                              playerId == "38" ~ "Jorge",
                              playerId == "39" ~ "Zak",
                              playerId == "40" ~ "Brandon"))

optadata <- optadata %>%
  rename(finalX = endX,
         finalY = endY)
passnetPlot <- plot_passnet(data = optadata, data_type = "opta", 
                            team_name = "Team 1", scale_stat = "EPV", 
                            scale_color = "#fec44f")

passnetPlot
