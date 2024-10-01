library(tidyverse)
library(nflfastR)
library(ggimage)
library(gt)
library(ggthemes)
library(dplyr)
#library(ggradar)

pbp <- load_pbp(2023)

sf_rp <- pbp %>%
  filter(pass == 1 | rush == 1) %>%
  filter(!is.na(epa)) %>% 
  filter(posteam == "SF" | defteam == "SF") %>%
  mutate(winstreak = ifelse(week <=5, "1-5","6-8"))
  

results <- sf_rp %>%
  filter(posteam == "SF") %>% 
  group_by(winstreak) %>% 
  summarize(epa_play = mean(epa), success = mean(success, na.rm = TRUE),adot = mean(air_yards, na.rm = TRUE), pass_rate = mean(pass, na.rm = TRUE), cpoe = mean(cpoe, na.rm = TRUE), shotgun_rate = mean(shotgun, na.rm = TRUE))


niners<-results %>%
  mutate(pass_rate = 100*round(pass_rate,3), epa_play = round(epa_play,2), success = 100*round(success,3), shotgun_rate = 100*round(shotgun_rate,3), cpoe = round(cpoe,2), adot = round(adot,2) ) %>%
  gt() %>%
  tab_header(
    title = md("49ers Offensive Performance"),
    subtitle = md("By Callan Capitolo | @CapAnalytics7 | nflfastR"),
  ) %>% 
  cols_align(align = "center") %>%
  cols_align(align = "center") %>%
  cols_label(epa_play = "EPA/Play",
             success = "Success Rate %",
              adot = "ADoT",
              cpoe = "CPOE",
              pass_rate = "Pass Rate %",
              shotgun_rate = "Shotgun Rate %",
              winstreak = "Weeks") %>%
  gtExtras::gt_theme_538()
gtsave(niners, "49ersOffense.png") 

# Reshape data for ggplot
results_long <- gather(results, key = "Metric", value = "Value", -winstreak)

# Radar plot
ggplot(results_long, aes(x = Metric, y = Value, group = winstreak, color = winstreak)) +
  geom_point() +
  geom_line() +
  coord_polar(start = 0) +
  theme_minimal() +
  labs(title = "Radar Chart Example", subtitle = "Values for Different Categories")

pass_results <- sf_rp %>% 
  filter(posteam == "SF",pass == 1) %>% 
  group_by(winstreak) %>% 
  summarize(pass_epa_play = mean(epa), pass_success = mean(success, na.rm = T), pass_cpoe = mean(cpoe, na.rm = T), shotgun_rate = mean(shotgun, na.rm = T), qb_scramble = mean(qb_scramble, na.rm = T), qb_hit = mean(qb_hit, na.rm = T),
            sack = mean(sack,na.rm = T))

sf_rp %>% 
  filter(posteam == "SF",pass == 1) %>% 
  group_by(winstreak, shotgun) %>% 
  summarize(pass_epa_play = mean(epa), pass_success = mean(success, na.rm = T), pass_cpoe = mean(cpoe, na.rm = T), shotgun_rate = mean(shotgun, na.rm = T), qb_scramble = mean(qb_scramble, na.rm = T), qb_hit = mean(qb_hit, na.rm = T),
            sack = mean(sack,na.rm = T))

rush_results <- sf_rp %>% 
  filter(posteam == "SF",rush == 1) %>% 
  group_by(winstreak) %>% 
  summarize(run_epa_play = mean(epa), run_success = mean(success, na.rm = T), shotgun_rate = mean(shotgun, na.rm = T))

sf_rp %>% 
  filter(posteam == "SF",rush == 1, shotgun == 1) %>% 
  group_by(winstreak, run_location) %>% 
  summarize(count = n(), run_epa_play = mean(epa), run_success = mean(success, na.rm = T), shotgun_rate = mean(shotgun, na.rm = T))

sf_rp %>% 
  filter(posteam == "SF",rush == 1) %>% 
  group_by(winstreak, rusher_player_name) %>% 
  summarize(count = n(), run_epa_play = mean(epa), run_success = mean(success, na.rm = T), shotgun_rate = mean(shotgun, na.rm = T)) %>% 
  filter(count >5)


test <- pbp %>% 
  filter((pass == 1 | rush ==1) & week>5) %>% 
  group_by(posteam) %>% 
  summarize(epa_play = mean(epa), success = mean(success, na.rm =T))


#Follow Up Analysis


sf_rp_10 <- pbp %>%
  filter(pass == 1 | rush == 1) %>%
  filter(!is.na(epa)) %>% 
  filter(posteam == "SF" & week == 10)

week10_left<-sf_rp_10 %>% 
  filter(posteam == "SF",rush == 1, run_location == "left") %>% 
  group_by(run_location) %>% 
  summarize(runs = n(), run_epa_play = mean(epa), run_success = mean(success, na.rm = T)) %>% 
  mutate(weeks = "10")



comparison <- rbind(week10_left,losing_left) %>% 
  select(weeks,runs, run_epa_play, run_success)

trent_effect<-comparison %>%
  mutate(run_epa_play = round(run_epa_play,2), run_success = 100*round(run_success,3)) %>%
  gt() %>%
  tab_header(
    title = md("The Trent Williams Effect"),
    subtitle = md("49ers Efficiency Running to the Left"),
  ) %>% 
  cols_align(align = "center") %>%
  cols_align(align = "center") %>%
  cols_label(runs = "# of Runs",
            run_epa_play = "EPA/Rush",
             run_success = "Success Rate %",
             weeks = "Week") %>%
  gtExtras::gt_theme_538()
gtsave(trent_effect, "TrentWilliams.png") 

losing_left<- pbp %>% 
  filter(posteam == "SF",rush == 1, week %in% c(6,7,8), run_location == "left") %>% 
  group_by(run_location) %>% 
  summarize(runs = n(), run_epa_play = mean(epa), run_success = mean(success, na.rm = T)) %>% 
  mutate(weeks = "6-8")


