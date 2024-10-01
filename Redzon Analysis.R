library(tidyverse)
library(nflfastR)
library(vip)
library(ggimage)

pbp_2023 <- load_pbp(2023)
pbp_2023<-pbp_2023 %>% 
  mutate(success_rate = ifelse(
    (down == 1 & yards_gained/ydstogo >= 0.5) |
      (down == 2 & yards_gained/ydstogo >= 0.7) |
      (down %in% c(3, 4) & yards_gained/ydstogo >= 1),
    1,
    0
  ))
pbp_rp_2023 <- pbp_2023 %>%
  filter(pass == 1 | rush == 1) %>%
  filter(!is.na(epa))
pbp_red <- pbp_rp_2023 %>% 
  filter(yardline_100 <=20)
pbp_not_red <- pbp_rp_2023 %>% 
  filter(yardline_100 >20)

red_eff <- pbp_red %>% 
  group_by(posteam) %>% 
  summarize(count = n(), red_epa_play = mean(epa), red_adot = mean(air_yards), touchdown = mean(touchdown, na.rm = T))

not_red_eff <- pbp_not_red %>% 
  group_by(posteam) %>% 
  summarize(count = n(), epa_play = mean(epa), adot = mean(air_yards))

total_eff <- not_red_eff %>% 
  left_join(red_eff,"posteam")

total_eff <- total_eff %>%
  left_join(teams_colors_logos, by = c("posteam" = "team_abbr"))

total_eff %>%
  ggplot(aes(x = epa_play, y = red_epa_play))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_image(aes(image = team_logo_espn), size = 0.05, asp = 16/9)+
  theme_bw()+
  labs(x = "Non Red Zone EPA/Play",
       y = "Red Zone EPA/Play", title = "Red Zone Offensive Performance vs Non Red Zone Offensive Performance Following Week 9 TNF",
       caption = "Callan Capitolo | @CapAnalytics7 | nflfastR")
ggsave("Week9RedZoneEfficiency.png", width = 14, height =10, dpi = "retina")

phi_red_eff <- pbp_red %>% 
  filter(posteam == "PHI") %>% 
  group_by(posteam) %>% 
  summarize(count = n(), epa_play = mean(epa), adot = mean(air_yards, na.rm = T), success = mean(success, na.rm = T), tds = sum(touchdown), pass_rate = mean(pass,na.rm = T),
            run_rate = mean(rush, na.rm = T), cpoe = mean(cpoe, na.rm = T), shotgun_rate = mean(shotgun, na.rm = T), qb_scramble = mean(qb_scramble, na.rm = T)) %>% 
  mutate(red = "Yes")

pbp_red %>% 
  filter(posteam == "PHI", rush ==1) %>% 
  group_by(as.factor(shotgun), rusher_player_name) %>% 
  summarize(epa_play = mean(epa), success = mean(success), count = n())

pbp_red %>% 
  filter(posteam == "PHI", pass ==1) %>% 
  group_by(as.factor(shotgun)) %>% 
  summarize(epa_play = mean(epa), success = mean(success), count = n())

phi_not_red_eff <- pbp_not_red %>% 
  filter(posteam == "PHI") %>% 
  group_by(posteam) %>% 
  summarize(count = n(), epa_play = mean(epa), adot = mean(air_yards, na.rm = T), success = mean(success, na.rm = T), tds = sum(touchdown, na.rm = T), pass_rate = mean(pass,na.rm = T),
            run_rate = mean(rush, na.rm = T), cpoe = mean(cpoe, na.rm = T), shotgun_rate = mean(shotgun, na.rm = T), qb_scramble = mean(qb_scramble, na.rm = T)) %>% 
  mutate(red = "No")

pbp_not_red %>% 
  filter(posteam == "PHI", rush ==1) %>% 
  group_by(posteam) %>% 
  summarize(epa_play = mean(epa), success = mean(success))

pbp_not_red %>% 
  filter(posteam == "PHI", pass ==1) %>% 
  group_by(posteam) %>% 
  summarize(epa_play = mean(epa), success = mean(success))
rbind(phi_red_eff,phi_not_red_eff)


defense = pbp_red %>% 
  group_by(defteam) %>% 
  summarize(count = n(), red_epa_play = mean(epa), red_adot = mean(air_yards), touchdown = mean(touchdown, na.rm = T))


install.packages("ggplot2")
# Load necessary packages
library(ggplot2)

# Create a sample dataset with adjustments
set.seed(123)
n_points <- 10
time_steps <- 50
data <- data.frame(
  ID = rep(1:n_points, each = time_steps),
  Time = rep(1:time_steps, times = n_points),
  X = rnorm(n_points * time_steps),
  Y = rnorm(n_points * time_steps),
  Adjustment_X = rnorm(n_points * time_steps, mean = 2, sd = 0.5),  # Adjustments in X
  Adjustment_Y = rnorm(n_points * time_steps, mean = -2, sd = 0.5)  # Adjustments in Y
)

# Apply adjustments to the X and Y coordinates
data$X_adjusted <- data$X + data$Adjustment_X
data$Y_adjusted <- data$Y + data$Adjustment_Y

# Create a scatter plot with tails
scatter_plot_with_tails <- ggplot(data, aes(x = X_adjusted, y = Y_adjusted)) +
  geom_path(aes(group = ID, color = ID), alpha = 0.5) +  # Connect points to create tails
  geom_point(aes(color = ID), size = 3) +  # Plot points
  labs(title = "Scatter Plot with Tails") +
  theme_minimal()

# Display the scatter plot with tails
print(scatter_plot_with_tails)


# Render and save the animation as a GIF

