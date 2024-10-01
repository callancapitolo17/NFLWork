library(ggimage)
library(gt)
library(nflfastR)
library(tidyverse)
library(png)
library(nflreadr)
library(gtExtras)
library(chromote)
library(glue)
library(magick)
library(grid)
library(ggrepel)
library(cowplot)
library(patchwork)
ftn_data <- nflreadr::load_ftn_charting(2022:2023) %>%
  select(-week, -season)
pbp <- load_pbp(2022:2023)
participation <- load_participation(2022:2023) %>% 
  select(-old_game_id)
pbp <- pbp %>%
  left_join(ftn_data, by = c("game_id" = "nflverse_game_id",
                             "play_id" = "nflverse_play_id")) %>% 
  left_join(participation,by = c("game_id" = "nflverse_game_id",
                                 "play_id" = "play_id"))

nfl22 <- pbp %>% 
  filter(pass == 1 | rush == 1) %>% 
  filter(qb_kneel != 1, qb_spike !=1) %>% 
  mutate(explosive = ifelse((yards_gained>20 & pass_attempt == 1) | (yards_gained >12 & (qb_scramble == 1 | rush == 1)),1,0),
         negative = ifelse(yards_gained < 0, 1,0))

#Cluster Work----

library(cluster) # For hierarchical clustering and k-means
library(factoextra) # For visualizing clustering results
library(mclust) # For Gaussian Mixture Models

nfl99all <- load_pbp(2006:2023)
nfl99 <- nfl99all %>% 
  filter(pass == 1 | rush == 1) %>% 
  mutate(explosive = ifelse((yards_gained>20 & pass_attempt == 1) | (yards_gained >12 & (qb_scramble == 1 | rush == 1)),1,0),
         negative = ifelse(yards_gained < 0, 1,0))

season_qb_stats <- nfl99 %>% 
  filter(season>=2006) %>% 
  mutate(short_throw = ifelse(air_yards<=10,1,0),
         medium_throw = ifelse(air_yards>10&air_yards<=20,1,0),
         long_throw = ifelse(air_yards>20,1,0)) %>% 
  mutate(qb_fumbler = ifelse(is.na(fumbled_1_player_id) | fumbled_1_player_id != id,0,1)) %>% 
  group_by(id,season) %>% 
  summarize(name = max(name),dropbacks = sum(pass, na.rm = T),`EPA/Pass` = mean(epa[pass == 1 & qb_scramble!= 1],na.rm = T), aDoT = mean(air_yards,na.rm = T), `EPA/Rush` = mean(epa[rush == 1 | qb_scramble == 1]), `Success Rate/DB` = mean(success[pass == 1]), CPOE = mean(cpoe,na.rm = T),dropbacks = sum(pass),
            `Sack Rate` = mean(sack[pass==1], na.rm = T), `Explosive DB Rate` = mean(explosive[pass==1], na.rm = T),
            `Int Rate` = mean(interception[!is.na(air_yards)], na.rm = T), `% of Yards From YAC` = sum(yards_after_catch,na.rm = T)/sum(yards_gained[complete_pass == 1],na.rm = T),
            `Comp Rate` = sum(complete_pass,na.rm = T)/sum(!is.na(air_yards)), `Negative DB Rate` = mean(negative[pass == 1], na.rm = T), `Early Down EPA/DB` = mean(epa[pass == 1 & down %in% c(1,2)],na.rm = T),
            `EPA Outside Num Throws` = mean(epa[pass_location %in% c("left","right")],na.rm = T),
            `EPA Btw Num Throws` = mean(epa[pass_location == "middle"], na.rm = T), `Scramble Rate` = mean(qb_scramble[pass == 1]),`EPA/xPass DB` = mean(epa[xpass>0.85 & pass==1], na.rm = T),
            `EPA/Late Down DB` = mean(epa[down %in% c(3,4)], na.rm = T), `EPA/Trailing DB` = mean(epa[posteam_score<defteam_score],na.rm = T),
            `EPA/Tied/Winning DB` = mean(epa[posteam_score>=defteam_score],na.rm = T), `EPA/Garbage Time DB` = mean(epa[(def_wp>0.9 | def_wp < 0.1) & pass == 1],na.rm = T), `EPA/Non Garbage Time DB` = mean(epa[(def_wp <= 0.9 | def_wp >= 0.1) & pass == 1],na.rm = T), 
            `EPA/Shotgun DB` = mean(epa[pass == 1 & shotgun == 1], na.rm = T), `EPA/Under Center DB` = mean(epa[pass == 1 & shotgun == 0],na.rm =T),
            `EPA/Clutch DB` = mean(epa[abs(posteam_score-defteam_score)<=8 & game_seconds_remaining <= 300],na.rm = T),
            `WPA/Play` = mean(wpa,na.rm = T), `YAC Over Expected/Catch` = mean(yards_after_catch-xyac_mean_yardage,na.rm = T), 
            `Fumble Rate` = mean(qb_fumbler, na.rm = T), `EPA/Red Zone DB` = mean(epa[yardline_100 <= 200], na.rm = T),
            `EPA Lost/Turnover` = mean(epa[(qb_fumbler == 1 & fumble_lost == 1) | interception == 1],na.rm = T), `EPA/Scramble` = mean(epa[qb_scramble == 1], na.rm = T), `EPA/Deep Pass` = mean(epa[long_throw == 1],na.rm = T),
            `EPA/Medium Pass` = mean(epa[medium_throw == 1],na.rm = T),
            `EPA/Short Pass` = mean(epa[short_throw == 1],na.rm = T)) %>% 
  filter(dropbacks>175) %>% 
  drop_na()


season_qb_stats_numeric <- season_qb_stats %>% 
  ungroup() %>% 
  select(-id,-season,-dropbacks,-name) %>% 
  mutate(`Sack Rate` = `Sack Rate`*-1,`Int Rate` = `Int Rate`*-1,
         `Negative DB Rate` = `Negative DB Rate`*-1,
         `Fumble Rate` = `Fumble Rate`*-1,
         `EPA Lost/Turnover` = `EPA Lost/Turnover`*-1) %>% 
  mutate(across(everything(), scale))
# Apply Gaussian Mixture Model clustering
gmm_result <- Mclust(season_qb_stats_numeric)

# Add the GMM cluster results back to the original data frame
season_qb_stats$cluster_gmm <- as.factor(gmm_result$classification)

# Visualize the clusters
fviz_mclust(gmm_result, "classification", geom = "point", pointsize = 1.5)


# Summarize the statistics by clusters (e.g., for k-means)
cluster_summary <- as.data.frame(season_qb_stats) %>%
  group_by(cluster_gmm) %>%
  summarize(across(where(is.numeric), median, na.rm = TRUE))

print(cluster_summary)


test <- season_qb_stats %>% 
  select(cluster_gmm, everything()) %>% 
  mutate_if(is.numeric, ~round(.,3))

#4 Explosive passers who push down field and scrabmlers, 3 best qbs

# Visualize the means across clusters
# Extract the means for each cluster
gmm_means <- gmm_result$parameters$mean

# Convert to a data frame for easier viewing
gmm_means_df <- as.data.frame(t(gmm_means))
colnames(gmm_means_df) <- paste("Cluster", 1:ncol(gmm_means_df))

# Get the posterior probabilities for each cluster
posterior_probs <- gmm_result$z

# Assign each observation to the cluster with the highest posterior probability
season_qb_stats$gmm_posterior_cluster <- apply(posterior_probs, 1, which.max)

# Calculate the mean of each variable weighted by posterior probabilities
posterior_weighted_means <- apply(posterior_probs, 2, function(w) colMeans(season_qb_stats_numeric * w))

# Convert to data frame
posterior_weighted_means_df <- as.data.frame(t(posterior_weighted_means))
colnames(posterior_weighted_means_df) <- paste("Cluster", 1:ncol(posterior_weighted_means_df))

# View the weighted means of each variable across clusters
print(posterior_weighted_means_df)

# Calculate variance within each cluster
gmm_vars <- gmm_result$parameters$variance$sigma

# Compute the ratio of each variable's variance in each cluster to the total variance
explained_variance <- apply(gmm_vars, 1:2, function(x) sum(x) / sum(gmm_vars))

# Convert to data frame
explained_variance_df <- as.data.frame(explained_variance)
colnames(explained_variance_df) <- paste("Cluster", 1:ncol(explained_variance_df))

# View explained variance by cluster
print(explained_variance_df)

# Perform PCA on the scaled data
pca_result <- prcomp(season_qb_stats_numeric, scale. = TRUE)

# Plot the PCA results colored by GMM cluster
fviz_pca_ind(pca_result, geom.ind = "point", habillage = season_qb_stats$gmm_posterior_cluster, 
             addEllipses = TRUE, ellipse.level = 0.95, palette = "jco")

# View the PCA loadings to see which variables contribute most to each principal component
print(pca_result$rotation)

# test <- season_qb_stats %>% 
#   ungroup() %>% 
#   mutate(cluster = ifelse(cluster_gmm == 1, "Below Average Scrambler QB", ifelse(cluster_gmm == 2,"Below Average Pocket Passer",
#                                                                                  ifelse(cluster_gmm == 3, "Elite Pocket Production QB","Gunslinger Mobile QB"))))

# Comparison Building -----
qb_all_stats <- nfl22 %>% 
  mutate(qb_fumbler = ifelse(is.na(fumbled_1_player_id) | fumbled_1_player_id != id,0,1)) %>% 
  group_by(season,id) %>% 
  summarize(name = max(name),`EPA/Pass` = mean(epa[pass == 1 & qb_scramble!= 1],na.rm = T), aDoT = mean(air_yards,na.rm = T), `EPA/Rush` = mean(epa[rush == 1 | qb_scramble == 1]), `Success Rate/DB` = mean(success[pass == 1]), CPOE = mean(cpoe,na.rm = T),dropbacks = sum(pass),
          `Pressure Rate` = (sum(was_pressure,na.rm = T)+sum(sack,na.rm = T))/sum(pass,na.rm=T), `Explosive DB Rate` = mean(explosive[pass==1], na.rm = T),
          `Int Rate` = mean(interception[!is.na(air_yards)], na.rm = T), `% of Yards From YAC` = sum(yards_after_catch,na.rm = T)/sum(yards_gained[complete_pass == 1],na.rm = T),
          `Comp Rate` = sum(complete_pass,na.rm = T)/sum(!is.na(air_yards)), `Negative DB Rate` = mean(negative[pass == 1], na.rm = T), `Early Down EPA/DB` = mean(epa[pass == 1 & down %in% c(1,2)],na.rm = T),
          `EPA/Outside Num Pass` = mean(epa[pass_location %in% c("left","right")],na.rm = T),
          `EPA/Btw Num Pass` = mean(epa[pass_location == "middle"], na.rm = T), `Scramble Rate` = mean(qb_scramble[pass == 1]),`EPA/xPass DB` = mean(epa[xpass>0.85 & pass==1], na.rm = T),
          `EPA/Late Down DB` = mean(epa[down %in% c(3,4)], na.rm = T), `EPA/4Q DB` = mean(epa[qtr>=4 & pass == 1]), `EPA/Trailing DB` = mean(epa[posteam_score<defteam_score],na.rm = T),
          `EPA/Tied/Winning DB` = mean(epa[posteam_score>=defteam_score],na.rm = T), `EPA/Garbage Time DB` = mean(epa[(def_wp>0.9 | def_wp < 0.1) & pass == 1],na.rm = T), `EPA/Non Garbage Time DB` = mean(epa[(def_wp <= 0.9 | def_wp >= 0.1) & pass == 1],na.rm = T), 
          `EPA/Shotgun DB` = mean(epa[pass == 1 & shotgun == 1], na.rm = T), `EPA/Under Center DB` = mean(epa[pass == 1 & shotgun == 0],na.rm =T),`EPA/Motion DB` = mean(epa[is_motion == 1 & pass == 1],na.rm =T), `EPA/No Motion DB` = mean(epa[is_motion == 0 & pass == 1],na.rm = T),
            , `EPA/Out of Pocket DB` = mean(epa[is_qb_out_of_pocket == 1 & pass == 1], na.rm = T),
            `EPA/In Pocket DB` = mean(epa[is_qb_out_of_pocket == 0 & pass == 1], na.rm = T), `Int Worthy Throw Rate` = mean(is_interception_worthy[!is.na(air_yards)] , na.rm = T),
            `Catchable Ball Rate` = mean(is_catchable_ball[!is.na(air_yards) & is_throw_away == 0], na.rm = T),`Catchable Ball Drop Rate` = mean(is_drop[!is.na(air_yards) & is_catchable_ball == 1],na.rm = T),
            `Contested Throw Rate` = mean(is_contested_ball[!is.na(air_yards)]), `EPA/Blitz DB` = mean(epa[n_blitzers>0],na.rm = T), `EPA/Non Blitz DB` = mean(epa[n_blitzers<=0],na.rm = T), 
            `EPA/1st Read DB` = mean(epa[read_thrown %in% c(1,"DES")],na.rm = T), `EPA/Non 1st Read DB` = mean(epa[read_thrown %in% c(2,"SD","CHK",0)],na.rm = T),
          `Shotgun Rate` = mean(qb_location == "S",na.rm = T), `EPA/Shotgun DB` = mean(epa[qb_location == "S" & pass == 1],na.rm = T),
          `Under Center Rate` = mean(qb_location == "U",na.rm = T),`EPA/Under Center DB` = mean(epa[qb_location == "U" & pass == 1],na.rm = T),
          `QB Sack Fault Rate` = mean(is_qb_fault_sack[sack == 1],na.rm = T), `EPA/Clutch DB` = mean(epa[abs(posteam_score-defteam_score)<=8 & game_seconds_remaining <= 300],na.rm = T),
          `EPA/PA DB` = mean(epa[pass == 1 & is_play_action == 1],na.rm = T), `EPA/ Non PA DB` = mean(epa[pass == 1 & is_play_action == 0],na.rm = T),
          `WPA/Play` = mean(wpa,na.rm = T), `YAC Over Expected/Catch` = mean(yards_after_catch-xyac_mean_yardage,na.rm = T), 
          `EPA/Light Box DB` = mean(epa[n_defense_box <7],na.rm = T), `Fumble Rate` = mean(qb_fumbler, na.rm = T), `EPA/Red Zone DB` = mean(epa[yardline_100 <= 200], na.rm = T),
          `EPA Lost/Turnover` = mean(epa[(qb_fumbler == 1 & fumble_lost == 1) | interception == 1],na.rm = T), `EPA/Scramble` = mean(epa[qb_scramble == 1], na.rm = T), `EPA/Deep Pass` = mean(epa[air_yards>=20],na.rm = T),
          `Pressure to Sack Rate` = sum(sack,na.rm = T)/(sum(was_pressure,na.rm = T)+sum(sack,na.rm = T)), `EPA/Pressure` = mean(epa[was_pressure == TRUE | sack == 1],na.rm = T),
          `Avg Time to Throw` = mean(time_to_throw,na.rm = T)) %>% 
  filter(dropbacks>60) %>% 
  mutate(name = paste (name,season, sep = " ")) %>% 
  ungroup %>% 
  select(-id,-dropbacks,-season)


qb_all_stats_clustered <- season_qb_stats %>% 
  mutate(name = paste(name,season, sep = " ")) %>% 
  ungroup() %>% 
  mutate(cluster = ifelse(cluster_gmm == 1, "Below Average Scrambler", ifelse(cluster_gmm == 2,"Below Average Pocket Passer",
                                                                              ifelse(cluster_gmm == 3, "Elite Pocket Production","Gunslinger & Mobile")))) %>% 
  select(name,cluster) %>% 
  right_join(qb_all_stats,by = c("name"))

replace_with_values_and_ranks <- function(column) {
  values <- column
  ranks <- rank(column*-1,ties.method = "max")
}



qb_rank <- as.data.frame(cbind(qb_all_stats_clustered$name,(apply(qb_all_stats_clustered %>% 
                                                 select(-name,-cluster) %>% 
                                                 mutate(`Pressure Rate` = `Pressure Rate`*-1,`Int Rate` = `Int Rate`*-1,
                                                        `Negative DB Rate` = `Negative DB Rate`*-1, `Int Worthy Throw Rate` = `Int Worthy Throw Rate` * -1, `QB Sack Fault Rate` = `QB Sack Fault Rate`*-1,
                                                        `Fumble Rate` = `Fumble Rate`*-1,`Pressure to Sack Rate` = `Pressure to Sack Rate`*-1, `Contested Throw Rate` = `Contested Throw Rate`*-1), 2, replace_with_values_and_ranks)))) %>% 
  pivot_longer(cols = -`V1`, names_to = "statistic", values_to = "rank") %>% 
  mutate(rank = as.numeric(rank))


qbs_all_comb<- qb_all_stats_clustered %>% pivot_longer(cols = c(-name,-cluster), names_to = "statistic", values_to = "value") %>% 
  inner_join(qb_rank, by = c("name" = "V1", "statistic"))





          
          
           
          
#deep throw rate, epa 4th, light box epa
          

production_stats <- c("EPA/Pass", "EPA/Rush", "aDoT", "Avg Time to Throw", "Success Rate/DB", "CPOE","Explosive DB Rate","Negative DB Rate","Pressure to Sack Rate","EPA/xPass DB","WPA/Play","% of Yards From YAC")
advanced_prod_stats = c("EPA/Outside Num Pass", "EPA/Btw Num Pass","Comp Rate","Pressure Rate", "Int Rate", "Int Worthy Throw Rate",
                        "Fumble Rate", "EPA Lost/Turnover","YAC Over Expected/Catch", "Catchable Ball Rate", "Catchable Ball Drop Rate", "Contested Throw Rate")
situational_stats <- c("Early Down EPA/DB","EPA/Late Down DB", "EPA/Trailing DB", "EPA/Tied/Winning DB", "EPA/Garbage Time DB", "EPA/Non Garbage Time DB", "EPA/Clutch DB", "EPA/Red Zone DB",
                       "Scramble Rate", "EPA/Scramble","EPA/1st Read DB", "EPA/Non 1st Read DB")
playcall_stats <- c("EPA/Motion DB", "EPA/No Motion DB", "Shotgun Rate", "EPA/Shotgun DB","EPA/Under Center DB",
                    "EPA/PA DB", "EPA/ Non PA DB", "EPA/Pressure" ,"EPA/Blitz DB", "EPA/Non Blitz DB","EPA/Out of Pocket DB", "EPA/In Pocket DB")

player1 <- "A.Dalton 2023" #Rework to be more flexible and incorporate more players
player2 <- "B.Young 2023"
#Production Stats ----
prod_qbs <- qbs_all_comb %>% 
  filter(statistic %in% production_stats) %>% 
  mutate(statistic = factor(statistic, levels = production_stats))

temp <- (360/n_distinct(prod_qbs$statistic))/2

myAng <- seq(-temp, -360+temp, length.out = n_distinct(prod_qbs$statistic))

ang <- ifelse(myAng < -90, myAng+180,myAng)

ang <- ifelse(ang < -90, ang+180, ang)


pizza_prod<- prod_qbs %>% 
  mutate(rank = max(rank) + 1 - rank) %>% 
  filter(name %in% c(player1,player2)) %>% 
  mutate(name = factor(name, levels = c(player1, player2))) %>% 
  arrange(name) %>% 
  ggplot(aes(x = statistic, y = rank, fill = paste(name,cluster, sep = "\n"), label = value))+
  geom_bar(stat = "identity", position = position_identity(), alpha = 0.6)+
  # geom_label(color = "white", size=2.5, fontface="bold", show.legend = FALSE, position = position_jitterdodge())+
  coord_polar()+
  geom_bar(aes(y = max(qbs_all_comb$rank)/n_distinct(name)),stat = "identity", width =1, alpha = 0.1, fill = "grey")+
  geom_hline(yintercept = seq(1, max(qbs_all_comb$rank), by = max(qbs_all_comb$rank)),
             color = "white",
             size = 1)+
  geom_vline(xintercept = seq(.5, n_distinct(prod_qbs$statistic), by = 1),
             color = "white",
             size = .5)+
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color="white"),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "bold"),
        plot.title = element_text(hjust = .5, colour = "white", face = "bold", size = 16),
        plot.subtitle = element_text(hjust = .5, colour = "white", size = 8),
        plot.background = element_rect(fill = "black", color="black"),
        panel.background = element_rect(fill = "black", color="black"),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(face = "bold", size = 8, colour = "white"),
        axis.title = element_blank(),
        axis.text.x = element_text(face = "bold", size = 12, angle = ang)) +
  labs(x = NULL, y = NULL)+
  scale_fill_brewer(palette = "Set1")+
  annotate("text", x = (pi/12) * 2, y = seq(10,n_distinct(qbs_all_comb$name), by = 10), label = seq(10*n_distinct(qbs_all_comb$name)%/%10, 10 ,by = -10), hjust = 1.15, 
           color = "White", size = 5)


tab_prod <- prod_qbs %>%
  filter(name %in% c(player1,player2)) %>%
  select(-cluster) %>% 
  mutate(value = round(value,3)) %>% 
  mutate(name = factor(name, levels = c(player1, player2))) %>% 
  pivot_wider(names_from = name, values_from = c(value,rank)) %>% 
  rename("Value1" = paste("value_",player1,sep = ""), "Rank1" = paste("rank_",player1,sep = ""), "Value2" = paste("value_",player2,sep = ""),"Rank2" = paste("rank_",player2,sep = "")) %>%
  select(statistic,Value1,Rank1,Value2,Rank2) %>% 
  mutate(statistic = factor(statistic, levels = production_stats)) %>% 
  arrange(statistic) %>% 
  gt() %>% 
  cols_align(align = "center") %>% 
  tab_spanner(label = player1, columns = c(Value1,Rank1)) %>% 
  tab_spanner(label = player2, columns = c(Value2,Rank2)) %>% 
  cols_label(Value1 = "Value", Rank1 = "Rank",Value2 = "Value",Rank2 = "Rank") %>% 
  gtExtras::gt_theme_538() %>% 
  tab_options(
    table.background.color = "black",        # Set the entire table background to black
    heading.background.color = "black",      # Set the header background to black
    column_labels.background.color = "black", # Set the column label background to black
    row_group.background.color = "black",    # Set row group background to black (if any)
    summary_row.background.color = "black",  # Set summary row background to black (if any)
    grand_summary_row.background.color = "black", # Set grand summary row background to black (if any)
    footnotes.background.color = "black",    # Set footnotes background to black
    source_notes.background.color = "black", # Set source notes background to black
    table.border.top.color = "black",        # Set table top border to black
    table.border.bottom.color = "black",     # Set table bottom border to black
    heading.border.bottom.color = "black",   # Set header bottom border to black
    column_labels.border.top.color = "black",# Set column label top border to black
    column_labels.border.bottom.color = "black" # Set column label bottom border to black
  ) %>% 
  gt_hulk_col_numeric(columns = c(Rank1,Rank2),reverse = TRUE) %>%
  tab_style(
    style = cell_text(size = px(16), weight = "bold", color = "white"),  # Change font and size for column labels
    locations = cells_column_labels(columns = everything())
  ) %>% 
  tab_style(
    style = cell_text(color = "white", size = px(16)),  # Change font and size for the body text
    locations = cells_body(columns = c(statistic, Value1, Value2))
  )

if (exists("f") && inherits(f, "ChromoteSession")) {
  try(f$shutdown(), silent = TRUE)
}

# Start a new session
f <- ChromoteSession$new()

gtsave(tab_prod, "prod_temp_table.png")

table_image <- image_read("prod_temp_table.png")
table_image_transparent <- image_transparent(table_image, "white")
table_grob <- rasterGrob(table_image, interpolate = TRUE)
spacer <- ggplot() + theme_void() + theme(panel.background = element_rect(fill = "black"))

pizza_prod + table_grob+spacer + plot_layout(ncol = 3, widths = c(6,3,.1))& 
  theme(
    plot.background = element_rect(fill = "black", color = "black"),
    panel.background = element_rect(fill = "black", color = "black"),
    plot.caption = element_text(size = 10, hjust = 0.5, face = "italic", color = "white"),
    plot.title = element_text(size =20,hjust = 0.5, face = "italic", color = "white"),
    plot.subtitle = element_text(hjust = 0.5, face = "bold", color = "white"))&
  plot_annotation(
    caption = glue("Compared to {n_distinct(qb_all_stats$name)} QB seasons with 60+ dropbacks in 2022-2023"),
    title = "2022-2023 QB Production Comparison",
    subtitle =  "@CapAnalytics7 | nflfastR"
  )
ggsave("QB_Comp_Prod.png", bg = "black", ,width = 14, height =10, dpi = "retina") 

#Advanced Prod Stats ----
advanced_prod <- qbs_all_comb %>% 
  filter(statistic %in% advanced_prod_stats) %>% 
  mutate(statistic = factor(statistic, levels = advanced_prod_stats))

temp <- (360/n_distinct(prod_qbs$statistic))/2

myAng <- seq(-temp, -360+temp, length.out = n_distinct(advanced_prod$statistic))

ang <- ifelse(myAng < -90, myAng+180,myAng)

ang <- ifelse(ang < -90, ang+180, ang)

pizza_prod_advanced<- advanced_prod %>% 
  mutate(rank = max(rank) + 1 - rank) %>% 
  filter(name %in% c(player1,player2)) %>% 
  mutate(name = factor(name, levels = c(player1, player2))) %>% 
  ggplot(aes(x = statistic, y = rank, fill = paste(name,cluster, sep = "\n"), label = value))+
  geom_bar(stat = "identity", position = position_identity(), alpha = 0.6)+
  # geom_label(color = "white", size=2.5, fontface="bold", show.legend = FALSE, position = position_jitterdodge())+
  coord_polar()+
  geom_bar(aes(y = max(qbs_all_comb$rank)/n_distinct(name)),stat = "identity", width =1, alpha = 0.1, fill = "grey")+
  geom_hline(yintercept = seq(1, max(qbs_all_comb$rank), by = max(qbs_all_comb$rank)),
             color = "white",
             size = 1)+
  geom_vline(xintercept = seq(.5, n_distinct(prod_qbs$statistic), by = 1),
             color = "white",
             size = .5)+
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color="white"),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "bold"),
        plot.title = element_text(hjust = .5, colour = "white", face = "bold", size = 16),
        plot.subtitle = element_text(hjust = .5, colour = "white", size = 8),
        plot.background = element_rect(fill = "black", color="black"),
        panel.background = element_rect(fill = "black", color="black"),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(face = "bold", size = 8, colour = "white"),
        axis.title = element_blank(),
        axis.text.x = element_text(face = "bold", size = 12, angle = ang)) +
  labs(x = NULL, y = NULL)+
  scale_fill_brewer(palette = "Set1")+
  annotate("text", x = (pi/12) * 2, y = seq(10,n_distinct(qbs_all_comb$name), by = 10), label = seq(10*n_distinct(qbs_all_comb$name)%/%10, 10 ,by = -10), hjust = 1.15, 
           color = "White", size = 5)


tab_prod_adv <- advanced_prod %>%
  filter(name %in% c(player1,player2)) %>% 
  select(-cluster) %>% 
  mutate(value = round(value,3)) %>% 
  mutate(name = factor(name, levels = c(player1, player2))) %>% 
  pivot_wider(names_from = name, values_from = c(value,rank)) %>% 
  rename("Value1" = paste("value_",player1,sep = ""), "Rank1" = paste("rank_",player1,sep = ""), "Value2" = paste("value_",player2,sep = ""),"Rank2" = paste("rank_",player2,sep = "")) %>%
  select(statistic,Value1,Rank1,Value2,Rank2) %>% 
  mutate(statistic = factor(statistic, levels = advanced_prod_stats)) %>% 
  arrange(statistic) %>% 
  gt() %>% 
  cols_align(align = "center") %>% 
  tab_spanner(label = player1, columns = c(Value1,Rank1)) %>% 
  tab_spanner(label = player2, columns = c(Value2,Rank2)) %>% 
  cols_label(Value1 = "Value", Rank1 = "Rank",Value2 = "Value",Rank2 = "Rank") %>% 
  gtExtras::gt_theme_538() %>% 
  tab_options(
    table.background.color = "black",        # Set the entire table background to black
    heading.background.color = "black",      # Set the header background to black
    column_labels.background.color = "black", # Set the column label background to black
    row_group.background.color = "black",    # Set row group background to black (if any)
    summary_row.background.color = "black",  # Set summary row background to black (if any)
    grand_summary_row.background.color = "black", # Set grand summary row background to black (if any)
    footnotes.background.color = "black",    # Set footnotes background to black
    source_notes.background.color = "black", # Set source notes background to black
    table.border.top.color = "black",        # Set table top border to black
    table.border.bottom.color = "black",     # Set table bottom border to black
    heading.border.bottom.color = "black",   # Set header bottom border to black
    column_labels.border.top.color = "black",# Set column label top border to black
    column_labels.border.bottom.color = "black" # Set column label bottom border to black
  ) %>% 
  gt_hulk_col_numeric(columns = c(Rank1,Rank2),reverse = TRUE) %>%
  tab_style(
    style = cell_text(size = px(16), weight = "bold", color = "white"),  # Change font and size for column labels
    locations = cells_column_labels(columns = everything())
  ) %>% 
  tab_style(
    style = cell_text(color = "white", size = px(16)),  # Change font and size for the body text
    locations = cells_body(columns = c(statistic, Value1, Value2))
  )

if (exists("f") && inherits(f, "ChromoteSession")) {
  try(f$shutdown(), silent = TRUE)
}

# Start a new session
f <- ChromoteSession$new()

gtsave(tab_prod_adv, "prod_adv_temp_table.png")

table_image <- image_read("prod_adv_temp_table.png")
table_image_transparent <- image_transparent(table_image, "white")
table_grob_adv <- rasterGrob(table_image, interpolate = TRUE)
spacer <- ggplot() + theme_void() + theme(panel.background = element_rect(fill = "black"))

pizza_prod_advanced + table_grob_adv+spacer + plot_layout(ncol = 3, widths = c(6,3,.1))& 
  theme(
    plot.background = element_rect(fill = "black", color = "black"),
    panel.background = element_rect(fill = "black", color = "black"),
    plot.caption = element_text(size = 10, hjust = 0.5, face = "italic", color = "white"),
    plot.title = element_text(size =20,hjust = 0.5, face = "italic", color = "white"),
    plot.subtitle = element_text(hjust = 0.5, face = "bold", color = "white"))&
  plot_annotation(
    caption = glue("Compared to {n_distinct(qb_all_stats$name)} QB seasons with 60+ dropbacks in 2022-2023"),
    title = "2022-2023 QB Advanced Production Comparison",
    subtitle =  "@CapAnalytics7 | nflfastR"
  )
ggsave("QB_Comp_Prod_Adv.png", bg = "black", ,width = 14, height =10, dpi = "retina") 

#Situation Radar -----
sit_stat <- qbs_all_comb %>% 
  filter(statistic %in% situational_stats) %>% 
  mutate(statistic = factor(statistic, levels = situational_stats))

temp <- (360/n_distinct(sit_stat$statistic))/2

myAng <- seq(-temp, -360+temp, length.out = n_distinct(sit_stat$statistic))

ang <- ifelse(myAng < -90, myAng+180,myAng)

ang <- ifelse(ang < -90, ang+180, ang)
radial_labels <- seq(n_distinct(qb_all_stats$name),1,by = -5)

label_data <- data.frame(
  x = rep(0, length(radial_labels)),  # Center at x = 0
  y = radial_labels,  # Position along the radial axis (y)
  label = radial_labels  # Labels for the radial positions
)

pizza_prod_advanced<- sit_stat %>% 
  mutate(rank = max(rank) + 1 - rank) %>% 
  filter(name %in% c(player1,player2)) %>% 
  mutate(name = factor(name, levels = c(player1, player2))) %>% 
  ggplot(aes(x = statistic, y = rank, fill = paste(name,cluster, sep = "\n"), label = value))+
  geom_bar(stat = "identity", position = position_identity(), alpha = 0.6)+
  # geom_label(color = "white", size=2.5, fontface="bold", show.legend = FALSE, position = position_jitterdodge())+
  coord_polar()+
  geom_bar(aes(y = max(qbs_all_comb$rank)/n_distinct(name)),stat = "identity", width =1, alpha = 0.1, fill = "grey")+
  geom_hline(yintercept = seq(1, max(qbs_all_comb$rank), by = max(qbs_all_comb$rank)),
             color = "white",
             size = 1)+
  geom_vline(xintercept = seq(.5, n_distinct(sit_stat$statistic), by = 1),
             color = "white",
             size = .5)+
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color="white"),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "bold"),
        plot.title = element_text(hjust = .5, colour = "white", face = "bold", size = 16),
        plot.subtitle = element_text(hjust = .5, colour = "white", size = 8),
        plot.background = element_rect(fill = "black", color="black"),
        panel.background = element_rect(fill = "black", color="black"),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(face = "bold", size = 8, colour = "white"),
        axis.title = element_blank(),
        axis.text.x = element_text(face = "bold", size = 12, angle = ang)) +
  labs(x = NULL, y = "Rank")+
  scale_fill_brewer(palette = "Set1")+
  # geom_text(aes(y = rank, label = round(rank, 1)), color = "white", size = 5, fontface = "bold", position = position_identity())
  annotate("text", x = (pi/12) * 2, y = seq(10,n_distinct(qbs_all_comb$name), by = 10), label = seq(10*n_distinct(qbs_all_comb$name)%/%10, 10 ,by = -10), hjust = 1.15, 
           color = "White", size = 5)

sit_stat_tab <- sit_stat %>%
  filter(name %in% c(player1,player2)) %>% 
  mutate(value = round(value,3)) %>% 
  select(-cluster) %>% 
  mutate(name = factor(name, levels = c(player1, player2))) %>% 
  pivot_wider(names_from = name, values_from = c(value,rank)) %>% 
  rename("Value1" = paste("value_",player1,sep = ""), "Rank1" = paste("rank_",player1,sep = ""), "Value2" = paste("value_",player2,sep = ""),"Rank2" = paste("rank_",player2,sep = "")) %>%
  select(statistic,Value1,Rank1,Value2,Rank2) %>% 
  arrange(statistic) %>% 
  gt() %>% 
  cols_align(align = "center") %>% 
  tab_spanner(label = player1, columns = c(Value1,Rank1)) %>% 
  tab_spanner(label = player2, columns = c(Value2,Rank2)) %>% 
  cols_label(Value1 = "Value", Rank1 = "Rank",Value2 = "Value",Rank2 = "Rank") %>% 
  gtExtras::gt_theme_538() %>% 
  tab_options(
    table.background.color = "black",        # Set the entire table background to black
    heading.background.color = "black",      # Set the header background to black
    column_labels.background.color = "black", # Set the column label background to black
    row_group.background.color = "black",    # Set row group background to black (if any)
    summary_row.background.color = "black",  # Set summary row background to black (if any)
    grand_summary_row.background.color = "black", # Set grand summary row background to black (if any)
    footnotes.background.color = "black",    # Set footnotes background to black
    source_notes.background.color = "black", # Set source notes background to black
    table.border.top.color = "black",        # Set table top border to black
    table.border.bottom.color = "black",     # Set table bottom border to black
    heading.border.bottom.color = "black",   # Set header bottom border to black
    column_labels.border.top.color = "black",# Set column label top border to black
    column_labels.border.bottom.color = "black" # Set column label bottom border to black
  ) %>% 
  gt_hulk_col_numeric(columns = c(Rank1,Rank2),reverse = TRUE) %>%
  tab_style(
    style = cell_text(size = px(16), weight = "bold", color = "white"),  # Change font and size for column labels
    locations = cells_column_labels(columns = everything())
  ) %>% 
  tab_style(
    style = cell_text(color = "white", size = px(16)),  # Change font and size for the body text
    locations = cells_body(columns = c(statistic, Value1, Value2))
  )

if (exists("f") && inherits(f, "ChromoteSession")) {
  try(f$shutdown(), silent = TRUE)
}

# Start a new session
f <- ChromoteSession$new()

gtsave(sit_stat_tab, "sit_table.png")

table_image <- image_read("sit_table.png")
table_image_transparent <- image_transparent(table_image, "white")
table_grob_adv <- rasterGrob(table_image, interpolate = TRUE)
spacer <- ggplot() + theme_void() + theme(panel.background = element_rect(fill = "black"))

pizza_prod_advanced + table_grob_adv+spacer + plot_layout(ncol = 3, widths = c(6,3,.1))& 
  theme(
    plot.background = element_rect(fill = "black", color = "black"),
    panel.background = element_rect(fill = "black", color = "black"),
    plot.caption = element_text(size = 10, hjust = 0.5, face = "italic", color = "white"),
    plot.title = element_text(size =20,hjust = 0.5, face = "italic", color = "white"),
    plot.subtitle = element_text(hjust = 0.5, face = "bold", color = "white"))&
  plot_annotation(
    caption = glue("Compared to {n_distinct(qb_all_stats$name)} QB seasons with 60+ dropbacks in 2022-2023"),
    title = "2022-2023 QB Situational Comparison",
    subtitle =  "@CapAnalytics7 | nflfastR"
  )
ggsave("QB_Sit_Stat.png", bg = "black", ,width = 14, height =10, dpi = "retina") 


# Playcall Stats----
play_stat <- qbs_all_comb %>% 
  filter(statistic %in% playcall_stats) %>% 
  mutate(statistic = factor(statistic, levels = playcall_stats)) #gets stats to be laid out in order

#angles of text
temp <- (360/n_distinct(play_stat$statistic))/2 

myAng <- seq(-temp, -360+temp, length.out = n_distinct(play_stat$statistic))

ang <- ifelse(myAng < -90, myAng+180,myAng)

ang <- ifelse(ang < -90, ang+180, ang)

pizza_prod_play<- play_stat %>% 
  mutate(rank = max(rank) + 1 - rank) %>% 
  filter(name %in% c(player1,player2)) %>% 
  mutate(name = factor(name, levels = c(player1, player2))) %>% 
  ggplot(aes(x = statistic, y = rank, fill = paste(name,cluster, sep = "\n"), label = value))+
  # geom_label(color = "white", size=2.5, fontface="bold", show.legend = FALSE, position = position_jitterdodge())+
  coord_polar()+
  geom_bar(aes(y = max(qbs_all_comb$rank)/n_distinct(name)),stat = "identity", width =1, alpha = 0.1, fill = "grey")+
  geom_bar(stat = "identity", position = position_identity(), alpha = 0.6)+#makes bars correctly
  geom_hline(yintercept = seq(1, max(qbs_all_comb$rank), by = max(qbs_all_comb$rank)),
             color = "white",
             size = 1)+
  geom_vline(xintercept = seq(.5, n_distinct(play_stat$statistic), by = 1),
             color = "white",
             size = .5)+
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color="white"),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "bold"),
        plot.title = element_text(hjust = .5, colour = "white", face = "bold", size = 16),
        plot.subtitle = element_text(hjust = .5, colour = "white", size = 8),
        plot.background = element_rect(fill = "black", color="black"),
        panel.background = element_rect(fill = "black", color="black"),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(face = "bold", size = 8, colour = "white"),
        axis.title = element_blank(),
        axis.text.x = element_text(face = "bold", size = 12, angle = ang)) +
  labs(x = NULL, y = NULL)+
  scale_fill_brewer(palette = "Set1")+
  annotate("text", x = (pi/12) * 2, y = seq(10,n_distinct(qbs_all_comb$name), by = 10), label = seq(10*n_distinct(qbs_all_comb$name)%/%10, 10 ,by = -10), hjust = 1.15, 
           color = "White", size = 5)


play_tab <- play_stat %>%
  filter(name %in% c(player1,player2)) %>% 
  mutate(value = round(value,3)) %>% 
  mutate(name = factor(name, levels = c(player1, player2))) %>% 
  select(-cluster) %>% 
  pivot_wider(names_from = name, values_from = c(value,rank)) %>% 
  rename("Value1" = paste("value_",player1,sep = ""), "Rank1" = paste("rank_",player1,sep = ""), "Value2" = paste("value_",player2,sep = ""),"Rank2" = paste("rank_",player2,sep = "")) %>%
  select(statistic,Value1,Rank1,Value2,Rank2) %>% 
  arrange(statistic) %>% 
  gt() %>% 
  cols_align(align = "center") %>% 
  tab_spanner(label = player1, columns = c(Value1,Rank1)) %>% 
  tab_spanner(label = player2, columns = c(Value2,Rank2)) %>% 
  cols_label(Value1 = "Value", Rank1 = "Rank",Value2 = "Value",Rank2 = "Rank") %>% 
  gtExtras::gt_theme_538() %>% 
  tab_options(
    table.background.color = "black",        # Set the entire table background to black
    heading.background.color = "black",      # Set the header background to black
    column_labels.background.color = "black", # Set the column label background to black
    row_group.background.color = "black",    # Set row group background to black (if any)
    summary_row.background.color = "black",  # Set summary row background to black (if any)
    grand_summary_row.background.color = "black", # Set grand summary row background to black (if any)
    footnotes.background.color = "black",    # Set footnotes background to black
    source_notes.background.color = "black", # Set source notes background to black
    table.border.top.color = "black",        # Set table top border to black
    table.border.bottom.color = "black",     # Set table bottom border to black
    heading.border.bottom.color = "black",   # Set header bottom border to black
    column_labels.border.top.color = "black",# Set column label top border to black
    column_labels.border.bottom.color = "black" # Set column label bottom border to black
  ) %>% 
  gt_hulk_col_numeric(columns = c(Rank1,Rank2),reverse = TRUE) %>%
  tab_style(
    style = cell_text(size = px(16), weight = "bold", color = "white"),  # Change font and size for column labels
    locations = cells_column_labels(columns = everything())
  ) %>% 
  tab_style(
    style = cell_text(color = "white", size = px(16)),  # Change font and size for the body text
    locations = cells_body(columns = c(statistic, Value1, Value2))
  )

if (exists("f") && inherits(f, "ChromoteSession")) {
  try(f$shutdown(), silent = TRUE)
}

# Start a new session
f <- ChromoteSession$new()

gtsave(play_tab, "play_table.png")

table_image <- image_read("play_table.png")
table_image_transparent <- image_transparent(table_image, "white")
table_grob_play <- rasterGrob(table_image, interpolate = TRUE)
spacer <- ggplot() + theme_void() + theme(panel.background = element_rect(fill = "black"))

pizza_prod_play + table_grob_play+spacer + plot_layout(ncol = 3, widths = c(6,3,.1))& 
  theme(
    plot.background = element_rect(fill = "black", color = "black"),
    panel.background = element_rect(fill = "black", color = "black"),
    plot.caption = element_text(size = 10, hjust = 0.5, face = "italic", color = "white"),
    plot.title = element_text(size =20,hjust = 0.5, face = "italic", color = "white"),
    plot.subtitle = element_text(hjust = 0.5, face = "bold", color = "white"))&
  plot_annotation(
    caption = glue("Compared to {n_distinct(qb_all_stats$name)} QB seasons with 60+ dropbacks in 2022-2023"),
    title = "2022-2023 QB Play Type Production Comparison",
    subtitle =  "@CapAnalytics7 | nflfastR"
  )
ggsave("QB_Playcall_Stat.png", bg = "black", ,width = 14, height =10, dpi = "retina") 







