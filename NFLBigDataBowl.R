# the tidyverse and dplyr!!
suppressWarnings(library(tidyverse))
suppressWarnings(library(dplyr))
suppressWarnings(library(stringr))
suppressWarnings(library(gganimate))

# Input data files are available in the read-only "../input/" directory
# big data bowl files are in ../input/nfl-big-data-bowl-2025 so set to variable
data_path <- '../input/nfl-big-data-bowl-2025'

# Function definition! 
#' Load base NFL data files (excluding tracking data)
#' 
#' @param data_path Character string specifying the base path where data files are stored
#' @return List containing data frames for games, plays, players, and player_play data
#' @export
load_base_data <- function(data_path) {
  tryCatch({
    games <- read_csv(file.path(data_path, "games.csv"), show_col_types = FALSE)
    plays <- read_csv(file.path(data_path, "plays.csv"), show_col_types = FALSE)
    players <- read_csv(file.path(data_path, "players.csv"), show_col_types = FALSE)
    player_play <- read_csv(file.path(data_path, "player_play.csv"), show_col_types = FALSE)
    
    list(
      games = games,
      plays = plays,
      players = players,
      player_play = player_play
    )
  }, error = function(e) {
    message("Error loading data: ", e$message)
    return(NULL)
  })
}

# Load our base data real quick.
base_data <- 
  suppressWarnings(load_base_data(data_path))

# Function definition!
#' Load tracking data for specified weeks
#' 
#' @param data_path Character string specifying the base path where tracking data files are stored
#' @param weeks Numeric vector of weeks to load
#' @return Combined tracking data for specified weeks
#' @export
load_tracking_data <- function(data_path, weeks) {
  # Ensure weeks is a vector
  if (!is.vector(weeks)) weeks <- c(weeks)
  
  # Try to load each week's data
  tracking_data <- map_df(weeks, function(week) {
    file_path <- file.path(data_path, sprintf("tracking_week_%d.csv", week))
    
    tryCatch({
      read_csv(file_path, show_col_types = FALSE)
    }, error = function(e) {
      warning(sprintf("Could not load tracking data for week %d", week))
      return(NULL)
    })
  })
  
  if (nrow(tracking_data) == 0) return(NULL)
  return(tracking_data)
}

# not going to load in any tracking info just yet. 
# but we'll show names of base_data just for clarity.
names(base_data)

# Function definition!
#' Join the base_data sets together
#' 
#' @param base_data list of data frame objects: `plays`, `games`, `player_play`, `players`
#' @export
create_pbp_full <- function(base_data){
  x <- 
    base_data$plays %>%
    dplyr::left_join(., base_data$games, by='gameId') %>%
    ## join player_plays on gameId and playId
    dplyr::left_join(., base_data$player_play, by=c('gameId', 'playId')) %>%
    ## join players on nflId. i don't think we really need this but sure.
    dplyr::left_join(., base_data$players, by='nflId') %>%
    dplyr::mutate(team=teamAbbr)
  return(x)
  
}

# Function definition!
#' gets a singular play or frame from the full dataset 
#' 
#' @param base_data list of data frame objects: `plays`, `games`, `player_play`, `players`
#' @param week - week of season 1-9
#' @param gameid - gameId desired
#' @param playid - play desired
#' @param frameid - frame desired. if 0 then entire play data is returned. 
#' @export
get_pre_snap_moment <- function(base_data, week, gameid, playid, frameid){
  ## join df_base and filter!!
  pbp_full <-
    create_pbp_full(base_data) %>%
    dplyr::filter(gameId==gameid, playId==playid)
  
  tracking_df <- 
    load_tracking_data(data_path, week) %>%
    ## Add columns to data
    ## going to take advantage of right triangles and trig here
    ## -derive constants we can use to flip coordinates - delta x/y sign
    ## -force new direction (dir_new) to be between 0-90 degrees
    ## -use new direction and constants to derive endpoints for arrows 
    ##    assuming length of arrow should be proportional to acceleration
    ##    recall: sohcahtoa! sin(theta) = opposite_length/hypotenuse_length
    ##                       cos(theta) = adjacent_length/hypotenuse_length
    ## there's probably a better way to do this directly in ggplot w/ the dir
    ## variable but here we are.
    ## does this belong here or in plotting function? for now we keep here. 
    mutate(deltaxsign=if_else(dir<=180,1,-1),
           deltaysign=if_else(dir <= 90 | dir >= 270, 1, -1),
           dir_new = if_else(dir <=90,
                             dir, 
                             if_else(dir<=180, 
                                     180-dir, 
                                     if_else(dir<=270, 270-dir, 360-dir)))
    )
  
  if(frameid > 0) {
    tracking_df <-
      tracking_df %>%
      dplyr::filter(frameId==frameid)
  }
  
  
  return(pbp_full %>% 
           dplyr::inner_join(., tracking_df, by=c('gameId', 'playId', 'nflId')))
}

# Function definition!
#' Create our directional data viz!
#' 
#' @param playid - play desired. 
#' @param frameid - frame desired. if 0 then entire play data is returned. 
#' @param df_plot - dataframe for plots - singular frame 
#' @param position_groups_to_show - 'all', 'offense', 'defense'. who to show vectors for.
#' @param scale_arrows_by - increase or decrease length of vectors
#' @param movement_vector - use acceleration or speed as vector basis. 'a' or 's'
#' @param max_play_desc_len - maximum length of the play description to show in subtitle.
#' @export
build_directional_plot <- function(playid,
                                   frameid,                                
                                   df_plot, 
                                   position_groups_to_show='all', 
                                   scale_arrows_by = 1, 
                                   movement_vector='a',
                                   max_play_desc_len=80){
  
  
  if(movement_vector=='a'){
    df_plot$scalar <-  df_plot$a
    title <- 'Pre-Snap Acceleration Vectors for'
  }
  if(movement_vector=='s'){
    df_plot$scalar <-  df_plot$s
    title <- 'Pre-Snap Speed Vectors for'
  }
  
  sample_play <- 
    df_plot %>% 
    dplyr::filter(playId==playid, frameId==frameid) %>%
    mutate(xnew = 2*scale_arrows_by*scalar*cos(dir_new*pi/180) * deltaxsign + x,
           ynew= 2*scale_arrows_by*scalar*sin(dir_new*pi/180) * deltaysign + y,
           playdescriptiontruncated=str_sub(playDescription, 1, max_play_desc_len)
    )
  
  names(sample_play) <-
    tolower(names(sample_play))
  
  title <-
    paste(title,
          sample_play %>% distinct(hometeamabbr) %>% pull(.), 
          'v', 
          sample_play %>% distinct(visitorteamabbr) %>% pull(.), 
          sample_play %>% distinct(gamedate) %>% pull(.), 
          sep = ' ')
  
  subtitle <-
    paste0(
      paste0(' FrameId: ', 
             sample_play %>% distinct(frameid) %>% pull(.)),
      paste0(' Event: ', 
             sample_play %>% distinct(event) %>% pull(.)),
      #paste0(' PlayId: ', 
      #       sample_play %>% distinct(playid) %>% pull(.)),
      paste0(' PlayDesc: ', 
             sample_play %>% distinct(playdescriptiontruncated) %>% pull(.)),
      sep='\n'
    )
  
  ## ideally config
  def_positions <-
    c('SS','DE','ILB','FS','CB','DT', 'NT', 'OLB')
  
  ## ideally config
  off_positions <-
    c('WR','TE', 'T', 'QB', 'RB', 'C', 'G', 'FB')
  
  positions_to_show <- 
    c(off_positions, def_positions)
  
  if(position_groups_to_show %in% c('all', 'ALL', 'All')){
    positions_to_show <- 
      c(off_positions, def_positions)
  }
  
  if(position_groups_to_show %in% c('Offense', 'O', 'offense', 'o')){
    positions_to_show <- 
      off_positions
  }
  
  if(position_groups_to_show %in% c('Defense', 'D', 'defense','d')){
    positions_to_show <- 
      def_positions
  }
  
  ## lets build the plot!! taking code from another notebook to help build the field 
  ## https://www.kaggle.com/jesucristo/animated-visualization
  
  # General field boundaries
  xmin <- 0
  xmax <- 160/3
  hash.right <- 38.35
  hash.left <- 12
  hash.width <- 3.3
  
  
  ## Specific boundaries for a given play
  ymin <- max(round(min(sample_play$x, na.rm = TRUE) - 10, -1), 0)
  ymax <- min(round(max(sample_play$x, na.rm = TRUE) + 10, -1), 120)
  df.hash <- expand.grid(x = c(0, 23.36667, 29.96667, xmax), y = (10:110))
  df.hash <- df.hash %>% filter(!(floor(y %% 5) == 0))
  df.hash <- df.hash %>% filter(y < ymax, y > ymin)
  
  ## build plot
  ggplot() +
    scale_size_manual(values = c(6, 4, 6), guide = FALSE) + 
    scale_shape_manual(values = c(21, 16, 21), guide = FALSE) +
    scale_fill_manual(values = c("#e31837", "#654321", "#002244"), guide = FALSE) + 
    scale_colour_manual(values = c("black", "#654321", "#c60c30"), guide = FALSE) + 
    annotate("text", x = df.hash$x[df.hash$x < 55/2], 
             y = df.hash$y[df.hash$x < 55/2], label = "_", hjust = 0, vjust = -0.2, na.rm=TRUE) + 
    annotate("text", x = df.hash$x[df.hash$x > 55/2], 
             y = df.hash$y[df.hash$x > 55/2], label = "_", hjust = 1, vjust = -0.2,na.rm=TRUE) + 
    annotate("segment", x = xmin, 
             y = seq(max(10, ymin), min(ymax, 110), by = 5), 
             xend =  xmax, 
             yend = seq(max(10, ymin), min(ymax, 110), by = 5),na.rm=TRUE) + 
    annotate("text", x = rep(hash.left, 11), y = seq(10, 110, by = 10), 
             label = c("G   ", seq(10, 50, by = 10), rev(seq(10, 40, by = 10)), "   G"), 
             angle = 270, size = 4,na.rm=TRUE) + 
    annotate("text", x = rep((xmax - hash.left), 11), y = seq(10, 110, by = 10), 
             label = c("   G", seq(10, 50, by = 10), rev(seq(10, 40, by = 10)), "G   "), 
             angle = 90, size = 4,na.rm=TRUE) + 
    annotate("segment", x = c(xmin, xmin, xmax, xmax), 
             y = c(ymin, ymax, ymax, ymin), 
             xend = c(xmin, xmax, xmax, xmin), 
             yend = c(ymax, ymax, ymin, ymin), colour = "black",na.rm=TRUE) + 
    geom_point(data = sample_play, 
               aes(x = (xmax-y), y = x, 
                   shape = team, fill = team, 
                   group = nflid, size = team, 
                   color = team), 
               alpha = 1, 
               size=8, 
               na.rm=TRUE) + 
    geom_segment(data = sample_play %>% filter(position %in% positions_to_show),
                 aes(x = (xmax-y), y = x, 
                     xend = (xmax-ynew), yend = xnew),
                 arrow = arrow(angle=20, type='closed', length=unit(0.2, 'cm')),na.rm=TRUE) +
    geom_text(data = sample_play, 
              aes(x = (xmax-y), y = x, 
                  label = jerseynumber), 
              color = "white", 
              vjust = 0.3, 
              size = 5, na.rm=TRUE) + 
    labs(
      x=NULL,
      y=NULL,
      title = title,
      subtitle = str_wrap(subtitle,80),
      caption = paste0("Data: 2025 Big Data Bowl Pre-Snap Tracking Data\n", 
                       "GameId: ", sample_play %>% distinct(gameid) %>% pull(.), 
                       '\n',
                       "PlayId: ", sample_play %>% distinct(playid) %>% pull(.),
                       '\n',
                       "Viz: Zeno Muscarella @ZenoIsMyName_ on Twitter/X"
      )
    ) +
    ylim(ymin, ymax) +
    theme_bw() 
  
}

# set the gamid, playid, and week. frameid set to 0.
gameid <- 2022091808 #2022100901
playid <- 3360 #3981    
week   <- 2 #5    
frameid<- 0 ## set to zero so we get the full dataset.

## get our full dataset for the play
df_play_full <-
  get_pre_snap_moment(base_data, week, gameid, playid, 0)

## choose a fame and create plot
suppressWarnings(build_directional_plot(frameid=1, 
                                        playid=playid,
                                        df_plot=df_play_full, 
                                        position_groups_to_show='all', 
                                        scale_arrows_by = 2, 
                                        movement_vector='a',
                                        max_play_desc_len=80))
build_directional_plot(frameid=53, ## frame of the snap event
                       playid=playid,
                       df_plot=df_play_full, 
                       position_groups_to_show='all', 
                       scale_arrows_by = 2,
                       movement_vector='a',
                       max_play_desc_len=80)

build_directional_gif <- function(playid, 
                                  df_plot, 
                                  position_groups_to_show='all', 
                                  scale_arrows_by = 1, 
                                  movement_vector = 'a' ,
                                  max_play_desc_len=80){
  
  if(movement_vector=='a'){
    df_plot$scalar <-  df_plot$a
    title <- 'Pre-Snap Acceleration Vectors for'
  }
  if(movement_vector=='s'){
    df_plot$scalar <-  df_plot$s
    title <- 'Pre-Snap Speed Vectors for' # Lazy? Sure am!
  }
  
  sample_play <- 
    df_plot %>% 
    dplyr::filter(playId==playid) %>%
    mutate(xnew = 2*scale_arrows_by*scalar*cos(dir_new*pi/180) * deltaxsign + x,
           ynew= 2*scale_arrows_by*scalar*sin(dir_new*pi/180) * deltaysign + y,
           playdescriptiontruncated=str_sub(playDescription, 1, max_play_desc_len)
    )
  
  names(sample_play) <-
    tolower(names(sample_play))
  
  title <-
    paste(title,
          sample_play %>% distinct(hometeamabbr) %>% pull(.), 
          'v', 
          sample_play %>% distinct(visitorteamabbr) %>% pull(.), 
          sample_play %>% distinct(gamedate) %>% pull(.), 
          sep = ' ')
  
  
  sample_play$playframeid <- 
    as.factor(sample_play$frameid)
  
  title <-
    paste(title,
          sample_play %>% distinct(hometeamabbr) %>% pull(.), 
          'v', 
          sample_play %>% distinct(visitorteamabbr) %>% pull(.), 
          sample_play %>% distinct(gamedate) %>% pull(.), 
          sep = ' ')
  
  subtitle <-
    paste0(
      #paste0(' FrameId: ', 
      #       sample_play %>% distinct(frameid) %>% pull(.)),
      #paste0(' Event: ', 
      #       sample_play %>% distinct(event) %>% pull(.)),
      #paste0(' PlayId: ', 
      #       sample_play %>% distinct(playid) %>% pull(.)),
      paste0(' PlayDesc: ', 
             sample_play %>% distinct(playdescriptiontruncated) %>% pull(.)),
      sep='\n'
    )
  
  ## ideally config
  def_positions <-
    c('SS','DE','ILB','FS','CB','DT', 'NT', 'OLB')
  
  ## ideally config
  off_positions <-
    c('WR','TE', 'T', 'QB', 'RB', 'C', 'G', 'FB')
  
  positions_to_show <- 
    c(off_positions, def_positions)
  
  if(position_groups_to_show %in% c('all', 'ALL', 'All')){
    positions_to_show <- 
      c(off_positions, def_positions)
  }
  
  if(position_groups_to_show %in% c('Offense', 'O', 'offense', 'o')){
    positions_to_show <- 
      off_positions
  }
  
  if(position_groups_to_show %in% c('Defense', 'D', 'defense','d')){
    positions_to_show <- 
      def_positions
  }
  
  ## lets build the plot!! taking code from another notebook to help build the field 
  ## https://www.kaggle.com/jesucristo/animated-visualization
  
  # General field boundaries
  xmin <- 0
  xmax <- 160/3
  hash.right <- 38.35
  hash.left <- 12
  hash.width <- 3.3
  
  
  ## Specific boundaries for a given play
  ymin <- max(round(min(sample_play$x, na.rm = TRUE) - 10, -1), 0)
  ymax <- min(round(max(sample_play$x, na.rm = TRUE) + 10, -1), 120)
  df.hash <- expand.grid(x = c(0, 23.36667, 29.96667, xmax), y = (10:110))
  df.hash <- df.hash %>% filter(!(floor(y %% 5) == 0))
  df.hash <- df.hash %>% filter(y < ymax, y > ymin)
  
  ## build plot
  p <-
    ggplot() +
    scale_size_manual(values = c(6, 4, 6), guide = FALSE) + 
    scale_shape_manual(values = c(21, 16, 21), guide = FALSE) +
    scale_fill_manual(values = c("#e31837", "#654321", "#002244"), guide = FALSE) + 
    scale_colour_manual(values = c("black", "#654321", "#c60c30"), guide = FALSE) + 
    annotate("text", x = df.hash$x[df.hash$x < 55/2], 
             y = df.hash$y[df.hash$x < 55/2], label = "_", hjust = 0, vjust = -0.2, na.rm=TRUE) + 
    annotate("text", x = df.hash$x[df.hash$x > 55/2], 
             y = df.hash$y[df.hash$x > 55/2], label = "_", hjust = 1, vjust = -0.2,na.rm=TRUE) + 
    annotate("segment", x = xmin, 
             y = seq(max(10, ymin), min(ymax, 110), by = 5), 
             xend =  xmax, 
             yend = seq(max(10, ymin), min(ymax, 110), by = 5),na.rm=TRUE) + 
    annotate("text", x = rep(hash.left, 11), y = seq(10, 110, by = 10), 
             label = c("G   ", seq(10, 50, by = 10), rev(seq(10, 40, by = 10)), "   G"), 
             angle = 270, size = 4,na.rm=TRUE) + 
    annotate("text", x = rep((xmax - hash.left), 11), y = seq(10, 110, by = 10), 
             label = c("   G", seq(10, 50, by = 10), rev(seq(10, 40, by = 10)), "G   "), 
             angle = 90, size = 4,na.rm=TRUE) + 
    annotate("segment", x = c(xmin, xmin, xmax, xmax), 
             y = c(ymin, ymax, ymax, ymin), 
             xend = c(xmin, xmax, xmax, xmin), 
             yend = c(ymax, ymax, ymin, ymin), colour = "black",na.rm=TRUE) + 
    geom_point(data = sample_play, 
               aes(x = (xmax-y), y = x, 
                   shape = team, fill = team, 
                   group = nflid, size = team, 
                   color = team), 
               alpha = 1, 
               size=8, 
               na.rm=TRUE) + 
    geom_segment(data = sample_play %>% filter(position %in% positions_to_show),
                 aes(x = (xmax-y), y = x, 
                     xend = (xmax-ynew), yend = xnew),
                 arrow = arrow(angle=20, type='closed', length=unit(0.2, 'cm')),na.rm=TRUE) +
    geom_text(data = sample_play, 
              aes(x = (xmax-y), y = x, 
                  label = jerseynumber), 
              color = "white", 
              vjust = 0.3, 
              size = 5, na.rm=TRUE) + 
    labs(
      x=NULL,
      y=NULL,
      title = title,
      subtitle = str_wrap(subtitle,80),
      caption = paste0("Data: 2025 Big Data Bowl Pre-Snap Tracking Data\n", 
                       "GameId: ", 
                       sample_play %>% distinct(gameid) %>% pull(.), 
                       '\n',
                       "PlayId: ", 
                       sample_play %>% distinct(playid) %>% pull(.),
                       '\n',
                       "Viz: Zeno Muscarella @ZenoIsMyName_ on Twitter/X"
      )
    ) +
    ylim(ymin, ymax) +
    theme_bw()
  
  # now for animation!
  p + 
    gganimate::transition_states(playframeid,
                                 transition_length = 2,
                                 state_length = 2)
  
}

anim <- build_directional_gif(playid=playid, 
                              df_plot=df_play_full, 
                              position_groups_to_show='all', 
                              scale_arrows_by = 1.5, 
                              movement_vector = 'a',
                              max_play_desc_len=80)
anim
gganimate::anim_save('play.gif', anim)