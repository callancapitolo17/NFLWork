library(plotly)
library(nflfastR)
library(tidyverse)
# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("NFL 2023-24 Landscape by Week"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          # selectInput("year", "NFL Season Since 1999:", choices = 1999:2023, selected = 2023),
          selectInput("xside", "X Axis Side of Ball:", choices = c("Offense","Defense"), selected = "Defense"),
          selectInput("xstat", "X Axis Stat:", choices = c("yards_gained","success","epa"), selected = "epa"),
          selectInput("yside", "Y Axis Side of Ball:", choices = c("Offense","Defense"), selected = "Offense"),
          selectInput("ystat", "Y Axis Stat:", choices = c("yards_gained","success","epa"), selected = "epa")),

        # Show a plot of the generated distribution
        mainPanel(
           plotlyOutput("plot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$plot <- renderPlotly({
    pbp_rp <- load_pbp(2023) %>%
      select(week,yards_gained,success,epa,posteam,defteam,pass,rush) %>%
      filter(pass == 1 | rush == 1,!is.na(epa))
    # pbp_rp <- pbp_rp %>% 
    #   filter(season == as.numeric(input$year))
    x_side <- ifelse(input$xside == "Defense","defteam","posteam")
    y_side <- ifelse(input$yside == "Defense","defteam","posteam")
    x_df <- data.frame()  # Initialize an empty list to store results
    
    for (i in c(1:max(pbp_rp$week))) {
      result <- pbp_rp %>%
        filter(week %in% c(1:i)) %>% 
        group_by(!!sym(x_side)) %>%
        summarize(cumulative_x_stat = round(mean(!!sym(input$xstat),na.rm =T),3)) %>% 
        mutate(week = i) %>% 
        do(data.frame(.))
      x_df <- bind_rows(x_df, result)
    }
    x_df <- x_df %>% 
      rename("team" = x_side)
    
    y_df <- data.frame()  # Initialize an empty list to store results
    
    for (i in c(1:max(pbp_rp$week))) {
      result <- pbp_rp %>%
        filter(week %in% c(1:i)) %>% 
        group_by(!!sym(y_side)) %>%
        summarize(cumulative_y_stat = round(mean(!!sym(input$ystat),na.rm =T),3)) %>% 
        mutate(week = i) %>% 
        do(data.frame(.))
      y_df <- bind_rows(y_df, result)
    }
    y_df %>% 
      rename("team" = y_side) %>%
      inner_join(x_df, by = c("team","week")) %>% 
      left_join(teams_colors_logos[c("team_abbr","team_color")], by = c("team" = "team_abbr")) %>% 
    plot_ly(y = ~cumulative_y_stat, x = ~cumulative_x_stat, frame = ~week,
            showlegend = F) %>% 
      add_text(
        text = ~team,  # Specify the text labels
        x = ~cumulative_x_stat, y = ~cumulative_y_stat,  # Coordinates for text annotations
        textposition = "top center",  # Adjust text position as needed
        showlegend = FALSE,  # Hide legend for text annotations
        textfont = list(color = ~team_color)
      ) %>%
      layout(xaxis = list(title = paste("Cumulative",input$xside,ifelse(input$xstat == "yards_gained","Yards",
                                                           ifelse(input$xstat == "success","Success Rate", input$xstat)),"Per Play", sep =" ")), 
             yaxis = list(title = paste("Cumulative",input$yside,ifelse(input$ystat == "yards_gained","Yards",
                                                           ifelse(input$ystat == "success","Success Rate", input$ystat)),"Per Play", sep =" ")),
             title = list(text = "NFL Landscape"),
             margin = list(l = 50, r = 10, b = 50, t = 50)
             )  %>% 
      animation_opts(
        1250, easing = "cubic-in-out", redraw = FALSE
      )
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
