#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(plotly)
library(nflfastR)
library(dplyr)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("NFL Efficiency Interactive Landscape"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            numericInput(
                        "NFL_Season",
                        "NFL Season:",
                        min = 2000,
                        max = 2025,
                        value = 2020)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotlyOutput("plot1")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

  output$plot1 <- renderPlotly({
    pbp <- load_pbp(as.integer(input$NFL_Season))
    
    pbp_rp <- pbp %>%
      filter(pass == 1 | rush == 1) %>%
      filter(!is.na(epa))
    total_offensive_efficiency <- pbp_rp %>%
      group_by(posteam) %>%
      summarise(offensive_epa = mean(epa))
    
    total_defensive_efficiency <- pbp_rp %>%
      group_by(defteam) %>%
      summarise(defensive_epa = mean(epa))
    
    total_efficiency_both <- total_offensive_efficiency %>%
      left_join(total_defensive_efficiency, by = c("posteam" = "defteam"))
    
    total_efficiency_both <- total_efficiency_both %>%
      left_join(teams_colors_logos, by = c("posteam" = "team_abbr"))
    
    p<-total_efficiency_both %>%
      ggplot(aes(x = defensive_epa, y = offensive_epa,label = team_nick))+
      #geom_image(aes(image = team_logo_espn), size = 0.05, asp = 16/9)+
      geom_hline(yintercept = 0)+
      geom_vline(xintercept = 0)+
      #geom_image(aes(image = team_logo_espn), size = 0.05, asp = 16/9)+
      theme_minimal()+
      scale_x_reverse()+
      labs(x = "Defensive Efficiency",
           y = "Offensive Efficiency", title = "NFL Offensive & Defensive Efficiency", subtitle = "Efficiency Represented by EPA/Play",
           caption = "Callan Capitolo | @CapAnalytics7 | nflfastR")+
      theme(axis.text.x = element_blank(), axis.text.y = element_blank())+
      geom_text(aes(x = 0.1, y = 0.25, label = "Bad Defense, Good Offense"), size = 3) +
      geom_text(aes(x = -0.15, y = 0.25, label = "Good Defense, Good Offense"), size = 3) +
      geom_text(aes(x = 0.1, y = -0.25, label = "Bad Defense, Bad Offense"), size = 3) +
      geom_text(aes(x = -0.15, y = -0.25, label = "Good Defenese, Bad Offense"), size = 3)
    p_plotly <- ggplotly(p,tooltip = "text")
    return(p_plotly)
  })
}
#focus on just getting the graph to work
# Run the application 
shinyApp(ui = ui, server = server)
