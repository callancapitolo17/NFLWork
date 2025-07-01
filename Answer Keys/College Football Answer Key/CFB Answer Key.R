library(cfbfastR)
#https://www.sportsdataverse.org/packages
test <- cfbfastR::load_cfb_schedules(year = 2024)
cfbd_betting_lines(year = 2018, week = 12, team = "Florida State")
Sys.setenv(CFBD_API_KEY = "UvbKzWwAzL8wiIg7kW8CLJRHa0ztKtrPErp4IgaBfO3TcP/4718HWzPCz/p0GhYq")
