bets <- read.csv("transactions.csv")
bets %>% 
  mutate(date_placed = as.Date(time_placed, format = "%m/%d/%Y %H:%M:%S", tz = "GMT")) %>% 
  filter(date_placed > as.Date(("2025-03-19")))
