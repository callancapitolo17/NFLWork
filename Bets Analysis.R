bets <- read.csv("transactions.csv")
bets %>% 
  mutate(date_placed = as.Date(time_placed, format = "%m/%d/%Y %H:%M:%S", tz = "GMT")) %>% 
  filter(date_placed > as.Date(("2025-03-19"))) %>% 
  filter(tags != "NFL draft ", tags != "March ",sports != "American Football") %>% 
  group_by(is.na(ev)) %>% 
  summarize(sum(profit), count = n(), ep = sum(amount * ev,na.rm = T))

bets %>% 
  mutate(date_placed = as.Date(time_placed, format = "%m/%d/%Y %H:%M:%S", tz = "GMT")) %>% 
  filter(date_placed > as.Date(("2025-03-19"))) %>% 
  filter(tags != "NFL draft ", tags != "March ",sports != "American Football") %>% 
  filter(is.na(ev)) %>% 
  group_by(sportsbook) %>% 
  summarize(sum(profit), count = n(), ep = sum(amount * ev,na.rm = T))

bets %>% 
  mutate(date_placed = as.Date(time_placed, format = "%m/%d/%Y %H:%M:%S", tz = "GMT")) %>% 
  filter(date_placed > as.Date(("2025-03-19"))) %>% 
  filter(tags != "NFL draft ", tags != "March ",sports != "American Football") %>% 
  filter(!is.na(ev)) %>% 
  group_by(sportsbook) %>% 
  summarize(sum(profit), count = n(), wagered = sum(amount),ep = sum(amount * ev,na.rm = T), xRoi = ep/wagered)
