library(dplyr)
data <- readxl::read_xlsx("AnalystDataUseCase.xlsx")
str(data)
summary(data) 
unique(data$event_name) #multiple events, possible game comparisons?
unique(data$ticket_type) #type of buyer, possible exploration by game?
unique(data$ticket_status) #relevant not sure? I think it's not that relevant but should explore
unique(data$full_price) #price, look at spend by other features
unique(data$zip) #locations, distance to stadium, fan segmentation --> group buyers where do they come from?
unique(data$num_seats) #number of seats purchased, relationship between seats and price paid?
unique(data$row_name) #not really important
unique(data$add_datetime) #not sure if relevant without game date
unique(data$section_name) #this could be interesting for other exploration, pricing? - not sure if relevant
unique(data$event_name) #does game have date in it? not sure, but I think so
unique(data)

local_zips <- c(
  # 919xx
  91902, 91910, 91911, 91913, 91914, 91915, 91932, 91941, 91942, 91945, 91950, 91977, 91978,
  # 920xx
  92007, 92014, 92019, 92020, 92021, 92024, 92037, 92040, 92064, 92067, 92071, 92075, 92091,
  # 921xx
  92101, 92102, 92103, 92104, 92105, 92106, 92107, 92108, 92109, 92110, 92111,
  92113, 92114, 92115, 92116, 92117, 92118, 92119, 92120,
  92121, 92122, 92123, 92124,
  92126, 92127, 92128, 92129, 92130, 92131,
  92136, 92139, 92140, 92145, 92152, 92154, 92155, 92163, 92173
)

enhanced_data <- data %>% 
  mutate(local = ifelse(zip %in% local_zips, "Local",ifelse(is.na(zip), "Unknown" ,"Visitor"))) %>% 
  mutate(enhanced_type = ifelse(grepl("Group", ticket_type), "Group", "Other"))

enhanced_data %>% 
  group_by(local) %>% 
  summarize(count = n(), avg_seats_per_transaction = mean(num_seats), avg_price = mean(full_price), weighted_avg_spend = sum(num_seats*full_price)/sum(num_seats),group = mean(enhanced_type == "Group"),
            non_group_avg_seats = mean(num_seats[enhanced_type != "Group"]), avg_non_group_spend_per_transaction= mean(full_price[enhanced_type != "Group"])) #visitors spend more, but less groups
enhanced_data %>% 
  group_by(enhanced_type) %>% 
  summarize(mean(full_price), weighted_avg_spend = sum(num_seats*full_price)/sum(num_seats), mean(num_seats), n())
