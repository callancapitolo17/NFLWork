library(dplyr)
library(ggplot2)
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
            non_group_avg_seats = mean(num_seats[enhanced_type != "Group"]), avg_non_group_spend_per_transaction= mean(full_price[enhanced_type != "Group"]),
            non_group_weighted_avg_spend = sum(num_seats[enhanced_type!= "Group"]*full_price[enhanced_type!= "Group"])/sum(num_seats[enhanced_type!= "Group"])) #visitors spend more, but less groups
enhanced_data %>% 
  group_by(enhanced_type) %>% 
  summarize(mean(full_price), weighted_avg_spend = sum(num_seats*full_price)/sum(num_seats), mean(num_seats), n())

enhanced_data %>% 
  group_by(event_name) %>% 
  summarize(total_tickets_local = sum(num_seats[local == "Local"]), total_ticket_visit = sum(num_seats[local == "Visitor"])) %>% 
  mutate(pct_local = total_tickets_local/(total_tickets_local+ total_ticket_visit))

enhanced_data %>% 
  mutate(purchase_date = as.Date(add_datetime, origin = "1899-12-30")) %>% 
  mutate(
    mmdd = substr(event_name, nchar(event_name) - 3, nchar(event_name)),        # Extract last 4 chars
    event_date = as.Date(paste0("2025", mmdd), format = "%Y%m%d")               # Assuming events are in 2025
  ) %>% 
  mutate(days_before_event = as.integer(event_date - purchase_date)) %>% 
  group_by(local) %>% 
  summarize(mean(days_before_event))

enhanced_data %>% 
  mutate(purchase_date = as.Date(add_datetime, origin = "1899-12-30")) %>% 
  mutate(
    mmdd = substr(event_name, nchar(event_name) - 3, nchar(event_name)),        # Extract last 4 chars
    event_date = as.Date(paste0("2025", mmdd), format = "%Y%m%d")               # Assuming events are in 2025
  ) %>% 
  mutate(days_before_event = as.integer(event_date - purchase_date)) %>% 
  filter(local != "Unknown") %>% 
  group_by(days_before_event,local) %>% 
  summarize(transactions = n()) %>% 
  ungroup() %>% 
  group_by(local) %>% 
  mutate(pct_transactions = transactions/sum(transactions)) %>% 
  ggplot(aes(x = days_before_event, y = pct_transactions, colour = local))+
  geom_smooth(se = F, lwd = 3)+
  labs(x = "Days Before Event of Transaction", y = "% of Transactions",
       title = "When Should SDFC Target Different Fan Groups?",
       subtitle = "Visiting fans purchase earlier than local fans",
       color = "Fan Type")+
  theme_bw()+
  theme(
    panel.grid = element_blank(),              # removes gridlines
    plot.title = element_text(hjust = 0.5, face = "bold", size = 36),    # centers the title
    plot.subtitle = element_text(hjust = 0.5,size = 30),
    axis.text = element_text(face = "bold",size = 24),
    axis.title = element_text(face = "bold", size = 24),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(colour = "black", face = "bold", size = 18)
  )+
  scale_y_continuous(labels = percent_format())
ggsave("TicketPurchasersTimePurchase.png", width = 14, height = 10, dpi = "retina", bg = "white")
  
#local fans and vistiing fans behave differently. Local fans buy earlier, spend less per ticket and buy less tickets
library(scales)
enhanced_data %>% filter(local != "Unknown") %>% 
  group_by(local) %>% 
  summarize(purchasers = n()) %>% 
  ungroup() %>% 
  mutate(pct_purchasers = purchasers/sum(purchasers)) %>% 
  ggplot(aes(x= local, y = pct_purchasers))+
  geom_bar(stat = "identity", show.legend = FALSE) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "Where Do SDFC's Ticket Purchasers Come From?",
    x = NULL,
    y = "Percent of Purchasers",
    subtitle = "Approximately 2/3 of Ticket Purchasers are Local"
  ) +
  theme_bw()+
  theme(
    panel.grid = element_blank(),              # removes gridlines
    plot.title = element_text(hjust = 0.5, face = "bold", size = 36),    # centers the title
    plot.subtitle = element_text(hjust = 0.5,size = 30),
    axis.text = element_text(face = "bold",size = 24),
    axis.title = element_text(face = "bold", size = 24)
  )
ggsave("TicketPurchasers.png", width = 14, height = 10, dpi = "retina", bg = "white")

