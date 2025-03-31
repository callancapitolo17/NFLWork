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
