library(tidyverse)

################################################
## Generate Sample file 
################################################

### STUDY SPECIFIC, THE CSV FILE SHOULD CONTAIN 3 COLUMN, ONSTUDY, EVENT_TIME, AND EVENT INDICATOR
### THE ONSTUDY AND EVENT TIME WILL NEED TO BE CONVERTED SO FIRST PATINET ON STUDY IS DAY 0 (NOT CALENDAR DATE )


sample <- data.frame(
  onstudy = c(floor(runif(500, 1, 365)), rep(NA,500)),
  event_time = floor(c(rep(0, 500), runif(500, 1, 365))),
  event_indicator = as.integer(c(rep(0,500), rep(1,500)))
) %>% 
  mutate(onstudy = ifelse(event_time > 0 , event_time - floor(log(event_time)), onstudy)) %>% 
  write_csv(here::here("sample.csv"))



