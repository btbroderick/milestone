##############################################
### This is for NOD stuff
##############################################


labs <-  data.frame(clinic = c(1,2,2,2,3,3,3,3,3,4,4,4,4),
                    type = c(1,1,2,1,3,3,3,2,2,3,1,2,1),
                    date = c(10,15,15,16,10,10,10,11,11,12,12,15,17),
                    ele = c(1,1,0,1,1,1,1,1,1,1,1,0,0))


x <- labs %>% 
  group_by(clinic, date) %>% 
  mutate(count = n()) %>% 
  filter(count < 3) %>% 
  select(-count) %>% 
  ungroup() %>%  
  group_by(clinic, date, type) %>% 
  mutate(count = n()) %>% 
  filter(count < 2) %>%
  ungroup() %>% 
  group_by(clinic, date) %>% 
  arrange(desc(ele)) %>% 
  filter(row_number() == 1)

