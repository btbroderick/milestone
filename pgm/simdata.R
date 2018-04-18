### FAKE DATA FOR BRENDAN TO TEST ###
### BASED OF ATOMIC SPECS         ###
library(tidyverse)

set.seed(1234)
n <- 700
lambda_c <-  0.09589402    # 75%    , 3-YEAR DFS RATE
lambda_t <- 0.05753497     # 84.147%, 3-YEAR DFS RATE

col1 <- runif(n=700)*1184  # 18 PER MONTH ENROLLMENT RATE SO 1184 DAYS TO ENROLL 700 PATIENTS
col2 <- col1+c(365.25*rexp(n/2, rate=lambda_c), 365.25*rexp(n/2, rate=lambda_t))
col3 <- I(col2 < 1280)    # 1280 DAYS SHOULD HAVE 85 EVENTS AND 1983 DAYS SHOULD AHVE 165 EVENTS

mat <- cbind(col1, col2, col3) %>% 
  as.data.frame()

write_csv(mat, here::here("sample2.csv"))

# USE THE CURRENT DATA (mat) TO PREDICT 165 EVENTS.  THE PREDICTED DATE SHOULD BE 1983 DAYS

