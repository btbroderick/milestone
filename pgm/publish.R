library(rsconnect)
library(rlang)
rsconnect::deployApp(appDir ="/people/biostat6/m146014/consult/milestone",
                     appFiles = c("app.r",
                                  "/www/", 
                                  "/screenshots/",
                                  "/pgm"),
                     server = "rcfshinypoc.mayo.edu",
                     account="M146014",
                     launch.browser = F,
                     appName = "Milestone Prediction")