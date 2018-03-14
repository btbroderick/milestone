###########################################
### Setup
###########################################
library(shiny)
library(shinythemes)
library(flexsurv)
library(quadprog)
library(Hmisc)
library(msm) 
source(here::here("pgm","utilsBayes1.r"))
source(here::here("pgm","utilsFreq.r"))
source(here::here("pgm","utilsWts.r"))

###########################################
### User Interface
###########################################
ui <- fluidPage(
  theme = shinytheme("united"),
  navbarPage(
    title = "Milestone"
  )
)
##########################################
### Server
##########################################
server <- function(input, output, session) {
  
}
##########################################
### Knit app together 
##########################################
shinyApp(ui, server)