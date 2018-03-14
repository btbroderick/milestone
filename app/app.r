library(shiny)
library(shinythemes)
library(flexsurv)
library(quadprog)
library(Hmisc)
library(msm) 

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