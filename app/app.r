###########################################
### Setup
###########################################
library(shiny)
library(shinythemes)
library(tidyverse)
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
    title = "Milestone",
    tabPanel("Main",
             sidebarPanel(
               numericInput("nE", label = "Landmark Event Number", value = NA),
               numericInput("N", )
               fileInput("inputfile", NULL, buttonLabel = "Upload", multiple = FALSE),
               tags$h6("*File upload format can be found in the 'About' tab")
             ),
             mainPanel(
               tableOutput("data_view")
             )
    ),
    tabPanel("About")
  )
)

##########################################
### Server
##########################################

server <- function(input, output, session) {
  output$data_view <- renderTable({
    read <- input$inputfile
    if (is.null(read)){
      return()
    }
    read_csv(read$datapath)
  })
}

##########################################
### Knit app together 
##########################################
shinyApp(ui, server)