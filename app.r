###########################################
### Setup
###########################################
library(shiny)
library(shinythemes)
library(DT)
library(tidyverse)
library(flexsurv)
library(quadprog)
library(Hmisc)
library(msm) 
library(VGAM)
source(here::here("pgm","utilsBayes1.r"))
source(here::here("pgm","utilsFreq.r"))
source(here::here("pgm","utilsWts.r"))
source(here::here("pgm", "helper.R"))
###########################################
### User Interface
###########################################

ui <- fluidPage(
  theme = shinytheme("united"),
  navbarPage(
    title = "Milestone",
    tabPanel("Main",
             sidebarPanel(
               numericInput("nE", label = "Landmark Event Number", value = NA, width = 300),
               numericInput("N", label = "Number of simultions", value = NA, width = 300),
               dateInput("study_date",label = "Study start date", value = NA, width = 300,format = "mm-dd-yyyy"),
               tags$h6("Date format: mm-dd-yyyy"),
               HTML("<br/>"),
               fileInput("inputfile", NULL, buttonLabel = "Upload", multiple = FALSE, width = 300),
               tags$h6("*File upload format can be found in the 'About' tab")
             ),
             mainPanel(
               tabsetPanel(
                 tabPanel("Data View",
                  dataTableOutput("data_view"),
                  tags$h3(textOutput("data_checks"))
                 )
               ),
               tabPanel("Calculate Milestone")
             )
    ),
    tabPanel("About")
  )
)

##########################################
### Server
##########################################

server <- function(input, output, session) {
  
  inputData <- eventReactive(input$inputfile, {
    read <- input$inputfile
    if (is.null(read)){
      return()
    }
    read_any_file(read$datapath)
  })
  
  data_check_text <- eventReactive(input$inputfile,{
    if (is.numeric(inputData()[[1]]) & !anyNA(inputData()[[1]])) {
      check1 <- "Column 1 is OK"
      }
    else {
      "Bad"
    }
  })
  
  output$data_checks <- renderText({
    data_check_text()
  })
  
  output$data_view <- renderDataTable({
    inputData()
  }, rownames= FALSE)
}

##########################################
### Knit app together 
##########################################
shinyApp(ui, server)