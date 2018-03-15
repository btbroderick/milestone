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
library(rmeta)
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
               numericInput("nE", label = "Landmark Event Number", value = 1000, width = 300),
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
                 ),
                 tabPanel("Calculate Milestone",
                          actionButton("calculate", label = "Run Milestone Prediction"),
                          plotOutput("forestPlot",width = "200%")
                 )
               )
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
  
  forestPlot <- eventReactive(input$calculate,{
    withProgress(value = 0, message =  "Calculating", {
      nE <- input$nE # landmark event number
      tempdat <- inputData()
      dat <- cbind(tempdat[[1]], tempdat[[2]]*tempdat[[3]])
      lambda <- 0.0003255076
      
      #Priors
      # Weibull prior, mean and varaince for lambda and k
      wP< - c(lambda, 50, 1, 50)
      
      # Gompertz prior, mean and variance for eta and b
      b <- lambda*log(log(2)+1)/log(2)
      gP <- c(1, 50, b, 50)
      
      # Lon-logistic prior, mean and variance for alpha and beta
      llP <- c(1/lambda, 50, 1, 50)
      
      # Log-normal prior, mean and varaince for mu and sigma
      mu <- -1*log(lambda)-log(2)/2
      lnP <- c(mu, 50, sqrt(log(2)), 50)
      
      cTime <- max(dat)
      
      incProgress(amount= .1, message = "Initialized values")
      
      #Frequentist Predictions
      freqRes <- getFreqInts(dat,nE,MM=200)
      
      incProgress(amount= .65, message = "Frequentist Predictions")
      
      #Bayes predictions
      BayesRes <- getBayesInt(dat,nE,wP,lnP,gP,llP,MM=800)
      mean <- c(freqRes[[1]], BayesRes[[1]])
      lower <- c(freqRes[[2]][,1], BayesRes[[2]][,1])
      upper <- c(freqRes[[2]][,2], BayesRes[[2]][,2])
      methodText <- cbind(c("Freq-Weibull", "Freq-LogNormal", "Freq-Gompertz", "Freq-LogLogistic",
                            "Freq-PredSyn(Avg)", "Freq-PredSyn(MSPE)", "Freq-PredSyn(Vote)",
                            "Bayes-Weibull", "Bayes-LogNormal", "Bayes-Gompertz", "Bayes-LogLogistic",
                            "Bayes-PredSyn(Avg)", "Bayes-PredSyn(MSPE)", "Bayes-PredSyn(Vote)"), rep(" ", 14))
      
      xmin <- floor(min(lower)/50)*50
      xmax <- ceil(max(upper)/50)*50
      
      incProgress(amount= .95, message = "Bayes Predictions")
      
      
      forestplot(methodText, mean, lower, upper, clip = c(xmin, xmax), zero=xmin,
                 xlab=c("Days since first pt on-study"), xticks=seq(xmin, xmax, by=100), boxsize=0.3)
      
      
    })
    
  })
  
  output$data_checks <- renderText({
    data_check_text()
  })
  
  output$data_view <- renderDataTable({
    inputData()
  }, rownames= FALSE)
  
  output$forestPlot <- renderPlot({
    forestPlot()
  })
}

##########################################
### Knit app together 
##########################################
shinyApp(ui, server)