###########################################
### Setup
###########################################
library(tidyverse)
library(markdown)
library(shiny)
library(shinythemes)
library(DT)
library(flexsurv)
library(quadprog)
library(Hmisc)
library(msm) 
library(VGAM)
library(rmeta)
library(ggplot2)
library(plotly)
theme_set(theme_classic())
source(here::here("pgm","utilsBayes1.r"))
source(here::here("pgm","utilsFreq.r"))
source(here::here("pgm","utilsWts.r"))
source(here::here("pgm", "helper.R"))
###########################################
### User Interface
###########################################

ui <- fluidPage(
  theme = "mayo_theme.css",
  navbarPage(
    title = "Milestone prediction",
    tabPanel("Main",
             sidebarPanel(
               numericInput("nE", label = "Landmark Event Number", value = 1000, width = 300),
               dateInput("study_date",label = "Study start date", value = "2018-01-01", width = 300,format = "yyyy-dd-mm"),
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
                          plotOutput("forestPlot", width = "100%"),
                          downloadButton("report", "Generate report")
                 )
               )
             )
    ),
    tabPanel("About",
             includeMarkdown(here::here("About.md"))
    )
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
  
  predictions <- eventReactive(input$calculate,{
    withProgress(value = 0, message =  "Calculating", {
      nE <- input$nE # landmark event number
      tempdat <- inputData()
      dat <- cbind(tempdat[[1]], tempdat[[2]] * tempdat[[3]])
      lambda <- 0.0003255076
      
      #Priors
      # Weibull prior, mean and varaince for lambda and k
      wP <- c(lambda, 50, 1, 50)
      
      # Gompertz prior, mean and variance for eta and b
      b <- lambda * log(log(2) + 1) / log(2)
      gP <- c(1, 50, b, 50)
      
      # Lon-logistic prior, mean and variance for alpha and beta
      llP <- c(1 / lambda, 50, 1, 50)
      
      # Log-normal prior, mean and varaince for mu and sigma
      mu <- -1 * log(lambda) - log(2) / 2
      lnP <- c(mu, 50, sqrt(log(2)), 50)
      
      cTime <- max(dat)
      
      incProgress(amount= .1, message = "Initialized values")
      
      #Frequentist Predictions
      freqRes <- getFreqInts(dat, nE, MM=200)
      
      incProgress(amount= .65, message = "Frequentist Predictions")
      
      #Bayes predictions
      BayesRes <- getBayesInt(dat, nE, wP, lnP, gP, llP, MM = 800)
      mean <- c(freqRes[[1]], BayesRes[[1]])
      lower <- c(freqRes[[2]][,1], BayesRes[[2]][,1])
      upper <- c(freqRes[[2]][,2], BayesRes[[2]][,2])
      methodText <- c("Freq-Weibull", "Freq-LogNormal", "Freq-Gompertz", "Freq-LogLogistic",
                            "Freq-PredSyn(Avg)", "Freq-PredSyn(MSPE)", "Freq-PredSyn(Vote)",
                            "Bayes-Weibull", "Bayes-LogNormal", "Bayes-Gompertz", "Bayes-LogLogistic",
                            "Bayes-PredSyn(Avg)", "Bayes-PredSyn(MSPE)", "Bayes-PredSyn(Vote)")
      
      xmin <- floor(min(lower) / 50) * 50
      xmax <- ceil(max(upper) / 50) * 50
      
      incProgress(amount = .95, message = "Bayes Predictions")
      
      plotdata <- data.frame(method = methodText,
                              mean = as.Date(mean, origin = input$study_date) , 
                              lower = as.Date(lower, origin = input$study_date), 
                              upper = as.Date(upper, origin = input$study_date))
      
      plotdata
      
    })
  })
  
  
  output$report <- downloadHandler(
    filename = "report.html",
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "report.Rmd")
      file.copy("report.Rmd", tempReport, overwrite = TRUE)
      
      # Set up parameters to pass to Rmd document
      params <- list(data = predictions())
      
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  )
  
  
  output$data_checks <- renderText({
    data_check_text()
  })
  
  output$data_view <- renderDataTable({
    inputData()
  }, rownames= FALSE)
  
  output$forestPlot <- renderPlot({
    p <- ggplot(predictions(), aes(x = method, y = mean, ymin = lower, ymax = upper)) +
      geom_pointrange() +
      coord_flip() + 
      scale_y_date(labels = scales::date_format("%d/%m/%Y")) +
      labs(y = "Predicted milestone date", x = "") 
    p
  })
}

##########################################
### Knit app together 
##########################################
shinyApp(ui, server)