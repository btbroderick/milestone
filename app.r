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
library(scales)
library(ggstance, lib.loc = here::here("rpkgs"))
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
               textInput("study_title", label = "Study Name", value = "Enter Study Name", width = 300),
               numericInput("nE", label = "Mileston (number of events)", value = 1000, width = 300),
               dateInput("study_date",label = "First patient enrollment date", value = "2018-01-01", width = 300,format = "yyyy-mm-dd"),
               tags$h6("Date format: yyyy-mm-dd"),
               HTML("<br/>"),
               fileInput("inputfile", NULL, buttonLabel = "Upload", multiple = FALSE, width = 300),
               tags$h6("*File upload format can be found in the 'About' tab"),
               radioButtons(inputId="calculation", label="How are you calculating lambda?", 
                            choices=c("option1","option2","option3"), selected = "option3"),
               conditionalPanel(
                   condition = "input.calculation == 'option1'",
                   numericInput("method1param1", label = "Parameter 1", value = 0),
                   numericInput("method2param2", label = "Parameter 2", value = 0)
                 ),
                 conditionalPanel(
                   condition = "input.calculation == 'option2'",
                   numericInput("method2param1", label = "Parameter 1", value = 0),
                   numericInput("method2param2", label = "Parameter 2", value = 0)
                 ),
                 conditionalPanel(
                   condition = "input.calculation == 'option3'",
                   numericInput("lambda", label = "Lambda", value = 0.0003255076)
                 )
               
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
                          tableOutput("table1"),
                          tableOutput("table2"),
                          downloadButton("report", "Generate report")
                 )
               )
             )
    ),
    tabPanel("About",
             includeMarkdown(here::here("markdown", "About.md"))
    ),
    tabPanel("Bayesian Prior Tab",
             includeMarkdown(here::here("markdown", "Bayesian.md")))
  )
)

##########################################
### Server
##########################################

server <- function(input, output, session) {
  
  lambda <- reactive({
    if (input$calculation == 'option1'){
      0.0003255076
    } else if (input$calculation == 'option2'){
      0.0003255076
    } else input$lambda
  })
  
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
      #lambda <- 0.0003255076
      
      #Priors
      # Weibull prior, mean and varaince for lambda and k
      wP <- c(lambda(), 50, 1, 50)
      
      # Gompertz prior, mean and variance for eta and b
      b <- lambda() * log(log(2) + 1) / log(2)
      gP <- c(1, 50, b, 50)
      
      # Lon-logistic prior, mean and variance for alpha and beta
      llP <- c(1 / lambda(), 50, 1, 50)
      
      # Log-normal prior, mean and varaince for mu and sigma
      mu <- -1 * log(lambda()) - log(2) / 2
      lnP <- c(mu, 50, sqrt(log(2)), 50)
      
      cTime <- max(dat)
      
      incProgress(amount= .1, message = "Initialized values")
      
      #Frequentist Predictions
      set.seed(7)
      freqRes <- getFreqInts(dat, nE, MM=200)
      
      incProgress(amount= .65, message = "Frequentist Predictions")
      
      #Bayes predictions
      set.seed(7)
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
                              upper = as.Date(upper, origin = input$study_date)) %>% 
        mutate(type = case_when(
          str_detect(method, pattern = "Freq") ~ "Frequentist",
          str_detect(method, pattern = "Bayes") ~ "Bayesian"),
          label = c("Weibull", "Log-Normal", "Gompertz", "Log-Logistic","Predictive Synthesis (Average)",
                    "Predictive Synthesis (MSPE)", "Predictive Synthesis (Vote)","Weibull", "Log-Normal", 
                    "Gompertz", "Log-Logistic","Predictive Synthesis (Average)",
                    "Predictive Synthesis (MSPE)", "Predictive Synthesis (Vote)"))
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
      params <- list(data = predictions(), 
                     study = input$study_title,
                     first_date = input$study_date,
                     number_events = sum(inputData()[[3]]),
                     milestone = input$nE)
      
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
    p <- ggplot(predictions(), aes(x = mean, y = label, xmin = lower, xmax = upper)) +
      geom_pointrangeh() +
      facet_grid(type ~ ., scale = "free", switch="both") + 
      scale_x_date(labels = date_format("%Y-%m-%d")) +
      labs(y = "Days since first patient enrolled", x = "")
    p
  })
  
  output$table1 <- function(){
    all <- predictions() %>%
      mutate(label2 = case_when(
        method == "Bayes-Gompertz" ~ "Gompertz",
        method == "Bayes-LogLogistic" ~ "Log-Logistic",
        method == "Bayes-LogNormal" ~ "Log-Normal",
        method == "Bayes-PredSyn(Avg)" ~ "Predicitive Synthesis (Average)",
        method == "Bayes-PredSyn(MSPE)" ~ "Predicitive Synthesis (MSPE)",
        method == "Bayes-PredSyn(Vote)" ~ "Predicitive Synthesis (Vote)",
        method == "Bayes-Weibull" ~ "Weibull",
        method == "Freq-Gompertz" ~ "Gompertz",
        method == "Freq-LogLogistic" ~ "Log-Logistic",
        method == "Freq-LogNormal" ~ "Log-Normal",
        method == "Freq-PredSyn(Avg)" ~ "Predicitive Synthesis (Average)",
        method == "Freq-PredSyn(MSPE)" ~ "Predicitive Synthesis (MSPE)",
        method == "Freq-PredSyn(Vote)" ~ "Predicitive Synthesis (Vote)",
        method == "Freq-Weibull" ~ "Weibull")) %>%
      select(label2, lower, mean, upper) %>%
      setNames(c("Method", "Lower bound", "Prediction", "Upper bound"))
    
    freq <- filter(all, row_number() <=7)
    kable(freq, "html") %>%
      kable_styling(bootstrap_options = c("striped", "hover"))
  }
  
  output$table2 <- function(){
    all <- predictions() %>%
      mutate(label2 = case_when(
        method == "Bayes-Gompertz" ~ "Gompertz",
        method == "Bayes-LogLogistic" ~ "Log-Logistic",
        method == "Bayes-LogNormal" ~ "Log-Normal",
        method == "Bayes-PredSyn(Avg)" ~ "Predicitive Synthesis (Average)",
        method == "Bayes-PredSyn(MSPE)" ~ "Predicitive Synthesis (MSPE)",
        method == "Bayes-PredSyn(Vote)" ~ "Predicitive Synthesis (Vote)",
        method == "Bayes-Weibull" ~ "Weibull",
        method == "Freq-Gompertz" ~ "Gompertz",
        method == "Freq-LogLogistic" ~ "Log-Logistic",
        method == "Freq-LogNormal" ~ "Log-Normal",
        method == "Freq-PredSyn(Avg)" ~ "Predicitive Synthesis (Average)",
        method == "Freq-PredSyn(MSPE)" ~ "Predicitive Synthesis (MSPE)",
        method == "Freq-PredSyn(Vote)" ~ "Predicitive Synthesis (Vote)",
        method == "Freq-Weibull" ~ "Weibull")) %>%
      select(label2, lower, mean, upper) %>%
      setNames(c("Method", "Lower bound", "Prediction", "Upper bound"))
    
    bayes <- filter(all, row_number() > 7)
    
    kable(bayes, "html") %>%
      kable_styling(bootstrap_options = c("striped", "hover"))
  }
  
}

##########################################
### Knit app together 
##########################################
shinyApp(ui, server)