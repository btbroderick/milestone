###########################################
### Setup
###########################################
library(tidyverse)
library(rlang)
library(markdown)
library(shiny)
library(shinythemes)
library(DT)
library(flexsurv)
library(quadprog)
library(Hmisc)
library(VGAM)
library(msm)
library(rmeta)
library(ggplot2)
library(scales)
library(ggstance)
library(plotly)
library(knitr)
library(kableExtra)
library(shinyjs)
theme_set(theme_classic())
source("pgm/utilsBayes1.r")
source("pgm/utilsFreq.r")
source("pgm/utilsWts.r")
source("pgm/helper.R")
###########################################
### User Interface
###########################################

ui <- fluidPage(
  theme = "mayo_theme.css",
  withMathJax(),
  useShinyjs(),
  navbarPage(
    title = "Milestone prediction",
    tabPanel("Main",
             sidebarPanel(
               textInput("study_title", label = "Study Name", value = "Enter Study Name", width = 300),
               numericInput("nE", label = "Milestone (number of events)", value = 165, width = 300),
               numericInput("totAcc", label = "totAcc", value = 701, width = 300),
               dateInput("study_date",label = "First patient enrollment date", value = "2018-01-01", width = 300,format = "yyyy-mm-dd"),
               tags$h6("Date format: yyyy-mm-dd"),
               HTML("<br/>"),
               fileInput("inputfile", NULL, buttonLabel = "Upload", multiple = FALSE, width = 300),
               tags$h6("*File upload format can be found in the 'About' tab"),
               radioButtons(inputId="calculation", label="If you wish to use defult priors for baysian estimates please check on of the follow options to provide historic event rate; otherwise see the 'Custom Prior Distributions' tab.", 
                            choices=c("Cumulative Survial Percentage",
                                      "Median Survival Time",
                                      "Hazard Rate"), selected = "Hazard Rate"),
               conditionalPanel(
                   condition = "input.calculation == 'Cumulative Survial Percentage'",
                   numericInput("survProp", label = "Survival Percentage (0-1)", value = .2, min = 0, max = 1, 100),
                   numericInput("cutoff", label = "Number of Days", value = 30, width = 300)
                 ),
                 conditionalPanel(
                   condition = "input.calculation == 'Median Survival Time'",
                   numericInput("medianTime", label = "Days", value = 1, width = 300)
                 ),
                 conditionalPanel(
                   condition = "input.calculation == 'Hazard Rate'",
                   numericInput("lambda", label = "Lambda", value = 0.0003255076, width = 300)
                 ),
               numericInput("seed", label = "Set Random Seed", value = 7 , width = 300)
               
             ),
             mainPanel(
               tabsetPanel(
                 
                 tabPanel("Data View",
                          dataTableOutput("data_view"),
                          tags$h3(textOutput("data_checks"))
                 ),
                 tabPanel("Customize Prior Distributions",
                          checkboxInput("customdist", "Customize the distributions?", FALSE),
                          div(
                          id = "reset",
                          tags$h3("Weibullprior"),
                          div(style="display:inline-block",numericInput(inputId="meanlambda", label="Mean of Lambda", value = 0.0003255076, width = 150)),
                          div(style="display:inline-block",numericInput(inputId="varlambda", label="Variance of Lambda", value = 10, width = 150)),
                          div(style="display:inline-block",numericInput(inputId="meank", label="Mean of k", value = 1, width = 150)),
                          div(style="display:inline-block",numericInput(inputId="vark", label="Variance of k", value = 10, width = 150)),
                          tags$h3("Gompertz"),
                          div(style="display:inline-block",numericInput(inputId="meaneta", label="Mean of ETA", value = 1, width = 150)),
                          div(style="display:inline-block",numericInput(inputId="vareta", label="Variance of ETA", value = 10, width = 150)),
                          div(style="display:inline-block",numericInput(inputId="meanb", label="Mean of b", value = 0.0002472905, width = 150)),
                          div(style="display:inline-block",numericInput(inputId="varb", label="Variance of b", value = 10, width = 150)),
                          tags$h3("Log-logistic prior"),
                          div(style="display:inline-block",numericInput(inputId="meanalpha", label="Mean of Alpha", value = 3072.125, width = 150)),
                          div(style="display:inline-block",numericInput(inputId="varalpha", label="Variance of Alpha", value = 30721.25, width = 150)),
                          div(style="display:inline-block",numericInput(inputId="meanbeta", label="Mean of Beta", value = 1, width = 150)),
                          div(style="display:inline-block",numericInput(inputId="varbeta", label="Variance of Beta", value = 10, width = 150)),
                          tags$h3("Log-normal"),
                          div(style="display:inline-block",numericInput(inputId="meanmu", label="Mean of Mu", value = 7.683551, width = 150)),
                          div(style="display:inline-block",numericInput(inputId="varmu", label="Variance of Mu", value =  76.83551, width = 150)),
                          div(style="display:inline-block",numericInput(inputId="meansigma", label="Mean of Sigma", value = 0.8325546, width = 150)),
                          div(style="display:inline-block",numericInput(inputId="varsigma", label="Variance of Sigma", value = 10, width = 150))),
                          tags$hr(),
                          actionButton("reset_input", "Reset All Priors to Default")
                          ),
                 tabPanel("Calculate Milestone",
                          actionButton("calculate", label = "Run Milestone Prediction"),
                          plotOutput("forestPlot", width = "100%"),
                          tags$h3("Bayesian"),
                          tableOutput("table2"),
                          tags$h3("Frequentist"),
                          tableOutput("table1"),
                          downloadButton("report", "Generate report")
                 )
               )
             )
    ),
    tabPanel("About",
             includeMarkdown("markdown/About.md")
    ),
    tabPanel("Bayesian Prior Tab",
             includeMarkdown("markdown/Bayesian.md"))
  )
)

##########################################
### Server
##########################################

server <- function(input, output, session) {
  
  lambda <- reactive({
    if (input$calculation == 'Cumulative Survial Percentage'){
      (-log(input$survProp))/input$cutoff
    } else if (input$calculation == 'Median Survival Time'){
      -log(0.5)/input$medianTime
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
      totAcc <- input$totAcc 
      tempdat <- inputData()
      dat <- cbind(tempdat[[1]], tempdat[[2]] * tempdat[[3]])
      #lambda <- 0.0003255076
      
      if(input$customdist)
      {
        lambda <- input$meanlambda
        varlambda <- input$varlambda
        meank <- input$meank
        vark <- input$vark
        meaneta <- input$meaneta
        vareta <- input$vareta
        meanb <- input$meanb
        varb <- input$varb
        meanalpha <- input$meanalpha
        varalpha <- input$varalpha
        meanbeta <- input$meanbeta
        varbeta <- input$varbeta
        meanmu <- input$meanmu
        varmu <- input$varmu
        meansigma <- input$meansigma
        varsigma <- input$varsigma
      } 
      else 
      {
        lambda <- lambda()
        varlambda <- 10*max(isolate(lambda()),1)
        meank <- 1
        vark <- 10
        meaneta <- 1
        vareta <- 10
        meanb <- isolate(lambda())*log(log(2)+1)/log(2)
        varb <- 10*max(meanb,1)
        meanalpha <- 1/isolate(lambda())
        varalpha <- 10*max(meanalpha,1)
        meanbeta <- 1
        varbeta <- 10
        meanmu <- -1*log(isolate(isolate(lambda())))-log(2)/2
        varmu <- 10*max(meanmu,1)
        meansigma <- sqrt(log(2))
        varsigma <- 10
      }
      
      #Priors
      # Weibull prior, mean and varaince for lambda and k
      wP <- c(lambda, varlambda, meank, vark)
      
      # Gompertz prior, mean and variance for eta and b
      #b <- lambda * log(log(2) + 1) / log(2)
      gP <- c(meaneta, vareta, meanb, varb)
      
      # Lon-logistic prior, mean and variance for alpha and beta
      llP <- c(meanalpha, varalpha, meanbeta, varbeta)
      
      # Log-normal prior, mean and varaince for mu and sigma
      #mu <- -1 * log(isolate(isolate(lambda()))) - log(2) / 2
      lnP <- c(meanmu, varmu, meansigma, varsigma)
      
      cTime <- max(dat)
      
      incProgress(amount= .1, message = "Initialized values")
      
      #Frequentist Predictions
      set.seed(input$seed)
      freqRes <- getFreqInts(dat, nE, MM=200, totAcc = totAcc)
      
      incProgress(amount= .65, message = "Frequentist Predictions")
      
      #Bayes predictions
      set.seed(input$seed)
      BayesRes <- getBayesInt(dat, nE, wP, lnP, gP, llP, MM = 800, totAcc = totAcc)
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
                    "Predictive Synthesis (MSPE)", "Predictive Synthesis (Vote)"),
          label = factor(label,levels = c(
            "Gompertz",
            "Log-Logistic",
            "Log-Normal",
            "Weibull",
            "Predictive Synthesis (Vote)",
            "Predictive Synthesis (Average)",
            "Predictive Synthesis (MSPE)"
          )),
          color = case_when(
            type == "Bayesian" & label == "Predictive Synthesis (MSPE)" ~ "A",
            TRUE ~ "B"
          ))

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
    group.colors <- c(A = "#2ca25f", B = "#636363")
    p <- ggplot(predictions(), aes(x = mean, y = label, xmin = lower, xmax = upper, color = color)) +
      geom_pointrangeh() +
      theme(legend.position="none") +
      scale_color_manual(values=group.colors) +
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
        method == "Freq-Weibull" ~ "Weibull"),
        order = c(4,5,7,6,2,1,3,8,10,14,12,4,2,6)) %>%
      arrange(order) %>% 
      select(label, lower, mean, upper, type) %>%
      setNames(c("Method", "Lower bound", "Prediction", "Upper bound","type"))
    
    freq <- filter(all, type == "Frequentist") %>% 
      select(-type)
    
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
        method == "Freq-Weibull" ~ "Weibull"),
        order = c(4,5,7,6,2,1,3,8,10,14,12,4,2,6)) %>%
      arrange(order) %>% 
      select(label2, lower, mean, upper, type) %>%
      setNames(c("Method", "Lower bound", "Prediction", "Upper bound","type"))
    
    bayes <- filter(all, type == "Bayesian") %>% 
      select(-type)
    
    kable(bayes, "html") %>%
      kable_styling(bootstrap_options = c("striped", "hover"))
  }
  
  ##### Reset button ########
  observeEvent(input$reset_input, {
    #reset("reset")
    ## Weibull
    ## chuh
    updateml <- isolate(lambda())
    updatevl <- 10*max(isolate(lambda()),1)
    updateNumericInput(session, "meanlambda", value = updateml)
    updateNumericInput(session, "varlambda", value = updatevl)
    updateNumericInput(session, "meank", value = 1)
    updateNumericInput(session, "vark", value = 10)
    
    ## Gompertz
    updatemb <- isolate(lambda())*log(log(2)+1)/log(2)
    updatevb <- 10*max(updatemb,1)
    updateNumericInput(session, "meaneta", value = 1)
    updateNumericInput(session, "vareta", value = 10)
    updateNumericInput(session, "meanb", value = updatemb)
    updateNumericInput(session, "varb", value = updatevb)
    
    ## Log-logistic
    updatema <- 1/isolate(lambda())
    updateva <- 10*max(updatema,1)
    updateNumericInput(session, "meanalpha", value = updatema)
    updateNumericInput(session, "varalpha", value = updateva)
    updateNumericInput(session, "meanbeta", value = 1)
    updateNumericInput(session, "varbeta", value = 10)
    
    ## Log-logistic
    updatemu <- -1*log(isolate(isolate(lambda())))-log(2)/2
    updatevmu <- 10*max(updatemu,1)
    updatems <- sqrt(log(2))
    updatevs <- 10
    updateNumericInput(session, "meanmu", value = updatemu)
    updateNumericInput(session, "varmu", value = updatevmu)
    updateNumericInput(session, "meansigma", value = updatems)
    updateNumericInput(session, "varsigma", value = updatevs)
  })
}

##########################################
### Knit app together 
##########################################
shinyApp(ui, server)