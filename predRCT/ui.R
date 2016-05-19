#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(survival)
library(devtools)
install_github("gcastella/predRCT")
library(predRCT)
library(MASS)

# Define UI for application that draws a histogram
shinyUI(navbarPage("predRCT",
  tabPanel("Data", 
           sidebarLayout(
             sidebarPanel(
               p("Upload your data:"),
               fileInput('file1', 'Choose file to upload',
                         accept = c(
                           'text/csv',
                           'text/comma-separated-values',
                           'text/tab-separated-values',
                           'text/plain',
                           '.csv',
                           '.tsv'
                         )
               ),
               tags$hr(),
               checkboxInput('header', 'Header', TRUE),
               radioButtons('sep', 'Separator',
                            c(Comma=',',
                              Semicolon=';',
                              Tab='\t',
                              Space=''),
                            ','),
               tags$hr(),
               shiny::actionButton("pujar", label = "Visualize data"),
               tags$hr(),
               shiny::uiOutput("noms_inc"), 
               shiny::uiOutput("noms_delta"), 
               shiny::uiOutput("noms_y"), 
               shiny::uiOutput("noms_ab"),
               tags$hr(),
               shiny::uiOutput("boto_noms")
               ),
             mainPanel(
               dataTableOutput('taula')
             )
          )
  ),
  
  tabPanel("Partition",
    sidebarLayout(
     sidebarPanel(
      shiny::sliderInput(inputId = "particio", label = "Choose a length for the periods:",
                         min = 1, max = 500, value = 30)
     ),
     mainPanel(shiny::plotOutput("particioplot"))
    )
  ),
  tabPanel("Inclusions",
           sidebarLayout(
             sidebarPanel(
               p("Chose your model here:"),
               radioButtons('model', 'Choose the underlaying distribution of the data:',
                            c("Poisson distr"="glm",
                              "Negative binomial"="glm.nb"),
                            "glm"),
               tags$hr(),
               shiny::numericInput('grau', 'Choose the order of the polinomial:', value = 1, min = 1, max = 5),
               tags$hr(),
               radioButtons('link', 'Choose the link function:',
                            c("Logarithmic"="log",
                              "Identity"="identity"),
                            "log"),
               tags$hr(),
               radioButtons('weights', 'Choose the weights of the model:',
                            c("No weights"="rep(1, length(counts$time.periods))",
                              "Linear weights"="counts %>% with(time.periods/sum(time.periods))",
                              "Quadratic weights"="counts %>% with(time.periods^2/sum(time.periods^2))")),
               tags$hr(),
               shiny::numericInput('partition', 'Length of the periods: ', value = 90, min = 1),
               tags$hr(),
               shiny::numericInput('nextp', "Future time periods to predict:", value = 15, min = 1, max = 30),
               tags$hr(),
               shiny::numericInput('totaln', "Total number of patients to achieve:", value = 1264, min = 1, max = 30),
               tags$hr(),
               shiny::actionButton("fit", label = "Fit")
             ),
             mainPanel(
               shiny::plotOutput("inclusionsplot")
             )
          )
  ),
  tabPanel("Events",
           sidebarLayout(
             sidebarPanel(
               radioButtons('survf', 'Chose a survival function:', 
                            c(Weibul = "weibull", 
                              Lognormal = "lognormal", 
                              "Log-logistic" = "loglogistic",
                              Exponential = "exponential", 
                              Gaussian = "gaussian", 
                              Logistic = "logistic",
                              "Prolonged KM" = "km")),
               tags$hr(),
               shiny::actionButton("survb", label = "Plot"),
               tags$hr(),
               shiny::numericInput('pred', "Future time periods to predict:", value = 5, min = 1, max = 30),
               tags$hr(),
               shiny::numericInput('totale', "Total number of events to achieve:", value = 221, min = 1, max = 30),
               tags$hr(),
               shiny::actionButton("eventsb", label = "Calculate")
             ),
             mainPanel(
               shiny::plotOutput("survsplot"),
               shiny::plotOutput("eventsplot")
             )
           )
  )
))
