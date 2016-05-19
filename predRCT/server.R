#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(survival)
library(devtools)
install_github("gcastella/predRCT", force = TRUE)
library(predRCT)
library(MASS)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  aux <- eventReactive(input$pujar, {
    read.table(file = input$file1$datapath, 
               header = input$header, 
               sep = input$sep)
  })
  
  
  output$noms_inc <- renderUI({
    if (is.null(aux())){
      return(NULL)
    } else {
      shiny::selectInput("inc", "Variable with the time of inclusions:", colnames(aux()))
    }
  })
  output$noms_y <- renderUI({
    if (is.null(aux())){
      return(NULL)
    } else {
      shiny::selectInput("y", "Variable with the time of follow-up:", colnames(aux()))
    }
  })
  output$noms_ab <- renderUI({
    if (is.null(aux())){
      return(NULL)
    } else {
      shiny::selectInput("ab", "Abandon indicator:", colnames(aux()))
    }
  })
  output$noms_delta <- renderUI({
    if (is.null(aux())){
      return(NULL)
    } else {
      shiny::selectInput("delta", "Censoring indicator:", colnames(aux()))
    }
  })
  output$boto_noms <- renderUI({
    if (is.null(aux())){
      return(NULL)
    } else {
      shiny::actionButton("botoopcions", "Read variables")
    }
  })
  output$noms_delta <- renderUI({
    if (is.null(aux())){
      return(NULL)
    } else {
      shiny::selectInput("delta", "Censoring indicator:", colnames(aux()))
    }
  })

  
  dades <- eventReactive(input$botoopcions, {
    auxvect <- c(inclusion = input$inc,
                 y = input$y, 
                 delta = input$delta,
                 abandon = input$ab)
    dd <- aux()[auxvect]
    names(dd) <- names(auxvect)
    return(dd)
  })
  
  output$particioplot <- renderPlot(expr = periods_plot(period.length = input$particio, 
                                                        time.inclusion = dades()$inclusion, 
                                                        use.letters = FALSE, 
                                                        call.mfrow = TRUE))

  aux_fit <- eventReactive(input$fit, {
    isaacc <- dades()
    isaacc <- period_adj(data = isaacc, period.length = input$partition)
    last.period <- as.period(max(isaacc$y + isaacc$inclusion, na.rm = TRUE), period.length = input$partition)
    time.periods <- isaacc$inclusion %>%
      as.period(period.length = input$partition) %>%
      factor(levels = 1:last.period)

    counts <- time.periods %>%
      table %>%
      cumsum %>%
      data.frame %>%
      add_rownames(var = "time.periods")

    counts <- transmute(counts,
                        time.periods = as.numeric(time.periods),
                        counts = diff(c(0, .)))
    formu <- as.formula(paste0("counts ~ poly(time.periods, degree = ", input$grau, ", raw = TRUE)"))
    aux <- if(input$model == "glm") eval(parse(text = paste0("poisson(link = ", input$link, ")"))) else input$link
    out <- list(input$model,
                formula = formu,
                weights = eval(parse(text = input$weights)),
                family = aux
    )
    if(input$model != "glm")
    names(out)[4] <- "link"
    return(inclusionsCount(time.inclusion = isaacc$inclusion, 
                           model.args = out, 
                           predict.next = input$nextp, 
                           boot = TRUE, 
                           period.length = input$partition, 
                           boot.samples = 500))
  })
  
  aux_surv <- eventReactive(input$survb, {
    isaacc <- dades()
    isaacc <- period_adj(data = isaacc, period.length = input$partition)
    last.period <- as.period(max(isaacc$y + isaacc$inclusion, na.rm = TRUE), period.length = input$partition)
    
    isaacc.surv <- Surv(subset(isaacc$y, isaacc$y > 0), 
                        subset(isaacc$delta, isaacc$y > 0))
    isaacc.km <- survfit(isaacc.surv ~ 1)
    if(input$survf != "km"){
      out <- survreg(isaacc.surv ~ 1, dist = input$survf)
    } else {
      out <- isaacc.km
    }
    res <- list(survFUN(how.many = 1000, object = out, period.length = input$partition))
    names(res) <- input$survf
    return(list(b = isaacc.km, a = res))
  })
  aux_ev <- eventReactive(input$eventsb, {
    isaacc <- dades()
    isaacc <- period_adj(data = isaacc, period.length = input$partition)
    last.period <- as.period(max(isaacc$y + isaacc$inclusion, na.rm = TRUE), period.length = input$partition)
    
    return(eventsCount(data = isaacc, new.inclusions = aux_fit()$counts[1:input$pred], period.length = input$partition, past = TRUE, FUN = aux_surv()$a[[1]]))
  })
  
  output$survsplot <- renderPlot(survs_plot(surv.funs = aux_surv()$a, survfit.object = aux_surv()$b, period.length = input$partition, call.mfrow = FALSE, use.letters = FALSE, mar = c(3, 3, 3, 2)))
  output$eventsplot <- renderPlot(expr = events_plot(list(aux_ev()), E = input$totale, period.length = input$partition, use.letters = FALSE))
  output$inclusionsplot <- renderPlot(expr = inclusions_plot(list(aux_fit()), period.length = input$partition, N = input$totaln, use.letters = FALSE, call.mfrow = TRUE))
  output$taula <- renderDataTable(aux())
})
