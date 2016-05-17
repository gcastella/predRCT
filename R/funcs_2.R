# OPTIONS -----------------------------------------------------------------

# .onLoad <- function(libname = find.package("predRCT"), pkgname = "predRCT") {
#   # CARE: side effects
#   options(RCT_names_match = c(inclusion = "inclusion",
#                               y = "y",
#                               delta = "delta",
#                               abandon = "abandon"))
#   
# }

options(RCT_names_match = c(inclusion = "inclusion",
                            y = "y",
                            delta = "delta",
                            abandon = "abandon"))


# GENERAL FUNCTIONS -------------------------------------------------------

#' Time Period
#'
#' Transform time unit to period of time
#'
#' @param times Vector of calendar times (elapsed from the begining of the study) to transform to periods. Dates are also allowed but origin argument must be given then.
#' @param period.length Length of each time period.
#' @param numeric Logical. If \code{TRUE}, returns the periods as a numeric. Defaults to \code{TRUE}
#' @param origin Date at which the trial began. Ignored unless times are dates. 
#'
#' @details
#' Time periods are closed on the right, and 0 is included in the first one.
#' Period labels, when \code{numeric} is \code{FALSE}, come from \code{cut}.
#'
#' @return Returns a character (or numeric when \code{numeric = TRUE}) vector of length values with
#' the corresponding time periods for each time.
#' 
#' @examples  
#' as.period(c(0, 10, 89, 90, 91, 200), period.length = 90, numeric = TRUE)
#' as.period(c(0, 10, 89, 90, 91, 200), period.length = 90, numeric = FALSE) 
#' as.period(as.Date("10-2-2000", "%d-%m-%Y"), period.length = 90, origin = as.Date("01-1-2000", "%d-%m-%Y")) 
#' 
#' @export
as.period <- function(times, period.length, numeric = TRUE, origin = NULL){
  stopifnot(all(times >= 0))
  if(!is.numeric(times)){
    if(is.null(origin)) stop("If the vector 'times' has dates, origin must be given.")
    times <- as.numeric(times - origin)
  }
  until <- max(times, na.rm = T) %/% period.length + 2
  out <- cut(times, 
             breaks = seq(0, 
                          by = period.length, 
                          length.out = until), 
             include.lowest = TRUE)
  out <- if(numeric) as.numeric(out) else as.character(out)
  return(out)
}


#' Removing last time period
#'
#' Removes last time period from the data, and changes those patients that had event in this period. Useful in case that
#' the last time period has not finished.
#'
#' @param data Data to remove last time period. A data frame with variables inclusion, y, abandon and delta at least. 
#' Other variables won't be modified.
#' @param period.length Length of each period.
#' @param last.period Last period to take into account for the adjustment.
#' @param origin Passed to \code{\link{as.period}}.
#' @param numeric Passed to \code{\link{as.period}}.
#' @param varnames By default \code{getOption("RCT_names_match")}. This is a character vector with names "inclusion", 
#' "delta" and "y" (in any order), and its corresponding name in the argument data (gives the matching for accessing
#'  to inclusion times, censore indicator, and time of follow-up).
#' @details 
#' Usually the calendar time partition does not match exactly with the current follow-up time of the RCT. In these cases, the last
#' time period is shorter than the others and one solution is to ignore the data in this time period, removing those patients and
#' correcting censores and events from others.
#' 
#' @return The same data frame but without the patients included in the last period of time (removed), and some 
#' modified delta's.
#' 
#' @export
#' 
period_adj <- function(data, period.length, 
                       last.period = NULL, 
                       varnames = getOption("RCT_names_match"), 
                       origin = NULL,
                       numeric = TRUE){
  
  if(is.null(last.period)) last.period <- as.period(max(data[, varnames["y"]] + data[, varnames["inclusion"]], na.rm = TRUE), 
                                                    period.length = period.length,
                                                    numeric = numeric,
                                                    origin = origin)

  start_lastp <- which(as.period(times = data[, varnames["inclusion"]], 
                                 period.length = period.length, 
                                 origin = origin,
                                 numeric = numeric) >= last.period)
  if(length(start_lastp) > 0)  data <- data[- start_lastp, ]

  end_lastp <- which(as.period(data[, varnames["y"]] + data[, varnames["inclusion"]], 
                               period.length = period.length, 
                               origin = origin,
                               numeric = numeric) >= last.period)
  if(length(end_lastp) > 0){
    data[end_lastp, varnames["y"]] <- (last.period - 1) * period.length - data[end_lastp, varnames["inclusion"]]
    data[end_lastp, varnames["delta"]] <- 0
  }

  return(data)
}


# COUNTING MODELS ---------------------------------------------------------

# limit <- function(periods, trim = 0.001, maxobs = max(periods), ...){
#   alpha1 <- -2 * log(trim / (2 - trim)) / maxobs
#   alpha0 <- -0.5 * maxobs * alpha1
#   num <- exp(periods * alpha1 + alpha0)
#   return(num / (1 + num))
# }
# 
# pred_asymptote <- function(periods, model, asymptote = min(model$y, na.rm = T), ...){
#   model <- eval(model, parent.frame())
#   if(missing(periods)){
#     predict(model, type = "resp") + asymptote * limit(periods, ...)
#   } else {
#     df <- data.frame(periods)
#     colnames(df) <- colnames(get_all_vars(model, data = model$model))[2]
#     predict(model, newdata = df, type = "resp") + asymptote * limit(periods, ...)
#   }
# }



#' Predicting the number of inclusions in a time period
#' 
#' Predict the number of total inclusions in future time periods using the passed model for the counts.
#' 
#' @param time.inclusion Inclusion times for the patients in the study.
#' @param period.inclusion Period of inclusion of the patients in the study. Either this or time.inlcusion has to be passed.
#' @param model.args List with the arguments needed for the model, called via do.call(model.args[[1]], model.args[-1]).
#' @param do.before, do.after Expressions to be executed before and after fitting the model.
#' @param predict.next Number of periods of which to predict the inclusions.
#' @param boot Calculate or not bootstrap CI intervals.
#' @param boot.samples Number of samples for the CI intervals.
#' @param period.length, origin, numeric Arguments passed to as.period if time.inclusion is used instead of period.inclusion.
#' 
#' @details 
#' Data can be passed using time.inclusion with the days elapsed from the beginning of the RCT until the inclusion of each patient, or period.inclusion
#' with the period of the study that the patient was recruited. In the first case peroid.length, origin and/or numeric can be passed also.
#' For fitting a poisson glm for instance, set model.args as list("glm", formula, family = poisson). Then do.before and do.after can be used to modify
#' or change the data as needed. substitute(df2$counts) can be used to get the data.
#' @return  A list with the counts predicted in each period using the model given, the cummulative counts, bootstrap CI for the counts and cummulative counts
#' if asked, the data used and the last period.
#' 
#' @export
inclusionsCount <- function(time.inclusion = NULL,
                            period.inclusion = NULL,
                            model.args = list("glm",
                                              formula = counts ~ time.periods + I(time.periods^2),
                                              family = poisson),
                            do.before = NULL,
                            do.after = NULL,
                            predict.next = 10,
                            boot = TRUE,
                            boot.samples = 500,
                            period.length,
                            origin = NULL,
                            numeric = TRUE){

  stopifnot(!is.null(time.inclusion) | !is.null(period.inclusion))
  if(!is.null(time.inclusion)) period.inclusion <- as.period(time.inclusion, period.length = period.length, numeric = TRUE)
  dat <- data.frame(inc_period = period.inclusion)

  last.period <- max(period.inclusion, na.rm = TRUE)

  bootstrap_fun <- function(dat, indices){
    df <- dat[indices, , drop = FALSE]
    df <- data.frame(table(factor(df[[1]], levels = 1:last.period)))
    names(df) <- c("time.periods", "counts")
    df$time.periods <- as.numeric(df$time.periods)
    df2 <- df
    if(!is.null(do.before)) df$counts <- do.call(do.before[[1]], c(list(df$counts), do.before[-1]))
    mod <- do.call(model.args[[1]], c(model.args[-1], list(data = df)))
    newdata <- data.frame(time.periods = min(df$time.periods, na.rm = TRUE):(last.period + predict.next))
    res <- predict(mod, newdata = newdata, type = "response")
    if(!is.null(do.after)) res <- do.call(do.after[[1]], c(list(res), do.after[-1]))
    return(res)
  }

  counts_out <- bootstrap_fun(dat = dat, indices = TRUE)

  bts_ci <- NULL
  cum_bts_ci <- NULL
  if(boot){
    bts <- replicate(boot.samples, expr = bootstrap_fun(dat = dat, sample(nrow(dat), replace = TRUE)))
    bts_ci <- apply(bts, 1, quantile, probs = c(0.025, 0.975))
    cum_bts <- apply(bts[-(1:last.period), ], 2, cumsum) + nrow(dat)
    cum_bts_ci <- apply(cum_bts, 1, quantile, probs = c(0.025, 0.975))
  }

  out <- list()
  out$counts <- counts_out
  out$bootstrap.ci <- bts_ci
  out$cummulative.counts <- cumsum(out$counts[-(1:last.period)]) + nrow(dat)
  out$bootstrap.cummulative.ci <- cum_bts_ci
  out$last.period <- last.period
  out$data <- data.frame(table(factor(period.inclusion, levels = 1:last.period)))
  names(out$data) <- c("time.periods", "counts")
  out$data$time.periods <- as.numeric(out$data$time.periods)
  return(invisible(out))
}


# SURVIVAL FUNCTIONS ------------------------------------------------------

#' Generate Survival Function
#'
#' Taking survfit or survreg output to create the corresponding survival function.
#'
#' @param object survfit or survreg object.
#' @param period.length Length of each time period.
#' @param how.many How many quantiles you want for survreg case.
#'
#' @return A function with one parameter: the vector of times to calculate the survival probabilities.
#'
#' @seealso \code{\link{survs_plot}}, \code{\link{eventsCount}}.
#' @import survival
#' @export
survFUN <- function(object, period.length, how.many = 1000){

  if("survfit" %in% class(object)){
    last <- length(summary(object)$n.risk)
    svv <- c(1, summary(object)$surv)
    t <- c(0, summary(object)$time)
    svv_last <- summary(object)$surv[last]
    lambda <- summary(object)$n.event[last] / summary(object)$n.risk[last]
    t_last <- summary(object)$time[last]

    return(function(times){
      quins <- findInterval(times, t)
      quins[quins == 0] <- 1
      ifelse(times <= 0,
             1,
             ifelse(times >= t_last,
                    svv_last * (1 - lambda) ^ ((times - t_last) / period.length),
                    svv[quins]
             )
      )
    })
  } else if("survreg" %in% class(object)){
    preds <- data.frame(valor = predict(object,
                                        type = "quantile",
                                        p = seq(0, 1, length = how.many))[1, ],
                        surv = 1 - seq(0, 1, length = how.many))
    return(function(times){
      quins <- findInterval(x = times, vec = preds$valor)
      quins[quins == 0] <- 1
      preds[quins, "surv"]
    })
  } else {
    stop("object must have 'survfit' or 'survreg' class.")
  }
}

# TIME-TO-EVENT FUNCTIONS -------------------------------------------------


#' Calculate proabilities of event
#' 
#' Calculate proabilities of event for a set/subset of patients in a time interval, using a given survival function.
#' 
#' @param patients Subset of indices for which to calculate the probabilities
#' @param data Dataset with the times of inclussion, withdrawals, censoring indicator, time of follow-up, ...
#' @param start, end Start or end of the interval of which to calculate the probabilities.
#' @param last.time Number of days elapsed since the beginning of the trial.
#' @param FUN survival function used for calculating the probabilities. Usually output from survFUN().
#' 
#' @return A vector with the probabilities with the same length as nrow(data[patients, , drop = FALSE]).
#' 
#' @export
#' 
probAB <- function(patients = TRUE, data, start, end, last.time, FUN){

  datapat <- data[patients, ]

  ifelse(datapat$delta == 1 &
           datapat$y + datapat$inclusion <= end &
           datapat$y + datapat$inclusion > start,
         1,
         ifelse(datapat$inclusion > last.time,
                FUN(start - datapat$inclusion) - FUN(end - datapat$inclusion),
                ifelse(datapat$delta == 0 & datapat$inclusion <= last.time &
                         datapat$abandon == 0 &
                         datapat$y + datapat$inclusion <= end &
                         datapat$y + datapat$inclusion > start,
                       1 - FUN(end - datapat$inclusion) / FUN(datapat$y),
                       ifelse(datapat$delta == 0 &
                                datapat$inclusion <= last.time &
                                datapat$abandon == 0 &
                                datapat$y + datapat$inclusion <= start,
                              (FUN(start - datapat$inclusion) - FUN(end - datapat$inclusion)) / FUN(datapat$y),
                              0
                       )
                )
         )
  )
}


#' Total number of expected events
#'
#' Calculates the total number of expected events (past and/or future) in a RCT in each time period.
#'
#' @param data Defaults to NULL. Data frame with the inclusion times (named as inclusion), follow-up times (y) or ending times (end),
#' withdrawal indicator (abandon) and censoring indicator (delta). You can provide time.inclusion, delta,
#' abandon and time.end or time.followup instead.
#' @param time.inclusion Ignored if data is not NULL. Vector with times of inclusion (numeric or date).
#' @param time.end Ignored if data is not NULL. Vector with the times of end of follow-up (numeric or date). One
#' of time.followup or time.end must be supplied (if data is not null)
#' @param time.followup Ignored if data is not NULL. Numeric vector with the total time of follow-up. One
#' of time.followup or time.end must be supplied (if data is not null)
#' @param delta Ignored if data is not NULL. Censoring indicator: 1 if the event is observed, 0 otherwise
#' @param abandon Ignored if data is not NULL. Withdrawal indicator: 1 if the patient is a withdrawal, abandons,
#' or dies for other causes (etc). 0 otherwise.
#' @param new.inclusions Numeric vector with the new patients included in each future time period.
#' @param origin Date that marks the first included patient in the RCT.
#' @param period.length Length of each time period.
#' @param past Logical. Calculate the expected number of events in each past time periods in the data provided?
#' @param FUN Survival function for calculating probabilities of event in each period.
#'
#' @return
#' A list with
#' \item{events}{Numeric with the cummulative number of events in each period (past and future).}
#' \item{percent.events}{Numeric with the percentage of events in the sample (sample size is modified in every period).}
#' \item{new.inclusions}{Argument new.inclusions.}
#' \item{last.time}{Highest time in the data: max(inclusion + followup).}
#' \item{data}{Data frame constructed from the arguments provided (or just the data argument with some modifications, if provided).}
#' \item{survFun}{Argument FUN.}
#' \item{origin}{Argument origin or the origin assumed for the data (first inclusion, when the RCT started).}
#'
#' @seealso \code{\link{survs_plot}}, \code{\link{events_plot}}, \code{\link{inclusionsCount}}
#'
#' @export
eventsCount <- function(data = NULL,
                        time.inclusion = NULL,
                        time.end = NULL, time.followup = NULL,
                        delta, abandon,
                        new.inclusions = NULL,
                        origin = NULL,
                        period.length,
                        past = FALSE,
                        boot = TRUE,
                        boot.samples = 500,
                        varnames = getOption("RCT_names_match"), 
                        surv.object,
                        how.many = 1000,
                        FUN = NULL){

  stopifnot(!is.null(time.end) | !is.null(time.followup) | !is.null(data))
  if(is.null(data)){
    if(is.null(time.followup)) time.followup <- time.end - time.inclusion
    data <- data.frame(inclusion = time.inclusion,
                       abandon = abandon,
                       delta = delta,
                       y = time.followup)

  } else {
    if(is.null(data[, varnames["y"]])){
      data[, varnames["y"]] <- as.numeric(data[, varnames["end"]] - data[, varnames["inclusion"]])
    }
    data <- data[, varnames]
    names(data) <- names(varnames)
  }
  data2 <- data
  if("Date" %in% class(data[, varnames["inclusion"]])){
    if(is.null(origin)) origin <- min(data[, varnames["inclusion"]], na.rm = T)
    data[, varnames["inclusion"]] <- as.numeric(data[, varnames["inclusion"]] - origin) + 1
  } else if(is.null(origin)){
    origin <- 0
    warning("0 was set as the origin.")
  }

  bootstrap_fun <- function(dat, indices){
    data <- dat[indices, , drop = FALSE]
    data2 <- data
    if(is.null(FUN)){
      FUN <- survFUN(object = update(surv.object, data = data), 
                     period.length = period.length, 
                     how.many = how.many)
    }
    
    stopifnot(ncol(data) == 4)
    last.period <- as.period(max(data[, varnames["inclusion"]] + data[, varnames["y"]], na.rm = TRUE), 
                             period.length = period.length)
    last.time <- last.period * period.length

    if(is.null(new.inclusions) | past){
      events_past <- NULL
      lep <- length(events_past)
      while(lep < last.period){
        events_past <- append(x = events_past,
                              values = sum(probAB(data = data,
                                                  start = 0,
                                                  end = period.length * (lep + 1),
                                                  last.time = last.time,
                                                  FUN = FUN)))
        lep <- length(events_past)
      }
      names(events_past) <- as.period(seq(along.with = events_past) * period.length, period.length = period.length, numeric = FALSE)
      n_past <- cumsum(table(as.period(data[, varnames["inclusion"]], period.length = period.length)))[seq_along(events_past)]
    } else {
      events_past <- sum(probAB(data = data,
                                start = 0,
                                end = last.time,
                                last.time = last.time,
                                FUN = FUN))
      names(events_past) <- paste0("[0,", last.time, "]")
      n_past <- nrow(data2)
    }

    if(!is.null(new.inclusions)){
      events_future <- NULL
      lep <- length(events_future)
      max.period <- length(new.inclusions)
      while(lep < max.period){
        if(new.inclusions[lep + 1] > 0){
          newdata <- data.frame(inclusion = runif(new.inclusions[lep + 1],
                                                  min = last.time + period.length * lep,
                                                  max = last.time + period.length * (lep + 1)),
                                y = 0,
                                abandon = 0,
                                delta = 0)
          data <- rbind(data, newdata[names(data)])
        }
        events_future <- append(x = events_future,
                                values = sum(probAB(data = data,
                                                    start = 0,
                                                    end =  last.time + period.length * (lep + 1),
                                                    last.time = last.time,
                                                    FUN = FUN)))
        lep <- length(events_future)
      }
      names(events_future) <- as.period(last.time + seq(along.with = events_future) * period.length, period.length = period.length, numeric = FALSE)
      n_future <- nrow(data2) + cumsum(new.inclusions)[seq_along(events_future)]
    } else {
      events_future <- NULL
      n_future <- NULL
    }
    
    ev_obs <- cumsum(
      table(as.period(times = data2$inclusion + data2$y,
                      period.length = period.length),
            data2$delta)[, 2]
    )
    
    out <- list()
    out$events <- c(events_past, events_future)
    out$events_observed <- ev_obs
    out$new.inclusions <- new.inclusions
    out$last.time <- last.time
    out$data <- data2
    out$survFUN <- FUN
    out$origin <- origin
    return(out)
  }
  out <- bootstrap_fun(data, TRUE)

  if(boot){
    bts <- replicate(simplify = FALSE, 
                     boot.samples, 
                     expr = bootstrap_fun(dat = data, 
                                          sample(nrow(data), replace = TRUE)))
    ev_mat <- do.call(rbind, lapply(bts, getElement, name = "events"))
    out$events_ci <- apply(ev_mat, 2, quantile, probs = c(0.025, 0.975))
  }
  
  return(invisible(out))
}


# PLOTS -------------------------------------------------------------------

#' Plot the output of inclusionsCount
#' 
#' Produces two plots: one with the counts of each time period and the predicted ones, and the other with the cummulative counts.
#'
#' @param inclusionsCount.list A list of outputs of inclusionsCount().
#' @param origin Date of the beginning of the trial.
#' @param period.length Length of the time periods.
#' @param N Total number of included patients to achieve.
#' @param use.letters Index all elements in the list with a letter when plotting.
#' @param call.mfrow Wether to display all plots in the same window.
#' @param ... Other arguments passed to par().
#' 
#' @return Just produces one or more plots.
#'  
#' @export
inclusions_plot <- function(inclusionsCount.list, 
                            origin = NULL, period.length = NULL, N, 
                            use.letters = TRUE, call.mfrow = TRUE, ...){
  par_opts <- par()
  on.exit(suppressWarnings(par(par_opts)))
  if(all(c("counts", "cummulative.counts") %in% names(inclusionsCount.list))) inclusionsCount.list <- list(inclusionsCount.list)
  if(! all(sapply(inclusionsCount.list, is.list))) stop("Argument 'inclusionsCount.list' must be a list of different outputs from 'inclusionsCount' function.")
 
  if(call.mfrow) par(mfrow = c(length(inclusionsCount.list), 2))
  par(...)
  
  for(i in seq_along(inclusionsCount.list)){
    aux <- inclusionsCount.list[[i]]
    bool <- 0
    if(!is.null(origin) & !is.null(period.length)){
      bool <- 1
      x_axis <- seq.Date(from = origin, by = period.length, along.with = aux$counts)
      x_axis <- format(x_axis, "%m'%y")
    } else if(!is.null(names(aux$counts))){
      x_axis <- names(aux$counts)
    } else {
      x_axis <- seq_along(aux$counts)
    }
    plot(seq_along(aux$data$counts), aux$data$counts,
         type = "b", xaxt = "n", pch = 16,
         xlab = "Time Periods", ylab = "New inclusions",
         xlim = range(seq_along(x_axis)), ylim = c(0, max(aux$data$counts, na.rm = TRUE)))
    if(bool){
      axis(side = 1, at = seq_along(x_axis) - 0.5, labels = x_axis, las = 2)
    } else {
      axis(side = 1, at = seq_along(x_axis), labels = x_axis)
    }
    lines(seq_along(x_axis), aux$counts, col = 2, lwd = 1)
    if(!is.null(aux$bootstrap.ci)){
      lines(seq_along(x_axis), aux$bootstrap.ci[1, ], col = 2, lty = 2)
      lines(seq_along(x_axis), aux$bootstrap.ci[2, ], col = 2, lty = 2)
    }
    if(use.letters & i <= length(LETTERS)) mtext(text = LETTERS[i], cex = 1.1, side = 3, at = -3, line = 1)
    plot(seq_along(aux$data$counts), cumsum(aux$data$counts),
         type = "b", xaxt = "n", pch = 16,
         xlab = "Time Periods", ylab = "Total inclusions",
         xlim = range(seq_along(x_axis)), ylim = c(0, if(max(c(N, aux$cummulative.counts), na.rm = TRUE) > 1.5 * N) 1.5 * N else max(c(N, aux$cummulative.counts), na.rm = TRUE)))
    if(bool){
      axis(side = 1, at = seq_along(x_axis) - 0.5, labels = x_axis, las = 2)
    } else {
      axis(side = 1, at = seq_along(x_axis), labels = x_axis)
    }
    lines(seq_along(x_axis)[-(1:nrow(aux$data))],
          aux$cummulative.counts, col = 2, lwd = 1.5)
    if(!missing(N)) abline(h = N, lty = 2)
    if(!is.null(aux$bootstrap.ci)){
      lines(seq_along(x_axis)[-(1:nrow(aux$data))],
            aux$bootstrap.cummulative.ci[1, ], col = 2, lty = 2)
      lines(seq_along(x_axis)[-(1:nrow(aux$data))],
            aux$bootstrap.cummulative.ci[2, ], col = 2, lty = 2)
    }
    legend("bottomright", legend = c("Observed", "Predicted"), lty = c(1, 1), col = 1:2, bty = "n")
  }
}


#' Plot the output of eventsCount
#' 
#' Plots the cummulative number of events observed and expected using the survival functions given to eventsCount().
#'
#' @param eventsCount.list A list of outputs of eventsCount().
#' @param origin Date of the beginning of the trial.
#' @param period.length Length of the time periods.
#' @param E Total number of events to achieve.
#' @param use.letters Index all elements in the list with a letter when plotting.
#' @param call.mfrow Wether to display all plots in the same window.
#' @param ... Other arguments passed to par().
#' 
#' @return Just produces one or more plots.
#'  
#' @export
events_plot <- function(eventsCount.list, 
                        origin = NULL, period.length = NULL, E, 
                        use.letters = TRUE, call.mfrow = TRUE, ...){
  par_opts <- par()
  on.exit(suppressWarnings(par(par_opts)))
  if(all(c("events", "percent.events") %in% names(eventsCount.list))) eventsCount.list <- list(eventsCount.list)
  if(! all(sapply(eventsCount.list, is.list))) stop("Argument 'eventsCount.list' must be a list of different outputs from 'eventsCount' function.")
  
  if(call.mfrow) par(mfrow = c(length(eventsCount.list), 1))
  par(...)
  
  for(i in seq_along(eventsCount.list)){
    aux <- eventsCount.list[[i]]
    bool <- 0
    if(!is.null(origin) & !is.null(period.length)){
      bool <- 1
      x_axis <- seq.Date(from = origin, by = period.length, along.with = aux$events)
      x_axis <- format(x_axis, "%m'%y")
    } else if(!is.null(names(aux$events))){
      x_axis <- names(aux$events)
    } else {
      x_axis <- seq_along(aux$events)
    }

    plot(seq_along(aux$events_observed), aux$events_observed,
         type = "b", xaxt = "n", pch = 16,
         xlab = "Time Periods", ylab = "Number of Events",
         xlim = range(seq_along(x_axis)), ylim = c(0, max(aux$events, na.rm = TRUE)))
    if(bool){
      axis(side = 1, at = seq_along(x_axis) - 0.5, labels = x_axis, las = 2)
    } else {
      axis(side = 1, at = seq_along(x_axis), labels = x_axis)
    }
    if(!missing(E)) abline(h = E, lty = 2)
    lines(seq_along(x_axis), aux$events, col = 2, lwd = 2)
    if(!is.null(aux$events_ci)){
      lines(seq_along(x_axis), aux$events_ci[1, ], col = 2, lwd = 1, lty = 2)
      lines(seq_along(x_axis), aux$events_ci[2, ], col = 2, lwd = 1, lty = 2)
    }
    if(use.letters & i <= length(LETTERS)) mtext(text = LETTERS[i], cex = 1.1, side = 3, at = -2, line = 1)
    legend("bottomright", legend = c("Observed", "Expected"), lty = c(1, 1), col = 1:2, bty = "n")
  }
}

#' Plot histogram for different lengths of the time periods
#'
#' Decide which calendar partition suits better the data by a visual inspection
#'
#' @param period.length Vector with the lengths for the periods. One plot is produced for each.
#' @param time.inclusion Vector with the inclusion times.
#' @param letters Whether to plot alphabet letters top-left corner of each plot.
#' @param ... Further arguments passed to par.
#'
#' @return No output is produced. Only a plot.
#'
#' @export
periods_plot <- function(period.length, time.inclusion, 
                         use.letters = TRUE, 
                         call.mfrow = TRUE, ...){
  mf <- statTools::descomp(length(period.length))
  if(call.mfrow) par(mfrow = mf)
  par(...)
  quin <- seq_along(period.length)
  names(quin) <- period.length
  for(i in period.length){
    talls <- seq(min(time.inclusion, na.rm = TRUE),
                 max(time.inclusion, na.rm = TRUE) + i,
                 by = i)
    prds <- as.period(times = time.inclusion, period.length = i, numeric = TRUE)

    hist(time.inclusion[-which(max(prds, na.rm = TRUE) == prds)],
         breaks = talls,
         include.lowest = TRUE, right = TRUE,
         col = "grey", xlab = "Calendar time", main = paste0("Period length of ", i))
    if(use.letters) mtext(at = -max(time.inclusion, na.rm = TRUE) / 20, line = 1, text = paste0("(", LETTERS[quin[as.character(i)]], ")"), side = 3)
  }
}

#' Plot survival functions
#'
#' Plot survival functions returned by \code{survFUN()}. KM estimation of the survival function is also plotted (survfit.object).
#'
#' @param surv.funs List of functions returned by \code{survFUN}.
#' @param survfit.object Object returned by survfit.
#' @param period.length Length of each period.
#' @param xlim Maximum xlim for the plot.
#' @param call.par Call \code{par()} with mfrow for plotting all graphics in the same figure.
#' @param letters Whether to plot alphabet letters top-left corner of each plot.
#' @param rev Defaults to TRUE. See argument rev in \code{statTools::descomp}. It is reverses the vector passed to \code{mfrow} argument.
#' @param ... Further arguments passed to par.
#'
#' @return No output is produced. Only a plot.
#'
#' @export
survs_plot <- function(surv.funs,
                       survfit.object,
                       period.length,
                       xlim = 2*max(summary(survfit.object)$time),
                       call.mfrow = TRUE,
                       use.letters = call.mfrow,
                       rev = TRUE, ...){

  mf <- statTools::descomp(length(surv.funs), sum = TRUE, rev = rev)
  if(call.mfrow) par(mfrow = mf)
  par(...)
  
  for(i in seq_along(surv.funs)){
    xs <- seq(0, xlim, length.out = 100)
    plot(survfit.object, conf.int = FALSE, xlim = c(0, xlim))
    lines(xs, surv.funs[[i]](xs), col = 2, type = "l", lwd = 2)
    if(use.letters) mtext(at = - xlim / 20, line = 1, text = LETTERS[i], side = 3, cex = 1.1)
    legend("topright", legend = names(surv.funs)[i], bty = "n", cex = 1.1)
  }
}

