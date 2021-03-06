% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funcs_2.R
\name{period_adj}
\alias{period_adj}
\title{Removing last time period}
\usage{
period_adj(data, period.length, last.period = NULL,
  varnames = getOption("RCT_names_match"), origin = NULL, numeric = TRUE)
}
\arguments{
\item{data}{Data to remove last time period. A data frame with variables inclusion, y, abandon and delta at least.
Other variables won't be modified.}

\item{period.length}{Length of each period.}

\item{last.period}{Last period to take into account for the adjustment.}

\item{varnames}{By default \code{getOption("RCT_names_match")}. This is a character vector with names "inclusion",
"delta" and "y" (in any order), and its corresponding name in the argument data (gives the matching for accessing
 to inclusion times, censore indicator, and time of follow-up).}

\item{origin}{Passed to \code{\link{as.period}}.}

\item{numeric}{Passed to \code{\link{as.period}}.}
}
\value{
The same data frame but without the patients included in the last period of time (removed), and some
modified delta's.
}
\description{
Removes last time period from the data, and changes those patients that had event in this period. Useful in case that
the last time period has not finished.
}
\details{
Usually the calendar time partition does not match exactly with the current follow-up time of the RCT. In these cases, the last
time period is shorter than the others and one solution is to ignore the data in this time period, removing those patients and
correcting censores and events from others.
}

