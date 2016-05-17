# DATA --------------------------------------------------------------------

#' Inclussions and ACS events of a RCT.
#'
#' Times of inclusion and time of follow-up until event or censorship for 881 patients.
#'
#' @format A dataset with 881 rows and 6 variables:
#' \describe{
#' \item{id}{Id of the patient.}
#' \item{inclusion}{Days elapsed between the start of the RCT and the inclusion of the patient.}
#' \item{abandon}{Binary indicator: 1 withdrawal, 0 still in the trial.}
#' \item{delta}{Censoring indicator (1 for events).}
#' \item{y}{Follow-up time until the event or a censorship.}
#' \item{end}{Days elapsed between the inclusion and the last follow-up visit or event.}
#' }
#'
#' @docType data
#' @name RCTdata
#' @usage data(RCTdata)
#' @export
"RCTdata"
