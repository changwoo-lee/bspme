#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL


#' @name NO2_Jan2012
#' @title Daily average NO2 concentrations in and around the Harris County, Texas, in Jan 2012
#'
#' @description
#' This dataset contains daily average NO2 (nitrogen dioxide) concentrations obtained from 21 monitoring stations in and around Harris County, Texas, in January 2012.
#'
#' @usage data(NO2_Jan2012)
#'
#' @format A data frame with 651 (21 sites x 31 days) rows and 5 variables:
#' \describe{
#'   \item{date}{date in POSIXct format}
#'   \item{site_name}{monitoring station name}
#'   \item{lon}{monitoring station longitude}
#'   \item{lat}{monitoring station latitude}
#'   \item{lnNO2}{natural logarithm of daily average NO2 concentrations, measured in parts per billion by volume (ppbv)}
#' }
NULL


#' @name health_sim
#' @title Simulated health data
#'
#' @description
#' Simulated health data associated with ln(NO2) concentration on Jan 10, 2012. For details, see \code{health_sim.R}.
#'
#' @usage data(health_sim)
#'
#' @format A data frame with n = 2000 rows and 6 variables:
#' \describe{
#'   \item{Y}{simulated continuous health outcome}
#'   \item{Ybinary}{simulated binary health outcome}
#'   \item{lon}{simulated health subject longitude}
#'   \item{lat}{simulated health subject latitude}
#'   \item{Z}{simulated covariate (p=1) that is not subject to measurement error}
#'   \item{X_true}{true ln(NO2) exposure used for simulating health outcome}
#' }
NULL
