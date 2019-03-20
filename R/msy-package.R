#' @docType package
#'
#' @name msy-package
#'
#' @aliases msy
#'
#' @title Estimation of equilibrium reference points for fisheries
#'
#' @description
#' Methods to estimate equilibrium reference points for fisheries data.
#' Currently data must be converted into FLStock objects of the FLR (Fisheries
#' Library in R) style, defined in the R package FLCore
#'
#' @details
#' \emph{Model fitting and simulation:}
#' \tabular{ll}{
#'   \code{\link{eqsr_fit}} \tab fitting stock recruit models to data\cr
#'   \code{\link{eqsim_run}} \tab simulation am 'equilibrium' population state\cr
#' }
#' \emph{Plotting:}
#' \tabular{ll}{
#'   \code{\link{eqsr_plot}} \tab plot stock recuitment fit\cr
#'   \code{\link{eqsim_plot}} \tab plot summary of simulation showing reference points\cr
#'   \code{\link{eqsim_plot_range}} \tab plot summary of MSY ranges reference points\cr
#' }
#' \emph{Example data:}
#' \tabular{ll}{
#'   \code{\link{icesStocks}} \tab A list of various stocks\cr
#' }
#'
#' @author John Simmonds, Einar Hjorleifsson, Carmen Fernandez and Colin Millar.
#'
#' @references
#' ICES (2015) Report of the Workshop to consider F MSY ranges for stocks in
#' ICES categories 1 and 2 in Western Waters (WKMSYREF4).
#' \href{http://ices.dk/sites/pub/Publication\%20Reports/Expert\%20Group\%20Report/acom/2015/WKMSYREF4/01\%20WKMSYREF4\%20Report.pdf}{01
#' WKMSYREF4 Report.pdf}
#'
#' ICES (2017) ICES fisheries management reference points for category 1 and 2
#' stocks.
#' DOI: \href{https://doi.org/10.17895/ices.pub.3036}{10.17895/ices.pub.3036}
#'
#' Buckland, S.T., K.P. Burnham & N.H. Augustin (1997). Model selection: An integral part of inference.
#' Biometrics 53, 603-618.
#' DOI: \href{https://doi.org/10.2307/2533961}{10.2307/2533961}
#'
#' To explore the code of the package see the GitHub repo:
#' \href{https://github.com/ices-tools-prod/msy}{ices-tools-prod/msy}
#'
#' @import graphics
#' @import mgcv

NULL
