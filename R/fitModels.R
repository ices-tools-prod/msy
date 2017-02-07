#' Stock recruitment fit
#'
#' Unmaintained version of eqsr_fit.
#'
#' @param stk FLStock object
#' @param nsamp Number of samples
#' @param models A character vector containing sr-models to use. User can set
#' any combination of "ricker","segreg","bevholt".
#' @param method A character vector. Currently only "Buckland" is implemented.
#' @param runid A character vector specifying run name
#' @param remove.years A vector specifying the years to remove
#' @param delta A value, used in method "Simmonds" (not implemented)
#' @param nburn An integer, used in method Simmonds (not implemented)
#' @return A list containing the following objects:
#' \itemize{
#' \item fit data.frame containing the alpha (a), beta (b), cv and model names.
#' The number of rows correspond to the value set in nsamp in the function call.
#' \item pred A vector of predicted recruitment values. The length of the vector
#' corresponds to the value set in nsamp in the function call.
#' \item fits The parameters in the stock recruitment model corresponding to the
#' "best fit" of any given model.
#' \item data data.frame containing the recruitment (rec), spawning stock
#' biomass (ssb) and year used in the fitting of the data.
#' \item stknam A character vector containing stock name
#' \item stk FLStock object, same as provided as input by the user.
#' }
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}

fitModels <- function(stk, nsamp = 5000, models = c("ricker","segreg","bevholt"),
               method = "Buckland",
               runid = NULL, remove.years = NULL, delta = 1.3, nburn = 10000)
{
  message("NOTE: THIS FUNCTION IS NO LONGER MAINTAINED, USE FUNCTION eqsr_fit")
  dms <- FLCore::dims(stk)
  rage <- dms $ min
  if (rage == 0)
  {
    data <- data.frame(rec = FLCore::stock.n(stk)[1,drop=TRUE],
                       ssb = FLCore::ssb(stk)[drop=TRUE],
                       year = with(dms, 1:year + minyear - 1))
  } else
  {
    data <- data.frame(rec = FLCore::stock.n(stk)[1,-seq(rage),drop=TRUE],
                       ssb = FLCore::ssb(stk)[1,seq(dms$year - rage),drop=TRUE],
                       year = with(dms, (rage+1):year + minyear - 1))
  }

  if (!is.null(remove.years)) {
    data $ ssb[data $ year %in% remove.years] <- NA
  }
  #--------------------------------------------------------
  # tidy data - remove nas
  #--------------------------------------------------------
  data <- data[stats::complete.cases(data),]

  if (is.null(runid)) runid <- FLCore::name(stk)

  method <- match.arg(method, c("Buckland","Simmonds","King","Cadigan"))
  if (!is.character(models)) stop("models arg should be character vector giving names of stock recruit models")

  if (method == "Buckland") {
    c(fitModelsBuck(data, runid, nsamp, models), list(stk = stk))
  } else
  {
    cat("The", method, "is not ready yet!  Working on it!\n")
  }
}
