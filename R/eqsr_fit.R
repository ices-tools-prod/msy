#' Stock recruitment fitting
#'
#' Fits one or more stock recruitment relationship to data containted in an
#' FLStock object. If more than one stock recruit relationship is provided, the
#' models are weighted based on smooth AIC weighting (See Buckland et al., 1997).
#'
#' @param stk FLStock object
#' @param nsamp Number of samples (iterations) to take from the stock recruitment
#'              fit.
#' @param models A character vector containing stock recruitment models to use
#'               in the model averaging. User can set any combination of
#'               "Ricker", "Segreg", "Bevholt", "Smooth_hockey".
#' @param id.sr A character vector specifying an id or name for the stock
#'              recruitment fit being run. The default is to use the slot "name"
#'              in the stk parameter is provided
#' @param remove.years A vector specifying the years to remove from the model
#'                     fitting.
#' @param rshift lag ssb by aditional years (default = 0).  As an example, for
#'   some herring stocks, age 1 (1 winter ring) fish were spawned 2 years
#'   previously, in this case, rshift = 1.
#' @return A list containing the following objects:
#' \itemize{
#'   \item `sr.sto` data.frame containing the alpha (a), beta (b), cv and model
#'         names. The number of rows correspond to the value set of `nsamp` in
#'         the function call.
#'   \item `sr.det` The parameters in the stock recruitment model corresponding
#'         to the "best fit" of any given model.
#'   \item `pRec` A matrix of predicted recruitment values. The number of rows
#'         corresponds to the value set in `nsamp` in the function call. The
#'         number of columns matches with the number of years used in the model
#'         fitting.
#'   \item `stk` An FLStock object, same as provided as input by the user.
#'   \item `rby` A data.frame containing the recruitment (rec), spawning stock
#'         biomass (ssb) and year used in the fitting of the data.
#'   \item `id.sr` A string containing run name (taken from the `id.sr` argument)
#'
#' }
#'
#' @references
#'
#' Buckland, S.T., K.P. Burnham & N.H. Augustin (1997). Model selection: An integral part of inference.
#' Biometrics 53, 603-618.
#' DOI: \href{https://doi.org/10.2307/2533961}{10.2307/2533961}
#'
#' @seealso
#' \code{\link{eqsr_plot}} plots a simulation of predictive recruitment
#' from the fit, and shows a summary of the contributions of each stock
#' recruitment model to the model average fit.
#'
#' @examples
#' \dontrun{
#' data(icesStocks)
#' FIT <- eqsr_fit(icesStocks$saiNS,
#'                 nsamp = 1000,
#'                 models = c("Ricker", "Segreg"))
#'
#' # summary of individual fits
#' FIT$sr.det
#' }
#'
#' @export
eqsr_fit <- function(stk, nsamp = 1000, models = c("Ricker","Segreg","Bevholt"),
                     id.sr = FLCore::name(stk), remove.years = NULL, rshift = 0)
{
  # some checks on the model argument
  if (!is.character(models)) stop("models arg should be character vector giving names of stock recruit models")
  if(any(models %in% c("ricker","segreg","bevholt"))) {
    stop("Please note that the msy stock-recruitment functions have been renamed:
   ricker -> Ricker
   bevholt -> Bevholt
   segreg ->  Segreg
   smooth_hockey -> Smooth_hockey
   This was done to resolve conflicts with same named functions in the FLCore-package.
   ERGO: use a capital in the first letter if you want to call these functions")
  }

  # get correct recruitment vector for each SSB
  # dims$min is the minimum age => recruitment age
  dms <- FLCore::dims(stk)
  rec <- c(FLCore::rec(stk))
  if (dms$min > 0 | rshift > 0)
  {
    ssb_lag <- dms$min + rshift
    rec <- c(rec[-seq(ssb_lag)], rep(NA, ssb_lag))
  }

  # combine all required data together
  # note year is the year that SSB had that value
  data <-
    data.frame(year = with(dms, minyear:maxyear),
               rec = rec,
               ssb = c(FLCore::ssb(stk)),
               fbar = c(FLCore::fbar(stk)),
               landings = c(FLCore::landings(stk)),
               catch = c(FLCore::catch(stk)))

  # remove years with NA recruitment
  data <- data[stats::complete.cases(data),]

  # which years to use in the fit
  data$remove.years <- FALSE
  if (!is.null(remove.years)) {
    remove.years <- data$year[data$year %in% remove.years]
    message("removing (ssb) years:\n\t",
            paste(remove.years, collapse = ", "),
            "\n  from the recruitment fitting procedure.")
    data$remove.years <- data$year %in% remove.years
  }

  # run model averaging
  srfit <- eqsr_Buckland(data[!data$remove.years,c("year", "rec", "ssb")],
                         nsamp,
                         models)

  # create output object
  out <- c(srfit, list(stk = stk, rby = data, id.sr = id.sr))
  class(out) <- c("eqsr_fit", c("list", "vector"))

  out
}
