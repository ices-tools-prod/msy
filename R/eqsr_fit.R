#' Stock recruitment fit
#'
#' Fits a stock recruitment relationship to data containted in an FLStock object.
#' If more than one stock recruit relationship is provided, the models are weighted
#' based on smooth AIC weighting (See Buckland et al.).
#'
#' @param stk FLStock object
#' @param nsamp Number of samples (iterations)
#' @param models A character vector containing sr-models to use. User can set
#' any combination of "Ricker", "Segreg", "Bevholt", "Smooth_hockey".
#' @param method A character vector. Currently only "Buckland" is implemented.
#' @param id.sr A character vector specifying an id for the stock recruitment
#' model. If not specified (default) the slot "name" in the FLStock is used.
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
#' @author Colin Millar \email{colin.millar@@ices.dk}
#' @export
eqsr_fit <- function(stk, nsamp = 5000, models = c("Ricker","Segreg","Bevholt"),
                     method = "Buckland",
                     id.sr = NULL, remove.years = NULL, delta = 1.3, nburn = 10000)
{

  if(any(models %in% c("ricker","segreg","bevholt"))) {
    return(cat("Please note that the msy stock-recruitment functions have been renamed:
   ricker -> Ricker
   bevholt -> Bevholt
   segreg ->  Segreg
   smooth_hockey -> Smooth_hockey
   This was done to resolve conflicts with same named functions in the FLCore-package.
   ERGO: use a capital in the first letter if you want to call these functions"))
  }

  dms <- FLCore::dims(stk)
  rage <- dms $ min
  if (rage == 0)
  { x = FLCore::stock.n(stk)[1,drop=TRUE]
  } else {
    x = c(FLCore::stock.n(stk)[1,-seq(rage),drop=TRUE],rep(NA,rage))
  }

  rby <- data.frame(year = with(dms, minyear:maxyear),
                      rec = x,
                      ssb = FLCore::ssb(stk)[drop=TRUE],
                      fbar = FLCore::fbar(stk)[drop=TRUE],
                      landings=FLCore::landings(stk)[drop=TRUE],
                      catch=FLCore::catch(stk)[drop=TRUE])

  row.names(rby) <- NULL
  rby <- rby[!is.na(rby$rec),]
  # This is for stuff later down the pipes
  data <- rby[,1:3]

  # EINAR: strange that here only the ssb is set to as NA
  #        question how this affect what happens further down the line
  if (!is.null(remove.years)) {
    data $ ssb[data $ year %in% remove.years] <- NA
  }
  #--------------------------------------------------------
  # tidy data - remove nas
  #--------------------------------------------------------
  data <- data[stats::complete.cases(data),]

  if (is.null(id.sr)) id.sr <- FLCore::name(stk)

  method <- match.arg(method, c("Buckland","Simmonds","King","Cadigan"))
  if (!is.character(models)) stop("models arg should be character vector giving names of stock recruit models")

  if (method == "Buckland") {
    return(c(eqsr_Buckland(data, nsamp, models), list(stk = stk,rby=rby, id.sr=id.sr)))
  } else
  {
    cat("The", method, "is not ready yet!  Working on it!\n")
  }
}
