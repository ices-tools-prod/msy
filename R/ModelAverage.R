
#' @title Stock recruitment fit
#'
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
  data <- data[complete.cases(data),]

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



#' stock recruitment function
#'
#'
#' @param data data.frame containing stock recruitment data
#' @param runid A character vector
#' @param nsamp Number of samples
#' @param models A character vector
#' @param ... Additional arguements
#' @return log recruitment according to model
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
fitModelsBuck <- function(data, runid, nsamp = 5000, models = c("ricker","segreg","bevholt"), ...)
{

#--------------------------------------------------------
# Fit models
#--------------------------------------------------------

  nllik <- function(param, ...) -1 * llik(param, ...)
  ndat <- nrow(data)
  fit <- lapply(1:nsamp, function(i)
    {
      sdat <- data[sample(1:ndat, replace = TRUE),]
       
      fits <- lapply(models, function(mod) nlminb(initial(mod, sdat), nllik, data = sdat, model = mod, logpar = TRUE))
    
      best <- which.min(sapply(fits, "[[", "objective"))
    
      with(fits[[best]], c(a = exp(par[1]), b = exp(par[2]), cv = exp(par[3]), model = best))
    })
    
  fit <- as.data.frame(do.call(rbind, fit))
  fit $ model <- models[fit $ model]
  
#--------------------------------------------------------
# get posterior distribution of estimated recruitment
#--------------------------------------------------------
  pred <- t(sapply(seq(nsamp), function(j) exp(match.fun(fit $ model[j]) (fit[j,], sort(data $ ssb))) ))


#--------------------------------------------------------
# get best fit for each model
#--------------------------------------------------------
  fits <- 
    do.call(rbind,
      lapply(models, 
           function(mod) 
               with(nlminb(initial(mod, data), nllik, data = data, model = mod, logpar = TRUE), 
                 data.frame(a = exp(par[1]), b = exp(par[2]), cv = exp(par[3]), model = mod))))

  
  list(fit = fit, pred = pred, fits = fits, data = data, stknam = runid)
}


#' stock recruitment function
#'
#'
#' @param data data.frame containing stock recruitment data
#' @param runid A character vector
#' @param delta A value
#' @param nburn An integer, specifying the burn-in period
#' @param nsamp An integer, specifying the number of samples
#' @param models A character vector specifying stock-recruitment models
#' @param ... Additional arguments
#' @return log recruitment according to model
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
fitModelsSimmonds <- function(data, runid, delta = 1.3, nburn = 10000, nsamp = 5000, models = c("ricker","segreg","bevholt"), ...)
{


#--------------------------------------------------------
# Fit models
#--------------------------------------------------------

  if (!is.character(models)) stop("models arg should be character vector giving names of stock recruit models")

  # fit individual stock recruit relationships using Metropolis hasings MCMC algorithm
  fits <- lapply(models, function(x) cbind.data.frame(mod = x, MH(nburn + nsamp, nburn, data, delta = delta[1], model = x), stringsAsFactors = FALSE))
  names(fits) = models

#--------------------------------------------------------
# get posterior summaries on log scale for RJMCMC (NOT YET!)
#--------------------------------------------------------

  lfits <- lapply(fits, function(x) {x[-(1:2)] <- lapply(x[-(1:2)], log); x})

  mu <- lapply(fits, function(x) colMeans(x[-(1:2)]))
  Sig <- lapply(fits, function(x) var(x[-(1:2)]))


#--------------------------------------------------------
# make stock recruit object with correct probabilities of each model
#--------------------------------------------------------
 
  # we can remove unlikely models here....  probably not nessisary with so few
  

  #TODO this will be the BMA
  # Johns way for now
  post.mod.prob <- sapply(fits, function(x) 1/mean(exp(-x $ llik)))  
  mod.samp <- sample(models, nsamp, replace = TRUE, prob = post.mod.prob / sum(post.mod.prob))
  
  fit <- do.call(rbind, lapply(seq(mod.samp), function(i) fits[[mod.samp[i]]][i,]))

#--------------------------------------------------------
# get posterior distribution of estimated recruitment
#--------------------------------------------------------
  preds <- t(sapply(seq(nsamp), function(j) exp(match.fun(fit $ model[j]) (fit[j,], sort(data $ ssb))) ))
  
  list(fit = fit, preds = preds, data = data, stknam = runid)
}



#' plot simulated predictive distribution of recruitment
#'
#'
#' @param fit an fitted MCMC returned from \code{fitModels}
#' @param n Number of random recruitment draws to plot
#' @return NULL produces a plot
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
SRplot <- function (fit, n = 5000) 
{
  modset <- fit $ fit
  data <- fit $ data

  ssb <- data $ ssb
  rec <- data $ rec

  #TODO Blim <- res $ Blim
  #TODO PBlim <- res $ PBlim

  mn <- length(ssb)
  minSSB <- min(ssb, max(ssb)*0.05)
  maxSSB <- max(ssb)*1.1
  maxrec <- max(rec*1.5)
  
  x2 <- reshape2::melt(tapply(fit$fit$model,fit$fit$model,length))
  x2$lab <- paste(x2$Var1,round(x2$value/sum(x2$value),2))

  plot(ssb, rec, xlim = c(0, maxSSB), ylim = c(0, maxrec), type = "n", 
       xlab = "SSB ('000 t)", ylab="Recruits", main = paste("Predictive distribution of recruitment\nfor", fit $ stknam))

  out <-
  do.call(rbind, lapply(sample(1:nrow(modset), 500), 
    function(i)
    {
      fssb <- runif(500, minSSB, maxSSB)
      FUN <-  match.fun(modset $ model[i])
      frec <- exp( FUN(modset[i,], fssb) + rnorm(500, sd = modset $ cv[i]) )
      #points(fssb, frec, pch = 20, col = paste0(grey(0), "05"), cex = 0.0625)

      data.frame(ssb = fssb, rec = frec)
    }))
  points(out$ssb[1:n], out$rec[1:n], pch = 20, col = paste0(grey(0), "05"), cex = 1)
  out $ grp <- with(out, floor(10 * (ssb - min(ssb)) / (max(ssb) - min(ssb) + 0.001)))
  out $ mid.grp <- with(out, (grp + 0.5) / 10 * (max(ssb) - min(ssb)) + min(ssb))

  #TODO use 
  summ <- with(out, 
    t(simplify2array( tapply(rec, grp, quantile, c(0.5, .05, .95)) )))

  mid.grp <- sort(unique(out $ mid.grp))

  lines(mid.grp, summ[,1], col = 7, lwd = 3)
  lines(mid.grp, summ[,2], col = 4, lwd = 3)
  lines(mid.grp, summ[,3], col = 4, lwd = 3)

  x <- fit $ fits
  y <- seq(1, round(max(ssb)), length = 100)
  sapply(1:nrow(x), function(i) lines(y, exp(match.fun(as.character(x$model[i])) (x[i,], y)), col = "black", lwd = 2, lty = i))

  lines(ssb, rec, col = 10)
  points(ssb, rec, pch = 19, col = 10, cex = 1.25)
  
  for (i in 1:nrow(x2)) {
    text(0.2*maxSSB,maxrec*(1-i/10),x2$lab[i],cex=0.9)
  }

  
  #TODO plot Blim
  #lines(PBlim$mids,0.1*maxrec/max(PBlim$counts)*PBlim$counts,col=5,lwd=2)
  #lines(c(Blim,Blim),c(0,maxrec),col=5,lwd=2)
}


