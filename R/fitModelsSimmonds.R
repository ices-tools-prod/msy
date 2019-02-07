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
#' @author Colin Millar \email{colin.millar@@ices.dk}
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
  Sig <- lapply(fits, function(x) stats::var(x[-(1:2)]))


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
