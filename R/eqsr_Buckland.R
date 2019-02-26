
eqsr_Buckland <- function(data, nsamp = 5000, models = c("Ricker","Segreg","Bevholt"), ...)
{
  # useful objects
  nllik <- function(param, ...) -1 * llik(param, ...)
  ndat <- nrow(data)

  #--------------------------------------------------------
  # get best fit for each model
  #--------------------------------------------------------
  fits <-
    do.call(rbind,
            lapply(models,
                   function(mod)
                     with(stats::nlminb(initial(mod, data), nllik, data = data, model = mod, logpar = TRUE),
                          data.frame(a = exp(par[1]), b = exp(par[2]), cv = exp(par[3]), model = mod))))
  row.names(fits) <- NULL

  if (nsamp > 0) {
    #--------------------------------------------------------
    # Fit models on bootstrap resamples
    #--------------------------------------------------------
    fit <- lapply(1:nsamp, function(i)
    {
      sdat <- data[sample(1:ndat, replace = TRUE),]

      fits <- lapply(models, function(mod) stats::nlminb(initial(mod, sdat), nllik, data = sdat, model = mod, logpar = TRUE))

      best <- which.min(sapply(fits, "[[", "objective"))

      with(fits[[best]], c(a = exp(par[1]), b = exp(par[2]), cv = exp(par[3]), model = best))
    })

    fit <- as.data.frame(do.call(rbind, fit))
    fit$model <- models[fit $ model]

    # summarise and join to deterministic fit
    tmp <- plyr::ddply(fit, "model", plyr::summarise, n = length(model))
    tmp$prop <- tmp$n / sum(tmp$n)
    fits <- plyr::join(fits, tmp, by = "model")

    #--------------------------------------------------------
    # get posterior distribution of estimated recruitment
    #--------------------------------------------------------
    pred <- t(sapply(seq(nsamp), function(j) exp(match.fun(fit$model[j]) (fit[j,], sort(data$ssb))) ))
    dimnames(pred) <- list(model = fit$model, ssb = data$ssb)
  } else {
    fit <- pred <- NULL
    fits$n <- 0
    fits$prop <- 0
  }

  list(sr.sto = fit, sr.det = fits, pRec = pred)
}
