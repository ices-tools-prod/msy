
eqsr_Buckland <- function(data, nsamp = 5000, models = c("Ricker","Segreg","Bevholt"), ...)
{
  # useful objects
  nllik <- function(param, ...) -1 * llik(param, ...)
  ndat <- nrow(data)

  #--------------------------------------------------------
  # get best fit for each model
  #--------------------------------------------------------
  onefit <- function(mod) {
    fit <-
      stats::nlminb(
        initial(mod, data),
        nllik, data = data,
        model = mod, logpar = TRUE,
        control = list(iter.max = 500, eval.max = 500)
      )
    out <-
      data.frame(
        a = exp(fit$par[1]),
        b = exp(fit$par[2]),
        cv = exp(fit$par[3]),
        llik = -1 * fit$objective,
        model = mod,
        stringsAsFactors = FALSE
      )
    out
  }
  sr.det <- do.call(rbind, lapply(models, onefit))
  row.names(sr.det) <- NULL

  if (nsamp > 0) {
    #--------------------------------------------------------
    # Fit models on bootstrap resamples
    #--------------------------------------------------------
    sr.sto <- lapply(1:nsamp, function(i)
    {
      sdat <- data[sample(1:ndat, replace = TRUE),]

      fits <- lapply(models, function(mod) stats::nlminb(initial(mod, sdat), nllik, data = sdat, model = mod, logpar = TRUE))

      best <- which.min(sapply(fits, "[[", "objective"))

      with(fits[[best]], c(a = exp(par[1]), b = exp(par[2]), cv = exp(par[3]), model = best))
    })

    sr.sto <- as.data.frame(do.call(rbind, sr.sto))
    sr.sto$model <- models[sr.sto$model]

    # summarise and join to deterministic fit
    tmp <- table(sr.sto$model)
    sr.det$n <- unname(tmp[sr.det$model])
    sr.det$prop <- sr.det$n / sum(sr.det$n)
  } else {
    sr.sto <- NULL
    sr.det$n <- 0
    sr.det$prop <- 0
  }

  #list(sr.sto = fit, sr.det = fits, pRec = pred)
  list(sr.sto = sr.sto, sr.det = sr.det)
}
