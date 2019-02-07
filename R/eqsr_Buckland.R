
eqsr_Buckland <- function(data, nsamp = 5000, models = c("Ricker","Segreg","Bevholt"), ...)
{

  ## dummy
  model <- 0
  #--------------------------------------------------------
  # Fit models
  #--------------------------------------------------------

  nllik <- function(param, ...) -1 * llik(param, ...)
  ndat <- nrow(data)
  fit <- lapply(1:nsamp, function(i)
  {
    sdat <- data[sample(1:ndat, replace = TRUE),]

    fits <- lapply(models, function(mod) stats::nlminb(initial(mod, sdat), nllik, data = sdat, model = mod, logpar = TRUE))

    best <- which.min(sapply(fits, "[[", "objective"))

    with(fits[[best]], c(a = exp(par[1]), b = exp(par[2]), cv = exp(par[3]), model = best))
  })

  fit <- as.data.frame(do.call(rbind, fit))
  fit $ model <- models[fit $ model]

  #--------------------------------------------------------
  # get posterior distribution of estimated recruitment
  #--------------------------------------------------------
  pred <- t(sapply(seq(nsamp), function(j) exp(match.fun(fit $ model[j]) (fit[j,], sort(data $ ssb))) ))
  #pred <- t(sapply(seq(nsamp), function(j) exp(get(fit $ model[j], , pos = "package:msy", mode = "function") (fit[j,], sort(data $ ssb))) ))


  #--------------------------------------------------------
  # get best fit for each model
  #--------------------------------------------------------
  fits <-
    do.call(rbind,
            lapply(models,
                   function(mod)
                     with(stats::nlminb(initial(mod, data), nllik, data = data, model = mod, logpar = TRUE),
                          data.frame(a = exp(par[1]), b = exp(par[2]), cv = exp(par[3]), model = mod))))

  tmp <- plyr::ddply(fit,c("model"), plyr::summarise, n=length(model))
  tmp$prop <- tmp$n/sum(tmp$n)
  fits <- plyr::join(fits,tmp,by="model")

  dimnames(pred) <- list(model=fit$model,ssb=data$ssb)

  return(list(sr.sto = fit, sr.det = fits, pRec = pred))
}
