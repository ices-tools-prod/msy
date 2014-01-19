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
#' @export
eqsr_fit <- function(stk, nsamp = 5000, models = c("ricker","segreg","bevholt"), 
                      method = "Buckland",
                      runid = NULL, remove.years = NULL, delta = 1.3, nburn = 10000) 
{
  dms <- dims(stk)
  rage <- dms $ min
  if (rage == 0)
  {
    data <- data.frame(rec = stock.n(stk)[1,drop=TRUE],
                       ssb = ssb(stk)[drop=TRUE],
                       year = with(dms, 1:year + minyear - 1)) 
  } else
  {
    data <- data.frame(rec = stock.n(stk)[1,-seq(rage),drop=TRUE],
                       ssb = ssb(stk)[1,seq(dms$year - rage),drop=TRUE],
                       year = with(dms, (rage+1):year + minyear - 1)) 
  }
  
  if (!is.null(remove.years)) {
    data $ ssb[data $ year %in% remove.years] <- NA
  }
  #--------------------------------------------------------
  # tidy data - remove nas
  #--------------------------------------------------------
  data <- data[complete.cases(data),]
  
  if (is.null(runid)) runid <- name(stk)
  
  method <- match.arg(method, c("Buckland","Simmonds","King","Cadigan"))
  if (!is.character(models)) stop("models arg should be character vector giving names of stock recruit models")
  
  if (method == "Buckland") {
    c(eqsr_Buckland(data, runid, nsamp, models), list(stk = stk))
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
eqsr_Buckland <- function(data, runid, nsamp = 5000, models = c("ricker","segreg","bevholt"), ...)
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




#' plot simulated predictive distribution of recruitment
#'
#'
#' @param fit an fitted MCMC returned from \code{eqsr_fit}
#' @param n Number of random recruitment draws to plot
#' @param ggPlot Flag, if FALSE (default) plot base graphics, if true
#' do a ggplot
#' @return NULL produces a plot
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
eqsr_plot <- function (fit, n = 5000, ggPlot=FALSE) 
{
  modset <- fit $ fit
  data <- fit $ data
  
  #ssb <- data $ ssb
  #rec <- data $ rec
  
  #TODO Blim <- res $ Blim
  #TODO PBlim <- res $ PBlim
  
  #mn <- length(data$ssb)
  minSSB <- min(data$ssb, max(data$ssb)*0.0125)
  maxSSB <- max(data$ssb)*1.1
  maxrec <- max(data$rec*1.5)
  
  x2 <- melt(tapply(fit$fit$model,fit$fit$model,length))
  x2$lab <- paste(x2$Var1,round(x2$value/sum(x2$value),2))
  names(x2) <- c("variable","value","Model")
  x2$variable <- as.character(x2$variable)
  
  x3 <- melt(tapply(fit$fit$model,fit$fit$model,length))
  x3$prop <- round(x2$value/sum(x2$value),2)
  # 1) Take 500 samples from nrow(modset) - default value is 5000
  # 2) Make 500 uniformly distributed SSBs (fssb)
  # 3) Generate recruitment value
  # A total of 250000 (500*500) recruitment values are thus generated. The
  # weight that each model gets when generating the recruits is equivalent
  # to the number of occurrences in the fit$fit object, i.e table(fit$fit$models)
  
  out <-
    do.call(rbind, lapply(sample(1:nrow(modset), 500), 
                          function(i)
                          {
                            fssb <- runif(500, minSSB, maxSSB)
                            FUN <-  match.fun(modset $ model[i])
                            frec <- exp( FUN(modset[i,], fssb) + rnorm(500, sd = modset $ cv[i]) )
                            srModel <- modset$model[i]
                            #points(fssb, frec, pch = 20, col = paste0(grey(0), "05"), cex = 0.0625)
                            data.frame(ssb = fssb, rec = frec, variable=srModel)
                          }))
  # group the ssbs into 10 bins
  out $ grp <- with(out, floor(10 * (ssb - min(ssb)) / (max(ssb) - min(ssb) + 0.001)))
  # find the midvalue of ssb within each group
  out $ mid.grp <- with(out, (grp + 0.5) / 10 * (max(ssb) - min(ssb)) + min(ssb))
  # calculate the recruitment median and 5th and 95th percentile within each
  # ssb group and then plot the distribution
  summ <- with(out, 
               t(simplify2array( tapply(rec, grp, quantile, c(0.5, .05, .95)) )))
  mid.grp <- sort(unique(out $ mid.grp))
  d1 <- data.frame(ssb=mid.grp,p50=summ[,1],p05=summ[,2],p95=summ[,3])
  
  if(!ggPlot) {
    plot(data$ssb, data$rec, xlim = c(0, maxSSB), ylim = c(0, maxrec), type = "n", 
         xlab = "SSB ('000 t)", ylab="Recruits", main = paste("Predictive distribution of recruitment\nfor", fit $ stknam))
    points(out$ssb[1:n], out$rec[1:n], pch = 20, col = paste0(grey(0), "05"), cex = 1)
    lines(mid.grp, summ[,1], col = 7, lwd = 3)
    lines(mid.grp, summ[,2], col = 4, lwd = 3)
    lines(mid.grp, summ[,3], col = 4, lwd = 3)
    # plot the best fit for each model as a line
    x <- fit $ fits
    y <- seq(1, round(max(data$ssb)), length = 100)
    sapply(1:nrow(x), function(i) lines(y, exp(match.fun(as.character(x$model[i])) (x[i,], y)), col = "black", lwd = 2, lty = i))
    # plot the observation points
    lines(data$ssb, data$rec, col = "red")
    points(data$ssb, data$rec, pch = 19, col = "red", cex = 1.25)
    # plot the model weights on the graph
    for (i in 1:nrow(x2)) {
      text(0.2*maxSSB,maxrec*(1-i/10),x2$lab[i],cex=0.9)
    }
  } else { # ggplot
    
    # best model fit
    x <- fit $ fits
    ssb <- seq(1, round(max(data$ssb)), length = 100)
    z <- sapply(1:nrow(x), function(i) rec <- exp(match.fun(as.character(x$model[i])) (x[i,], ssb)))
    d2 <- as.data.frame(cbind(ssb,z))
    names(d2) <- c("ssb",x$model)
    d2 <- melt(d2,id.vars="ssb")
    d2 <- join(d2,x2[,c("variable","Model")])
    out <- join(out,x2[,c("variable","Model")])
    ggplot(out[1:n,]) + 
      theme_bw() +
      geom_point(aes(x=ssb,y=rec,colour=variable),size=1) +
      geom_line(data=d1,aes(x=ssb,y=p05),colour="yellow") +
      geom_line(data=d1,aes(x=ssb,y=p95),colour="yellow") +
      geom_line(data=d1,aes(ssb,p50),col="yellow",lwd=2) +
      geom_line(data=d2,aes(ssb,value,colour=variable),lwd=1) +
      #scale_colour_brewer(palette="Set1") +
      scale_colour_manual(values=c("bevholt"="#F8766D","ricker"="#00BA38","segreg"="#619CFF"),
                          labels=paste(x3$Var1,x3$prop)) +
      coord_cartesian(ylim=c(0,quantile(out$rec[1:n],0.99))) +
      geom_path(data=fit$data,aes(ssb,rec),col="grey") +
      geom_text(data=fit$data,aes(ssb,rec,label=substr(year,3,4)),angle=45,size=3) +
      theme(legend.position = c(0.15,0.80)) +
      labs(x="Spawning stock biomass",y="Recruitment",colour="Model")
    
  }
  
  #"#F8766D" "#00BA38" "#619CFF"
  #TODO plot Blim
  #lines(PBlim$mids,0.1*maxrec/max(PBlim$counts)*PBlim$counts,col=5,lwd=2)
  #lines(c(Blim,Blim),c(0,maxrec),col=5,lwd=2)
}