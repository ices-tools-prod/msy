

#' @title Stock recruitment fit
#'
#'
#' @param stk FLStock object
#' @param nsamp Number of samples (iterations)
#' @param models A character vector containing sr-models to use. User can set
#' any combination of "ricker","segreg","bevholt".
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
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
eqsr_fit <- function(stk, nsamp = 5000, models = c("ricker","segreg","bevholt"), 
                     method = "Buckland",
                     id.sr = NULL, remove.years = NULL, delta = 1.3, nburn = 10000) 
{
  dms <- dims(stk)
  rage <- dms $ min
  if (rage == 0)
  { x = stock.n(stk)[1,drop=TRUE]
  } else {
    x = c(stock.n(stk)[1,-seq(rage),drop=TRUE],rep(NA,rage))
  }
  
  rby <- data.frame(year = with(dms, minyear:maxyear),
                      rec = x,
                      ssb = ssb(stk)[drop=TRUE],
                      fbar = fbar(stk)[drop=TRUE],
                      landings=landings(stk)[drop=TRUE],
                      catch=catch(stk)[drop=TRUE])

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
  data <- data[complete.cases(data),]
  
  if (is.null(id.sr)) id.sr <- name(stk)
  
  method <- match.arg(method, c("Buckland","Simmonds","King","Cadigan"))
  if (!is.character(models)) stop("models arg should be character vector giving names of stock recruit models")
  
  if (method == "Buckland") {
    return(c(eqsr_Buckland(data, nsamp, models), list(stk = stk,rby=rby, id.sr=id.sr)))
  } else
  {
    cat("The", method, "is not ready yet!  Working on it!\n")
  }
}



#' stock recruitment function
#'
#'
#' @param data data.frame containing stock recruitment data
#' @param nsamp Number of samples
#' @param models A character vector
#' @param ... Additional arguements
#' @return log recruitment according to model
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
eqsr_Buckland <- function(data, nsamp = 5000, models = c("ricker","segreg","bevholt"), ...)
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
  pred <- t(sapply(seq(nsamp), function(j) exp(get(fit $ model[j], , pos = "package:msy", mode = "function") (fit[j,], sort(data $ ssb))) ))
  
  
  #--------------------------------------------------------
  # get best fit for each model
  #--------------------------------------------------------
  fits <- 
    do.call(rbind,
            lapply(models, 
                   function(mod) 
                     with(nlminb(initial(mod, data), nllik, data = data, model = mod, logpar = TRUE), 
                          data.frame(a = exp(par[1]), b = exp(par[2]), cv = exp(par[3]), model = mod))))
  
  tmp <- plyr::ddply(fit,c("model"), plyr::summarise, n=length(model))
  tmp$prop <- tmp$n/sum(tmp$n)
  fits <- plyr::join(fits,tmp,by="model")
  
  dimnames(pred) <- list(model=fit$model,ssb=data$ssb)
  
  return(list(sr.sto = fit, sr.det = fits, pRec = pred))
}




#' plot simulated predictive distribution of recruitment
#'
#'
#' @param fit an fitted MCMC returned from \code{eqsr_fit}
#' @param n Number of random recruitment draws to plot
#' @param x.mult max.y (ssb) as a multiplier of maximum observed ssb
#' @param y.mult max.x (rec) as a multiplier of maismum observed rec
#' @param ggPlot Flag, if FALSE (default) plot base graphics, if true
#' do a ggplot
#' @return NULL produces a plot
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
eqsr_plot <- function (fit, n = 5000, x.mult=1.1, y.mult=1.4, ggPlot=FALSE, Scale=1) 
{
  #x.mult <- 1.1
  #y.mult <- 1.4
  # dummy stuff
  Ftarget <- p05 <- p95 <- p50 <- variable <- value <- year <- Model <- NULL
  
  modset <- fit$sr.sto
  data <- fit$rby[,1:3]
  
  minSSB <- min(data$ssb, max(data$ssb)*0.0125)
  maxSSB <- max(data$ssb)*x.mult
  maxrec <- max(data$rec* y.mult)
  ##############################################################################
  # very strange way to do things
  out <-
    do.call(rbind, lapply(sample(1:nrow(modset), 500), 
                          function(i)
                          {
                            fssb <- runif(500, minSSB, maxSSB)
                            FUN <-  match.fun(modset $ model[i])
                            frec <- exp( FUN(modset[i,], fssb) + rnorm(500, sd = modset $ cv[i]) )
                            srModel <- modset$model[i]
                            #points(fssb, frec, pch = 20, col = paste0(grey(0), "05"), cex = 0.0625)
                            data.frame(ssb = fssb, rec = frec, model=srModel)
                          }))
  # group the ssbs into 10 bins
  out $ grp <- with(out, floor(10 * (ssb - min(ssb)) / (max(ssb) - min(ssb) + 0.001)))
  # find the midvalue of ssb within each group
  out $ mid.grp <- with(out, (grp + 0.5) / 10 * (max(ssb) - min(ssb)) + min(ssb))
  tmp <- fit$sr.det
  tmp$Model <- paste(tmp$model,tmp$prop)
  out <- join(out,tmp[,c("model","Model")],by="model")
  # calculate the recruitment median and 5th and 95th percentile within each
  # ssb group and then plot the distribution
  summ <- with(out, 
               t(simplify2array( tapply(rec, grp, quantile, c(0.5, .05, .95)) )))
  mid.grp <- sort(unique(out $ mid.grp))
  
  # For ggplot2
  Percentiles <- data.frame(ssb=mid.grp,p50=summ[,1],p05=summ[,2],p95=summ[,3])
  # end of very strange things
  #############################################################################
  
  if(!ggPlot) {
    
    plot(data$ssb, data$rec, xlim = c(0, maxSSB), ylim = c(0, maxrec), type = "n", 
         xlab = "SSB ('000 t)", ylab="Recruits", main = paste("Predictive distribution of recruitment\nfor", fit $id.sr))
    points(out$ssb[1:n], out$rec[1:n], pch = 20, col = paste0(grey(0), "05"), cex = 1)
    lines(mid.grp, summ[,1], col = 7, lwd = 3)
    lines(mid.grp, summ[,2], col = 4, lwd = 3)
    lines(mid.grp, summ[,3], col = 4, lwd = 3)
    # plot the best fit for each model as a line
    x <- fit $ sr.det[,1:4]
    y <- seq(1, round(maxSSB), length = 100)
    sapply(1:nrow(x), function(i) lines(y, exp(match.fun(as.character(x$model[i])) (x[i,], y)), col = "black", lwd = 2, lty = i))
    # plot the observation points
    lines(data$ssb, data$rec, col = "red")
    points(data$ssb, data$rec, pch = 19, col = "red", cex = 1.25)
    # plot the model weights on the graph
    for (i in 1:nrow(fit$sr.det)) {
      text(0.2*maxSSB,maxrec*(1-i/10),paste(fit$sr.det$model[i],round(fit$sr.det$prop[i],2)),cex=0.9)
    }
    
  } else { # ggplot
    
    #
    x <- fit$sr.det
    ssb <- seq(1,round(max(maxSSB)),length=100)
    z <- sapply(1:nrow(x), function(i) rec <- exp(match.fun(as.character(x$model[i])) (x[i,], ssb)))
    modelLines <- as.data.frame(cbind(ssb,z))
    names(modelLines) <- c("ssb",paste(x$model,x$prop))
    modelLines <- melt(modelLines,id.var="ssb",variable.name="Model",value.name="rec")
    
    out$ssb <- out$ssb/Scale
    out$rec <- out$rec/Scale
    out$mid.grp <- out$mid.grp/Scale
    Percentiles$ssb <- Percentiles$ssb/Scale
    Percentiles$p50 <- Percentiles$p50
    Percentiles$p05 <- Percentiles$p05
    Percentiles$p95 <- Percentiles$p95
    
    modelLines$ssb <- modelLines$ssb/Scale
    modelLines$rec <- modelLines$rec/Scale
    
    fit$rby$ssb <- fit$rby$ssb/Scale
    fit$rby$rec <- fit$rby$rec/Scale
    i <- sample(nrow(out),n)
    
    ggplot(out[i,]) + 
      theme_bw() +
      geom_point(aes(x=ssb,y=rec,colour=Model),size=1) +
      geom_line(data=Percentiles,aes(x=ssb,y=p05),colour="yellow") +
      geom_line(data=Percentiles,aes(x=ssb,y=p95),colour="yellow") +
      geom_line(data=Percentiles,aes(ssb,p50),col="yellow",lwd=2) +
      geom_line(data=modelLines,aes(ssb,rec,colour=Model),lwd=1) +
      coord_cartesian(ylim=c(0,quantile(out$rec[i],0.99))) +
      geom_path(data=fit$rby,aes(ssb,rec),col="black",linetype=2) +
      geom_text(data=fit$rby,aes(ssb,rec,label=substr(year,3,4)),size=4,col="black",angle=45) +
      theme(legend.position = c(0.20,0.85)) +
      labs(x="Spawning stock biomass",y="Recruitment",colour="Model")
    
  }
}