#' @title Summarise predictive distribution of recruitment
#' 
#' @description The input is a list object that is returned from \code{eqsr_fit}.
#' The return object is a list of data.frames (see below) that can be used to
#' get e.g. a summary plot using the \code{eqsr_plot2} function.
#' 
#' The function is based on the same algorithm already in SRplot, the difference
#' is that
#' here the objects used in the plot are returned to the user. The philosophy is
#' to separate analysis from plot. Output hopefully facilitates user in
#' generating alternative plot to that provided by the \code{SRplot} function.
#' 
#' @export
#' 
#' @param fit A list containing fitted MCMC. Returned from \code{eqsr_fit}
#' @return A list containing the following data.frames:
#' \itemize{
#' \item data The original ssb-rec estimations
#' \item Parameters The estimated alpha (a), beta (b) and cv estimates.
#' \item Models The best parameter fit of each stock recruitment model and the
#' number (n) and proportion (prop) of stochastic recruitment from each model.
#' \item srBestFit ssb and rec values for each best fit model, facilitates plotting a
#' lines for each model.
#' \item recStochastic Contains stochastic recruiment for each model.
#' \item recPercentiles The 0.50, 0.05 and 0.95 percentiles of the stochastic
#' recruitment, grouped along 10 ssb-groups.
#' }
#' 
#' @author Einar Hjorleifsson \email{einar.hjorleifsson@@gmail.com}
#' 

eqsr_summary <- function (fit) 
  {
  
  Parameters <- fit $ fit
  data <- fit $ data
  
  minSSB <- min(data$ssb, max(data$ssb)*0.0125)
  maxSSB <- max(data$ssb)*1.1
  maxrec <- max(data$rec*1.5)
  
  # Contribution of the different stock recruitment models
  Models <- melt(tapply(fit$fit$model,fit$fit$model,length))
  names(Models) <- c("model","n")
  Models$model <- as.character(Models$model)
  Models$prop <- round(Models$n/sum(Models$n),2)
  Models <- join(Models,fit$fits,by="model")
  
  
  # 1) Take 500 samples from nrow(Parameters) - default value is 5000
  # 2) Make 500 uniformly distributed SSBs (fssb)
  # 3) Generate recruitment value
  # A total of 250000 (500*500) recruitment values are thus generated. The
  # weight that each model gets when generating the recruits is equivalent
  # to the number of occurrences in the fit$fit object, i.e table(fit$fit$models)
  
  recStochastic <-
    do.call(rbind, lapply(sample(1:nrow(Parameters), 500), 
                          function(i)
                          {
                            fssb <- runif(500, minSSB, maxSSB)
                            FUN <-  match.fun(Parameters $ model[i])
                            frec <- exp( FUN(Parameters[i,], fssb) + rnorm(500, sd = Parameters $ cv[i]) )
                            srModel <- Parameters$model[i]
                            data.frame(ssb = fssb, rec = frec, model=srModel)
                          }))
  
  # Calculate the recruitement percentile distribution over a certain ssb-range
  # group the ssbs into 10 bins
  recStochastic $ ssb_grp <- with(recStochastic, floor(10 * (ssb - min(ssb)) / (max(ssb) - min(ssb) + 0.001)))
  # find the midvalue of ssb within each group
  recStochastic $ ssb_grp.mid <- with(recStochastic, (ssb_grp + 0.5) / 10 * (max(ssb) - min(ssb)) + min(ssb))
  # calculate the rec median and 5th and 95th percentile within each
  # ssb group and then plot the distribution
  recPercentiles <- with(recStochastic, 
               t(simplify2array( tapply(rec, ssb_grp, quantile, c(0.5, .05, .95)) )))
  mid.grp <- sort(unique(recStochastic $ ssb_grp.mid))
  recPercentiles <- data.frame(ssb=mid.grp,p50=recPercentiles[,1],p05=recPercentiles[,2],p95=recPercentiles[,3])
  
  # Line for best model fit
  ssb <- seq(1, round(max(data$ssb)), length = 100)
  x <- fit $ fits
  rec <- sapply(1:nrow(x), function(i) rec <- exp(match.fun(as.character(x$model[i])) (x[i,], ssb)))
  srBestFit <- as.data.frame(cbind(ssb,rec))
  names(srBestFit) <- c("ssb",x$model)
  srBestFit <- melt(srBestFit,id.vars="ssb")
  names(srBestFit) <- c("ssb","model","rec")

  return(list(data=data,Parameters=Parameters,Models=Models,srBestFit=srBestFit,
                recStochastic=recStochastic,recPercentiles=recPercentiles))
}

#' @title Summarise predictive distribution of recruitment
#'
#'
#' @param x A list from eqsr_summary
#' @param n An integer specifying the number of stochastic recruitment points
#' to include in the plot
#' @param ggPlot A flag, if TRUE (default) returns a ggplot. Base plot (FALSE)
#' not yet implemented.
#' @return A plot
#' @author Einar Hjorleifsson
#' @export
#' 
eqsr_plot2 <- function (x,n=5000,ggPlot=TRUE) 
  {
  if(ggPlot) {
    ggplot(x$recStochastic[1:n,]) + 
      theme_bw() +
      geom_line(data=x$recPercentiles,aes(x=ssb,y=p05),colour="yellow") +
      geom_line(data=x$recPercentiles,aes(x=ssb,y=p95),colour="yellow") +
      geom_line(data=x$recPercentiles,aes(ssb,p50),col="yellow",lwd=2) +
      geom_point(aes(ssb,rec,colour=model),size=1) +
      geom_line(data=x$srBestFit,aes(ssb,rec,colour=model),lwd=1) +
      #scale_colour_brewer(palette="Set1") +
      scale_colour_manual(values=c("bevholt"="#F8766D","ricker"="#00BA38","segreg"="#619CFF"),
                          labels=paste(x$Models$model,x$Models$prop)) +
      coord_cartesian(ylim=c(0,quantile(x$recStochastic$rec[1:n],0.99))) +
      geom_path(data=x$data,aes(ssb,rec),col="grey") +
      geom_text(data=x$data,aes(ssb,rec,label=substr(year,3,4)),angle=45,size=3) +
      theme(legend.position = c(0.15,0.80)) +
      labs(x="Spawning stock biomass",y="Recruitment",colour="Model")
  } else {
    stop("Base graph not implemented yet")
  }
}



