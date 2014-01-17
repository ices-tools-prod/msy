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

#' @title Plot the results from eqsim
#'
#' @description This is a modification of the \code{Eqplot}-function where the mean
#' catch can be turned on-off by users choise with some additional deletion
#' of vertical lines that were deemed unessessary at the wkmsyref2 meeting.
#' 
#' @param sim An object returned from the function EqSim 
#' @param fit An object returned from the function fitModels
#' @param Blim Value for the Blim
#' @param Bpa Value for the Bpa
#' @param ymax vector of three values, dictating maximum y-value on the recruitment
#' plot, the maximum value on the ssb-plot and the maximum value on the catch plot (in
#' that order)
#' @param xmax Value for the maximum F to plot, if missing will use the whole
#' F-range simulated.
#' @param plot Flag, if TRUE (default) provide a plot as part of the output
#' @param plotMeanCatch Flag, if FALSE (default) does not provide values for the
#' mean catch.
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
eqsim_plot <- 
  function (sim, fit, Blim, Bpa = 1.4 * Blim, ymax = c(NA,NA,NA), xmax, plot = TRUE,
            plotMeanCatch = FALSE)
  {
    
    stk <- fit $ stk
    
    Nmod <- dim(sim $ ssbsa)[3]
    Nyrs <- dim(sim $ ssbsa)[2]
    
    Fscan <- sim $ Fscan
    
    
    catm <- apply(sim $ catsa, 1, mean)
    maxcatm <- which.max(catm)
    catsam <- apply(sim $ catsa, c(1,3), mean)
    maxpf <- apply(catsam, 2, which.max)
    fmsy <- Fscan[maxpf]
    
    msym <- mean(fmsy)
    vcum <- median(fmsy)
    fmsy.dens <- density(fmsy)
    vmode <- fmsy.dens $ x[which.max(fmsy.dens $ y)]
    
    pssb1 <- apply(sim $ ssbsa > Blim, 1, mean)
    pssb2 <- apply(sim $ ssbsa > Bpa, 1, mean)
    
    pp1 <- max(which(pssb1>.95))
    grad <- diff(Fscan[pp1 + 0:1]) / diff(pssb1[pp1 + 0:1])
    flim <- Fscan[pp1] + grad * (0.95 - pssb1[pp1])  # linear interpolation i think..
    
    
    maint <- fit $ stknam  
    
    rec <- fit $ data $ rec
    ssb <- fit $ data $ ssb
    
    Catchs <- catch(stk)[, 1:length(ssb), drop = TRUE]
    FbarO <- fbar(stk)[, 1:length(ssb), drop = TRUE]
    
    recs <- sim $ recs
    ssbs <- sim $ ssbs
    cats <- sim $ cats
    ssbsa <- sim $ ssbsa
    
    NF <- length(Fscan)
    pp1 <- max(which(pssb1>.50))
    grad <- diff(Fscan[pp1 + 0:1]) / diff(pssb1[pp1 + 0:1])
    flim50 <- Fscan[pp1]+grad*(0.5-pssb1[pp1]) # linear interpolation i think..
    
    pp1 <- max(which(pssb1>.90))
    grad <- diff(Fscan[pp1 + 0:1]) / diff(pssb1[pp1 + 0:1])
    flim10 <- Fscan[pp1]+grad*(0.9-pssb1[pp1]) # linear interpolation i think..
    
    maxcatm <- which.max(catm)
    
    
    if (plot) {
      op <- par(mfrow = c(2, 2), mar = c(2.5, 4, 1.5, 1), oma = c(0, 0, 0, 0), 
                cex.axis = 0.75, tcl = 0.25, mgp = c(0, 0.25, 0), las = 1)
      
      # recruits versus Fbar
      if(missing(xmax)) xmax <- max(Fscan)
      y.max <- if (!is.na(ymax[1])) ymax[1] else max(1.25*rec)
      
      plot(Fscan, recs[7,], type = "l", lty = 4, 
           ylim = c(0, y.max), xlim = c(0, xmax), ylab = "", xlab = "")
      title(ylab = "Recruitment", xlab = "F bar", 
            cex.lab = 0.75, line = 2.5, cex.main = 0.75)
      mtext(text = paste(maint," a) Recruits"), cex = 0.75, side = 3, line = 0.5)
      lines(Fscan, recs[6,], lty = 3)
      lines(Fscan, recs[5,], lty = 2)
      lines(Fscan, recs[4,], lty = 1)
      lines(Fscan, recs[3,], lty = 2)
      lines(Fscan, recs[2,], lty = 3)
      lines(Fscan, recs[1,], lty = 4)
      points(FbarO, rec, pch = 21, cex = .75, bg = 1)
      lines(c(flim, flim), c(0, y.max), col = "red")
      text(x = flim, y = 0, "Flim", cex = 0.7, col="red")
      lines(rep(vcum,2), c(0,y.max), lty = 1, col = "blue")
      text(x = vcum, y = max(1.25*rec) * 0.9, "Fmsy",cex=0.7,col="blue")
      
      # recruits versus SSB
      y.max <- if (!is.na(ymax[2])) ymax[2] else max(1.25*ssb)
      plot(Fscan, ssbs[7,], type = "l", lty = 4, 
           ylim = c(0, y.max), xlim = c(0, xmax), ylab = "", xlab = "")
      title(ylab = "SSB", xlab = "F bar", 
            cex.lab = 0.75, line = 2.5, cex.main = 0.75)
      mtext(text = "b) Spawning Stock Biomass", cex = 0.75, side = 3, line = 0.5)
      lines(Fscan, ssbs[6,], lty=3)
      lines(Fscan, ssbs[5,], lty=2)
      lines(Fscan, ssbs[4,], lty=1)
      lines(Fscan, ssbs[3,], lty=2)
      lines(Fscan, ssbs[2,], lty=3)
      lines(Fscan, ssbs[1,], lty=4)
      lines(c(0,xmax), c(Blim, Blim),col="red")
      text(x = 0.1, y = Blim * 1.1, "Blim", cex = 0.7,col="red")
      points(FbarO, ssb, pch = 21, cex = .75, bg = 1)
      lines(c(flim, flim), c(0, Blim), col = "red")
      text(x = flim, y = 0, "Flim", cex = 0.7, col="red")
      lines(rep(vcum,2), c(0,y.max), lty = 1, col = "blue")
      text(x = vcum, y  = max(1.25*ssb) * 0.9, "Fmsy",cex=0.7,col="blue")
    }
    FCrash5  <- Fscan[which.max(cats[2,]):NF][ which(cats[2, which.max(cats[2,]):NF] < 0.05*max(cats[2,]) )[1] ]
    
    FCrash50 <- Fscan[which.max(cats[4,]):NF][ which(cats[4, which.max(cats[4,]):NF] < 0.05*max(cats[4,]) )[1] ]
    
    if (plot) {
      # catch versus Fbar
      y.max <- if (!is.na(ymax[3])) ymax[3] else max(1.25*Catchs)
      plot(Fscan, cats[7,], type = "l", lty = 4,
           ylim = c(0, y.max), xlim = c(0, xmax), ylab = "", xlab = "")
      title(ylab = "Catch", xlab = "F bar", 
            cex.lab = 0.75, line = 2.5, cex.main = 0.75)
      mtext(text = "c) Catch", cex = 0.75, side = 3, line = 0.5)
      #points(fbarsa[ssbsa<BLIM], catsa[ssbsa<BLIM],pch=21,col=6,bg=6,cex=0.125)
      #points(fbarsa[ssbsa>BLIM]+.0075,catsa[ssbsa>BLIM],pch=21,col=7,bg=7,cex=0.125)
      lines(Fscan, cats[7,], lty=4)
      lines(Fscan, cats[6,], lty=3)
      lines(Fscan, cats[5,], lty=2)
      lines(Fscan, cats[4,], lty=1)
      lines(Fscan, cats[3,], lty=2)
      lines(Fscan, cats[2,], lty=3)
      lines(Fscan, cats[1,], lty=4)
      points(FbarO, Catchs, pch = 21, cex = .75, bg = 1)
      lines(c(flim, flim), c(0, y.max), col = "red")
      lines(c(vcum,vcum),c(0,y.max),col="blue")
      #lines(c(FCrash5, FCrash5), c(0, y.max), col = 5)
      #lines(c(FCrash50, FCrash50), c(0, y.max), col = 5)
      if(plotMeanCatch) {
        lines(Fscan, catm, lty=1, col = "grey")
        lines(rep(Fscan[maxcatm], 2), c(0, y.max), lty = 4, col = "grey",lwd=1)
        text(x=Fscan[maxcatm],y=0,"mean Fmsy",cex=0.7,col="grey")
      }
      text(x = flim, y = 0, "Flim", cex = 0.7, col="red")
      lines(rep(vcum,2), c(0,y.max), lty = 1, col = "blue")
      text(x = vcum,  y = max(1.25*Catchs) * 0.9, "Fmsy",cex=0.7,col="blue")
      
      # F versus SSB probability profile
      plot(Fscan, 1-pssb1, type = "l", col="red",
           ylim = c(0,1), xlim = c(0,xmax), ylab = "", xlab = "")
      title(ylab = "Prob MSY, SSB<Bpa or Blim", xlab = "F bar", 
            cex.lab = 0.75, line = 2.5, cex.main = 0.75)
      mtext(text = "d) Prob MSY and Risk to SSB", cex = 0.75, side = 3, line = 0.5)
      lines(Fscan, 1 - pssb2, lty = 4, col="darkgreen")
      
      
      #text(x = max(Fscan[pssb2 > 0.5]) - .05, y = 0.5, "SSB<Bpa", cex = 0.75, col="darkgreen")
      text(x = 0.1, y = 0.8, "SSB<Bpa", cex = 0.75, col="darkgreen")
      text(x = 0.1, y = 0.7, "SSB<Blim", cex = 0.75, col="red")
      
      #lines(c(flim,flim), c(0,1), col = "red")
      lines(c(flim,flim), c(0,0.05), col = "red")
      lines(c(0,flim), c(0.05,0.05), lty = 2, col = "red")
      text(x = 0.1, y = 0.075, "5%", cex = 0.75, col = "red")
      #lines(c(flim10,flim10), c(0,1), col = "darkgreen")
      lines(c(flim10,flim10), c(0,0.10), col = "orange")
      lines(c(0,flim10), c(0.1,0.1), lty = 2, col = "orange")
      text(x = 0.05, y = 0.125, "10%", cex = 0.75, col = "orange")
      lines(fmsy.dens $ x, cumsum(fmsy.dens $ y * diff(fmsy.dens $ x)[1]), col = "blue",lwd=2)
      lines(c(0, vcum), c(0.5, 0.5), lty = 2, col = "blue")
      text(x = 0.1, y = 0.9, "p(Fmsy)", cex = 0.75, col = "blue")
      lines(rep(vcum,2), c(0,0.50), lty = 1, col = "blue")
      lines(c(0,vcum),c(0.5,0.5),lty=1,col="blue")
      text(x=vcum,y=0.05,"Fmsy",col="blue",cex = 0.75)
      lines(c(Fscan[maxcatm],Fscan[maxcatm]), c(0,1), col = "brown", lwd=1)
    }
    
    
    out <- c(Blim, Bpa)
    outF <- c(flim, flim10, flim50, vcum, Fscan[maxcatm], FCrash5, FCrash50)
    outB <- approx(Fscan, ssbs[4,], xout = outF) $ y
    outC <- approx(Fscan, cats[4,], xout = outF) $ y
    
    outTable <- rbind(outF, outB, outC)
    rownames(outTable) <- c("F","SSB","Catch") 
    colnames(outTable) <- c("Flim","Flim10","Flim50","MSY:median","Maxmeanland","FCrash5","FCrash50")
    
    list(Blim = Blim, Bpa = Bpa, Refs = outTable)
  }


