#' plot simulated predictive distribution of recruitment
#'
#'
#' @param fit an fitted MCMC returned from \code{eqsr_fit}
#' @param n Number of random recruitment draws to plot
#' @param x.mult max.y (ssb) as a multiplier of maximum observed ssb
#' @param y.mult max.x (rec) as a multiplier of maismum observed rec
#' @param ggPlot Flag, if FALSE (default) plot base graphics, if true
#' do a ggplot
#' @param Scale Numeric value for scaling varibles in plot.
#' @return NULL produces a plot
#' @author Colin Millar \email{colin.millar@@ices.dk}
#' @export
eqsr_plot <- function (fit, n = 5000, x.mult=1.1, y.mult=1.4, ggPlot=FALSE, Scale=1)
{
  #x.mult <- 1.1
  #y.mult <- 1.4
  ## dummy stuff
  Ftarget <- p05 <- p95 <- p50 <- variable <- value <- year <- Model <- rec <- 0

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
                            fssb <- stats::runif(500, minSSB, maxSSB)
                            FUN <-  match.fun(modset $ model[i])
                            frec <- exp( FUN(modset[i,], fssb) + stats::rnorm(500, sd = modset $ cv[i]) )
                            srModel <- modset$model[i]
                            #points(fssb, frec, pch = 20, col = paste0(grDevices::grey(0), "05"), cex = 0.0625)
                            data.frame(ssb = fssb, rec = frec, model=srModel)
                          }))
  # group the ssbs into 10 bins
  out $ grp <- with(out, floor(10 * (ssb - min(ssb)) / (max(ssb) - min(ssb) + 0.001)))
  # find the midvalue of ssb within each group
  out $ mid.grp <- with(out, (grp + 0.5) / 10 * (max(ssb) - min(ssb)) + min(ssb))
  tmp <- fit$sr.det
  tmp$Model <- paste(tmp$model,tmp$prop)
  out <- plyr::join(out,tmp[,c("model","Model")],by="model")
  # calculate the recruitment median and 5th and 95th percentile within each
  # ssb group and then plot the distribution
  summ <- with(out,
               t(simplify2array( tapply(rec, grp, stats::quantile, c(0.5, .05, .95)) )))
  mid.grp <- sort(unique(out $ mid.grp))

  # For ggplot2
  Percentiles <- data.frame(ssb=mid.grp,p50=summ[,1],p05=summ[,2],p95=summ[,3])
  # end of very strange things
  #############################################################################

  if(!ggPlot) {

    plot(data$ssb, data$rec, xlim = c(0, maxSSB), ylim = c(0, maxrec), type = "n",
         xlab = "SSB ('000 t)", ylab="Recruits", main = paste("Predictive distribution of recruitment\nfor", fit $id.sr))
    points(out$ssb[1:n], out$rec[1:n], pch = 20, col = paste0(grDevices::grey(0), "05"), cex = 1)
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
    modelLines <- reshape2::melt(modelLines,id.var="ssb",variable.name="Model",value.name="rec")

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

    ggplot2::ggplot(out[i,]) +
      ggplot2::theme_bw() +
      ggplot2::geom_point(ggplot2::aes(x=ssb,y=rec,colour=Model),size=1) +
      ggplot2::geom_line(data=Percentiles,ggplot2::aes(x=ssb,y=p05),colour="yellow") +
      ggplot2::geom_line(data=Percentiles,ggplot2::aes(x=ssb,y=p95),colour="yellow") +
      ggplot2::geom_line(data=Percentiles,ggplot2::aes(ssb,p50),col="yellow",lwd=2) +
      ggplot2::geom_line(data=modelLines,ggplot2::aes(ssb,rec,colour=Model),lwd=1) +
      ggplot2::coord_cartesian(ylim=c(0, stats::quantile(out$rec[i],0.99))) +
      ggplot2::geom_path(data=fit$rby,ggplot2::aes(ssb,rec),col="black",linetype=2) +
      ggplot2::geom_text(data=fit$rby,ggplot2::aes(ssb,rec,label=substr(year,3,4)),size=4,col="black",angle=45) +
      ggplot2::theme(legend.position = c(0.20,0.85)) +
      ggplot2::labs(x="Spawning stock biomass",y="Recruitment",colour="Model")

  }
}
