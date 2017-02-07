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
      fssb <- stats::runif(500, minSSB, maxSSB)
      FUN <-  match.fun(modset $ model[i])
      frec <- exp( FUN(modset[i,], fssb) + stats::rnorm(500, sd = modset $ cv[i]) )
      #points(fssb, frec, pch = 20, col = paste0(grey(0), "05"), cex = 0.0625)

      data.frame(ssb = fssb, rec = frec)
    }))
  points(out$ssb[1:n], out$rec[1:n], pch = 20, col = paste0(grDevices::grey(0), "05"), cex = 1)
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
