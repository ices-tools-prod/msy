#' Plot Simulated Predictive Distribution of Recruitment
#'
#'
#'
#' @param fit an fitted stock recruit model returned from \code{eqsr_fit}
#' @param n Number of random recruitment draws to plot
#' @param x.mult max value for the y axis (ssb) as a multiplier of maximum
#'               observed ssb
#' @param y.mult max value for the x axis (rec) as a multiplier of maismum
#'               observed rec
#' @param ggPlot Flag, if FALSE (default) plot using base graphics, if TRUE
#'               do a ggplot
#' @param Scale Numeric value for scaling varibles in plot.
#' @param filter.lag Integer (odd) value giving the size of the
#'                   window used (as a percentage of the full range) to average
#'                   the median and precentiles of the simulated
#'                   recruitment (default is 15).  Recruitment is calculated for
#'                   100 SSB points, so a max window of 99 is allowed, but is
#'                   pretty silly. A filter.lag of 1, gives no smoothing, but
#'                   is typically very noisy. All of this is purely for plotting
#'                   purposes.
#'
#' @return NULL produces a plot
#'
#' @seealso
#' \code{\link{eqsr_fit}} Fits several stock recruitment models to a data set
#' and calculates the proportion contribution of each model based on a bootstrap
#' model averaging procedure.
#'
#' @examples
#'
#' #' \dontrun{
#' data(icesStocks)
#' FIT <- eqsr_fit(icesStocks$saiNS,
#'                 nsamp = 1000,
#'                 models = c("Ricker", "Segreg"))
#'
#' eqsr_plot(FIT, n = 20000)
#'
#' # Scale argument only available for ggPlot = TRUE
#' eqsr_plot(FIT, n = 20000, ggPlot = TRUE, Scale = 1000)
#' }
#'
#' @export
eqsr_plot <- function (fit, n = 20000, x.mult = 1.1, y.mult = 1.4,
                       ggPlot = FALSE, Scale = 1,
                       filter.lag = 15)
{
  # set a range of values to zero
  Ftarget <- p05 <- p95 <- p50 <- variable <- value <- year <- Model <- rec <- 0

  # get the draws from the SR parameter simulations
  modset <- fit$sr.sto
  # get the full data set
  data <- fit$rby[,c("year", "rec", "ssb")]

  # set up ranges
  minSSB <- min(data$ssb, max(data$ssb)*0.0125)
  maxSSB <- max(data$ssb)*x.mult
  maxrec <- max(data$rec* y.mult)

  # evaluate recruitment at 100 ssb points
  ssb_eval <- seq(minSSB, maxSSB, length.out = 100)
  # reduce n to get n samples pairs
  n_mods <- floor(n / 100)

  # a function to sample from the predictive distribution
  # of recruitment given a bootstrap SR fit.  Each row in modset
  # has a (potrentially) different form, parameter estimates and residual CV.
  sample_rec <- function(i) {
    # what SR model are we simulting from:
    FUN <-  match.fun(modset $ model[i])
    # simulate from _predictive_ distribrution of recruitment
    exp( FUN(modset[i,], ssb_eval) + stats::rnorm(length(ssb_eval), sd = modset $ cv[i]) )
  }

  # (up)sample the model fits
  ids <- sample(1:nrow(modset), n_mods, replace = TRUE)
  rec_sim <- sapply(ids, sample_rec)

  # form into a big DF
  out <- data.frame(grp = rep(1:length(ssb_eval), n_mods),
                    mid.grp = rep(ssb_eval, n_mods),
                    ssb = jitter(rep(ssb_eval, n_mods), 2), # jitter for nices plotting
                    rec = c(rec_sim),
                    model = rep(modset[ids,"model"], each = length(ssb_eval)))

  tmp <- fit$sr.det
  tmp$Model <- paste(tmp$model, tmp$prop)
  out <- plyr::join(out,tmp[, c("model", "Model")], by = "model")

  # calculate the recruitment median and 5th and 95th percentile within each
  # ssb group and then plot the distribution
  summ <- t(apply(rec_sim, 1, stats::quantile, c(0.5, .05, .95)))

  # smooth for plotting
  f_lag <- rep(1/filter.lag, filter.lag)
  Percentiles <-
    data.frame(
      ssb = ssb_eval,
      p50 = as.vector(filter(summ[,1], f_lag, sides = 2)),
      p05 = as.vector(filter(summ[,2], f_lag, sides = 2)),
      p95 = as.vector(filter(summ[,3], f_lag, sides = 2)))
  Percentiles <- Percentiles[stats::complete.cases(Percentiles),]

  if(!ggPlot) {

    # set up plot
    plot(0, 0, type = "n",
         xlim = c(0, maxSSB), ylim = c(0, maxrec), las = 1,
         xlab = "SSB ('000 t)", ylab="Recruits",
         main = paste("Predictive distribution of recruitment\nfor", fit$id.sr))

    points(out$ssb, out$rec,
           pch = 20, cex = 1,
           col = grDevices::grey(0, alpha = 0.02))

    lines(p50 ~ ssb, col = 7, lwd = 3, data = Percentiles)
    lines(p05 ~ ssb, col = 4, lwd = 3, data = Percentiles)
    lines(p95 ~ ssb, col = 4, lwd = 3, data = Percentiles)

    # plot the best fit for each model as a line
    x <- fit$sr.det[,1:4] # a            b        cv  model
    y <- seq(1, round(maxSSB), length = 100)
    for (i in 1:nrow(x)) {
      lines(y, exp(match.fun(as.character(x$model[i])) (x[i,], y)), col = "black", lwd = 2, lty = i)
    }

    # plot the observation points
    lines(data$ssb, data$rec, col = "red")
    points(data$ssb, data$rec, pch = 19, col = "red", cex = 1.25)
    # plot the model weights on the graph
    for (i in 1:nrow(fit$sr.det)) {
      text(0.2 * maxSSB,
           maxrec*(1 - i/10),
           paste(fit$sr.det$model[i], round(fit$sr.det$prop[i], 2)),
           cex = 0.9)
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
      ggplot2::geom_point(ggplot2::aes(x=ssb,y=rec,colour=Model),size=1, alpha = 0.2) +
      ggplot2::geom_line(data=Percentiles,ggplot2::aes(x=ssb,y=p05),colour="blue", lwd = 1.5) +
      ggplot2::geom_line(data=Percentiles,ggplot2::aes(x=ssb,y=p95),colour="blue", lwd = 1.5) +
      ggplot2::geom_line(data=Percentiles,ggplot2::aes(ssb,p50),col="yellow", lwd = 1.5) +
      ggplot2::geom_line(data=modelLines,ggplot2::aes(ssb,rec,colour=Model),lwd=1.5) +
      ggplot2::coord_cartesian(ylim=c(0, stats::quantile(out$rec[i],0.99))) +
      ggplot2::geom_path(data=fit$rby, ggplot2::aes(ssb, rec), col="black",linetype=2, lwd = 1) +
      ggplot2::geom_text(data=fit$rby, ggplot2::aes(ssb, rec, label = substr(year,3,4)),size=4,col="black",angle=45) +
      ggplot2::theme(legend.position = c(0.20, 0.85)) +
      ggplot2::labs(x = "Spawning stock biomass",
                    y = "Recruitment",
                    colour="Model")

  }
}
