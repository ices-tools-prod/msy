eqsim_range <-function (sim, interval=0.95)
  {

  data.95 <- sim$rbp
  x.95 <- data.95[data.95$variable == "Landings",]$Ftarget
  y.95 <- data.95[data.95$variable == "Landings",]$Mean
  x.95 <- x.95[2:length(x.95)]
  y.95 <- y.95[2:length(y.95)]

  # Plot curve with 95% line
  #windows(width = 10, height = 7)
  #par(mfrow = c(1,1), mar = c(5,4,2,1), mgp = c(3,1,0))
  #plot(x.95, y.95, ylim = c(0, max(y.95, na.rm = TRUE)),
  #     xlab = "Total catch F", ylab = "Mean landings")
  yield.p95 <- interval * max(y.95, na.rm = TRUE)
  #abline(h = yield.p95, col = "blue", lty = 1)

  # Fit loess smoother to curve
  x.lm <- stats::loess(y.95 ~ x.95, span = 0.2)
  lm.pred <- data.frame(x = seq(min(x.95), max(x.95), length = 1000),
                        y = rep(NA, 1000))
  lm.pred$y <- stats::predict(x.lm, newdata = lm.pred$x)
  #lines(lm.pred$x, lm.pred$y, lty = 1, col = "red")
  #points(x = sim$Refs["lanF","meanMSY"],
  #      y = predict(x.lm, newdata = sim$Refs["lanF","meanMSY"]),
  #       pch = 16, col = "blue")

  # Limit fitted curve to values greater than the 95% cutoff
  lm.pred.95 <- lm.pred[lm.pred$y >= yield.p95,]
  fmsy.lower <- min(lm.pred.95$x)
  fmsy.upper <- max(lm.pred.95$x)
  #abline(v = c(fmsy.lower, fmsy.upper), lty = 8, col = "blue")
  #abline(v = sim$Refs["lanF","meanMSY"], lty = 1, col = "blue")
  #legend(x = "bottomright", bty = "n", cex = 1.0,
  #       title = "F(msy)", title.col = "blue",
  #       legend = c(paste0("lower = ", round(fmsy.lower,3)),
  #                  paste0("mean = ", round(sim$Refs["lanF","meanMSY"],3)),
  #                  paste0("upper = ", round(fmsy.upper,3))))

  fmsy.lower.mean <- fmsy.lower
  fmsy.upper.mean <- fmsy.upper
  landings.lower.mean <- lm.pred.95[lm.pred.95$x == fmsy.lower.mean,]$y
  landings.upper.mean <- lm.pred.95[lm.pred.95$x == fmsy.upper.mean,]$y

  # Repeat for 95% of yield at F(05):
  f05 <- sim$Refs["catF","F05"]
  yield.f05 <- stats::predict(x.lm, newdata = f05)
  #points(f05, yield.f05, pch = 16, col = "green")
  yield.f05.95 <- interval * yield.f05
  #abline(h = yield.f05.95, col = "green")
  lm.pred.f05.95 <- lm.pred[lm.pred$y >= yield.f05.95,]
  f05.lower <- min(lm.pred.f05.95$x)
  f05.upper <- max(lm.pred.f05.95$x)
  #abline(v = c(f05.lower,f05.upper), lty = 8, col = "green")
  #abline(v = f05, lty = 1, col = "green")
  #legend(x = "right", bty = "n", cex = 1.0,
  #       title = "F(5%)", title.col = "green",
  #       legend = c(paste0("lower = ", round(f05.lower,3)),
  #                  paste0("estimate = ", round(f05,3)),
  #                  paste0("upper = ", round(f05.upper,3))))

  ################################################
  # Extract yield data (landings) - median version

  data.95 <- sim$rbp
  x.95 <- data.95[data.95$variable == "Landings",]$Ftarget
  y.95 <- data.95[data.95$variable == "Landings",]$p50

  # Plot curve with 95% line
  #windows(width = 10, height = 7)
  #par(mfrow = c(1,1), mar = c(5,4,2,1), mgp = c(3,1,0))
  #plot(x.95, y.95, ylim = c(0, max(y.95, na.rm = TRUE)),
  #     xlab = "Total catch F", ylab = "Median landings")
  yield.p95 <- interval * max(y.95, na.rm = TRUE)
  #abline(h = yield.p95, col = "blue", lty = 1)

  # Fit loess smoother to curve
  x.lm <- stats::loess(y.95 ~ x.95, span = 0.2)
  lm.pred <- data.frame(x = seq(min(x.95), max(x.95), length = 1000),
                        y = rep(NA, 1000))
  lm.pred$y <- stats::predict(x.lm, newdata = lm.pred$x)
  #lines(lm.pred$x, lm.pred$y, lty = 1, col = "red")

  # Find maximum of fitted curve - this will be the new median (F(msy)
  Fmsymed <- lm.pred[which.max(lm.pred$y),]$x
  Fmsymed.landings <- lm.pred[which.max(lm.pred$y),]$y

  # Overwrite Refs table
  sim$Refs2 <- sim$Refs
  sim$Refs2[,"medianMSY"] <- NA
  sim$Refs2["lanF","medianMSY"] <- Fmsymed
  sim$Refs2["landings","medianMSY"] <- Fmsymed.landings

  # Add maximum of medians to plot
  #points(x = sim$Refs["lanF","medianMSY"],
  #       y = predict(x.lm, newdata = sim$Refs["lanF","medianMSY"]),
  #       pch = 16, col = "blue")

  # Limit fitted curve to values greater than the 95% cutoff
  lm.pred.95 <- lm.pred[lm.pred$y >= yield.p95,]
  fmsy.lower <- min(lm.pred.95$x)
  fmsy.upper <- max(lm.pred.95$x)
  #abline(v = c(fmsy.lower, fmsy.upper), lty = 8, col = "blue")
  #abline(v = sim$Refs["lanF","medianMSY"], lty = 1, col = "blue")
  #legend(x = "bottomright", bty = "n", cex = 1.0,
  #       title = "F(msy)", title.col = "blue",
  #       legend = c(paste0("lower = ", round(fmsy.lower,3)),
  #                  paste0("median = ", round(sim$Refs["lanF","medianMSY"],3)),
  #                  paste0("upper = ", round(fmsy.upper,3))))

  fmsy.lower.median <- fmsy.lower
  fmsy.upper.median <- fmsy.upper
  landings.lower.median <- lm.pred.95[lm.pred.95$x == fmsy.lower.median,]$y
  landings.upper.median <- lm.pred.95[lm.pred.95$x == fmsy.upper.median,]$y

  # Repeat for 95% of yield at F(05):
  f05 <- sim$Refs["catF","F05"]
  yield.f05 <- stats::predict(x.lm, newdata = f05)
  #points(f05, yield.f05, pch = 16, col = "green")
  yield.f05.95 <- interval * yield.f05
  #abline(h = yield.f05.95, col = "green")
  lm.pred.f05.95 <- lm.pred[lm.pred$y >= yield.f05.95,]
  f05.lower <- min(lm.pred.f05.95$x)
  f05.upper <- max(lm.pred.f05.95$x)
  #abline(v = c(f05.lower,f05.upper), lty = 8, col = "green")
  #abline(v = f05, lty = 1, col = "green")
  #legend(x = "right", bty = "n", cex = 1.0,
  #       title = "F(5%)", title.col = "green",
  #       legend = c(paste0("lower = ", round(f05.lower,3)),
  #                  paste0("estimate = ", round(f05,3)),
  #                  paste0("upper = ", round(f05.upper,3))))

  # Estimate implied SSB for each F output

  x.95 <- data.95[data.95$variable == "Spawning stock biomass",]$Ftarget
  b.95 <- data.95[data.95$variable == "Spawning stock biomass",]$p50

  # Plot curve with 95% line
  #windows(width = 10, height = 7)
  #par(mfrow = c(1,1), mar = c(5,4,2,1), mgp = c(3,1,0))
  #plot(x.95, b.95, ylim = c(0, max(b.95, na.rm = TRUE)),
  #     xlab = "Total catch F", ylab = "Median SSB")

  # Fit loess smoother to curve
  b.lm <- stats::loess(b.95 ~ x.95, span = 0.2)
  b.lm.pred <- data.frame(x = seq(min(x.95), max(x.95), length = 1000),
                          y = rep(NA, 1000))
  b.lm.pred$y <- stats::predict(b.lm, newdata = b.lm.pred$x)
  #lines(b.lm.pred$x, b.lm.pred$y, lty = 1, col = "red")

  # Estimate SSB for median F(msy) and range
  b.msymed <- stats::predict(b.lm, newdata = Fmsymed)
  b.medlower <- stats::predict(b.lm, newdata = fmsy.lower.median)
  b.medupper <- stats::predict(b.lm, newdata = fmsy.upper.median)
  #abline(v = c(fmsy.lower.median, Fmsymed, fmsy.upper.median), col = "blue", lty = c(8,1,8))
  #points(x = c(fmsy.lower.median, Fmsymed, fmsy.upper.median),
  #       y = c(b.medlower, b.msymed, b.medupper), col = "blue", pch = 16)
  #legend(x = "topright", bty = "n", cex = 1.0,
  #       title = "F(msy)", title.col = "blue",
  #       legend = c(paste0("lower = ", round(b.medlower,0)),
  #                  paste0("median = ", round(b.msymed,0)),
  #                  paste0("upper = ", round(b.medupper,0))))

  # Update summary table with John's format

  sim$Refs2 <- sim$Refs2[,!(colnames(sim$Refs) %in% c("FCrash05","FCrash50"))]
  sim$Refs2 <- cbind(sim$Refs2, Medlower = rep(NA,6), Meanlower = rep(NA,6),
                       Medupper = rep(NA,6), Meanupper = rep(NA,6))

  NAforNumeric0 <- function(x) if (length(x)) x else NA

  sim$Refs2["lanF","Medlower"] <- NAforNumeric0(fmsy.lower.median)
  sim$Refs2["lanF","Medupper"] <- NAforNumeric0(fmsy.upper.median)
  sim$Refs2["lanF","Meanlower"] <- NAforNumeric0(fmsy.lower.mean)
  sim$Refs2["lanF","Meanupper"] <- NAforNumeric0(fmsy.upper.mean)

  sim$Refs2["landings","Medlower"] <- NAforNumeric0(landings.lower.median)
  sim$Refs2["landings","Medupper"] <- NAforNumeric0(landings.upper.median)
  sim$Refs2["landings","Meanlower"] <- NAforNumeric0(landings.lower.mean)
  sim$Refs2["landings","Meanupper"] <- NAforNumeric0(landings.upper.mean)

  sim$Refs2["lanB","medianMSY"] <- NAforNumeric0(b.msymed)
  sim$Refs2["lanB","Medlower"] <- NAforNumeric0(b.medlower)
  sim$Refs2["lanB","Medupper"] <- NAforNumeric0(b.medupper)

  # Reference point estimates
  #cat("\nReference point estimates:\n")
  return (sim)
}
