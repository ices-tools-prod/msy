
##############################################
# utility functions
##############################################



#' Get starting values for models
#'
#' a quick fix!!
#'
#' @param ...
#' @return vector of starting values
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
initial <- function(model, data)
{
  if (model == "segreg") {
    c(0,log(mean(data $ ssb)),0)
  } else {
    c(0,0,0)
  }
}

#' A function only available for R 2.15.1
#'
#' for back compatibility
#'
#' @param ...
#' @return result of paste with sep = ""
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
paste0 <- function(...) paste(..., sep = "")


#' stock recruitment function
#'
#'
#' @param ab the model parameters
#' @param ssb a vector of ssb
#' @return log recruitment according to model
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
loader <- function(p)
{
  if (p==0) cat("0%                       50%                     100%\n")
 str <- paste0(rep(c("\r[", "=", ">", " ", "]"), c(1, floor(p*50), 1, 50 - floor(p*50), 1)), collapse = "")  
 cat(str); flush.console()
 if (floor(p) == 1) cat("\n")
}


##############################################
# Main processing subroutines: 
#   1) fitModels
#   2) EqSim
##############################################




##### simulates the equilibrium results for a population
#' stock recruitment function
#'
#'
#' @param ab the model parameters
#' @param ssb a vector of ssb
#' @return log recruitment according to model
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
#' @examples
#' data(codEB)
#' fit <- fitModels(codEB, nsamp = 100) # a few samples for example
#' sim <- EqSim(fit, wt.years = c(2007, 2011), Fphi = 0.5, Fcv = 0.1)
#' hist(c(apply(sim $ ferr, c(1,3), function(x) acf(x, plot = FALSE) $ acf[2])), main = "histogram of estimates of ACF lag 1 in F errors")
#' arfits <- apply(sim $ ferr, c(1,3), function(x) ar(x)[c("order", "ar", "var.pred")])
#' summary(sapply(arfits, function(x) x $ order))
#' summary(sapply(arfits, function(x) x $ ar[1]))
#' summary(sapply(arfits, function(x) x $ var.pred))
EqSim <- function(fit, 
                  Nrun = 200, # number of years to run in total
                  wt.years = c(2007, 2011), # years sample weights, sel from
                  Fscan = seq(0, 1, len = 20), # F values to scan over
                  process.error = TRUE, # use predictive recruitment or mean recruitment?  (TRUE = predictive)
                  verbose = TRUE,
                  Btrigger = 0,
                  Fphi = 0,
                  Fcv = 0) 
{

  if (abs(Fphi) >= 1) stop("Fphi, the autocorelation parameter for log F should be between (-1, 1)")

  btyr1 <- wt.years[1]
  btyr2 <- wt.years[2] 
  flgsel <- 0
  flgmatwt <- 0
  keep <- min(Nrun, 50)

  SR <- fit $ fit
  data <- fit $ data
  stk <- fit $ stk

  # forecast settings (mean wt etc)
  stk.win <- window(stk, start = btyr1, end = btyr2)

  west <- matrix(stock.wt(stk.win), ncol = btyr2 - btyr1 + 1)
  weca <- matrix(catch.wt(stk.win), ncol = btyr2 - btyr1 + 1)
  sel <- matrix(harvest(stk.win), ncol = btyr2 - btyr1 + 1)
  Fbar <- matrix(fbar(stk.win), ncol = btyr2 - btyr1 + 1)
  sel <- sweep(sel, 2, Fbar, "/")

  if (flgsel == 0) { # take means of selection
    sel[] <- apply(sel, 1, mean)
  }
  if (flgmatwt==0){ # take means of wts
    west[] <- apply(west, 1, mean)
    weca[] <- apply(weca, 1, mean)
  } 

  Mat <- apply(mat(stk.win), 1, mean)[drop=TRUE]
  M <- apply(m(stk.win), 1, mean)[drop=TRUE] #mean(m(stk.win))
  Fprop <- apply(harvest.spwn(stk.win), 1, mean)[drop=TRUE] # vmean(harvest.spwn(stk.win))
  Mprop <- apply(m.spwn(stk.win), 1, mean)[drop=TRUE] # mean(m.spwn(stk.win))

  # get ready for the simulations
  Nmod <- nrow(SR)
  NF <- length(Fscan)
  ages <- dims(stk) $ age

  ssby <- Ferr <- array(0, c(Nrun,Nmod))
  Ny <- Fy <- WSy <- WCy <- Cy <- Wy <- array(0, c(ages, Nrun, Nmod))
  rsam <- array(sample(1:ncol(weca), Nrun * Nmod, TRUE), c(Nrun, Nmod)) 
  Wy[] <- c(weca[, c(rsam)])

  # initial recruitment
  R <- mean( data $ rec)
  ssbs <- cats <- recs <- array(0, c(7, NF))
  pssb1 <- pssb2 <- array(0, NF)
  ferr <- ssbsa <- catsa <- recsa <- array(0, c(NF, keep, Nmod))
  begin <- Nrun - keep + 1

  if (verbose) loader(0)
  for (i in 1:NF) {

    # The F value to test
    Fbar <- Fscan[i]

    # the selection patterns for the first year
    Zpre <- (Fbar * Fprop * sel[,rsam[1,]] + M * Mprop)
    Zpos <- (Fbar * (1-Fprop) * sel[,rsam[1,]] + M * (1-Mprop))
    # run Z out to age 50 ...
    Zcum <- c(0, cumsum(Fbar * sel[c(1:ages, rep(ages, 49 - ages)), rsam[1,]] + M[c(1:ages, rep(ages, 49 - ages))]))
    N1 <- R * exp(- unname(Zcum))

    # set up age structure in first year for all simulations
    Ny[,1,] <- c(N1[1:(ages-1)], sum(N1[ages:50]))

    # calculate ssb in first year using a different stock.wt for each simulation
    ssby[1,] <- colSums(Mat * Ny[,1,] * west[,rsam[1,]] / exp(Zpre))

    # loop over years
    for (j in 2:Nrun) {
      # get ssb from previous year
      SSB <- ssby[j-1,]

      # predict recruitment using various models
      if (process.error) {
        allrecs <- sapply(unique(SR $ mod), function(mod) exp(match.fun(mod) (SR, SSB) + rnorm(Nmod, 0, SR $ cv)))
      } else {
        allrecs <- sapply(unique(SR $ mod), function(mod) exp(match.fun(mod) (SR, SSB)))
      }
      select <- cbind(seq(Nmod), as.numeric(factor(SR $ mod, levels = unique(SR $ mod))))
      Ny[1,j,] <- allrecs[select]

      # apply HCR
      Fnext <- Fbar * pmin(1, SSB/Btrigger)      

      # apply some noise to the F
      Ferr[j,] <- Fphi * Ferr[j-1,] + rnorm(Nmod, 0, Fcv) 
      Fnext <- exp(Ferr[j,]) * Fnext

      # get a selection pattern for each simulation and apply this to get N
      Zpre <- rep(Fnext, each = length(Fprop)) * Fprop * sel[, rsam[j,]] + M * Mprop

      # get Fy
      Fy[    , j-1, ] <- rep(Fnext, each = ages) * sel[, rsam[j-1,]]

      Ny[  -1,   j, ] <- Ny[1:(ages-1), j-1, ] * exp(-Fy[1:(ages-1), j-1, ] - M[1:(ages-1)])
      Ny[ages,   j, ] <- Ny[ages, j, ] + Ny[ages, j-1, ] * exp(-Fy[ages, j-1, ] - M[ages])
      # calculate ssb and catch.n 
      ssby[j, ] <- apply(array(Mat * Ny[,j,] * west[, rsam[j,]] / exp(Zpre), c(ages, Nmod)), 2, sum)
      Cy[, j, ] <- Ny[, j-1, ] * Fy[, j-1, ] / (Fy[, j-1, ] + M) * (1 - exp(-Fy[, j-1, ] - M))
    }

    # convert to catch weight
    Cw <- Cy * Wy
    Cat <- apply(Cw, 2:3, sum)

    # summarise everything and spit out!
    quants <- c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
    ssbs[, i]   <- quantile(ssby[begin:Nrun, ], quants)
    cats[, i]   <- quantile(Cat[begin:Nrun, ], quants)
    recs[, i]   <- quantile(Ny[1, begin:Nrun, ], quants)
    ferr[i, , ] <- Ferr[begin:Nrun, ]
    ssbsa[i, , ] <- ssby[begin:Nrun, ]
    catsa[i, , ] <- Cat[begin:Nrun, ]
    recsa[i, , ] <- Ny[1, begin:Nrun, ]

      if (verbose) loader(i/NF)
  }

  list(ssbs = ssbs, cats = cats, recs = recs, ferr = ferr, ssbsa = ssbsa, catsa = catsa, recsa = recsa,
       Mat = Mat, M = M, Fprop = Fprop, Mprop = Mprop, west = west, weca = weca, sel = sel,
       Fscan = Fscan)
}


######## Creates equilibrium plots for a given Blim

#' stock recruitment function
#'
#'
#' @param ab the model parameters
#' @param ssb a vector of ssb
#' @return log recruitment according to model
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
Eqplot <- function (sim, fit, Blim, Bpa = 1.4 * Blim, ymax = c(NA,NA,NA), plot = TRUE)
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
  xmax <- max(Fscan)
  y.max <- if (!is.na(ymax[1])) ymax[1] else max(recs[7,], rec)

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
  lines(c(flim, flim), c(0, y.max), col = 3)

# recruits versus SSB
  y.max <- if (!is.na(ymax[2])) ymax[2] else max(ssbs[7,])
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
  lines(c(0,xmax), c(Blim, Blim))
  text(x = 0.1, y = Blim * 1.1, "Blim", cex = 0.7)
  points(FbarO, ssb, pch = 21, cex = .75, bg = 1)
  lines(c(flim, flim), c(0, y.max), col = 3)
}
  FCrash5  <- Fscan[which.max(cats[2,]):NF][ which(cats[2, which.max(cats[2,]):NF] < 0.05*max(cats[2,]) )[1] ]

  FCrash50 <- Fscan[which.max(cats[4,]):NF][ which(cats[4, which.max(cats[4,]):NF] < 0.05*max(cats[4,]) )[1] ]

if (plot) {
# catch versus Fbar
  y.max <- if (!is.na(ymax[3])) ymax[3] else max(cats[7,])
  plot(Fscan, cats[7,], type = "l", lty = 4,
       ylim = c(0, y.max), xlim = c(0, max(Fscan)), ylab = "", xlab = "")
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
  lines(c(flim, flim), c(0, y.max), col = 3)
  lines(c(FCrash5, FCrash5), c(0, y.max), col = 5)
  lines(c(FCrash50, FCrash50), c(0, y.max), col = 5)
  lines(Fscan, catm, lty=1, col = 2)
  lines(rep(Fscan[maxcatm], 2), c(0, y.max), lty = 1, col = 5)

# F versus SSB
  plot(Fscan, 1-pssb1, type = "l", lty = 2,
       ylim = c(0,1), xlim = c(0,max(Fscan)), ylab = "", xlab = "")
  title(ylab = "Prob MSY, SSB<Bpa or Blim", xlab = "F bar", 
        cex.lab = 0.75, line = 2.5, cex.main = 0.75)
  mtext(text = "d) Prob MSY and Risk to SSB", cex = 0.75, side = 3, line = 0.5)
  lines(Fscan, 1 - pssb2, lty = 4)
  
  
  text(x = max(Fscan[pssb2 > 0.5]) - .05, y = 0.5, "SSB<Bpa", cex = 0.75)
  text(x = max(Fscan[pssb1 > 0.7]) + .1, y = 0.3, "SSB<Blim", cex = 0.75)

  lines(c(flim,flim), c(0,1), col = 3)
  lines(c(0,flim), c(0.05,0.05), lty = 2, col = 3)
  text(x = 0.1, y = 0.075, "5%", cex = 0.75, col = 3)
  lines(c(flim10,flim10), c(0,1), col = "darkgreen")
  lines(c(0,flim10), c(0.1,0.1), lty = 2, col = "darkgreen")
  text(x = 0.05, y = 0.125, "10%", cex = 0.75, col = "darkgreen")
  lines(fmsy.dens $ x, cumsum(fmsy.dens $ y * diff(fmsy.dens $ x)[1]), col = 4)
#  lines(c(0, vcum), c(0.5, 0.5), lty = 2, col = 4)
  text(x = 0.9, y = 0.8, "Prob of Fmsy", cex = 0.75, col = 4)
  lines(rep(vcum,2), c(0,1), lty = 1, col = 4)

  lines(c(Fscan[maxcatm],Fscan[maxcatm]), c(0,1), col = 5)
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




############################################
# more plotting routines
##############################################
#' plots the ordered posterior densities accross MCMC iterations for each model considered in \code{fitModels}
#'
#'
#' @param fit an fitted MCMC returned from \code{fitModels}
#' @return NULL produces a plot
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
LLplot <- function(fit) 
{
  lliks <- sapply(fit $ fits, function(x) sort(x $ llik))

  plot(0, 0, type = "n", 
       main = paste("LL of Bayes model:", fit $ stknam), 
       xlab = "model order", ylab = "log likelihood", 
       ylim = range(lliks), xlim = c(1, nrow(lliks)))
  for (i in 2:ncol(lliks)-1) 
  {
    lines(lliks[,i], lty = i, col = i)
  }
  lines(sort(fit $ fits $ BMA $ llik), lwd = 2)
  legend(x = "bottomright", legend = c(names(fit $ fits)), 
         lty = c(2:ncol(lliks)-1,1), col = c(2:ncol(lliks)-1,1), lwd = c(rep(1, ncol(lliks)-1),2)) 
}


#' plots the fits for each model considered in \code{fitModels}
#'
#'
#' @param fit an fitted MCMC returned from \code{fitModels}
#' @return NULL produces one or several plots
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
LLmodplot <- function(fit)
{
  rec <- fit $ data $ rec[order(fit $ data $ ssb)]
  ssb <- sort(fit $ data $ ssb)

  for (mod in seq(fit $ fits)[2:length(fit $ fits)-1]) 
  { 
    model <- paste(names(fit $ fits)[mod], "fit")
    Rsym <- fit $ preds[[mod]]

    plot(ssb, rec, xlab = 'SSB', ylab = 'Recruits', 
         main = paste(fit $ stknam, model), 
         ylim = c(0, max(rec)), xlim = c(0, max(ssb)), type = "n")
    for (i in sample(1:nrow(Rsym), 1000)) {
      lines(ssb, Rsym[i,], col = paste0(grey(0), "10") )
    }
    for (i in c(0.05, 0.95)) {
      lines(ssb, apply(Rsym, 2, quantile, i), col = 4, lwd = 3)
    }
    lines(ssb, apply(Rsym, 2, quantile, 0.5), col = 7, lwd = 3)
    #TODO lines(ssb, Rsym[which.max(res[[mod]]),], col=1, lwd=2)
    lines(fit $ data $ ssb, fit $ data $ rec, col = "red")
    points(ssb, rec, pch = 19, col = "red")  
  }
}







