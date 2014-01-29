
##############################################
# utility functions
##############################################



#' Get starting values for models
#'
#' a quick fix!!
#'
#' @param model XXX
#' @param data XXX
#' @return vector of starting values
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
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
#' @param ... whatever
#' @return result of paste with sep = ""
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
paste0 <- function(...) paste(..., sep = "")


#' @title Progress function
#'
#'
#' @param p Value
#' @return NULL
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




#' @title simulates the equilibrium results for a population
#'
#' @description XXX
#' 
#' @param fit A list returned from the function fitModels
#' @param bio.years The years to sample maturity, weights and M from
#' @param bio.const A flag, if FALSE mean of the biological values from the years selected are used
#' @param sel.years The years to sample the selection patterns from
#' @param sel.const A flag, if FALSE mean of the selection patterns from the years selected are used
#' @param Fscan F values to scan over
#' @param Fcv Assessment error in the advisory year
#' @param Fphi Autocorrelation in assessment error in the advisory year
#' @param recruitment.trim A numeric vector with two log-value clipping the extreme
#' recruitment values from a continuous lognormal distribution. The values must
#' be set as c("high","low").
#' @param Btrigger If other than 0 (default) the target F applied is reduced by
#' SSB/Btrigger
#' @param Nrun The number of years to run in total (last 50 years from that will be retained)
#' @param process.error Use stochastic recruitment or mean recruitment?  (TRUE = predictive)
#' @param verbose Flag, if TRUE (default) indication of the progress of the
#' simulation is provided in the console. Useful to turn to FALSE when 
#' knitting documents.
#' @return A list containing the following objects:
#' \itemize{
#' \item ssbs A matrix containing the 0.025, 0.050, 0.25, 0.50, 0.75, 0.95, 0.975
#' percentiles of ssb for the scanned fishing mortalities (see Fscan below).
#' \item cats A matrix containing the 0.025, 0.050, 0.25, 0.50, 0.75, 0.95, 0.975
#' percentiles of catches for the scanned fishing mortalities (see Fscan below).
#' \item recs A matrix containing the 0.025, 0.050, 0.25, 0.50, 0.75, 0.95, 0.975
#' percentiles of rec for the scanned fishing mortalities (see Fscan below).
#' \item ferr An array 
#' \item ssbsa
#' \item catsa
#' \item recsa
#' \item Mat
#' \item M
#' \item Fprop
#' \item Mprop
#' \item west
#' \item weca
#' \item sel
#' \item Fscan
#' }
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export

EqSim <- function(fit,
                  bio.years = c(2007, 2012), # years sample weights, M and mat
                  bio.const = 0,
                  sel.years= c(2005, 2010), # years sample sel and discard proportion by number from
                  sel.const = 0,
                  Fscan = seq(0, 1, len = 20), # F values to scan over
                  Fcv = 0,
                  Fphi = 0,
                  recruitment.trim = c(3, -3),
                  Btrigger = 0,
                  Nrun = 200, # number of years to run in total
                  process.error = TRUE, # use predictive recruitment or mean recruitment? (TRUE = predictive)
                  verbose = TRUE)
{
  
  if (abs(Fphi) >= 1) stop("Fphi, the autocorelation parameter for log F should be between (-1, 1)")
  if ((recruitment.trim[1] + recruitment.trim[2])> 0) stop("recruitment truncation must be between a high - low range")
  btyr1 <- bio.years[1]
  btyr2 <- bio.years[2]
  slyr1 <- sel.years[1]
  slyr2 <- sel.years[2]
  
  keep <- min(Nrun, 50)
  
  SR <- fit $ fit
  data <- fit $ data
  stk <- fit $ stk
  
  # forecast settings (mean wt etc)
  stk.win <- window(stk, start = btyr1, end = btyr2)
  stk.winsel <- window(stk, start = slyr1  , end = slyr2)
  
  west <- matrix(stock.wt(stk.win), ncol = btyr2 - btyr1 + 1)
  weca <- matrix(catch.wt(stk.win), ncol = btyr2 - btyr1 + 1)
  wela <- matrix(landings.wt(stk.win), ncol = btyr2 - btyr1 + 1)
  Mat <- matrix(mat(stk.win), ncol = btyr2 - btyr1 + 1)
  M <- matrix(m(stk.win), ncol = btyr2 - btyr1 + 1)
  landings <- matrix(landings.n(stk.winsel), ncol = slyr2 - slyr1 + 1)
  catch <- matrix(catch.n(stk.winsel), ncol = slyr2 - slyr1 + 1)
  sel <- matrix(harvest(stk.winsel), ncol = slyr2 - slyr1 + 1)
  Fbar <- matrix(fbar(stk.winsel), ncol = slyr2 - slyr1  + 1)
  sel <- sweep(sel, 2, Fbar, "/")
  
  if (sel.const == TRUE) { # take means of selection
    sel[] <- apply(sel, 1, mean)
    landings[]  <- apply(landings, 1, mean)
    catch[]  <- apply(catch, 1, mean)
  }
  if (bio.const== TRUE){ # take means of wts Mat and M and ratio of landings to catch
    west[] <- apply(west, 1, mean)
    weca[] <- apply(weca, 1, mean)
    Mat[] <- apply(Mat, 1, mean)
    M[] <- apply(M, 1, mean) #me
  }
  land.cat= landings / catch  # ratio of number of landings to catch
  
  Fprop <- apply(harvest.spwn(stk.winsel), 1, mean)[drop=TRUE] # vmean(harvest.spwn(stk.win))
  Mprop <- apply(m.spwn(stk.win), 1, mean)[drop=TRUE] # mean(m.spwn(stk.win))
  
  # get ready for the simulations
  Nmod <- nrow(SR)
  NF <- length(Fscan)
  ages <- dims(stk) $ age
  
  ssby <- Ferr <- array(0, c(Nrun,Nmod))
  Ny <- Fy <- WSy <- WCy <- Cy <- Wy <- Wl <- Ry <- array(0, c(ages, Nrun, Nmod))
  rsam <- array(sample(1:ncol(weca), Nrun * Nmod, TRUE), c(Nrun, Nmod))
  rsamsel <- array(sample(1:ncol(sel), Nrun * Nmod, TRUE), c(Nrun, Nmod))
  Wy[] <- c(weca[, c(rsam)])
  Wl[] <- c(wela[, c(rsam)])
  Ry[]  <- c(land.cat[, c(rsamsel)])
  # initial recruitment
  R <- mean( data $ rec)
  ssbs <- cats <- lans <- recs <- array(0, c(7, NF))
  pssb1 <- pssb2 <- array(0, NF)
  ferr <- ssbsa <- catsa <- lansa <- recsa <- array(0, c(NF, keep, Nmod))
  begin <- Nrun - keep + 1
  
  # New from Simmonds' 29.1.2014
  resids= array(rnorm(Nmod*(Nrun+1), 0, SR$cv),c(Nmod, Nrun+1))     # set up random
  lims = t(array(SR$cv,c(Nmod,2))) * recruitment.trim
  for (k in 1:Nmod) { resids[k,resids[k,]>lims[1,k]]=lims[1,k]}
  for (k in 1:Nmod) { resids[k,resids[k,]<lims[2,k]]=lims[2,k]}
  # end New from Simmonds 29.1.2014
  
  if (verbose) loader(0)
  for (i in 1:NF) {
    
    # The F value to test
    Fbar <- Fscan[i]
    
    # the selection patterns for the first year
    Zpre <- ( sel[,rsamsel[1,]]*Fbar * Fprop + M[,rsam[1,]] * Mprop)
    Zpos <- (Fbar * (1-Fprop) * sel[,rsamsel[1,]] + M[,rsam[1,]] * (1-Mprop))
    # run Z out to age 50 ...
    Zcum <- c(0, cumsum(Fbar * sel[c(1:ages, rep(ages, 49 - ages)), rsamsel[1,]] + M[c(1:ages, rep(ages, 49 - ages)), rsam[1,]]))
    
    N1 <- R * exp(- unname(Zcum))
    
    # set up age structure in first year for all simulations
    Ny[,1,] <- c(N1[1:(ages-1)], sum(N1[ages:50]))
    
    # calculate ssb in first year using a different stock.wt and Mat selection and M for each simulation
    ssby[1,] <- colSums(Mat[,rsam[1,]] * Ny[,1,] * west[,rsam[1,]] / exp(Zpre)[])
    
    # loop over years
    for (j in 2:Nrun) {
      # get ssb from previous year
      SSB <- ssby[j-1,]
      
      # predict recruitment using various models
      if (process.error) {
        # Changes 29.1.2014
        # new random draws each time
        # allrecs <- sapply(unique(SR $ mod), function(mod) exp(match.fun(mod) (SR, SSB) + rnorm(Nmod, 0, SR $ cv)))
        # same random draws used for each F
        allrecs <- sapply(unique(SR$mod), function(mod) exp(match.fun(mod)(SR, SSB) + resids[,j]))
        # end Changes 29.1.2014
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
      Zpre <- rep(Fnext, each = length(Fprop)) * Fprop * sel[, rsamsel[j,]] + M[, rsam[j,]] * Mprop
      
      # get Fy
      Fy[ , j-1, ] <- rep(Fnext, each = ages) * sel[, rsamsel[j-1,]]
      
      Ny[ -1, j, ] <- Ny[1:(ages-1), j-1, ] * exp(-Fy[1:(ages-1), j-1, ] - M[1:(ages-1), rsam[j-1,]])
      Ny[ages, j, ] <- Ny[ages, j, ] + Ny[ages, j-1, ] * exp(-Fy[ages, j-1, ] - M[ages, rsam[j-1,]])
      # calculate ssb and catch.n
      ssby[j, ] <- apply(array(Mat[, rsam[j,]] * Ny[,j,] * west[, rsam[j,]] / exp(Zpre), c(ages, Nmod)), 2, sum)
      Cy[, j, ] <- Ny[, j-1, ] * Fy[, j-1, ] / (Fy[, j-1, ] + M[, rsam[j-1,]]) * (1 - exp(-Fy[, j-1, ] - M[, rsam[j-1,]]))
    }
    
    # convert to catch weight
    Cw <- Cy * Wy   # catch Numbers *catch wts
    land <- Cy*Ry*Wl # catch Numbers * Fraction (in number) landed and landed wts
    Lan=apply(land,2:3,sum)
    Cat <- apply(Cw, 2:3, sum)
    
    # summarise everything and spit out!
    quants <- c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
    ssbs[, i] <- quantile(ssby[begin:Nrun, ], quants)
    cats[, i] <- quantile(Cat[begin:Nrun, ], quants)
    lans[, i] <- quantile(Lan[begin:Nrun, ], quants)
    recs[, i] <- quantile(Ny[1, begin:Nrun, ], quants)
    ferr[i, , ] <- Ferr[begin:Nrun, ]
    ssbsa[i, , ] <- ssby[begin:Nrun, ]
    catsa[i, , ] <- Cat[begin:Nrun, ]
    lansa[i, , ] <- Lan[begin:Nrun, ]
    recsa[i, , ] <- Ny[1, begin:Nrun, ]
    
    if (verbose) loader(i/NF)
  }
  
  list(ssbs = ssbs, cats = cats, lans = lans, recs = recs, ferr = ferr, ssbsa = ssbsa, catsa = catsa, lansa = lansa, recsa = recsa,
       Mat = Mat, M = M, Fprop = Fprop, Mprop = Mprop, west = west, weca = weca, sel = sel, land.cat = land.cat,
       Fscan = Fscan)
}


######## Creates equilibrium plots for a given Blim

#' @title Plot summary of the equilibrium simulation
#'
#'
#' @param sim A list obtained from function \code{EqSim}
#' @param fit A list obtained from function \code{fitModels}
#' @param Blim A value for Blim
#' @param Bpa A value for Bpa. If not specified Bpa is set to 1.4 Blim
#' @param ymax A vector of length three indicating the maximum values for
#' each of the first three panel plots
#' @param plot Flag, if TRUE (default) returns both values and plot
#' @param yield Character vector of length 1, either "landings" (default) or "catch"
#' @param extreme.trim Call John Simmonds :-)
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
Eqplot <- function (sim, fit, Blim, Bpa = 1.4 * Blim, ymax = c(NA,NA,NA), 
                    plot = TRUE, yield='landings',
                    extreme.trim)
{
  
  stk <- fit $ stk
  
  Nmod <- dim(sim $ ssbsa)[3]
  Nyrs <- dim(sim $ ssbsa)[2]
  
  Fscan <- sim $ Fscan
  
  # Einar added 29.1.2014
  if(missing(extreme.trim)) {
    catm <- apply(sim $ catsa, 1, mean)
    lanm <- apply(sim $ lansa, 1, mean)
  } else {
    i <- sim$catsa > quantile(sim$catsa,extreme.trim[2]) |
      sim$catsa < quantile(sim$catsa,extreme.trim[1])
    sim$catsa[i] <- NA
    catm <- apply(sim $ catsa, 1, mean, na.rm=TRUE)
    i <- sim$lansa > quantile(sim$lansa,extreme.trim[2]) |
      sim$lansa < quantile(sim$lansa,extreme.trim[1])
    sim$lansa[i] <- NA
    lanm <- apply(sim $ lansa, 1, mean, na.rm=TRUE)
  }
    
  # end Einar added 29.1.2014
  
  
  maxcatm <- which.max(catm)
  maxlanm <- which.max(lanm)
  catsam <- apply(sim $ catsa, c(1,3), mean)
  lansam <- apply(sim $ lansa, c(1,3), mean)
  maxpf <- apply(catsam, 2, which.max)
  maxpfl <- apply(lansam, 2, which.max)
  #  using catch as default or landings if required
  if (yield=='landings') { fmsy <- Fscan[maxpfl]  } else { fmsy <- Fscan[maxpf] }
  
  msym <- mean(fmsy)
  vcum <- median(fmsy)
  fmsy.dens <- density(fmsy)
  vmode <- fmsy.dens $ x[which.max(fmsy.dens $ y)]
  
  pssb1 <- apply(sim $ ssbsa > Blim, 1, mean)
  pssb2 <- apply(sim $ ssbsa > Bpa, 1, mean)
  
  pp1 <- max(which(pssb1>.95))
  grad <- diff(Fscan[pp1 + 0:1]) / diff(pssb1[pp1 + 0:1])
  flim <- Fscan[pp1] + grad * (0.95 - pssb1[pp1]) # linear interpolation i think..
  
  
  maint <- fit $ stknam
  
  rec <- fit $ data $ rec
  ssb <- fit $ data $ ssb
  
  Catchs <- catch(stk)[, 1:length(ssb), drop = TRUE]
  Landings <- landings(stk)[, 1:length(ssb), drop = TRUE]
  
  FbarO <- fbar(stk)[, 1:length(ssb), drop = TRUE]
  
  recs <- sim $ recs
  ssbs <- sim $ ssbs
  cats <- sim $ cats
  lans <- sim $ lans
  ssbsa <- sim $ ssbsa
  catsa <- sim $catsa
  lansa <- sim $lansa
  
  
  NF <- length(Fscan)
  pp1 <- max(which(pssb1>.50))
  grad <- diff(Fscan[pp1 + 0:1]) / diff(pssb1[pp1 + 0:1])
  flim50 <- Fscan[pp1]+grad*(0.5-pssb1[pp1]) # linear interpolation i think..
  
  pp1 <- max(which(pssb1>.90))
  grad <- diff(Fscan[pp1 + 0:1]) / diff(pssb1[pp1 + 0:1])
  flim10 <- Fscan[pp1]+grad*(0.9-pssb1[pp1]) # linear interpolation i think..
  
  
  if (plot) {
    op <- par(mfrow = c(2, 2), mar = c(2.5, 4, 1.5, 1), oma = c(0, 0, 0, 0),
              cex.axis = 0.75, tcl = 0.25, mgp = c(0, 0.25, 0), las = 1)
    
    # recruits versus Fbar
    xmax <- max(Fscan)
    y.max <- if (!is.na(ymax[1])) ymax[1] else max(recs[6,], rec)
    
    plot(Fscan, recs[6,], type = "l", lty = 4,
         ylim = c(0, y.max), xlim = c(0, xmax), ylab = "", xlab = "")
    title(ylab = "Recruitment", xlab = "F bar",
          cex.lab = 0.75, line = 2.5, cex.main = 0.75)
    mtext(text = paste(maint," a) Recruits"), cex = 0.75, side = 3, line = 0.5)
    lines(Fscan, recs[4,], lty = 1)
    lines(Fscan, recs[2,], lty = 3)
    points(FbarO, rec, pch = 21, cex = .75, bg = 1)
    lines(c(flim, flim), c(0, y.max), col = 3)
    
    # recruits versus SSB
    y.max <- if (!is.na(ymax[2])) ymax[2] else max(ssbs[6,])
    plot(Fscan, ssbs[6,], type = "l", lty = 4,
         ylim = c(0, y.max), xlim = c(0, xmax), ylab = "", xlab = "")
    title(ylab = "SSB", xlab = "F bar",
          cex.lab = 0.75, line = 2.5, cex.main = 0.75)
    mtext(text = "b) Spawning Stock Biomass", cex = 0.75, side = 3, line = 0.5)
    lines(Fscan, ssbs[4,], lty=1)
    lines(Fscan, ssbs[2,], lty=3)
    lines(c(0,xmax), c(Blim, Blim))
    text(x = 0.1, y = Blim * 1.1, "Blim", cex = 0.7)
    points(FbarO, ssb, pch = 21, cex = .75, bg = 1)
    lines(c(flim, flim), c(0, y.max), col = 3)
  }
  FCrash5 <- Fscan[which.max(cats[2,]):NF][ which(cats[2, which.max(cats[2,]):NF] < 0.05*max(cats[2,]) )[1] ]
  
  FCrash50 <- Fscan[which.max(cats[4,]):NF][ which(cats[4, which.max(cats[4,]):NF] < 0.05*max(cats[4,]) )[1] ]
  
  if (plot) {
    if ( yield=="catch" )  {
      # catch versus Fbar
      y.max <- if (!is.na(ymax[3])) ymax[3] else max(cats[6,])
      plot(Fscan, cats[6,], type = "l", lty = 4,
           ylim = c(0, y.max), xlim = c(0, max(Fscan)), ylab = "", xlab = "")
      title(ylab = "Catch", xlab = "F bar",
            cex.lab = 0.75, line = 2.5, cex.main = 0.75)
      mtext(text = "c) Catch", cex = 0.75, side = 3, line = 0.5)
      #points(fbarsa[ssbsa<BLIM], catsa[ssbsa<BLIM],pch=21,col=6,bg=6,cex=0.125)
      #points(fbarsa[ssbsa>BLIM]+.0075,catsa[ssbsa>BLIM],pch=21,col=7,bg=7,cex=0.125)
      lines(Fscan, cats[6,], lty=3)
      lines(Fscan, cats[4,], lty=1)
      lines(Fscan, cats[2,], lty=3)
      
      points(FbarO, Catchs, pch = 21, cex = .75, bg = 1)
      lines(c(flim, flim), c(0, y.max), col = 3)
      lines(c(FCrash5, FCrash5), c(0, y.max), col = 5)
      lines(c(FCrash50, FCrash50), c(0, y.max), col = 5)
      lines(Fscan, catm, lty=1, col = 2)
      lines(rep(Fscan[maxcatm], 2), c(0, y.max), lty = 1, col = 5)
    }
    
    if (yield ==  "landings") {
      # catch versus Fbar
      y.max <- if (!is.na(ymax[3])) ymax[3] else max(lans[6,])
      plot(Fscan, lans[6,], type = "l", lty = 4,
           ylim = c(0, y.max), xlim = c(0, max(Fscan)), ylab = "", xlab = "")
      title(ylab = "Landings", xlab = "F bar",
            cex.lab = 0.75, line = 2.5, cex.main = 0.75)
      mtext(text = "c) Landings", cex = 0.75, side = 3, line = 0.5)
      #points(fbarsa[ssbsa<BLIM], catsa[ssbsa<BLIM],pch=21,col=6,bg=6,cex=0.125)
      #points(fbarsa[ssbsa>BLIM]+.0075,catsa[ssbsa>BLIM],pch=21,col=7,bg=7,cex=0.125)
      #  set of quantile lines for landings
      lines(Fscan, lans[6,], lty=3,col=6)
      lines(Fscan, lans[4,], lty=1,col=6)
      lines(Fscan, lans[2,], lty=3,col=6)
      
      points(FbarO, Landings, pch = 21, cex = .75, bg = 1)
      lines(c(flim, flim), c(0, y.max), col = 3)
      lines(c(FCrash5, FCrash5), c(0, y.max), col = 5)
      lines(c(FCrash50, FCrash50), c(0, y.max), col = 5)
      lines(Fscan, lanm, lty=1, col = 2)
      lines(rep(Fscan[maxlanm], 2), c(0, y.max), lty = 1, col = 5)
      
    }
    
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
    #lines(c(flim10,flim10), c(0,1), col = "darkgreen")
    #lines(c(0,flim10), c(0.1,0.1), lty = 2, col = "darkgreen")
    #text(x = 0.05, y = 0.125, "10%", cex = 0.75, col = "darkgreen")
    lines(fmsy.dens $ x, cumsum(fmsy.dens $ y * diff(fmsy.dens $ x)[1]), col = 4)
    # lines(c(0, vcum), c(0.5, 0.5), lty = 2, col = 4)
    text(x = 0.9, y = 0.8, "Prob of Fmsy", cex = 0.75, col = 4)
    lines(rep(vcum,2), c(0,1), lty = 1, col = 4)
    
    lines(c(Fscan[maxcatm],Fscan[maxcatm]), c(0,1), col = 5)
  }
  
  
  out <- c(Blim, Bpa)
  if (yield =="landings") {
    outF <- c(flim, flim10, flim50, vcum, Fscan[maxlanm], FCrash5, FCrash50)
    outC <- approx(Fscan, lans[4,], xout = outF) $ y
  }
  if (yield =="catch") {
    outF <- c(flim, flim10, flim50, vcum, Fscan[maxcatm], FCrash5, FCrash50)
    outC <- approx(Fscan, cats[4,], xout = outF) $ y
  }
  outB <- approx(Fscan, ssbs[4,], xout = outF) $ y
  outTable <- rbind(outF, outB, outC)
  rownames(outTable) <- c("F","SSB",yield)
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







