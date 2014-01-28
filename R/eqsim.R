#' @title simulates the equilibrium results for a population
#'
#' @description XXX
#' 
#' @export
#' 
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' 
#' @param fit A list returned from the function fitModels
#' @param bio.years The years to sample maturity, weights and M from
#' @param bio.const A flag, if FALSE mean of the biological values are used
#' @param sel.years The years to sample sel and discard proportion by number from
#' @param sel.const A flag, if FALSE mean selection values used
#' @param Fscan F values to scan over
#' @param Fcv Assessment error in the advisory year
#' @param Fphi Autocorrelation in assessment error in the advisory year
#' @param Blim This we know
#' @param Bpa This we know
#' @param Btrigger If other than 0 (default) the target F applied is reduced by
#' SSB/Btrigger 
#' @param Nrun The number of years to run in total
#' @param process.error Use stochastic recruitment or mean recruitment?  (TRUE = predictive)
#' @param verbose Flag, if TRUE (default) indication of the progress of the
#' simulation is provided in the console. Useful to turn to FALSE when 
#' knitting documents.

eqsim_run <- function(fit,
                      bio.years = c(2008, 2012), # years sample weights, M and mat
                      bio.const = FALSE,
                      sel.years= c(2008, 2012), # years sample sel and discard proportion by number from
                      sel.const = FALSE,
                      Fscan = seq(0, 1, len = 20), # F values to scan over
                      Fcv = 0,
                      Fphi = 0,
                      Blim,
                      Bpa,                 
                      Btrigger = 0,
                      Nrun = 200, # number of years to run in total
                      process.error = TRUE, # use predictive recruitment or mean recruitment? (TRUE = predictive)
                      verbose = TRUE)
{
  
  if (abs(Fphi) >= 1) stop("Fphi, the autocorelation parameter for log F should be between (-1, 1)")
  
  btyr1 <- bio.years[1]
  btyr2 <- bio.years[2]
  slyr1 <- sel.years[1]
  slyr2 <- sel.years[2]
  
  keep <- min(Nrun, 50)
  
  SR <- fit $ sr.sto
  data <- fit $ rby[,c("rec","ssb","year")]
  stk <- fit $ stk
  
  # forecast settings (mean wt etc)
  stk.win <- window(stk, start = btyr1, end = btyr2)
  stk.winsel <- window(stk, start = slyr1  , end = slyr2)
  
  littleHelper <- function(x,i) {
    x2 <- x
    x2[i] <- NA
    x2[] <- apply(x2,1,mean,na.rm=TRUE)
    x[i] <- x2[i]
    return(x)
  }
  
  west <- matrix(stock.wt(stk.win), ncol = btyr2 - btyr1 + 1)
  i <- west == 0
  if(any(i)) west <- littleHelper(west,i)
  weca <- matrix(catch.wt(stk.win), ncol = btyr2 - btyr1 + 1)
  i <- weca == 0
  if(any(i)) weca <- littleHelper(weca,i)
  wela <- matrix(landings.wt(stk.win), ncol = btyr2 - btyr1 + 1)
  if(any(i)) wela <- littleHelper(wela,i)
  
  Mat <- matrix(mat(stk.win), ncol = btyr2 - btyr1 + 1)
  M <- matrix(m(stk.win), ncol = btyr2 - btyr1 + 1)
  landings <- matrix(landings.n(stk.winsel), ncol = slyr2 - slyr1 + 1)
  # if zero, use 0.10 of minimum value
  
  catch <- matrix(catch.n(stk.winsel), ncol = slyr2 - slyr1 + 1)
  sel <- matrix(harvest(stk.winsel), ncol = slyr2 - slyr1 + 1)
  Fbar <- matrix(fbar(stk.winsel), ncol = slyr2 - slyr1  + 1)
  sel <- sweep(sel, 2, Fbar, "/")
  
  if (sel.const == TRUE) { # take means of selection
    sel[] <- apply(sel, 1, mean)
    landings[]  <- apply(landings, 1, mean)
    catch[]  <- apply(catch, 1, mean)
  }
  if (bio.const==TRUE){ # take means of wts Mat and M and ratio of landings to catch
    west[] <- apply(west, 1, mean)
    weca[] <- apply(weca, 1, mean)
    Mat[] <- apply(Mat, 1, mean)
    M[] <- apply(M, 1, mean) #me
  }
  land.cat= landings / catch  # ratio of number of landings to catch
  
  # TODO: Check if this is sensible
  i <- is.na(land.cat)
  if(any(i)) land.cat[i] <- 1
  
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
  
  dimnames(ssbs) <- dimnames(cats) <- 
    dimnames(lans) <- dimnames(recs) <- 
    list(quants=c("p025","p05","p25","p50","p75","p95","p975"),
         fmort=Fscan)

  rbp2dataframe <- function(x,variable) {
    x <- data.frame(t(x))
    x$variable <- variable
    x$Ftarget <- as.numeric(row.names(x))
    rownames(x) <- NULL
    return(x)
  }
  rbp <- rbind(rbp2dataframe(recs,"Recruitment"),
               rbp2dataframe(ssbs,"Spawning stock biomass"),
               rbp2dataframe(cats,"Catch"),
               rbp2dataframe(lans,"Landings"))
  rbp <- rbp[,c(9,8,1:7)]
  
  # STOCK REFERENCE POINTS
  
  FCrash05 <- Fscan[which.max(cats[2,]):NF][ which(cats[2, which.max(cats[2,]):NF] < 0.05*max(cats[2,]) )[1] ]
  FCrash50 <- Fscan[which.max(cats[4,]):NF][ which(cats[4, which.max(cats[4,]):NF] < 0.05*max(cats[4,]) )[1] ]
  
  catm <- apply(catsa, 1, mean)
  lanm <- apply(lansa, 1, mean)
  maxcatm <- which.max(catm)
  maxlanm <- which.max(lanm)
  
  catsam <- apply(catsa, c(1,3), mean)
  lansam <- apply(lansa, c(1,3), mean)
  maxpf <- apply(catsam, 2, which.max)
  maxpfl <- apply(lansam, 2, which.max)

  FmsyLan <- Fscan[maxpfl]
  msymLan <- mean(FmsyLan)
  vcumLan <- median(FmsyLan)
  fmsy.densLan <- density(FmsyLan)
  vmodeLan <- fmsy.densLan$x[which.max(fmsy.densLan$y)]
  
  FmsyCat <- Fscan[maxpf]
  msymCat <- mean(FmsyCat)
  vcumCat <- median(FmsyCat)
  fmsy.densCat <- density(FmsyCat)
  vmodeCat <- fmsy.densCat$x[which.max(fmsy.densCat$y)]
  
  pFmsyCat  <- data.frame(Ftarget=fmsy.densCat$x,
                          value=cumsum(fmsy.densCat$y * diff(fmsy.densCat$x)[1]),
                          variable="pFmsyCatch")
  pFmsyLan  <- data.frame(Ftarget=fmsy.densLan$x,
                          value=cumsum(fmsy.densLan$y * diff(fmsy.densLan$x)[1]),
                          variable="pFmsyLandings")
  pProfile <- rbind(pFmsyCat,pFmsyLan)
  
  # PA REFERENCE POINTS
  if(!missing(Blim)) {
    pBlim <- apply(ssbsa > Blim, 1, mean)
    
    i <- max(which(pBlim > .95))
    grad <- diff(Fscan[i + 0:1]) / diff(pBlim[i + 0:1])
    flim <- Fscan[i] + grad * (0.95 - pBlim[i]) # linear interpolation i think..
    
    i <- max(which(pBlim > .90))
    grad <- diff(Fscan[i + 0:1]) / diff(pBlim[i + 0:1])
    flim10 <- Fscan[i]+grad*(0.9-pBlim[i]) # linear interpolation i think..
    
    i <- max(which(pBlim > .50))
    grad <- diff(Fscan[i + 0:1]) / diff(pBlim[i + 0:1])
    flim50 <- Fscan[i]+grad*(0.5-pBlim[i]) # linear interpolation i think..
    
    pBlim <- data.frame(Ftarget = Fscan,value = 1-pBlim,variable="Blim")
    pProfile <- rbind(pProfile,pBlim)
  } else {
    flim <- flim10 <- flim50 <- Blim <- NA
  }
  
  if(!missing(Bpa)) {
    pBpa <- apply(ssbsa > Bpa, 1, mean) 
    pBpa <- data.frame(Ftarget = Fscan,value = 1-pBpa,variable="Bpa")
    pProfile <- rbind(pProfile,pBpa)
  } else {
    Bpa <- NA
  }
  
  # GENERATE REF-TABLE 
  catF <- c(flim, flim10, flim50, vcumCat, Fscan[maxcatm], FCrash05, FCrash50)
  lanF <- c(   NA,    NA,     NA, vcumLan, Fscan[maxlanm],       NA,       NA)
  catC <- approx(Fscan, cats[4,], xout = catF)$y
  lanC <- approx(Fscan, lans[4,], xout = lanF)$y
  catB <- approx(Fscan, ssbs[4,], xout = catF)$y
  lanB <- approx(Fscan, ssbs[4,], xout = lanF)$y
  
  Refs <- rbind(catF, lanF, catC, lanC, catB, lanB)
  rownames(Refs) <- c("catF","lanF","catch","landings","catB","lanB")
  colnames(Refs) <- c("Flim","Flim10","Flim50","medianMSY","meanMSY","FCrash05","FCrash50")
  
  #TODO: id.sim - user specified.
  
  return(list(ibya=list(Mat = Mat, M = M, Fprop = Fprop, Mprop = Mprop, 
                        west = west, weca = weca, sel = sel),
              rby=fit$rby, rbp=rbp, Blim=Blim, Bpa=Bpa, Refs = Refs,
              pProfile=pProfile,id.sim=fit$id.sr))
  
}


#' @title Plots of the results from eqsim
#'
#' @description XXX
#' 
#' @author Einar Hjorleifsson \email{einar.hjorleifsson@@gmail.com}
#' 
#' @export
#' 
#' @param sim An object returned from the function eqsim_run
#' @param ymax.multiplier A value that acts as a multiplier of the maximum observed
#' variable being plotted. E.g. 1.2 means that for each of the three panels a, b and c
#' the ymax is set to 1.2 of the maximum observed recruitment, spawning stock biomass
#' and yield (catch or landings, depending on user input.
#' @param catch Boolean, if TRUE (default) returns a plot based on catch. If false
#' returns a plot based on landings.

eqsim_plot <- function(sim, ymax.multiplier=1.2, catch=TRUE) 
{
  
  # littleHelper function
  littleHelper <- function(rbp,dat,Flim) {
    ymax <- max(dat[,2]*ymax.multiplier)
    with(rbp,plot(Ftarget,p95,type="l",lty=4,ylab="", xlab="", ylim=c(0,ymax)))
    #title(ylab=Variable, xlab="F bar", cex.lab=0.75, line=2.5, cex.main=0.75)
    #mtext(text = paste("c)", Variable), cex = 0.75, side = 3, line = 0.5)
    with(rbp,lines(Ftarget,p50, lty = 1))
    with(rbp,lines(Ftarget,p05, lty = 4))
    points(dat[,1],dat[,2],pch=21,cex=0.75,bg=1)
    abline(v=Flim,col="red")
    text(0.98*Flim,0,"Flim",srt=90,pos=3,col="red",cex=0.75)
  }
                           
                           
  
  
  rby <- sim$rby
  #for (i in c(1,2,4,5)) rby[,i] <- rby[,i]/Scale
  
  rbp <- sim$rbp
  #for(i in 3:9) rbp[,i] <- rbp[,i]/Scale
  
  refs <- sim$Refs
  #refs[3:6,] <- refs[3:6,]/Scale
  
  pProfile <- sim$pProfile
  
  Flim <- sim$Refs[1,1]
  FCrash5 <- sim$Refs[1,6]
  FCrash50 <- sim$Refs[1,7]
  op <- par(mfrow = c(2, 2), mar = c(2.5, 4, 1.5, 1), oma = c(0, 0, 0, 0),
              cex.axis = 0.75, tcl = 0.25, mgp = c(0, 0.25, 0), las = 1)
    
  # A: Recruitment plot
  Variable <- "Recruitment"
  i <- rbp$variable == Variable
  littleHelper(rbp[i,],dat=rby[,c("fbar","rec")],Flim)
  title(ylab=Variable, xlab="F bar", cex.lab=0.75, line=2.5, cex.main=0.75)
  mtext(text = paste(sim$id.sim," a) Recruits"), cex = 0.75, side = 3, line = 0.5)

  # B: SSB plot  
  Variable <- "Spawning stock biomass"
  i <- rbp$variable == Variable
  littleHelper(rbp[i,],dat=rby[,c("fbar","ssb")],Flim)
  title(ylab=Variable, xlab="F bar", cex.lab=0.75, line=2.5, cex.main=0.75)
  mtext(text = paste("b)", Variable), cex = 0.75, side = 3, line = 0.5)
  
  # C: Yield plot
  if (catch)  {
    # catch versus Fbar
    Variable <- "Catch"
    i <- rbp$variable == Variable
    littleHelper(rbp[i,],dat=rby[,c("fbar","catch")],Flim)
    title(ylab=Variable, xlab="F bar", cex.lab=0.75, line=2.5, cex.main=0.75)
    mtext(text = paste("c)", Variable), cex = 0.75, side = 3, line = 0.5)
    #add landings
    j <- rbp$variable == "Landings"
    with(rbp[j,],lines(Ftarget,p50, lty = 2))
    
    medianFmsy <- sim$Refs[1,4]
    catch_at_medianFmsy <- sim$Refs[3,4]
    meanFmsy <- sim$Refs[1,5]
    catch_at_meanFmsy <- sim$Refs[3,5]
    
    lines(rep(medianFmsy,2),c(0,catch_at_medianFmsy),col="brown")
    text(0.98*medianFmsy,0,"median Fmsy",srt=90,pos=4,col="brown",cex=0.75)
    lines(rep(meanFmsy,2),c(0,catch_at_meanFmsy),col="brown",lty=2)
    text(0.98*meanFmsy,0,"mean Fmsy",srt=90,pos=4,col="brown",cex=0.75)
    #lines(Fscan, catm, lty=1, col = 2)
    #lines(rep(Fscan[maxcatm], 2), c(0, y.max), lty = 1, col = 5)
  } else {
    # landings versus Fbar
    Variable <- "Landings"
    i <- rbp$variable == Variable
    littleHelper(rbp[i,],dat=rby[,c("fbar","landings")],Flim)
    title(ylab=Variable, xlab="F bar", cex.lab=0.75, line=2.5, cex.main=0.75)
    mtext(text = paste("c)", Variable), cex = 0.75, side = 3, line = 0.5)
    #add catch
    j <- rbp$variable == "Catch"
    with(rbp[j,],lines(Ftarget,p50, lty = 2))
    
    medianFmsy <- sim$Refs[2,4]
    landings_at_medianFmsy <- sim$Refs[4,4]
    meanFmsy <- sim$Refs[2,5]
    landings_at_meanFmsy <- sim$Refs[4,5]
    lines(rep(medianFmsy,2),c(0,landings_at_medianFmsy),col="brown")
    lines(rep(meanFmsy,2),c(0,landings_at_meanFmsy),col="brown",lty=2)
    #lines(Fscan, lanm, lty=1, col = 2)
    #lines(rep(Fscan[maxlanm], 2), c(0, y.max), lty = 1, col = 5)
  }
    
  # D: F versus SSB probability 
  
  Variable = "Blim"
  i <- pProfile$variable == Variable
  with(pProfile[i,],plot(Ftarget,value,type="l",lty=1,ylab="", xlab="",col="red"))
  
  title(ylab="Prob MSY, SSB<Bpa or Blim", xlab="F bar", cex.lab=0.75, line=2.5, cex.main=0.75)
  mtext(text = "d) Prob MSY and Risk to SSB", cex = 0.75, side = 3, line = 0.5)
  
  Variable = "Bpa"
  i <- pProfile$variable == Variable
  with(pProfile[i,],lines(Ftarget,value,type="l",lty=1,ylab="", xlab="",col="darkgreen"))
  
  Variable = "pFmsyCatch"
  i <- pProfile$variable == Variable
  with(pProfile[i,],lines(Ftarget,value,type="l",lty=1,ylab="", xlab="",col="cyan"))
  
  Variable = "pFmsyLandings"
  i <- pProfile$variable == Variable
  with(pProfile[i,],lines(Ftarget,value,type="l",lty=1,ylab="", xlab="",col="brown"))
    
  text(0.01,0.8, "SSB<Blim", cex = 0.75, col="red",pos=4)
  text(0.01,0.85, "SSB<Bpa", cex = 0.75, col="darkgreen",pos=4)
  text(0.01,0.9, "Prob of cFmsy", cex = 0.75, col="cyan",pos=4)
  text(0.01,0.95, "Prob of lFmsy", cex = 0.75, col="brown",pos=4)
  
  lines(c(Flim,Flim), c(0.0,0.05), lty = 2, col = "red")
  lines(c(0,Flim), c(0.05,0.05), lty = 2, col = "red")
  text(x = 0.1, y = 0.075, "5%", cex = 0.75, col = "red")
  
  #lines(rep(vcum,2), c(0,1), lty = 1, col = 4)
  #lines(c(Fscan[maxcatm],Fscan[maxcatm]), c(0,1), col = 5)
  
  
  #out <- c(Blim, Bpa)
  #if (yield =="landings") {
  #  outF <- c(flim, flim10, flim50, vcum, Fscan[maxlanm], FCrash5, FCrash50)
  #  outC <- approx(Fscan, lans[4,], xout = outF) $ y
  #}
  #if (yield =="catch") {
  #  outF <- c(flim, flim10, flim50, vcum, Fscan[maxcatm], FCrash5, FCrash50)
  #  outC <- approx(Fscan, cats[4,], xout = outF) $ y
  #}
  #outB <- approx(Fscan, ssbs[4,], xout = outF) $ y
  #outTable <- rbind(outF, outB, outC)
  #rownames(outTable) <- c("F","SSB",yield)
  #colnames(outTable) <- c("Flim","Flim10","Flim50","MSY:median","Maxmeanland","FCrash5","FCrash50")
  
  #list(Blim = Blim, Bpa = Bpa, Refs = outTable)
}  # end of eqsim_plot





#' @title Plots of the results from eqsim
#'
#' @description XXX
#' 
#' @author Einar Hjorleifsson \email{einar.hjorleifsson@@gmail.com}
#' 
#' @export
#' 
#' @param sim An object returned from the function eqsim_run
#' @param Scale A value, the scaling on the yaxis
#' @param plotit Boolean, if TRUE (default) returns a plot

eqsim_ggplot <- function(sim, Scale=1e3, plotit=TRUE) 
{
  
  # dummy
  Ftarget <- p05 <- p95 <- p50 <- variable <- value <- year <- NULL
  
  rby <- sim$rby
  for (i in c(1,2,4,5)) rby[,i] <- rby[,i]/Scale
  
  rbp <- sim$rbp
  for(i in 3:9) rbp[,i] <- rbp[,i]/Scale
  
  refs <- sim$Refs
  refs[3:6,] <- refs[3:6,]/Scale
  
  pProfile <- sim$pProfile
  
  i <- rbp$variable %in% "Recruitment"
  plotR <- ggplot(rbp[i,],aes(Ftarget)) + 
    geom_ribbon(aes(ymin=p05,ymax=p95),fill="grey90") +
    geom_line(aes(y=p50)) + 
    geom_vline(xintercept=refs[1,4],col="darkgreen",lwd=1) +
    geom_vline(xintercept=refs[1,5],col="yellow",lwd=1) +
    geom_vline(xintercept=refs[1,1],col="red",lwd=1) +
    geom_vline(xintercept=refs[2,4],col="darkgreen",linetype=2) +
    geom_vline(xintercept=refs[2,5],col="yellow",linetype=2) +
    facet_wrap(~ variable) +
    labs(y = "",x="") +
    geom_point(data=rby,aes(fbar,rec)) +
    coord_cartesian(ylim=c(0,rby$rec * 1.2),xlim=c(0,rby$fbar * 1.2)) 
  
  
  i <- rbp$variable %in% "Spawning stock biomass"
  plotSSB <- ggplot(rbp[i,],aes(Ftarget)) + 
    geom_ribbon(aes(ymin=p05,ymax=p95),fill="grey90") +
    geom_line(aes(y=p50)) + 
    geom_vline(xintercept=refs[1,4],col="darkgreen",lwd=1) +
    geom_vline(xintercept=refs[1,5],col="yellow",lwd=1) +
    geom_vline(xintercept=refs[1,1],col="red",lwd=1) +
    geom_vline(xintercept=refs[2,4],col="darkgreen",linetype=2) +
    geom_vline(xintercept=refs[2,5],col="yellow",linetype=2) +
    geom_point(data=rby,aes(fbar,ssb)) +
    facet_wrap(~ variable) +
    coord_cartesian(ylim=c(0,rby$ssb * 1.2),xlim=c(0,rby$fbar * 1.2)) +
    labs(y = "",x="")
  
  i <- rbp$variable %in% "Catch"
  plotCatch <- ggplot(rbp[i,],aes(Ftarget)) + 
    geom_ribbon(aes(ymin=p05,ymax=p95),fill="grey90") +
    geom_line(aes(y=p50)) + 
    geom_vline(xintercept=refs[1,4],col="darkgreen",lwd=1) +
    geom_vline(xintercept=refs[1,5],col="yellow",lwd=1) +
    geom_vline(xintercept=refs[1,1],col="red",lwd=1) +
    geom_vline(xintercept=refs[2,4],col="darkgreen",linetype=2) +
    geom_vline(xintercept=refs[2,5],col="yellow",linetype=2) +
    geom_point(data=rby,aes(fbar,catch)) +
    facet_wrap(~ variable) +
    coord_cartesian(ylim=c(0,rby$catch * 1.2),xlim=c(0,rby$fbar * 1.2)) +
    labs(y = "",x="")
  
  i <- rbp$variable %in% "Landings"
  plotLandings <- ggplot(rbp[i,],aes(Ftarget)) + 
    geom_ribbon(aes(ymin=p05,ymax=p95),fill="grey90") +
    geom_line(aes(y=p50)) + 
    geom_vline(xintercept=refs[1,4],col="darkgreen",lwd=1) +
    geom_vline(xintercept=refs[1,5],col="yellow",lwd=1) +
    geom_vline(xintercept=refs[1,1],col="red",lwd=1) +
    geom_vline(xintercept=refs[2,4],col="darkgreen",linetype=2) +
    geom_vline(xintercept=refs[2,5],col="yellow",linetype=2) +
    geom_point(data=rby,aes(fbar,landings)) +
    facet_wrap(~ variable) +
    coord_cartesian(ylim=c(0,rby$landings * 1.2),xlim=c(0,rby$fbar * 1.2)) +
    labs(y = "",x="")
  
  
  d2 <- rby[,c("fbar","catch","landings")]
  names(d2) <- c("Ftarget","Catch","Landings")
  d2 <- melt(d2,id.vars="Ftarget")
  d2$dummy <- "Yield"
  d2$Ftarget <- as.numeric(d2$Ftarget)
  
  i <- rbp$variable %in% c("Catch","Landings")
  d1 <- rbp[i,]
  d1$dummy <- "Yield"
  plotYield <- ggplot(d1,aes(Ftarget)) + 
    geom_ribbon(aes(ymin=p05,ymax=p95,fill=variable),alpha=0.2) +
    geom_line(aes(y=p05,colour=variable)) + 
    geom_line(aes(y=p95,colour=variable)) + 
    geom_line(aes(y=p50,colour=variable)) + 
    geom_vline(xintercept=refs[1,4],col="darkgreen",lwd=1) +
    geom_vline(xintercept=refs[1,5],col="yellow",lwd=1) +
    geom_vline(xintercept=refs[1,1],col="red",lwd=1) +
    geom_vline(xintercept=refs[2,4],col="darkgreen",linetype=2) +
    geom_vline(xintercept=refs[2,5],col="yellow",linetype=2) +
    geom_point(data=d2,aes(Ftarget,value,colour=variable)) +
    facet_wrap(~ dummy) +
    coord_cartesian(ylim=c(0,max(rby$catch) * 1.2),xlim=c(0,max(rby$fbar) * 1.2)) +
    labs(y = "",x="") +
    theme(legend.position=c(0.9,0.85))
  
  
  
  pProfile$dummy <- "Probability plot"
  plotProbs <- 
    ggplot(pProfile,aes(Ftarget,value,colour=variable)) + 
    geom_line() + 
    geom_hline(yintercept=0.05,colour="red") +
    coord_cartesian(xlim=c(0,rby$fbar * 1.2)) +
    labs(x="",y="") + facet_wrap(~ dummy) +
    theme(legend.position=c(0.1,0.80),
          legend.text=element_text(size = rel(0.5)),
          legend.title=element_text(size=rel(0.5)))
  
  if(plotit) {
    grid.arrange(plotR,plotSSB,plotYield,plotProbs, ncol=2)
  }
  
  return(list(plotR=plotR,plotSSB=plotSSB,plotCatch=plotCatch,
              plotLandings=plotLandings,plotProbs=plotProbs))
}