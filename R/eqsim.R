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
#' @param bio.const A flag, if FALSE mean of the biological values from the years selected are used
#' @param sel.years The years to sample the selection patterns from
#' @param sel.const A flag, if FALSE mean of the selection patterns from the years selected are used
#' @param Fscan F values to scan over
#' @param Fcv Assessment error in the advisory year
#' @param Fphi Autocorrelation in assessment error in the advisory year
#' @param Blim This we know
#' @param Bpa This we know
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
#' @param extreme.trim Call John Simmonds :-)

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
                      recruitment.trim = c(3, -3),
                      Btrigger = 0,
                      Nrun = 200, # number of years to run in total
                      process.error = TRUE, # use predictive recruitment or mean recruitment? (TRUE = predictive)
                      verbose = TRUE,
                      extreme.trim)
{
  
  if (abs(Fphi) >= 1) stop("Fphi, the autocorelation parameter for log F should be between (-1, 1)")
  if ((recruitment.trim[1] + recruitment.trim[2])> 0) stop("recruitment truncation must be between a high - low range")
  
  btyr1 <- bio.years[1]
  btyr2 <- bio.years[2]
  slyr1 <- sel.years[1]
  slyr2 <- sel.years[2]
  # Keep at most 50 simulation years (which will be the last 50 of the Nrun 
  #  forward simulated years)
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
  
  # 22.2.2014 Added weight of landings per comment from Carmen
  if (bio.const==TRUE){ # take means of wts Mat and M and ratio of landings to catch
    west[] <- apply(west, 1, mean)
    weca[] <- apply(weca, 1, mean)
    wela[] <- apply(wela, 1, mean)
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
  
  ssby <- Ferr <- array(0, c(Nrun,Nmod),dimnames=list(year=1:Nrun,iter=1:Nmod))
  Ny <- Fy <- WSy <- WCy <- Cy <- Wy <- Wl <- Ry <- array(0, c(ages, Nrun, Nmod),
                                                          dimnames=list(age=(range(stk)[1]:range(stk)[2]),year=1:Nrun,iter=1:Nmod))
  # TODO per note from Carmen:
  #  NOTE: If we want Ferr to be a stationary AR(1) process, it would make
  #        more sense to initialise Ferr as a Normal dist with zero mean and
  #        standard deviation of AR(1) marginal distribution, i.e. standard 
  #        deviation of initial Ferr = Fcv/sqrt(1- Fphi^2), instead of just
  #        initialising Ferr=0
  
  rsam <- array(sample(1:ncol(weca), Nrun * Nmod, TRUE), c(Nrun, Nmod))
  rsamsel <- array(sample(1:ncol(sel), Nrun * Nmod, TRUE), c(Nrun, Nmod))
  Wy[] <- c(weca[, c(rsam)])
  Wl[] <- c(wela[, c(rsam)])
  Ry[]  <- c(land.cat[, c(rsamsel)])
  # initial recruitment
  R <- mean( data $ rec)
  ssbs <- cats <- lans <- recs <- array(0, c(7, NF))
  
  ferr <- ssbsa <- catsa <- lansa <- recsa <- array(0, c(NF, keep, Nmod))
  begin <- Nrun - keep + 1
  
  # New from Simmonds' 29.1.2014
  #   Residuals of SR fits (1 value per SR fit and per simulation year 
  #     but the same residual value for all Fscan values):
  resids= array(rnorm(Nmod*(Nrun+1), 0, SR$cv),c(Nmod, Nrun+1))
  
  # Limit how extreme the Rec residuals can get:
  lims = t(array(SR$cv,c(Nmod,2))) * recruitment.trim
  for (k in 1:Nmod) { resids[k,resids[k,]>lims[1,k]]=lims[1,k]}
  for (k in 1:Nmod) { resids[k,resids[k,]<lims[2,k]]=lims[2,k]}
  # end New from Simmonds 29.1.2014
  
  if (verbose) loader(0)
  
  # Looping over each F value in Fscan. For each of the Nmod SR fits 
  # (replicates), do a forward simulation during Nrun years
  # There are Rec residuals for each SR fit and year, which take the same
  # values for all Fscan 
  for (i in 1:NF) {
    
    # The F value to test
    Fbar <- Fscan[i]
    
    ############################################################################
    # Population in simulation year 1:
    
    # Zpre: Z that occurs before spawning
    Zpre <- ( sel[,rsamsel[1,]]*Fbar * Fprop + M[,rsam[1,]] * Mprop)
    
    # Zpos: Z that occurs after spawning
    # Zpos not used anywhere
    Zpos <- (Fbar * (1-Fprop) * sel[,rsamsel[1,]] + M[,rsam[1,]] * (1-Mprop))
    
    # run Z out to age 50 ...
    # TODO:
    # Comments from Carmen: Zcum is a cumulative sum, but it is done in a strange way:
    #  There is a matrix of F-at-age and a matrix of M-at-age (each has 49 ages, Nmod replicates)
    #  The F and M matrices are summed, giving Z-at-age (49 ages, Nmod replicates)
    #  But then a cumsum is taken considering the Z-at-age matrix as a vector (i.e. not column-wise) ????
    #  This is strange, by applying "cumsum" treating Z-at-age as a vector, really only the first 50 values of
    #  the resulting "Zcum" make sense (all other values seem "wrong", or at least, meaningless)
    Zcum <- c(0, cumsum(Fbar * sel[c(1:ages, rep(ages, 49 - ages)), rsamsel[1,]] + M[c(1:ages, rep(ages, 49 - ages)), rsam[1,]]))
    # Carmen: Following from "Zcum", only first 50 elements of N1 make sense ????
    N1 <- R * exp(- unname(Zcum))
    
    # set up age structure in first year for all simulations
    # Comments from Carmen:
    #   Ny has dimension = (no. ages, no. simulation yrs "Nrun", no. SR fits "Nmod")
    #   With this code, we seem to be getting always the same population-at-age value for year 1
    #   instead of Nmod different values, as might have been intended ????
    #   (the whole problem is coming from Zcum ==> N1 ==> Ny[,1,] )
    Ny[,1,] <- c(N1[1:(ages-1)], sum(N1[ages:50]))
    
    # calculate ssb in first year using a different stock.wt and Mat selection and M for each simulation
    # Comments from Carmen:
    #   ssby has dimension = (no. simul yrs "Nrun", no. SR fits "Nmod")
    #   SSB in year 1:
    #   although Ny[,1,] has dim no.ages x Nmod, all Nmod values of Ny[,1,] are
    #   the same (because of Zcum issue)
    ssby[1,] <- colSums(Mat[,rsam[1,]] * Ny[,1,] * west[,rsam[1,]] / exp(Zpre)[])
    
    # Years 2 to Nrun:
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
      
      # Comment from Carmen:
      #  For each of the Nmod replicates, this selects the appropriate SR model
      #   type to use in that replicate
      #  Note that the order of SR model types that comes out in "select" is
      #   not necessarily the same order in which the SR model types were
      #   entered as inputs -- I presume the **next 2 lines** of code have
      #   been checked to avoid potential bugs due to this reordering  ???? 
      select <- cbind(seq(Nmod), as.numeric(factor(SR $ mod, levels = unique(SR $ mod))))
      Ny[1,j,] <- allrecs[select]
      
      # Comment from Carmen:
      #   Note: it seems that Rec is coded as occurring always at age 1
      #   (i.e. based on SSB in previous year)
      #   Some stocks have Rec at ages other than 1 (e.g. age 0) 
      #    -- is this a problem ????
      
      # apply HCR
      # (intended) Fbar to be applied in year j-1 (depends on SSB in year j-1):
      Fnext <- Fbar * pmin(1, SSB/Btrigger)
      
      # apply some noise to the F
      # Notes from Carmen:
      #  Assessment and/or implementation error (modifies intended F to get
      #  realised F)
      #  Error: AR(1) process on log(F) with autocorrel = Fphi, and
      #  conditional stand deviation = Fcv
      #  Might make more sense to have the "Ferr" matrix calculated before
      #  the Fscan loop starts so that the same errors in F are applied to
      #  all Fscan values ???? (as for Rec residuals)
      Ferr[j,] <- Fphi * Ferr[j-1,] + rnorm(Nmod, 0, Fcv)
      
      # realised Fbar in year j-1:
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
  
  
  # Einar amended 30.1.2014
  if(missing(extreme.trim)) {
    catm <- apply(catsa, 1, mean)
    lanm <- apply(lansa, 1, mean)
  } else {
    x <- catsa
    i <- x > quantile(x,extreme.trim[2]) |
      x < quantile(x,extreme.trim[1])
    x[i] <- NA
    catm <- apply(x, 1, mean, na.rm=TRUE)
    
    x <- lansa
    i <- x > quantile(x,extreme.trim[2]) |
      x < quantile(x,extreme.trim[1])
    x[i] <- NA
    lanm <- apply(x, 1, mean, na.rm=TRUE)
  }
  
  # end Einar amended 30.1.2014
    
  maxcatm <- which.max(catm)
  maxlanm <- which.max(lanm)
  
  # Einar added 29.1.2014
  rbp$Mean <- NA
  rbp$Mean[rbp$variable == "Catch"] <- catm
  rbp$Mean[rbp$variable == "Landings"] <- lanm
  # end Einar added 29.1.2014
  
  
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
  colnames(Refs) <- c("F05","F10","F50","medianMSY","meanMSY","FCrash05","FCrash50")
  
  #TODO: id.sim - user specified.
  
  return(list(ibya=list(Mat = Mat, M = M, Fprop = Fprop, Mprop = Mprop, 
                        west = west, weca = weca, sel = sel),
              rbya=list(ferr=ferr),
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
    text(0.98*Flim,0,"F05",srt=90,pos=3,col="red",cex=0.75)
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
    # Einar added 29.1.2014
    with(rbp[i,],lines(Ftarget,Mean, lty = 1, lwd=4,col="red"))    
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
    # Einar added 29.1.2014
    with(rbp[i,],lines(Ftarget,Mean, lty = 1, lwd=4,col="red"))   
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

eqsim_ggplot <- function(sim, Scale=1, plotit=TRUE) 
{
  
  # dummy
  Ftarget <- p05 <- p95 <- p50 <- variable <- value <- year <- NULL
  
  rby <- sim$rby
  for (i in c(2,3,5,6)) rby[,i] <- rby[,i]/Scale
  
  rbp <- sim$rbp
  for(i in 3:9) rbp[,i] <- rbp[,i]/Scale
  
  refs <- sim$Refs
  refs[3:6,] <- refs[3:6,]/Scale
  
  sim$Blim <- sim$Blim/Scale
  sim$Bpa <- sim$Bpa/Scale
  pProfile <- sim$pProfile
  
  i <- rbp$variable %in% "Recruitment"
  plotR <- ggplot(rbp[i,],aes(Ftarget)) + 
    theme_bw() +
    geom_ribbon(aes(ymin=p05,ymax=p95),fill="grey90") +
    geom_line(aes(y=p50)) + 
    geom_line(aes(y=Mean),linetype=2) +
    geom_vline(xintercept=refs[1,1],col="red",lwd=1) +
    annotate("text",x=refs[1,1],y=0,label="F05",col="red",hjust=0,vjust=0,angle=90) +
    geom_vline(xintercept=refs[1,5],col="darkgreen",lwd=1) +
    annotate("text",x=refs[1,5],y=0,label="Fmsy",col="darkgreen",hjust=0,vjust=0,angle=90) +
    facet_wrap(~ variable) +
    labs(y = "",x="") +
    geom_point(data=rby,aes(fbar,rec)) +
    coord_cartesian(ylim=c(0,rby$rec * 1.2),xlim=c(0,rby$fbar * 1.2)) 
  
  
  i <- rbp$variable %in% "Spawning stock biomass"
  plotSSB <- ggplot(rbp[i,],aes(Ftarget)) + 
    theme_bw() +
    geom_ribbon(aes(ymin=p05,ymax=p95),fill="grey90") +
    geom_line(aes(y=p50)) + 
    geom_hline(yintercept=sim$Blim,col="red",lwd=1) +
    annotate("text",x=0,y=sim$Blim,label="Blim",col="red",hjust=0,vjust=0) +
    geom_vline(xintercept=refs[1,1],col="red",lwd=1) +
    annotate("text",x=refs[1,1],y=0,label="F05",col="red",hjust=0,vjust=0,angle=90) +
    geom_vline(xintercept=refs[1,5],col="darkgreen",lwd=1) +
    annotate("text",x=refs[1,5],y=0,label="Fmsy",col="darkgreen",hjust=0,vjust=0,angle=90) +
    geom_point(data=rby,aes(fbar,ssb)) +
    facet_wrap(~ variable) +
    coord_cartesian(ylim=c(0,rby$ssb * 1.2),xlim=c(0,rby$fbar * 1.2)) +
    labs(y = "",x="")
  
  i <- rbp$variable %in% "Catch"
  plotCatch <- ggplot(rbp[i,],aes(Ftarget)) + 
    theme_bw() +
    geom_ribbon(aes(ymin=p05,ymax=p95),fill="grey90") +
    geom_line(aes(y=p50)) + 
    geom_line(aes(y=Mean),linetype=2) +
    geom_vline(xintercept=refs[1,1],col="red",lwd=1) +
    annotate("text",x=refs[1,1],y=0,label="F05",col="red",hjust=0,vjust=0,angle=90) +
    geom_vline(xintercept=refs[1,5],col="darkgreen",lwd=1) +
    annotate("text",x=refs[1,5],y=0,label="Fmsy",col="darkgreen",hjust=0,vjust=0,angle=90) +
    geom_point(data=rby,aes(fbar,catch)) +
    facet_wrap(~ variable) +
    coord_cartesian(ylim=c(0,rby$catch * 1.2),xlim=c(0,rby$fbar * 1.2)) +
    labs(y = "",x="")
  
  i <- rbp$variable %in% "Landings"
  plotLandings <- ggplot(rbp[i,],aes(Ftarget)) + 
    theme_bw() +
    geom_ribbon(aes(ymin=p05,ymax=p95),fill="grey90") +
    geom_line(aes(y=p50)) + 
    geom_line(aes(y=Mean),linetype=2) +
    geom_vline(xintercept=refs[1,1],col="red",lwd=1) +
    annotate("text",x=refs[1,1],y=0,label="F05",col="red",hjust=0,vjust=0,angle=90) +
    geom_vline(xintercept=refs[2,5],col="darkgreen",lwd=1) +
    annotate("text",x=refs[2,5],y=0,label="Fmsl",col="darkgreen",hjust=0,vjust=0,angle=90) +
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
    theme_bw() +
    geom_ribbon(aes(ymin=p05,ymax=p95,fill=variable),alpha=0.15) +
    #geom_line(aes(y=p05,colour=variable)) + 
    #geom_line(aes(y=p95,colour=variable)) + 
    geom_line(aes(y=p50,colour=variable)) + 
    geom_vline(xintercept=refs[1,1],col="red",lwd=1) +
    annotate("text",x=refs[1,1],y=0,label="F05",col="red",hjust=0,vjust=0,angle=90) +
    geom_vline(xintercept=refs[1,5],col="darkgreen",lwd=1) +
    annotate("text",x=refs[1,5],y=0,label="Fmsy",col="darkgreen",hjust=0,vjust=0,angle=90) +
    geom_vline(xintercept=refs[2,5],col="blue",lwd=1,linetype=2) +
    annotate("text",x=refs[2,5],y=0,label="Fmsl",col="blue",hjust=0,vjust=0,angle=90) +
    geom_point(data=d2,aes(Ftarget,value,colour=variable)) +
    facet_wrap(~ dummy) +
    coord_cartesian(ylim=c(0,max(rby$catch) * 1.2),xlim=c(0,max(rby$fbar) * 1.2)) +
    labs(y = "",x="") +
    theme(legend.position="none") +
    scale_colour_manual(values=c("Catch"="darkgreen","Landings"="blue")) +
    scale_fill_manual(values=c("Catch"="darkgreen","Landings"="blue")) +
    annotate("text",x=max(rby$fbar) * 1.2,y=max(rby$catch) * 1.1,label="Catch",colour="darkgreen",hjust=1) +
    annotate("text",x=max(rby$fbar) * 1.2,y=max(rby$landings) * 1.1,label="Landings",colour="blue",hjust=1)
  
  
  
  pProfile$dummy <- "Probability plot"
  df <- data.frame(x=rep(max(rby$fbar),4),
                   y=c(0.80,0.75,0.70,0.65),
                   variable=c("p(SSB<Blim)","p(SSB<Bpa)","Fmsy","Fmsy - landings"))
  plotProbs <- 
    ggplot(pProfile,aes(Ftarget,value,colour=variable)) + 
    scale_colour_manual(values=c("pFmsyCatch"="darkgreen",
                                 "pFmsyLandings"="blue",
                                 "Blim"="red",
                                 "Bpa"="orange",
                                 "p(SSB<Blim)"="red",
                                 "p(SSB<Bpa)"="orange",
                                 "Fmsy"="darkgreen",
                                 "Fmsy - landings"="blue")) + 
    theme_bw() +
    geom_line(lwd=1) + 
    geom_text(data=df,aes(x,y,label=variable,colour=variable)) +
    geom_hline(yintercept=0.05,colour="black") +
    coord_cartesian(xlim=c(0,rby$fbar * 1.2)) +
    labs(x="",y="") + facet_wrap(~ dummy) +
    theme(legend.position="none") 

  if(plotit) {
    vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(2, 2)))
    print(plotSSB + theme(panel.margin = unit(c(0,0,0,0), "cm"),
                          plot.margin = unit(c(0,0.25,0,0), "cm")),
          vp = vplayout(1, 1))
    print(plotR + theme(panel.margin = unit(c(0,0,0,0), "cm"),
                        plot.margin = unit(c(0,0.25,0,0), "cm")),
          vp=vplayout(1,2))
    print(plotYield + theme(panel.margin = unit(c(0,0,0,0), "cm"),
                            plot.margin = unit(c(0,0.25,0,0), "cm")),
          vp=vplayout(2,1))
    print(plotProbs + theme(panel.margin = unit(c(0,0,0,0), "cm"),
                            plot.margin = unit(c(0,0.25,0,0), "cm")),
          vp=vplayout(2,2))
  } else {
    return(list(plotR=plotR,plotSSB=plotSSB,plotCatch=plotCatch,
              plotLandings=plotLandings,plotYield=plotYield,plotProbs=plotProbs))
  }
}