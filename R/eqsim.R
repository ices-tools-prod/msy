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
#' @param SSBcv Spawning stock biomass error in the advisory year
#' @param rhologRec A flag for recruitment autocorrelation. If FALSE (default) 
#' then not applied.
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
                      SSBcv = 0,
                      rhologRec = FALSE,
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
  stk.win <- FLCore::window(stk, start = btyr1, end = btyr2)
  stk.winsel <- FLCore::window(stk, start = slyr1  , end = slyr2)
  
  littleHelper <- function(x,i) {
    x2 <- x
    x2[i] <- NA
    x2[] <- apply(x2,1,mean,na.rm=TRUE)
    x[i] <- x2[i]
    return(x)
  }
  
  west <- matrix(FLCore::stock.wt(stk.win), ncol = btyr2 - btyr1 + 1)
  i <- west == 0
  if(any(i)) west <- littleHelper(west,i)
  weca <- matrix(FLCore::catch.wt(stk.win), ncol = btyr2 - btyr1 + 1)
  i <- weca == 0
  if(any(i)) weca <- littleHelper(weca,i)
  wela <- matrix(FLCore::landings.wt(stk.win), ncol = btyr2 - btyr1 + 1)
  if(any(i)) wela <- littleHelper(wela,i)
  
  Mat <- matrix(FLCore::mat(stk.win), ncol = btyr2 - btyr1 + 1)
  M <- matrix(FLCore::m(stk.win), ncol = btyr2 - btyr1 + 1)
  landings <- matrix(FLCore::landings.n(stk.winsel), ncol = slyr2 - slyr1 + 1)
  # if zero, use 0.10 of minimum value
  
  catch <- matrix(FLCore::catch.n(stk.winsel), ncol = slyr2 - slyr1 + 1)
  sel <- matrix(FLCore::harvest(stk.winsel), ncol = slyr2 - slyr1 + 1)
  Fbar <- matrix(FLCore::fbar(stk.winsel), ncol = slyr2 - slyr1  + 1)
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
  
  Fprop <- apply(FLCore::harvest.spwn(stk.winsel), 1, mean)[drop=TRUE] # vmean(harvest.spwn(stk.win))
  Mprop <- apply(FLCore::m.spwn(stk.win), 1, mean)[drop=TRUE] # mean(m.spwn(stk.win))
  
  # get ready for the simulations
  Nmod <- nrow(SR)
  NF <- length(Fscan)
  ages <- FLCore::dims(stk) $ age
  
  ssby <- Ferr <- array(0, c(Nrun,Nmod),dimnames=list(year=1:Nrun,iter=1:Nmod))
  Ny <- Fy <- WSy <- WCy <- Cy <- Wy <- Wl <- Ry <- array(0, c(ages, Nrun, Nmod),
                                                          dimnames=list(age=(range(stk)[1]:range(stk)[2]),year=1:Nrun,iter=1:Nmod))
  # TODO per note from Carmen:
  #  NOTE: If we want Ferr to be a stationary AR(1) process, it would make
  #        more sense to initialise Ferr as a Normal dist with zero mean and
  #        standard deviation of AR(1) marginal distribution, i.e. standard 
  #        deviation of initial Ferr = Fcv/sqrt(1- Fphi^2), instead of just
  #        initialising Ferr=0
  #  2014-03-12: Changed per note form Carmen/John
  Ferr[1,] <- rnorm(n=Nmod, mean=0, sd=1)*Fcv/sqrt(1-Fphi^2)
  for(j in 2:Nrun) { Ferr[j,] <- Fphi*Ferr[j-1,] + Fcv*rnorm(n=Nmod, mean=0, sd=1) }
  
  # 2014-03-12: Changed per note form Carmen/John
  #  Errors in SSB: this is used when the ICES MSY HCR is applied for F
  SSBerr <- matrix(rnorm(n=Nrun*Nmod, mean=0, sd=1), ncol=Nmod) * SSBcv
  
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
  
  # 2014-03-12: Changed per note form Carmen/John
  #  Autocorrelation in Recruitment Residuals:    
  if(rhologRec){
    fittedlogRec <-  do.call(cbind, lapply( c(1:nrow(fit$sr.sto)), function(i){     
      FUN <- match.fun(fit$sr.sto$model[i])
      FUN(fit$sr.sto[i, ], fit$rby$ssb) } )  ) 
    # Calculate lag 1 autocorrelation of residuals: 
    rhologRec <- apply(log(fit$rby$rec)-fittedlogRec, 2, function(x){cor(x[-length(x)],x[-1])})
    # Draw residuals according to AR(1) process:
    for(j in 2:(Nrun+1)){ resids[,j] <- rhologRec * resids[,j-1] + resids[,j]*sqrt(1 - rhologRec^2) }
  }    
  
  
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
        ###### 2014-03-13  TMP COMMENT - ERROR OCCURS HERE
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
      # 2014-03-12: Changed per note form Carmen/John
      # Fnext <- Fbar * pmin(1, SSB/Btrigger)
      Fnext <- Fbar * pmin(1, SSB * exp(SSBerr[j,]) / Btrigger)      
      
      # apply some noise to the F
      # Notes from Carmen:
      #  Assessment and/or implementation error (modifies intended F to get
      #  realised F)
      #  Error: AR(1) process on log(F) with autocorrel = Fphi, and
      #  conditional stand deviation = Fcv
      #  Might make more sense to have the "Ferr" matrix calculated before
      #  the Fscan loop starts so that the same errors in F are applied to
      #  all Fscan values ???? (as for Rec residuals)
      
      # Outcommented 2014-03-12 because F-error already been drawn outside the
      #   loop, so this line here is no longer needed:
      # Ferr[j,] <- Fphi * Ferr[j-1,] + rnorm(Nmod, 0, Fcv)
      
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
    
    # 2014-03-12 Outcommented per note from Carmen/John - see below
    #x <- catsa
    #i <- x > quantile(x,extreme.trim[2]) |
    #  x < quantile(x,extreme.trim[1])
    #x[i] <- NA
    #catm <- apply(x, 1, mean, na.rm=TRUE)
    #
    #x <- lansa
    #i <- x > quantile(x,extreme.trim[2]) |
    #  x < quantile(x,extreme.trim[1])
    #x[i] <- NA
    #lanm <- apply(x, 1, mean, na.rm=TRUE)
    
    # 2014-03-12: Above replaced with the following per note from Carmen/John
    #  If we want to remove whole SR models, we could use the following code. But it is too extreme, it ends up getting rid of most models:
    # auxi2 <- array( apply(catsa, 1, function(x){auxi<-rep(TRUE,Nmod); auxi[x > quantile(x, extreme.trim[2]) | x < quantile(x, extreme.trim[1])] <- FALSE; x <- auxi } ), dim=c(keep,Nmod,NF))
    # auxi2 <- (1:Nmod)[apply(auxi2, 2, function(x){length(unique(as.vector(x)))})==1]
    # apply(catsa[,,auxi2],1,mean)
    
    # So I think the alternative is not to get rid of whole SR models, but of different SR models depending on the value of F:
    catm <- apply(catsa, 1, function(x){mean(x[x <= quantile(x, extreme.trim[2]) & x >= quantile(x, extreme.trim[1])])})
    lanm <- apply(lansa, 1, function(x){mean(x[x <= quantile(x, extreme.trim[2]) & x >= quantile(x, extreme.trim[1])])})
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
  
  # 2014-03-12 Ammendments per note from Carmen/John
  # CALCULATIONS:
  
  # Fmsy: value that maximises median LT catch or median LT landings 
  auxi <- approx(Fscan, cats[4, ],xout=seq(min(Fscan),max(Fscan),length=200))
  FmsyMedianC <- auxi$x[which.max(auxi$y)]   
  MSYMedianC <- max(auxi$y)
  # Value of F that corresponds to 0.95*MSY:
  FmsylowerMedianC <- auxi$x[ min( (1:length(auxi$y))[auxi$y/MSYMedianC >= 0.95] ) ]
  FmsyupperMedianC <- auxi$x[ max( (1:length(auxi$y))[auxi$y/MSYMedianC >= 0.95] ) ]
  
  auxi <- approx(Fscan, lans[4, ],xout=seq(min(Fscan),max(Fscan),length=200))
  FmsyMedianL <- auxi$x[which.max(auxi$y)]
  MSYMedianL <- max(auxi$y)
  
  # Value of F that corresponds to 0.95*MSY:
  FmsylowerMedianL <- auxi$x[ min( (1:length(auxi$y))[auxi$y/MSYMedianL >= 0.95] ) ]
  FmsyupperMedianL <- auxi$x[ max( (1:length(auxi$y))[auxi$y/MSYMedianL >= 0.95] ) ]
  
  F5percRiskBlim <- flim
  
  refs_interval <- data.frame(FmsyMedianC = FmsyMedianC,
                             FmsylowerMedianC = FmsylowerMedianC,
                             FmsyupperMedianC = FmsyupperMedianC,
                             FmsyMedianL = FmsyMedianL,
                             FmsylowerMedianL = FmsylowerMedianL,
                             FmsyupperMedianL = FmsyupperMedianL,
                             F5percRiskBlim = F5percRiskBlim,
                             Btrigger = Btrigger)
  
  # END 2014-03-12 Ammendments per note from Carmen/John
  
  sim <- list(ibya=list(Mat = Mat, M = M, Fprop = Fprop, Mprop = Mprop, 
                        west = west, weca = weca, sel = sel),
              rbya=list(ferr=ferr),
              rby=fit$rby,
              rbp=rbp,
              Blim=Blim,
              Bpa=Bpa,
              Refs = Refs,
              pProfile=pProfile,
              id.sim=fit$id.sr,
              refs_interval=refs_interval)
  
  sim <- eqsim_range(sim)
  
  return(sim)
  
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
  Ftarget <- p05 <- p95 <- p50 <- variable <- value <- year <- 
    Mean <- fbar <- rec <- ssb <- catch <- landings <- x <- y <- 0
  
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
  plotR <- 
    ggplot2::ggplot(rbp[i,],ggplot2::aes(Ftarget)) + 
    ggplot2::theme_bw() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=p05,ymax=p95),fill="grey90") +
    ggplot2::geom_line(ggplot2::aes(y=p50)) + 
    ggplot2::geom_line(ggplot2::aes(y=Mean),linetype=2) +
    ggplot2::geom_vline(xintercept=refs[1,1],col="red",lwd=1) +
    ggplot2::annotate("text",x=refs[1,1],y=0,label="F05",col="red",hjust=0,vjust=0,angle=90) +
    ggplot2::geom_vline(xintercept=refs[1,5],col="darkgreen",lwd=1) +
    ggplot2::annotate("text",x=refs[1,5],y=0,label="Fmsy",col="darkgreen",hjust=0,vjust=0,angle=90) +
    ggplot2::facet_wrap(~ variable) +
    ggplot2::labs(y = "",x="") +
    ggplot2::geom_point(data=rby,ggplot2::aes(fbar,rec)) +
    ggplot2::coord_cartesian(ylim=c(0,rby$rec * 1.2),xlim=c(0,rby$fbar * 1.2)) 
  
  
  i <- rbp$variable %in% "Spawning stock biomass"
  plotSSB <- 
    ggplot2::ggplot(rbp[i,],ggplot2::aes(Ftarget)) + 
    ggplot2::theme_bw() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=p05,ymax=p95),fill="grey90") +
    ggplot2::geom_line(ggplot2::aes(y=p50)) + 
    ggplot2::geom_hline(yintercept=sim$Blim,col="red",lwd=1) +
    ggplot2::annotate("text",x=0,y=sim$Blim,label="Blim",col="red",hjust=0,vjust=0) +
    ggplot2::geom_vline(xintercept=refs[1,1],col="red",lwd=1) +
    ggplot2::annotate("text",x=refs[1,1],y=0,label="F05",col="red",hjust=0,vjust=0,angle=90) +
    ggplot2::geom_vline(xintercept=refs[1,5],col="darkgreen",lwd=1) +
    ggplot2::annotate("text",x=refs[1,5],y=0,label="Fmsy",col="darkgreen",hjust=0,vjust=0,angle=90) +
    ggplot2::geom_point(data=rby,ggplot2::aes(fbar,ssb)) +
    ggplot2::facet_wrap(~ variable) +
    ggplot2::coord_cartesian(ylim=c(0,rby$ssb * 1.2),xlim=c(0,rby$fbar * 1.2)) +
    ggplot2::labs(y = "",x="")
  
  i <- rbp$variable %in% "Catch"
  plotCatch <- 
    ggplot2::ggplot(rbp[i,],ggplot2::aes(Ftarget)) + 
    ggplot2::theme_bw() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=p05,ymax=p95),fill="grey90") +
    ggplot2::geom_line(ggplot2::aes(y=p50)) + 
    ggplot2::geom_line(ggplot2::aes(y=Mean),linetype=2) +
    ggplot2::geom_vline(xintercept=refs[1,1],col="red",lwd=1) +
    ggplot2::annotate("text",x=refs[1,1],y=0,label="F05",col="red",hjust=0,vjust=0,angle=90) +
    ggplot2::geom_vline(xintercept=refs[1,5],col="darkgreen",lwd=1) +
    ggplot2::annotate("text",x=refs[1,5],y=0,label="Fmsy",col="darkgreen",hjust=0,vjust=0,angle=90) +
    ggplot2::geom_point(data=rby,ggplot2::aes(fbar,catch)) +
    ggplot2::facet_wrap(~ variable) +
    ggplot2::coord_cartesian(ylim=c(0,rby$catch * 1.2),xlim=c(0,rby$fbar * 1.2)) +
    ggplot2::labs(y = "",x="")
  
  i <- rbp$variable %in% "Landings"
  plotLandings <- 
    ggplot2::ggplot(rbp[i,],ggplot2::aes(Ftarget)) + 
    ggplot2::theme_bw() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=p05,ymax=p95),fill="grey90") +
    ggplot2::geom_line(ggplot2::aes(y=p50)) + 
    ggplot2::geom_line(ggplot2::aes(y=Mean),linetype=2) +
    ggplot2::geom_vline(xintercept=refs[1,1],col="red",lwd=1) +
    ggplot2::annotate("text",x=refs[1,1],y=0,label="F05",col="red",hjust=0,vjust=0,angle=90) +
    ggplot2::geom_vline(xintercept=refs[2,5],col="darkgreen",lwd=1) +
    ggplot2::annotate("text",x=refs[2,5],y=0,label="Fmsl",col="darkgreen",hjust=0,vjust=0,angle=90) +
    ggplot2::geom_point(data=rby,ggplot2::aes(fbar,landings)) +
    ggplot2::facet_wrap(~ variable) +
    ggplot2::coord_cartesian(ylim=c(0,rby$landings * 1.2),xlim=c(0,rby$fbar * 1.2)) +
    ggplot2::labs(y = "",x="")
  
  
  d2 <- rby[,c("fbar","catch","landings")]
  names(d2) <- c("Ftarget","Catch","Landings")
  d2 <- reshape2::melt(d2,id.vars="Ftarget")
  d2$dummy <- "Yield"
  d2$Ftarget <- as.numeric(d2$Ftarget)
  
  i <- rbp$variable %in% c("Catch","Landings")
  d1 <- rbp[i,]
  d1$dummy <- "Yield"
  plotYield <- 
    ggplot2::ggplot(d1,ggplot2::aes(Ftarget)) + 
    ggplot2::theme_bw() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=p05,ymax=p95,fill=variable),alpha=0.15) +
    #geom_line(ggplot2::aes(y=p05,colour=variable)) + 
    #geom_line(ggplot2::aes(y=p95,colour=variable)) + 
    ggplot2::geom_line(ggplot2::aes(y=p50,colour=variable)) + 
    ggplot2::geom_vline(xintercept=refs[1,1],col="red",lwd=1) +
    ggplot2::annotate("text",x=refs[1,1],y=0,label="F05",col="red",hjust=0,vjust=0,angle=90) +
    ggplot2::geom_vline(xintercept=refs[1,5],col="darkgreen",lwd=1) +
    ggplot2::annotate("text",x=refs[1,5],y=0,label="Fmsy",col="darkgreen",hjust=0,vjust=0,angle=90) +
    ggplot2::geom_vline(xintercept=refs[2,5],col="blue",lwd=1,linetype=2) +
    ggplot2::annotate("text",x=refs[2,5],y=0,label="Fmsl",col="blue",hjust=0,vjust=0,angle=90) +
    ggplot2::geom_point(data=d2,ggplot2::aes(Ftarget,value,colour=variable)) +
    ggplot2::facet_wrap(~ dummy) +
    ggplot2::coord_cartesian(ylim=c(0,max(rby$catch) * 1.2),xlim=c(0,max(rby$fbar) * 1.2)) +
    ggplot2::labs(y = "",x="") +
    ggplot2::theme(legend.position="none") +
    ggplot2::scale_colour_manual(values=c("Catch"="darkgreen","Landings"="blue")) +
    ggplot2::scale_fill_manual(values=c("Catch"="darkgreen","Landings"="blue")) +
    ggplot2::annotate("text",x=max(rby$fbar) * 1.2,y=max(rby$catch) * 1.1,label="Catch",colour="darkgreen",hjust=1) +
    ggplot2::annotate("text",x=max(rby$fbar) * 1.2,y=max(rby$landings) * 1.1,label="Landings",colour="blue",hjust=1)
  
  
  
  pProfile$dummy <- "Probability plot"
  df <- data.frame(x=rep(max(rby$fbar),4),
                   y=c(0.80,0.75,0.70,0.65),
                   variable=c("p(SSB<Blim)","p(SSB<Bpa)","Fmsy","Fmsy - landings"))
  plotProbs <- 
    ggplot2::ggplot(pProfile,ggplot2::aes(Ftarget,value,colour=variable)) + 
    ggplot2::scale_colour_manual(values=c("pFmsyCatch"="darkgreen",
                                 "pFmsyLandings"="blue",
                                 "Blim"="red",
                                 "Bpa"="orange",
                                 "p(SSB<Blim)"="red",
                                 "p(SSB<Bpa)"="orange",
                                 "Fmsy"="darkgreen",
                                 "Fmsy - landings"="blue")) + 
    ggplot2::theme_bw() +
    ggplot2::geom_line(lwd=1) + 
    ggplot2::geom_text(data=df,ggplot2::aes(x,y,label=variable,colour=variable)) +
    ggplot2::geom_hline(yintercept=0.05,colour="black") +
    ggplot2::coord_cartesian(xlim=c(0,rby$fbar * 1.2)) +
    ggplot2::labs(x="",y="") + 
    ggplot2::facet_wrap(~ dummy) +
    ggplot2::theme(legend.position="none") 

  if(plotit) {
    vplayout <- function(x, y) grid::viewport(layout.pos.row = x, layout.pos.col = y)
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 2)))
    print(plotSSB + ggplot2::theme(panel.margin = grid::unit(c(0,0,0,0), "cm"),
                          plot.margin = grid::unit(c(0,0.25,0,0), "cm")),
          vp = vplayout(1, 1))
    print(plotR + ggplot2::theme(panel.margin = grid::unit(c(0,0,0,0), "cm"),
                        plot.margin = grid::unit(c(0,0.25,0,0), "cm")),
          vp=vplayout(1,2))
    print(plotYield + ggplot2::theme(panel.margin = grid::unit(c(0,0,0,0), "cm"),
                            plot.margin = grid::unit(c(0,0.25,0,0), "cm")),
          vp=vplayout(2,1))
    print(plotProbs + ggplot2::theme(panel.margin = grid::unit(c(0,0,0,0), "cm"),
                            plot.margin = grid::unit(c(0,0.25,0,0), "cm")),
          vp=vplayout(2,2))
  } else {
    return(list(plotR=plotR,plotSSB=plotSSB,plotCatch=plotCatch,
              plotLandings=plotLandings,plotYield=plotYield,plotProbs=plotProbs))
  }
}


#' @title Calculate Fmsy range
#' 
#' @description XXX
#' 
#' @export
#' 
#' @param sim XXX
#' @param interval XXX
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
  x.lm <- loess(y.95 ~ x.95, span = 0.2)
  lm.pred <- data.frame(x = seq(min(x.95), max(x.95), length = 1000),
                        y = rep(NA, 1000))
  lm.pred$y <- predict(x.lm, newdata = lm.pred$x) 
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
  yield.f05 <- predict(x.lm, newdata = f05)
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
  x.lm <- loess(y.95 ~ x.95, span = 0.2)
  lm.pred <- data.frame(x = seq(min(x.95), max(x.95), length = 1000),
                        y = rep(NA, 1000))
  lm.pred$y <- predict(x.lm, newdata = lm.pred$x) 
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
  yield.f05 <- predict(x.lm, newdata = f05)
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
  b.lm <- loess(b.95 ~ x.95, span = 0.2)
  b.lm.pred <- data.frame(x = seq(min(x.95), max(x.95), length = 1000),
                          y = rep(NA, 1000))
  b.lm.pred$y <- predict(b.lm, newdata = b.lm.pred$x) 
  #lines(b.lm.pred$x, b.lm.pred$y, lty = 1, col = "red")
  
  # Estimate SSB for median F(msy) and range
  b.msymed <- predict(b.lm, newdata = Fmsymed)
  b.medlower <- predict(b.lm, newdata = fmsy.lower.median)
  b.medupper <- predict(b.lm, newdata = fmsy.upper.median)
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
  
  sim$Refs2["lanF","Medlower"] <- fmsy.lower.median
  sim$Refs2["lanF","Medupper"] <- fmsy.upper.median
  sim$Refs2["lanF","Meanlower"] <- fmsy.lower.mean
  sim$Refs2["lanF","Meanupper"] <- fmsy.upper.mean
  
  sim$Refs2["landings","Medlower"] <- landings.lower.median
  sim$Refs2["landings","Medupper"] <- landings.upper.median
  sim$Refs2["landings","Meanlower"] <- landings.lower.mean
  sim$Refs2["landings","Meanupper"] <- landings.upper.mean
  
  sim$Refs2["lanB","medianMSY"] <- b.msymed
  sim$Refs2["lanB","Medlower"] <- b.medlower
  sim$Refs2["lanB","Medupper"] <- b.medupper
  
  # Reference point estimates
  #cat("\nReference point estimates:\n")
  return (sim)
}


#' @title Calculate Fmsy range
#' 
#' @description XXX
#' 
#' @export
#' 
#' @param sim XXX
#' @param interval XXX
#' @param type XXX

eqsim_plot_range <- function (sim, interval=0.95, type="median")

  {
  
  data.95 <- sim$rbp
  
  if(type == "mean") {
    
    x.95 <- data.95[data.95$variable == "Landings",]$Ftarget
    y.95 <- data.95[data.95$variable == "Landings",]$Mean
    
    #x.95 <- x.95[2:length(x.95)]
    #y.95 <- y.95[2:length(y.95)]
    
    # Plot curve with 95% line
    # windows(width = 10, height = 7)
    par(mfrow = c(1,1), mar = c(5,4,2,1), mgp = c(3,1,0))
    plot(x.95, y.95, ylim = c(0, max(y.95, na.rm = TRUE)),
         xlab = "Total catch F", ylab = "Mean landings")
    yield.p95 <- interval * max(y.95, na.rm = TRUE)
    abline(h = yield.p95, col = "blue", lty = 1)
    
    # Fit loess smoother to curve
    x.lm <- loess(y.95 ~ x.95, span = 0.2)
    lm.pred <- data.frame(x = seq(min(x.95), max(x.95), length = 1000),
                          y = rep(NA, 1000))
    lm.pred$y <- predict(x.lm, newdata = lm.pred$x) 
    lines(lm.pred$x, lm.pred$y, lty = 1, col = "red")
    points(x = sim$Refs["lanF","meanMSY"], 
           y = predict(x.lm, newdata = sim$Refs["lanF","meanMSY"]),
           pch = 16, col = "blue")
    
    # Limit fitted curve to values greater than the 95% cutoff
    lm.pred.95 <- lm.pred[lm.pred$y >= yield.p95,]
    fmsy.lower <- min(lm.pred.95$x)
    fmsy.upper <- max(lm.pred.95$x)
    abline(v = c(fmsy.lower, fmsy.upper), lty = 8, col = "blue")
    abline(v = sim$Refs["lanF","meanMSY"], lty = 1, col = "blue")
    legend(x = "bottomright", bty = "n", cex = 1.0, 
           title = "F(msy)", title.col = "blue",
           legend = c(paste0("lower = ", round(fmsy.lower,3)),
                      paste0("mean = ", round(sim$Refs["lanF","meanMSY"],3)), 
                      paste0("upper = ", round(fmsy.upper,3))))
    
    fmsy.lower.mean <- fmsy.lower
    fmsy.upper.mean <- fmsy.upper
    landings.lower.mean <- lm.pred.95[lm.pred.95$x == fmsy.lower.mean,]$y	
    landings.upper.mean <- lm.pred.95[lm.pred.95$x == fmsy.upper.mean,]$y
    
    # Repeat for 95% of yield at F(05):
    f05 <- sim$Refs["catF","F05"]
    yield.f05 <- predict(x.lm, newdata = f05)
    points(f05, yield.f05, pch = 16, col = "green")
    yield.f05.95 <- interval * yield.f05
    abline(h = yield.f05.95, col = "green")
    lm.pred.f05.95 <- lm.pred[lm.pred$y >= yield.f05.95,]
    f05.lower <- min(lm.pred.f05.95$x)
    f05.upper <- max(lm.pred.f05.95$x)
    abline(v = c(f05.lower,f05.upper), lty = 8, col = "green")
    abline(v = f05, lty = 1, col = "green")
    legend(x = "right", bty = "n", cex = 1.0, 
           title = "F(5%)", title.col = "green",
           legend = c(paste0("lower = ", round(f05.lower,3)),
                      paste0("estimate = ", round(f05,3)),
                      paste0("upper = ", round(f05.upper,3))))
    
    return(invisible(NULL))
    
  }
  
  
  if(type == "median") {
  ################################################
  # Extract yield data (landings) - median version
  
  data.95 <- sim$rbp
  x.95 <- data.95[data.95$variable == "Landings",]$Ftarget
  y.95 <- data.95[data.95$variable == "Landings",]$p50
  
  
  # Plot curve with 95% line
  #windows(width = 10, height = 7)
  par(mfrow = c(1,1), mar = c(5,4,2,1), mgp = c(3,1,0))
  plot(x.95, y.95, ylim = c(0, max(y.95, na.rm = TRUE)),
       xlab = "Total catch F", ylab = "Median landings")
  yield.p95 <- interval * max(y.95, na.rm = TRUE)
  abline(h = yield.p95, col = "blue", lty = 1)
  
  # Fit loess smoother to curve
  x.lm <- loess(y.95 ~ x.95, span = 0.2)
  lm.pred <- data.frame(x = seq(min(x.95), max(x.95), length = 1000),
                        y = rep(NA, 1000))
  lm.pred$y <- predict(x.lm, newdata = lm.pred$x) 
  lines(lm.pred$x, lm.pred$y, lty = 1, col = "red")
  
  # Find maximum of fitted curve - this will be the new median (F(msy)
  Fmsymed <- lm.pred[which.max(lm.pred$y),]$x
  Fmsymed.landings <- lm.pred[which.max(lm.pred$y),]$y
  
  # Overwrite Refs table
  sim$Refs[,"medianMSY"] <- NA
  sim$Refs["lanF","medianMSY"] <- Fmsymed
  sim$Refs["landings","medianMSY"] <- Fmsymed.landings
  
  # Add maximum of medians to plot
  points(x = sim$Refs["lanF","medianMSY"], 
         y = predict(x.lm, newdata = sim$Refs["lanF","medianMSY"]),
         pch = 16, col = "blue")
  
  # Limit fitted curve to values greater than the 95% cutoff
  lm.pred.95 <- lm.pred[lm.pred$y >= yield.p95,]
  fmsy.lower <- min(lm.pred.95$x)
  fmsy.upper <- max(lm.pred.95$x)
  abline(v = c(fmsy.lower, fmsy.upper), lty = 8, col = "blue")
  abline(v = sim$Refs["lanF","medianMSY"], lty = 1, col = "blue")
  legend(x = "bottomright", bty = "n", cex = 1.0, 
         title = "F(msy)", title.col = "blue",
         legend = c(paste0("lower = ", round(fmsy.lower,3)),
                    paste0("median = ", round(sim$Refs["lanF","medianMSY"],3)), 
                    paste0("upper = ", round(fmsy.upper,3))))
  
  fmsy.lower.median <- fmsy.lower
  fmsy.upper.median <- fmsy.upper
  landings.lower.median <- lm.pred.95[lm.pred.95$x == fmsy.lower.median,]$y
  landings.upper.median <- lm.pred.95[lm.pred.95$x == fmsy.upper.median,]$y
  
  # Repeat for 95% of yield at F(05):
  f05 <- sim$Refs["catF","F05"]
  yield.f05 <- predict(x.lm, newdata = f05)
  points(f05, yield.f05, pch = 16, col = "green")
  yield.f05.95 <- interval * yield.f05
  abline(h = yield.f05.95, col = "green")
  lm.pred.f05.95 <- lm.pred[lm.pred$y >= yield.f05.95,]
  f05.lower <- min(lm.pred.f05.95$x)
  f05.upper <- max(lm.pred.f05.95$x)
  abline(v = c(f05.lower,f05.upper), lty = 8, col = "green")
  abline(v = f05, lty = 1, col = "green")
  legend(x = "right", bty = "n", cex = 1.0, 
         title = "F(5%)", title.col = "green",
         legend = c(paste0("lower = ", round(f05.lower,3)),
                    paste0("estimate = ", round(f05,3)),
                    paste0("upper = ", round(f05.upper,3))))
  
  return(invisible(NULL))
  
  }
  
  if(type == "ssb") {
    # Estimate implied SSB for each F output
    
    x.95 <- data.95[data.95$variable == "Spawning stock biomass",]$Ftarget
    b.95 <- data.95[data.95$variable == "Spawning stock biomass",]$p50
    
    # Plot curve with 95% line
    #windows(width = 10, height = 7)
    par(mfrow = c(1,1), mar = c(5,4,2,1), mgp = c(3,1,0))
    plot(x.95, b.95, ylim = c(0, max(b.95, na.rm = TRUE)),
         xlab = "Total catch F", ylab = "Median SSB")
    
    # Fit loess smoother to curve
    b.lm <- loess(b.95 ~ x.95, span = 0.2)
    b.lm.pred <- data.frame(x = seq(min(x.95), max(x.95), length = 1000),
                            y = rep(NA, 1000))
    b.lm.pred$y <- predict(b.lm, newdata = b.lm.pred$x) 
    lines(b.lm.pred$x, b.lm.pred$y, lty = 1, col = "red")
    
    # Estimate SSB for median F(msy) and range
    Fmsymed <- sim$Refs["landings","medianMSY"]
    fmsy.lower.median <- sim$Refs2["lanF","Medlower"]
    fmsy.upper.median <- sim$Refs2["lanF","Medupper"]
    
    b.msymed <- predict(b.lm, newdata = Fmsymed)
    b.medlower <- predict(b.lm, newdata = fmsy.lower.median)
    b.medupper <- predict(b.lm, newdata = fmsy.upper.median)
    
    abline(v = c(fmsy.lower.median, Fmsymed, fmsy.upper.median), col = "blue", lty = c(8,1,8))
    points(x = c(fmsy.lower.median, Fmsymed, fmsy.upper.median), 
           y = c(b.medlower, b.msymed, b.medupper), col = "blue", pch = 16)
    legend(x = "topright", bty = "n", cex = 1.0, 
           title = "F(msy)", title.col = "blue",
           legend = c(paste0("lower = ", round(b.medlower,0)),
                      paste0("median = ", round(b.msymed,0)),
                      paste0("upper = ", round(b.medupper,0))))
    
    return(invisible(NULL))
    
  }
  
  # Update summary table with John's format
  
  #sim$Refs <- sim$Refs[,!(colnames(sim$Refs) %in% c("FCrash05","FCrash50"))]
  #sim$Refs <- cbind(sim$Refs, Medlower = rep(NA,6), Meanlower = rep(NA,6), 
  #                     Medupper = rep(NA,6), Meanupper = rep(NA,6))
  
  #sim$Refs["lanF","Medlower"] <- fmsy.lower.median
  #sim$Refs["lanF","Medupper"] <- fmsy.upper.median
  #sim$Refs["lanF","Meanlower"] <- fmsy.lower.mean
  #sim$Refs["lanF","Meanupper"] <- fmsy.upper.mean
  
  #sim$Refs["landings","Medlower"] <- landings.lower.median
  #sim$Refs["landings","Medupper"] <- landings.upper.median
  #sim$Refs["landings","Meanlower"] <- landings.lower.mean
  #sim$Refs["landings","Meanupper"] <- landings.upper.mean
  
  #sim$Refs["lanB","medianMSY"] <- b.msymed
  #sim$Refs["lanB","Medlower"] <- b.medlower
  #sim$Refs["lanB","Medupper"] <- b.medupper
  
  # Reference point estimates
  #cat("\nReference point estimates:\n")
  #return(invisible(NULL))
  
  if(type == "carmen") {
    par(mfrow=c(2,2))
    
    # original from Carmen
    #auxi <- sim$rbp$p50[sim$rbp$variable=="Catch"]
    #plot(Fscan, auxi, type="l", main=paste("Median long-term catch, Btrigger=",Btrigger,sep=""),lwd=2,xlab="F",ylab="")
    #abline(v=FmsyMedianC, col=4,lwd=2)
    #abline(v=FmsylowerMedianC, col=4,lwd=2,lty=2)
    #abline(v=FmsyupperMedianC, col=4,lwd=2,lty=2)
    #abline(h=max(auxi)*0.95, col=4, lty=2)
    #abline(v=F5percRiskBlim, col=2,lwd=2)
    
    i <- sim$rbp$variable %in% "Catch"
    plot(sim$rbp$Ftarget[i], sim$rbp$p50[i], type="l",
         main=paste("Median long-term catch, Btrigger=",sim$refs_interval$Btrigger,sep=""),lwd=2,xlab="F",ylab="")
    abline(v=sim$refs_interval$FmsyMedianC, col=4,lwd=2)
    abline(v=sim$refs_interval$FmsylowerMedianC, col=4,lwd=2,lty=2)
    abline(v=sim$refs_interval$FmsyupperMedianC, col=4,lwd=2,lty=2)
    abline(h=max(sim$rbp$p50[i])*0.95, col=4, lty=2)
    abline(v=sim$refs_interval$F5percRiskBlim, col=2,lwd=2)
    
    # original from Carmen
    #auxi <- rbp$p50[rbp$variable=="Landings"]
    #plot(Fscan, auxi, type="l", main=paste("Median long-term landings, Btrigger=",Btrigger,sep=""),lwd=2,xlab="F",ylab="")
    #abline(v=FmsyMedianL, col=3,lwd=2)
    #abline(v=FmsylowerMedianL, col=3,lwd=2,lty=2)
    #abline(v=FmsyupperMedianL, col=3,lwd=2,lty=2)
    #abline(h=max(auxi)*0.95, col=3, lty=2)
    #abline(v=F5percRiskBlim, col=2,lwd=2)
    
    i <- sim$rbp$variable %in% "Landings"
    plot(sim$rbp$Ftarget[i], sim$rbp$p50[i], type="l",
         main=paste("Median long-term landings, Btrigger=",sim$refs_interval$Btrigger,sep=""),lwd=2,xlab="F",ylab="")
    abline(v=sim$refs_interval$FmsyMedianL, col=4,lwd=2)
    abline(v=sim$refs_interval$FmsylowerMedianL, col=4,lwd=2,lty=2)
    abline(v=sim$refs_interval$FmsyupperMedianL, col=4,lwd=2,lty=2)
    abline(h=max(sim$rbp$p50[i])*0.95, col=4, lty=2)
    abline(v=sim$refs_interval$F5percRiskBlim, col=2,lwd=2)
    
    # Carmen original
    #auxi <- approx(Fscan, rbp$p50[rbp$variable=="Spawning stock biomass"],xout=seq(min(Fscan),max(Fscan),length=200))
    #plot(auxi$x, auxi$y, type="l", main=paste("SSB: Median and 5th percentile, Btrigger=",Btrigger,sep=""),lwd=2,xlab="F",ylab="",ylim=c(0,max(auxi$y)))
    #abline(v=FmsyMedianL, col=3,lwd=2)
    #abline(v=FmsylowerMedianL, col=3,lwd=2,lty=2)
    #abline(v=FmsyupperMedianL, col=3,lwd=2,lty=2)
    #abline(v=F5percRiskBlim, col=2,lwd=2)
    #abline(v=0)
    
    #auxi <- approx(Fscan, rbp$p05[rbp$variable=="Spawning stock biomass"],xout=seq(min(Fscan),max(Fscan),length=200))
    #lines(auxi$x, auxi$y, lwd=2, col=2)
    #if(!missing(Blim)){abline(h=Blim, col=2, lty=2, lwd=2)}
    #abline(h=Bpa, col=4, lty=2, lwd=2)
    
    i <- sim$rbp$variable %in% "Spawning stock biomass"
    auxi <- approx(sim$rbp$Ftarget[i], sim$rbp$p50[i],xout=seq(min(sim$rbp$Ftarget[i]),max(sim$rbp$Ftarget[i]),length=200))
    plot(auxi$x, auxi$y, type="l", main=paste("SSB: Median and 5th percentile, Btrigger=",sim$refs_interval$Btrigger,sep=""),lwd=2,xlab="F",ylab="",ylim=c(0,max(auxi$y)))
    abline(v=sim$refs_interval$FmsyMedianL, col=3,lwd=2)
    abline(v=sim$refs_interval$FmsylowerMedianL, col=3,lwd=2,lty=2)
    abline(v=sim$refs_interval$FmsyupperMedianL, col=3,lwd=2,lty=2)
    abline(v=sim$refs_interval$F5percRiskBlim, col=2,lwd=2)
    abline(v=0)
    lines(auxi$x, auxi$y, lwd=2, col=2)
    if(!is.null(sim$Blim)) {abline(h=sim$Blim, col=2, lty=2, lwd=2)}
    
    # Carmen original
    #auxi <- approx(Fscan, rbp$p50[rbp$variable=="Recruitment"],xout=seq(min(Fscan),max(Fscan),length=200))
    #plot(auxi$x, auxi$y, type="l", main=paste("Rec: Median and 5th percentile, Btrigger=",Btrigger,sep=""),lwd=2,xlab="F",ylab="",ylim=c(0,max(auxi$y)))
    #abline(v=FmsyMedianL, col=3,lwd=2)
    #abline(v=FmsylowerMedianL, col=3,lwd=2,lty=2)
    #abline(v=FmsyupperMedianL, col=3,lwd=2,lty=2)
    #abline(v=F5percRiskBlim, col=2,lwd=2)
    #abline(v=0)
    
    #auxi <- approx(Fscan, rbp$p05[rbp$variable=="Recruitment"],xout=seq(min(Fscan),max(Fscan),length=200))
    #lines(auxi$x, auxi$y, lwd=2, col=2)
    
    i <- sim$rbp$variable %in% "Recruitment"
    auxi <- approx(sim$rbp$Ftarget[i], sim$rbp$p50[i],xout=seq(min(sim$rbp$Ftarget[i]),max(sim$rbp$Ftarget[i]),length=200))
    plot(auxi$x, auxi$y, type="l", main=paste("Rec: Median and 5th percentile, Btrigger=",sim$refs_interval$Btrigger,sep=""),lwd=2,xlab="F",ylab="",ylim=c(0,max(auxi$y)))
    abline(v=sim$refs_interval$FmsyMedianL, col=3,lwd=2)
    abline(v=sim$refs_interval$FmsylowerMedianL, col=3,lwd=2,lty=2)
    abline(v=sim$refs_interval$FmsyupperMedianL, col=3,lwd=2,lty=2)
    abline(v=sim$refs_interval$F5percRiskBlim, col=2,lwd=2)
    abline(v=0)
    lines(auxi$x, auxi$y, lwd=2, col=2)
    return(invisble(NULL))
  }
  
  
}
  