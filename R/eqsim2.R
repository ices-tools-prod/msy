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

eqsim_run2 <- function(fit,
                      bio.years = c(2008, 2012), # years sample weights, M and mat
                      bio.const = FALSE,
                      sel.years= c(2008, 2012), # years sample sel and discard proportion by number from
                      sel.const = FALSE,
                      Fscan = c(0.1,0.2), # F values to scan over
                      Fcv = 0,
                      Fphi = 0,
                      Blim,
                      Bpa,
                      recruitment.trim = c(3, -3),
                      Btrigger = 0,
                      n_years = 100, # number of years to run in total
                      process.error = TRUE, # use predictive recruitment or mean recruitment? (TRUE = predictive)
                      verbose = TRUE,
                      extreme.trim)
{
  
  # some minor input checking
  if (abs(Fphi) >= 1) stop("Fphi, the autocorelation parameter for log F should be between (-1, 1)")
  if ((recruitment.trim[1] + recruitment.trim[2])> 0) stop("recruitment truncation must be between a high - low range")
  
  ##############################################################################
  # 1. SETTING DIMENSIONS
  #  age
  a1 <- range(fit$stk)[1]
  a2 <- range(fit$stk)[2]
  n_ages <- dims(fit$stk)$age
  #  time
  y1 <- range(fit$stk)[5]
  y2 <- y1 + n_years
  n_years_keep <- min(n_years, 50)
  #  target
  HRATE <- Fscan
  n_target <- length(HRATE)
  #  iterations
  iter <- nrow(fit$sr.sto)
  
  ##############################################################################
  # 2. SETTING UP OBJECTS
  
  # Array of n_years and n_iter
  x <- array(0, c(n_years,iter),dimnames=list(year=1:n_years,iter=1:iter))
  SSBy <- Fy_error <- x
  
  # Array of n_ages, n_years and n_iter
  x <- array(0, c(n_ages, n_years, iter),dimnames=list(age=a1:a2,year=1:n_years,iter=1:iter))
  Nay <- Fay  <- Cay <- cWay <- lWay <- lRay <- x
  
  # TODO per note from Carmen:
  #  NOTE: If we want Fy_error to be a stationary AR(1) process, it would make
  #        more sense to initialise Fy_error as a Normal dist with zero mean and
  #        standard deviation of AR(1) marginal distribution, i.e. standard 
  #        deviation of initial Fy_error = Fcv/sqrt(1- Fphi^2), instead of just
  #        initialising Fy_error=0
  

  # Note, we could equivalent for the other ay data
  
  # initial recruitment
  R <- mean(fit$rby$rec)
  
  x <- array(0,
             c(7, n_target),
             dimnames=list(quants=c("p025","p05","p25","p50","p75","p95","p975"),target=Fscan))
  ssbs <- cats <- lans <- recs <- x
  
  x  <- array(0,
              c(n_target, n_years_keep, iter),
              dimnames=list(fmort=Fscan,year=1:n_years_keep,iter=1:iter))
  ferr <- SSB <- CATCH <- LAND <- REC <- x
  begin <- n_years - n_years_keep + 1
  
  # New from Simmonds' 29.1.2014
  #   Residuals of SR fits (1 value per SR fit and per simulation year 
  #     but the same residual value for all Fscan values):
  SR <- fit$sr.sto
  #data <- fit $ rby[,c("rec","ssb","year")]
  resids= array(rnorm(iter*(n_years+1), 0, SR$cv),c(iter, n_years+1))
  # Limit how extreme the Rec residuals can get:
  lims = t(array(SR$cv,c(iter,2))) * recruitment.trim
  for (k in 1:iter) { resids[k,resids[k,]>lims[1,k]]=lims[1,k]}
  for (k in 1:iter) { resids[k,resids[k,]<lims[2,k]]=lims[2,k]}
  # end New from Simmonds 29.1.2014
  
  
  
  
  ##############################################################################
  # A. data input
  
  # NOTE: Should make this as a function call, i.e. independent of the rest
  
  ## A1: little helper functions
  littleHelper <- function(x,i) {
    x2 <- x
    x2[i] <- NA
    x2[] <- apply(x2,1,mean,na.rm=TRUE)
    x[i] <- x2[i]
    return(x)
  }
  
  ## A2: Biological data
  b1 <- bio.years[1]
  b2 <- bio.years[2]
  stk.winbio <- window(fit$stk, start = b1, end = b2)  
  
  west <- matrix(stock.wt(stk.winbio), ncol = b2 - b1 + 1)
  i <- west == 0
  if(any(i)) west <- littleHelper(west,i)
  weca <- matrix(catch.wt(stk.winbio), ncol = b2 - b1 + 1)
  i <- weca == 0
  if(any(i)) weca <- littleHelper(weca,i)
  wela <- matrix(landings.wt(stk.winbio), ncol = b2 - b1 + 1)
  if(any(i)) wela <- littleHelper(wela,i)
  Mat <- matrix(mat(stk.winbio), ncol = b2 - b1 + 1)
  M <- matrix(m(stk.winbio), ncol = b2 - b1 + 1)
  
  # 22.2.2014 Added weight of landings per comment from Carmen
  if (bio.const==TRUE){ # take means of wts Mat and M and ratio of landings to catch
    west[] <- apply(west, 1, mean)
    weca[] <- apply(weca, 1, mean)
    wela[] <- apply(wela, 1, mean)
    Mat[] <- apply(Mat, 1, mean)
    M[] <- apply(M, 1, mean) #me
  }
  
  rand_bio <- array(sample(1:ncol(weca), n_years * iter, TRUE), c(n_years, iter))
  cWay[] <- c(weca[, c(rand_bio)])
  lWay[] <- c(wela[, c(rand_bio)])

  
  ## B2: Fisheries data
  s1 <- sel.years[1]
  s2 <- sel.years[2]
  stk.winsel <- window(fit$stk, start = s1  , end = s2)
  
  landings <- matrix(landings.n(stk.winsel), ncol = s2 - s1 + 1)
  catch <- matrix(catch.n(stk.winsel), ncol = s2 - s1 + 1)
  sel <- matrix(harvest(stk.winsel), ncol = s2 - s1 + 1)
  Fbar <- matrix(fbar(stk.winsel), ncol = s2 - s1  + 1)
  sel <- sweep(sel, 2, Fbar, "/")
  
  if (sel.const == TRUE) { # take means of selection
    sel[] <- apply(sel, 1, mean)
    landings[]  <- apply(landings, 1, mean)
    catch[]  <- apply(catch, 1, mean)
  }

  land.cat= landings / catch  # ratio of number of landings to catch
  # TODO: Check if this is sensible
  i <- is.na(land.cat)
  if(any(i)) land.cat[i] <- 1
  
  rand_sel <- array(sample(1:ncol(sel), n_years * iter, TRUE), c(n_years, iter))
  lRay[]  <- c(land.cat[, c(rand_sel)])
  
  # No variability assumed in the pFa and pMa. Here generate an age vector
  pFa <- apply(harvest.spwn(stk.winsel), 1, mean)[drop=TRUE] # vmean(harvest.spwn(stk.win))
  pMa <- apply(m.spwn(stk.winsel), 1, mean)[drop=TRUE] # mean(m.spwn(stk.win))
  
  ##############################################################################
  # B. Simulation setup
  
  if (verbose) loader(0)
  
  # Looping over each F value in Fscan. For each of the iter SR fits 
  # (replicates), do a forward simulation during n_years years
  # There are Rec residuals for each SR fit and year, which take the same
  # values for all Fscan 
  for (i in 1:n_target) {
    
    # The F value to test
    Fbar <- Fscan[i]
    
    ############################################################################
    # Population in simulation year 1:
    
    # Zpre: Z that occurs before spawning
    Zpre <- ( sel[,rand_sel[1,]]*Fbar * pFa + M[,rand_bio[1,]] * pMa)
    
    # Zpos: Z that occurs after spawning
    # Zpos not used anywhere
    # Zpos <- (Fbar * (1-pFa) * sel[,rand_sel[1,]] + M[,rand_bio[1,]] * (1-pMa))
    
    # run Z out to age 50 ...
    # TODO:
    # Comments from Carmen: Zcum is a cumulative sum, but it is done in a strange way:
    #  There is a matrix of F-at-age and a matrix of M-at-age (each has 49 ages, iter replicates)
    #  The F and M matrices are summed, giving Z-at-age (49 ages, iter replicates)
    #  But then a cumsum is taken considering the Z-at-age matrix as a vector (i.e. not column-wise) ????
    #  This is strange, by applying "cumsum" treating Z-at-age as a vector, really only the first 50 values of
    #  the resulting "Zcum" make sense (all other values seem "wrong", or at least, meaningless)
    Zcum <- c(0, cumsum(Fbar * sel[c(1:n_ages, rep(n_ages, 49 - n_ages)), rand_sel[1,]] + M[c(1:n_ages, rep(n_ages, 49 - n_ages)), rand_bio[1,]]))
    # Carmen: Following from "Zcum", only first 50 elements of N1 make sense ????
    N1 <- R * exp(- unname(Zcum))
    
    # set up age structure in first year for all simulations
    # Comments from Carmen:
    #   Nay has dimension = (no. ages, no. simulation yrs "n_years", no. SR fits "iter")
    #   With this code, we seem to be getting always the same population-at-age value for year 1
    #   instead of iter different values, as might have been intended ????
    #   (the whole problem is coming from Zcum ==> N1 ==> Nay[,1,] )
    Nay[,1,] <- c(N1[1:(n_ages-1)], sum(N1[n_ages:50]))
    
    # calculate ssb in first year using a different stock.wt and Mat selection and M for each simulation
    # Comments from Carmen:
    #   SSBy has dimension = (no. simul yrs "n_years", no. SR fits "iter")
    #   SSB in year 1:
    #   although Nay[,1,] has dim no.ages x iter, all iter values of Nay[,1,] are
    #   the same (because of Zcum issue)
    SSBy[1,] <- colSums(Mat[,rand_bio[1,]] * Nay[,1,] * west[,rand_bio[1,]] / exp(Zpre)[])
    
    # Years 2 to n_years:
    for (j in 2:n_years) {
      # get ssb from previous year
      SSBy_1 <- SSBy[j-1,]
      
      # predict recruitment using various models
      if (process.error) {
        # Changes 29.1.2014
        # new random draws each time
        # allrecs <- sapply(unique(SR $ mod), function(mod) exp(match.fun(mod) (SR, SSBy_1) + rnorm(iter, 0, SR $ cv)))
        # same random draws used for each F
        allrecs <- sapply(unique(SR$mod), function(mod) exp(match.fun(mod)(SR, SSBy_1) + resids[,j]))
        # end Changes 29.1.2014
      } else {
        allrecs <- sapply(unique(SR $ mod), function(mod) exp(match.fun(mod) (SR, SSBy_1)))
      }
      
      # Comment from Carmen:
      #  For each of the iter replicates, this selects the appropriate SR model
      #   type to use in that replicate
      #  Note that the order of SR model types that comes out in "select" is
      #   not necessarily the same order in which the SR model types were
      #   entered as inputs -- I presume the **next 2 lines** of code have
      #   been checked to avoid potential bugs due to this reordering  ???? 
      select <- cbind(seq(iter), as.numeric(factor(SR $ mod, levels = unique(SR $ mod))))
      Nay[1,j,] <- allrecs[select]
      
      # Comment from Carmen:
      #   Note: it seems that Rec is coded as occurring always at age 1
      #   (i.e. based on SSBy_1 in previous year)
      #   Some stocks have Rec at ages other than 1 (e.g. age 0) 
      #    -- is this a problem ????
      
      # apply HCR
      # (intended) Fbar to be applied in year j-1 (depends on SSBy_1 in year j-1):
      Fnext <- Fbar * pmin(1, SSBy_1/Btrigger)
      
      # apply some noise to the F
      # Notes from Carmen:
      #  Assessment and/or implementation error (modifies intended F to get
      #  realised F)
      #  Error: AR(1) process on log(F) with autocorrel = Fphi, and
      #  conditional stand deviation = Fcv
      #  Might make more sense to have the "Fy_error" matrix calculated before
      #  the Fscan loop starts so that the same errors in F are applied to
      #  all Fscan values ???? (as for Rec residuals)
      Fy_error[j,] <- Fphi * Fy_error[j-1,] + rnorm(iter, 0, Fcv)
      
      # realised Fbar in year j-1:
      Fnext <- exp(Fy_error[j,]) * Fnext
      
      # get a selection pattern for each simulation and apply this to get N
      Zpre <- rep(Fnext, each = length(pFa)) * pFa * sel[, rand_sel[j,]] + M[, rand_bio[j,]] * pMa
      
      # get Fay
      Fay[ , j-1, ] <- rep(Fnext, each = n_ages) * sel[, rand_sel[j-1,]]
      
      Nay[ -1, j, ] <- Nay[1:(n_ages-1), j-1, ] * exp(-Fay[1:(n_ages-1), j-1, ] - M[1:(n_ages-1), rand_bio[j-1,]])
      Nay[n_ages, j, ] <- Nay[n_ages, j, ] + Nay[n_ages, j-1, ] * exp(-Fay[n_ages, j-1, ] - M[n_ages, rand_bio[j-1,]])
      # calculate ssb and catch.n
      SSBy[j, ] <- apply(array(Mat[, rand_bio[j,]] * Nay[,j,] * west[, rand_bio[j,]] / exp(Zpre), c(n_ages, iter)), 2, sum)
      Cay[, j, ] <- Nay[, j-1, ] * Fay[, j-1, ] / (Fay[, j-1, ] + M[, rand_bio[j-1,]]) * (1 - exp(-Fay[, j-1, ] - M[, rand_bio[j-1,]]))
    }
    
    # convert to catch weight
    Cw <- Cay * cWay   # catch Numbers *catch wts
    land <- Cay*lRay*lWay # catch Numbers * Fraction (in number) landed and landed wts
    Lan=apply(land,2:3,sum)
    Cat <- apply(Cw, 2:3, sum)
    
    # summarise everything and spit out!
    quants <- c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
    ssbs[, i] <- quantile(SSBy[begin:n_years, ], quants)
    cats[, i] <- quantile(Cat[begin:n_years, ], quants)
    lans[, i] <- quantile(Lan[begin:n_years, ], quants)
    recs[, i] <- quantile(Nay[1, begin:n_years, ], quants)
    
    
    ferr[i, , ] <- Fy_error[begin:n_years, ]
    SSB[i, , ] <- SSBy[begin:n_years, ]
    CATCH[i, , ] <- Cat[begin:n_years, ]
    LAND[i, , ] <- Lan[begin:n_years, ]
    REC[i, , ] <- Nay[1, begin:n_years, ]
    
    if (verbose) loader(i/n_target)
  }
  

  
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
  
  FCrash05 <- Fscan[which.max(cats[2,]):n_target][ which(cats[2, which.max(cats[2,]):n_target] < 0.05*max(cats[2,]) )[1] ]
  FCrash50 <- Fscan[which.max(cats[4,]):n_target][ which(cats[4, which.max(cats[4,]):n_target] < 0.05*max(cats[4,]) )[1] ]
  
  
  # Einar amended 30.1.2014
  if(missing(extreme.trim)) {
    catm <- apply(CATCH, 1, mean)
    lanm <- apply(LAND, 1, mean)
  } else {
    x <- CATCH
    i <- x > quantile(x,extreme.trim[2]) |
      x < quantile(x,extreme.trim[1])
    x[i] <- NA
    catm <- apply(x, 1, mean, na.rm=TRUE)
    
    x <- LAND
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
  
  
  catsam <- apply(CATCH, c(1,3), mean)
  lansam <- apply(LAND, c(1,3), mean)
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
    pBlim <- apply(SSB > Blim, 1, mean)
    
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
    pBpa <- apply(SSB > Bpa, 1, mean) 
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
  
  return(list(ibya=list(Mat = Mat, M = M, pFa = pFa, pMa = pMa, 
                        west = west, weca = weca, sel = sel),
              rbya=list(ferr=ferr),
              rby=fit$rby, rbp=rbp, Blim=Blim, Bpa=Bpa, Refs = Refs,
              pProfile=pProfile,id.sim=fit$id.sr))
  
}