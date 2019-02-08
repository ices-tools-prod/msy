#' Simulates the Equilibrium Results for a Population.
#'
#' Simulate a fish stock forward in time given biological parameters, fishery
#' parameters and advice parameters.
#'
#' @param fit A list returned from the function fitModels
#' @param bio.years The years to sample maturity, weights and M from, given as
#'                  a vector of length 2, i.e. c(2010, 2015) select from the
#'                  years 2010 to 2015 inclusive.
#' @param bio.const A flag (default FALSE), if TRUE mean of the biological values from the
#'                  years selected are used
#' @param sel.years The years to sample the selection patterns from, given as
#'                  a vector of length 2, i.e. c(2010, 2015) select from the
#'                  years 2010 to 2015 inclusive.
#' @param sel.const A flag (default FALSE), if TRUE mean of the selection patterns from the
#'                  years selected are used
#' @param Fscan F values to scan over, i.e. seq(0, 2, by = 0.05)
#' @param Fcv Assessment error in the advisory year
#' @param Fphi Autocorrelation in assessment error in the advisory year
#' @param SSBcv Spawning stock biomass error in the advisory year
#' @param rhologRec A flag for recruitment autocorrelation, default (TRUE).
#' @param Blim SSB limit reference point
#' @param Bpa SSB precuationary reference point
#' @param recruitment.trim A numeric vector with two log-value clipping the
#'        extreme recruitment values from a continuous lognormal distribution.
#'        The values must be set as c("high","low").
#' @param Btrigger If other than 0 (default) the target F applied is reduced by
#'                 SSB/Btrigger. This is the "ICES Advice Rule".
#' @param Nrun The number of years to run in total (the last 50 years from that
#'             will be retained to compute equilibrium values from)
#' @param process.error Use stochastic recruitment or mean recruitment?
#'                      (TRUE uses the predictive distribution of recruitment,
#'                      model estimate of recruitment + simulated observation
#'                      error)
#' @param verbose Flag, if TRUE (default) indication of the progress of the
#'        simulation is provided in the console. Useful to turn to FALSE when
#'        knitting documents.
#' @param extreme.trim a pair of quantiles (low, high) which are used to trim
#'                     the equilibrium catch values, across simulations within
#'                     an F scenario, when calculating the mean catch and
#'                     landings for that F scenario.  These mean values
#'                     calculated accross simulations within an F scenario
#'                     are used to find which F scenario gave the maximum catch.
#'                     \code{extreme.trim} can therefore be used to stablise the
#'                     estimate of mean equilibrium catch and landings by F
#'                     scenario.
#' @return
#' A list containing the results from the forward simulation and the reference
#' points calculated from it.
#'
#' @details
#' Details of the steps required to evaluate reference points are given in
#' ICES (2017).  WHile, details of the calculation of MSY ranges is given in
#' ICES (2015).
#'
#' @references
#' ICES (2015) Report of the Workshop to consider F MSY ranges for stocks in
#' ICES categories 1 and 2 in Western Waters (WKMSYREF4).
#' \href{http://ices.dk/sites/pub/Publication\%20Reports/Expert\%20Group\%20Report/acom/2015/WKMSYREF4/01\%20WKMSYREF4\%20Report.pdf}{01
#' WKMSYREF4 Report.pdf}
#'
#' ICES (2017) ICES fisheries management reference points for category 1 and 2
#' stocks.
#' DOI: \href{https://doi.org/10.17895/ices.pub.3036}{10.17895/ices.pub.3036}
#'
#' @seealso
#'
#' \code{\link{eqsr_fit}} fits multiple stock recruitment models to a data set.
#'
#' \code{\link{eqsr_plot}} plots the results from eqsr_fit.
#'
#' \code{\link{eqsim_plot}} summary plot of the forward simulation showing estimates
#'   of various reference points.
#'
#' \code{\link{eqsim_plot_range}} summary plots of the forward simulation showing
#'   the estimates of MSY ranges (ICES, 2015)
#'
#' \code{\link{msy-package}} gives an overview of the package.
#'
#' @examples
#' \dontrun{
#' data(icesStocks)
#' FIT <- eqsr_fit(icesStocks$saiNS,
#'                 nsamp = 1000,
#'                 models = c("Ricker", "Segreg"))
#' SIM <-
#'   eqsim_run(
#'     FIT,
#'     bio.years = c(2004, 2013),
#'     sel.years = c(2004, 2013),
#'     Fcv = 0.24,
#'     Fphi = 0.42,
#'     Blim = 106000,
#'     Bpa = 200000,
#'     Fscan = seq(0, 1.2, len = 40)
#'    )
#' }
#'
#' @export
eqsim_run <- function(fit,
                      bio.years = c(2008, 2012), # years sample weights, M and mat
                      bio.const = FALSE,
                      sel.years= c(2008, 2012), # years sample sel and discard proportion by number from
                      sel.const = FALSE,
                      Fscan = seq(0, 1, len = 20), # F values to scan over
                      Fcv = 0,
                      Fphi = 0,
                      SSBcv = 0,
                      rhologRec = TRUE,
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
  if (diff(recruitment.trim) > 0) stop("recruitment truncation must be given as c(high, low)")

  if (verbose) icesTAF::msg("Setting up...")

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
  Ferr[1,] <- stats::rnorm(n=Nmod, mean=0, sd=1)*Fcv/sqrt(1-Fphi^2)
  for(j in 2:Nrun) { Ferr[j,] <- Fphi*Ferr[j-1,] + Fcv*stats::rnorm(n=Nmod, mean=0, sd=1) }

  # 2014-03-12: Changed per note form Carmen/John
  #  Errors in SSB: this is used when the ICES MSY HCR is applied for F
  SSBerr <- matrix(stats::rnorm(n=Nrun*Nmod, mean=0, sd=1), ncol=Nmod) * SSBcv

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
  resids= array(stats::rnorm(Nmod*(Nrun+1), 0, SR$cv),c(Nmod, Nrun+1))

  # 2014-03-12: Changed per note form Carmen/John
  #  Autocorrelation in Recruitment Residuals:
  if(rhologRec){
    fittedlogRec <-  do.call(cbind, lapply( c(1:nrow(fit$sr.sto)), function(i){
      FUN <- match.fun(fit$sr.sto$model[i])
      FUN(fit$sr.sto[i, ], fit$rby$ssb) } )  )
    # Calculate lag 1 autocorrelation of residuals:
    rhologRec <- apply(log(fit$rby$rec)-fittedlogRec, 2, function(x){stats::cor(x[-length(x)],x[-1])})
    # Draw residuals according to AR(1) process:
    for(j in 2:(Nrun+1)){ resids[,j] <- rhologRec * resids[,j-1] + resids[,j]*sqrt(1 - rhologRec^2) }
  }


  # Limit how extreme the Rec residuals can get:
  lims = t(array(SR$cv,c(Nmod,2))) * recruitment.trim
  for (k in 1:Nmod) { resids[k,resids[k,]>lims[1,k]]=lims[1,k]}
  for (k in 1:Nmod) { resids[k,resids[k,]<lims[2,k]]=lims[2,k]}
  # end New from Simmonds 29.1.2014

  if (verbose) icesTAF::msg("Running forward simulations.")
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
    ssbs[, i] <- stats::quantile(ssby[begin:Nrun, ], quants)
    cats[, i] <- stats::quantile(Cat[begin:Nrun, ], quants)
    lans[, i] <- stats::quantile(Lan[begin:Nrun, ], quants)
    recs[, i] <- stats::quantile(Ny[1, begin:Nrun, ], quants)


    ferr[i, , ] <- Ferr[begin:Nrun, ]
    ssbsa[i, , ] <- ssby[begin:Nrun, ]
    catsa[i, , ] <- Cat[begin:Nrun, ]
    lansa[i, , ] <- Lan[begin:Nrun, ]
    recsa[i, , ] <- Ny[1, begin:Nrun, ]

    if (verbose) loader(i/NF)
  }

  if (verbose) icesTAF::msg("Summarising simulations")

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
    catm <- apply(catsa, 1, function(x){mean(x[x <= stats::quantile(x, extreme.trim[2]) & x >= stats::quantile(x, extreme.trim[1])])})
    lanm <- apply(lansa, 1, function(x){mean(x[x <= stats::quantile(x, extreme.trim[2]) & x >= stats::quantile(x, extreme.trim[1])])})
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
  vcumLan <- stats::median(FmsyLan)
  fmsy.densLan <- stats::density(FmsyLan)
  vmodeLan <- fmsy.densLan$x[which.max(fmsy.densLan$y)]

  FmsyCat <- Fscan[maxpf]
  msymCat <- mean(FmsyCat)
  vcumCat <- stats::median(FmsyCat)
  fmsy.densCat <- stats::density(FmsyCat)
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
  catC <- stats::approx(Fscan, cats[4,], xout = catF)$y
  lanC <- stats::approx(Fscan, lans[4,], xout = lanF)$y
  catB <- stats::approx(Fscan, ssbs[4,], xout = catF)$y
  lanB <- stats::approx(Fscan, ssbs[4,], xout = lanF)$y

  Refs <- rbind(catF, lanF, catC, lanC, catB, lanB)
  rownames(Refs) <- c("catF","lanF","catch","landings","catB","lanB")
  colnames(Refs) <- c("F05","F10","F50","medianMSY","meanMSY","FCrash05","FCrash50")

  #TODO: id.sim - user specified.

  # 2014-03-12 Ammendments per note from Carmen/John
  # CALCULATIONS:

  # Fmsy: value that maximises median LT catch or median LT landings
  auxi <- stats::approx(Fscan, cats[4, ],xout=seq(min(Fscan),max(Fscan),length=200))
  FmsyMedianC <- auxi$x[which.max(auxi$y)]
  MSYMedianC <- max(auxi$y)
  # Value of F that corresponds to 0.95*MSY:
  FmsylowerMedianC <- auxi$x[ min( (1:length(auxi$y))[auxi$y/MSYMedianC >= 0.95] ) ]
  FmsyupperMedianC <- auxi$x[ max( (1:length(auxi$y))[auxi$y/MSYMedianC >= 0.95] ) ]

  auxi <- stats::approx(Fscan, lans[4, ],xout=seq(min(Fscan),max(Fscan),length=200))
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

  if (verbose) icesTAF::msg("Calculating MSY range values")

  sim <- eqsim_range(sim)

  return(sim)

}
