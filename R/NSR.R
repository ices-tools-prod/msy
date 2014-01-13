
       
#' mydf
#'
#'
#' @param rho smoothing parameter
#' @param data XXX
#' @param nknots XXX
#' @param my.df XXX
#' @return something
#' @author Noel Cadigan \email{Noel.Cadigan@@mi.mun.ca}
#' @export
mydf <- function(rho, data, nknots, my.df) 
{
  temp <- scam(log.recruit ~ s(stock.size, k = nknots, bs = "mpd", m = 2) + offset(offset),
               family = gaussian(link="identity"),
               data = data, optimizer = "nlm",
               weights = data $ wt, 
               sp = exp(rho))

  (sum(temp$edf) - my.df)^2                                  
}

#' mygcv
#'
#'
#' @param dat XXX
#' @param rho a parameter
#' @return something
#' @author Noel Cadigan \email{Noel.Cadigan@@mi.mun.ca}
#' @export
mygcv <- function(dat, rho)
{
  temp <- scam(log.recruit ~ s(stock.size, k = nknots, bs = "mpd", m = 2) + offset(offset),
               family = gaussian(link = "identity"), 
               data = dat, 
               optimizer = "nlm",
               weights = dat $ wt,
               sp = exp(rho))

  c(temp $ dev / (n * (1 - sum(temp $ edf) / n)^2))                            
}


#' function name
#'
#'
#' @param dat XXX
#' @param y XXX
#' @param np XXX
#' @param scamfit XXX
#' @param sp XXX
#' @param nboot XXX
#' @return something
#' @author Noel Cadigan \email{Noel.Cadigan@@mi.mun.ca}
#' @export
bootscam <- function(dat, y, np, scamfit, sp, nboot)
{
  ret <- matrix(NA, length(dat $ stock.size), nboot)
  ind <- dat $ wt == 1
  raw.resid <- scamfit $ residuals[ind]
  for(i in 1:nboot) { 
    boot.dat <- dat
    boot.y <- exp(log(y) + sample(raw.resid, n, replace = T))
    boot.dat $ recruit <- c(boot.y, rep(1, np))
    boot.dat $ log.recruit <- log(boot.dat $ recruit)
    scamfiti <- scam(scamfit $ formula, 
                     family = scamfit $ family, 
                     optimizer = scamfit $ optimizer,
                     data = boot.dat,
                     weights = boot.dat $ wt,
                     sp = sp, 
                     start = scamfit$coefficients)
    pred.rec <- exp(scamfiti$fitted.values)      
    pred.prod <- pred.rec / boot.dat$stock.size
    max.pred.rec <- max(pred.rec)      
    max.pred.prod <- max(pred.prod)

    ret[,i] <- scamfiti $ fitted.values 
  }      
  return(ret)
}

