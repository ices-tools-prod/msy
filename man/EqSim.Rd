\name{EqSim}
\alias{EqSim}
\title{simulates the equilibrium results for a population}
\usage{
  EqSim(fit, Nrun = 200, wt.years = c(2007, 2011),
    Fscan = seq(0, 1, len = 20), process.error = TRUE,
    verbose = TRUE, Btrigger = 0, Fphi = 0, Fcv = 0,
    flgsel = 0, flgmatwt = 0)
}
\arguments{
  \item{fit}{A list returned from the function fitModels}

  \item{Nrun}{The number of years to run in total}

  \item{wt.years}{The years to sample weights and selection
  pattern from}

  \item{Fscan}{F values to scan over}

  \item{process.error}{Use stochastic recruitment or mean
  recruitment?  (TRUE = predictive)}

  \item{verbose}{Flag, if TRUE (default) indication of the
  progress of the simulation is provided in the console.
  Useful to turn to FALSE when knitting documents.}

  \item{Btrigger}{If other than 0 (default) the target F
  applied is reduced by SSB/Btrigger}

  \item{Fphi}{Autocorrelation in assessment error in the
  advisory year}

  \item{Fcv}{Assessment error in the advisory year}

  \item{flgsel}{A flag, if FALSE (default) use the mean
  selection value over the "wt.years"}

  \item{flgmatwt}{A flag, if FALSE (default) use the mean
  stock and catch values over the "wt.years"}
}
\value{
  A list containing the following objects: \itemize{ \item
  ssbs A matrix containing the 0.025, 0.050, 0.25, 0.50,
  0.75, 0.95, 0.975 percentiles of ssb for the scanned
  fishing mortalities (see Fscan below). \item cats A
  matrix containing the 0.025, 0.050, 0.25, 0.50, 0.75,
  0.95, 0.975 percentiles of catches for the scanned
  fishing mortalities (see Fscan below). \item recs A
  matrix containing the 0.025, 0.050, 0.25, 0.50, 0.75,
  0.95, 0.975 percentiles of rec for the scanned fishing
  mortalities (see Fscan below). \item ferr An array \item
  ssbsa \item catsa \item recsa \item Mat \item M \item
  Fprop \item Mprop \item west \item weca \item sel \item
  Fscan }
}
\description{
  XXX
}
\examples{
data(codEB)
fit <- fitModels(codEB, nsamp = 100) # a few samples for example
sim <- EqSim(fit, wt.years = c(2007, 2011), Fphi = 0.5, Fcv = 0.1)
hist(c(apply(sim $ ferr, c(1,3), function(x) acf(x, plot = FALSE) $ acf[2])), main = "histogram of estimates of ACF lag 1 in F errors")
arfits <- apply(sim $ ferr, c(1,3), function(x) ar(x)[c("order", "ar", "var.pred")])
summary(sapply(arfits, function(x) x $ order))
summary(sapply(arfits, function(x) x $ ar[1]))
summary(sapply(arfits, function(x) x $ var.pred))
}
\author{
  Colin Millar \email{colin.millar@jrc.ec.europa.eu}
}

