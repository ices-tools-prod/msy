\name{EqSim}
\alias{EqSim}
\title{simulates the equilibrium results for a population}
\usage{
EqSim(fit, bio.years = c(2007, 2012), bio.const = 0, sel.years = c(2005,
  2010), sel.const = 0, Fscan = seq(0, 1, len = 20), Fcv = 0, Fphi = 0,
  recruitment.trim = c(3, -3), Btrigger = 0, Nrun = 200,
  process.error = TRUE, verbose = TRUE)
}
\arguments{
  \item{fit}{A list returned from the function fitModels}

  \item{bio.years}{The years to sample maturity, weights
  and M from}

  \item{bio.const}{A flag, if FALSE mean of the biological
  values from the years selected are used}

  \item{sel.years}{The years to sample the selection
  patterns from}

  \item{sel.const}{A flag, if FALSE mean of the selection
  patterns from the years selected are used}

  \item{Fscan}{F values to scan over}

  \item{Fcv}{Assessment error in the advisory year}

  \item{Fphi}{Autocorrelation in assessment error in the
  advisory year}

  \item{recruitment.trim}{A numeric vector with two
  log-value clipping the extreme recruitment values from a
  continuous lognormal distribution. The values must be set
  as c("high","low").}

  \item{Btrigger}{If other than 0 (default) the target F
  applied is reduced by SSB/Btrigger}

  \item{Nrun}{The number of years to run in total (last 50
  years from that will be retained)}

  \item{process.error}{Use stochastic recruitment or mean
  recruitment?  (TRUE = predictive)}

  \item{verbose}{Flag, if TRUE (default) indication of the
  progress of the simulation is provided in the console.
  Useful to turn to FALSE when knitting documents.}
}
\value{
A list containing the following objects: \itemize{ \item
ssbs A matrix containing the 0.025, 0.050, 0.25, 0.50,
0.75, 0.95, 0.975 percentiles of ssb for the scanned
fishing mortalities (see Fscan below). \item cats A matrix
containing the 0.025, 0.050, 0.25, 0.50, 0.75, 0.95, 0.975
percentiles of catches for the scanned fishing mortalities
(see Fscan below). \item recs A matrix containing the
0.025, 0.050, 0.25, 0.50, 0.75, 0.95, 0.975 percentiles of
rec for the scanned fishing mortalities (see Fscan below).
\item ferr An array \item ssbsa \item catsa \item recsa
\item Mat \item M \item Fprop \item Mprop \item west \item
weca \item sel \item Fscan }
}
\description{
XXX
}
\author{
Colin Millar \email{colin.millar@jrc.ec.europa.eu}
}

