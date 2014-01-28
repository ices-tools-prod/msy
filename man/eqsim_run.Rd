\name{eqsim_run}
\alias{eqsim_run}
\title{simulates the equilibrium results for a population}
\usage{
  eqsim_run(fit, bio.years = c(2008, 2012),
    bio.const = FALSE, sel.years = c(2008, 2012),
    sel.const = FALSE, Fscan = seq(0, 1, len = 20),
    Fcv = 0, Fphi = 0, Blim, Bpa, Btrigger = 0, Nrun = 200,
    process.error = TRUE, verbose = TRUE)
}
\arguments{
  \item{fit}{A list returned from the function fitModels}

  \item{bio.years}{The years to sample maturity, weights
  and M from}

  \item{bio.const}{A flag, if FALSE mean of the biological
  values are used}

  \item{sel.years}{The years to sample sel and discard
  proportion by number from}

  \item{sel.const}{A flag, if FALSE mean selection values
  used}

  \item{Fscan}{F values to scan over}

  \item{Fcv}{Assessment error in the advisory year}

  \item{Fphi}{Autocorrelation in assessment error in the
  advisory year}

  \item{Blim}{This we know}

  \item{Bpa}{This we know}

  \item{Btrigger}{If other than 0 (default) the target F
  applied is reduced by SSB/Btrigger}

  \item{Nrun}{The number of years to run in total}

  \item{process.error}{Use stochastic recruitment or mean
  recruitment?  (TRUE = predictive)}

  \item{verbose}{Flag, if TRUE (default) indication of the
  progress of the simulation is provided in the console.
  Useful to turn to FALSE when knitting documents.}
}
\description{
  XXX
}
\author{
  Colin Millar \email{colin.millar@jrc.ec.europa.eu}
}

