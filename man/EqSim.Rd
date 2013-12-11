\name{EqSim}
\alias{EqSim}
\title{stock recruitment function}
\usage{
  EqSim(fit, Nrun = 200, wt.years = c(2007, 2011),
    Fscan = seq(0, 1, len = 20), process.error = TRUE,
    verbose = TRUE, Btrigger = 0, Fphi = 0, Fcv = 0)
}
\arguments{
  \item{ab}{the model parameters}

  \item{ssb}{a vector of ssb}
}
\value{
  log recruitment according to model
}
\description{
  stock recruitment function
}
\author{
  Colin Millar \email{colin.millar@jrc.ec.europa.eu}
}

