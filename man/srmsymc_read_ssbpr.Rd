\name{srmsymc_read_ssbpr}
\alias{srmsymc_read_ssbpr}
\title{Read the ssb per recruit estimates from the srmsymc program}
\usage{
srmsymc_read_ssbpr(file = "simparssbpr.dat", trim = TRUE,
  longformat = TRUE, doPlot = TRUE)
}
\arguments{
  \item{file}{the file name, should not be anything else
  than \code{simparssbpr.dat}}

  \item{trim}{boolean. If TRUE (default) come cleaning
  (NEED TO SPECIFY) is done.}

  \item{longformat}{boolean. If TRUE (default) return a
  long format.}

  \item{doPlot}{boolean. If TRUE (default) returns a
  ggplot}
}
\description{
Reads the file \code{simparssbpr.dat} generate by the
\code{srmsymc}. In the \code{simparssbpr.dat} the first row
is the fishing mortality, the second row is the
deterministic estimates of ssb and the remaining rows
represents the catch for each iterations. Returns a list
containing two data.frames, one containing the
deterministic estimate, the other containing the stochastic
simulations.
}

