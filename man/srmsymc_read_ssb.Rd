\name{srmsymc_read_ssb}
\alias{srmsymc_read_ssb}
\title{Read the ssb estimates the srmsymc program}
\usage{
  srmsymc_read_ssb(path, file = "simparssb.dat", srno = 1,
    trim = TRUE, longformat = TRUE, doPlot = TRUE)
}
\arguments{
  \item{path}{character containing the path to the file. If
  missing it is assumed that the file is in the current
  working directory.}

  \item{file}{the file name, should not be anything else
  than \code{simpary.dat}}

  \item{srno}{integer. 1=Ricker, 2=Beverton-Holt,
  3='Hockey-stick.}

  \item{trim}{boolean. If TRUE (default) come cleaning
  (NEED TO SPECIFY) is done.}

  \item{longformat}{boolean. If TRUE (default) return a
  long format.}

  \item{doPlot}{boolean. If TRUE (default) returns a
  ggplot}
}
\description{
  Reads the file \code{simpary.dat} generate by the
  \code{srmsymc}. In the \code{simpary.dat} the first row
  is the fishing mortality, the second row is the
  deterministic estimates of ssb and the remaining rows
  represents the catch for each iterations. Returns a list
  containing two data.frames, one containing the
  deterministic estimate, the other containing the
  stochastic simulations.
}

