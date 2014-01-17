\name{srmsymc_read_par}
\alias{srmsymc_read_par}
\title{Read key parameter estimates from the srmsymc program}
\usage{
  srmsymc_read_par(file = "simpar.dat", path, srno = 1,
    minSSB = 0, trim = TRUE, longformat = TRUE)
}
\arguments{
  \item{file}{the file name, should not be anything else
  than \code{simpar.dat}}

  \item{path}{character containing the path to the file. If
  missing it is assumed that the file is in the current
  working directory.}

  \item{srno}{integer. 1=Ricker, 2=Beverton-Holt,
  3=Hockey-stick.}

  \item{minSSB}{numerical. Minimum observed SSB value. Used
  for the hockey-stick, excluding values below the lowest
  observed. If set to 0 (default) do not exclude values
  below the minimum.}

  \item{trim}{boolean. If TRUE (default) come cleaning
  (NEED TO SPECIFY) is done.}

  \item{longformat}{boolean. If TRUE (default) return a
  long format.}
}
\description{
  Reads the fiel \code{simpar.dat} generate by the
  \code{srmsymc}. Returns a list containing the stochastic
  and the deterministic estimates.
}

