\name{MH}
\alias{MH}
\title{Runs a Metropolis Hasting MCMC on a stock recruit model}
\usage{
  MH(nt, nburn, data, model, delta = 1.3)
}
\arguments{
  \item{nt}{The total number of MCMC iterations to run}

  \item{nburn}{The number of initial iterations to remove}

  \item{data}{a list or data.frame with elements \code{ssb}
  and \code{rec}}

  \item{delta}{a scaling factor for the MH proposal}

  \item{model}{the recruitment model to use}
}
\value{
  a data.frame of nt - nburn rows and columns for the model
  parameters and the log-likelihood
}
\description{
  Runs a Metropolis Hasting MCMC on a stock recruit model
}
\author{
  Colin Millar \email{colin.millar@jrc.ec.europa.eu}
}

