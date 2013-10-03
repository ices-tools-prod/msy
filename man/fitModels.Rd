\name{fitModels}
\alias{fitModels}
\title{stock recruitment function}
\usage{
  fitModels(stk, nsamp = 5000,
    models = c("ricker", "segreg", "bevholt"),
    method = "Buckland", runid = NULL, remove.years = NULL,
    delta = 1.3, nburn = 10000)
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

