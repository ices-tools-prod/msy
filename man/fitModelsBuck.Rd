\name{fitModelsBuck}
\alias{fitModelsBuck}
\title{stock recruitment function}
\usage{
  fitModelsBuck(data, runid, nsamp = 5000,
    models = c("ricker", "segreg", "bevholt"), ...)
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

