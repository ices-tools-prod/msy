\name{fitModelsSimmonds}
\alias{fitModelsSimmonds}
\title{stock recruitment function}
\usage{
  fitModelsSimmonds(data, runid, delta = 1.3,
    nburn = 10000, nsamp = 5000,
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

