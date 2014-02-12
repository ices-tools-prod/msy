\name{fitModelsBuck}
\alias{fitModelsBuck}
\title{stock recruitment function}
\usage{
fitModelsBuck(data, runid, nsamp = 5000, models = c("ricker", "segreg",
  "bevholt"), ...)
}
\arguments{
  \item{data}{data.frame containing stock recruitment data}

  \item{runid}{A character vector}

  \item{nsamp}{Number of samples}

  \item{models}{A character vector}

  \item{...}{Additional arguements}
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

