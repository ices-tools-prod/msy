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
\examples{
data(codEB)
fit <- fitModels(codEB, nsamp = 100) # a few samples for example
sim <- EqSim(fit, wt.years = c(2007, 2011), Fphi = 0.5, Fcv = 0.1)
hist(c(apply(sim $ ferr, c(1,3), function(x) acf(x, plot = FALSE) $ acf[2])), main = "histogram of estimates of ACF lag 1 in F errors")
arfits <- apply(sim $ ferr, c(1,3), function(x) ar(x)[c("order", "ar", "var.pred")])
summary(sapply(arfits, function(x) x $ order))
summary(sapply(arfits, function(x) x $ ar[1]))
summary(sapply(arfits, function(x) x $ var.pred))
}
\author{
  Colin Millar \email{colin.millar@jrc.ec.europa.eu}
}

