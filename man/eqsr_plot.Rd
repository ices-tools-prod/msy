\name{eqsr_plot}
\alias{eqsr_plot}
\title{plot simulated predictive distribution of recruitment}
\usage{
  eqsr_plot(fit, n = 5000, x.mult = 1.1, y.mult = 1.4,
    ggPlot = FALSE)
}
\arguments{
  \item{fit}{an fitted MCMC returned from \code{eqsr_fit}}

  \item{n}{Number of random recruitment draws to plot}

  \item{x.mult}{max.y (ssb) as a multiplier of maximum
  observed ssb}

  \item{y.mult}{max.x (rec) as a multiplier of maismum
  observed rec}

  \item{ggPlot}{Flag, if FALSE (default) plot base
  graphics, if true do a ggplot}
}
\value{
  NULL produces a plot
}
\description{
  plot simulated predictive distribution of recruitment
}
\author{
  Colin Millar \email{colin.millar@jrc.ec.europa.eu}
}

