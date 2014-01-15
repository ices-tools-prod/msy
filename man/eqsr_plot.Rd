\name{eqsr_plot}
\alias{eqsr_plot}
\title{plot simulated predictive distribution of recruitment}
\usage{
  eqsr_plot(fit, n = 5000, ggPlot = FALSE)
}
\arguments{
  \item{fit}{an fitted MCMC returned from \code{eqsr_fit}}

  \item{n}{Number of random recruitment draws to plot}

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

