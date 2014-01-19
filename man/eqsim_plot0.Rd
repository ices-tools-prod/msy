\name{eqsim_plot0}
\alias{eqsim_plot0}
\title{Plot the results from eqsim}
\usage{
  eqsim_plot0(sim, Blim, Bpa = 1.4 * Blim,
    ymax = c(NA, NA, NA), plot = TRUE, yield = "landings")
}
\arguments{
  \item{sim}{An object returned from the function
  eqsim_run}

  \item{Blim}{Value for the Blim}

  \item{Bpa}{Value for the Bpa}

  \item{ymax}{vector of three values, dictating maximum
  y-value on the recruitment plot, the maximum value on the
  ssb-plot and the maximum value on the catch plot (in that
  order)}

  \item{plot}{Flag, if TRUE (default) provide a plot as
  part of the output}

  \item{yield}{Character}
}
\description{
  XXX
}
\author{
  Colin Millar \email{colin.millar@jrc.ec.europa.eu}
}

