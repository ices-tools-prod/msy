\name{eqsim_plot}
\alias{eqsim_plot}
\title{Plot the results from eqsim}
\usage{
  eqsim_plot(sim, fit, Blim, Bpa = 1.4 * Blim,
    ymax = c(NA, NA, NA), xmax, plot = TRUE,
    plotMeanCatch = FALSE)
}
\arguments{
  \item{sim}{An object returned from the function EqSim}

  \item{fit}{An object returned from the function
  fitModels}

  \item{Blim}{Value for the Blim}

  \item{Bpa}{Value for the Bpa}

  \item{ymax}{vector of three values, dictating maximum
  y-value on the recruitment plot, the maximum value on the
  ssb-plot and the maximum value on the catch plot (in that
  order)}

  \item{xmax}{Value for the maximum F to plot, if missing
  will use the whole F-range simulated.}

  \item{plot}{Flag, if TRUE (default) provide a plot as
  part of the output}

  \item{plotMeanCatch}{Flag, if FALSE (default) does not
  provide values for the mean catch.}
}
\description{
  This is a modification of the \code{Eqplot}-function
  where the mean catch can be turned on-off by users choise
  with some additional deletion of vertical lines that were
  deemed unessessary at the wkmsyref2 meeting.
}
\author{
  Colin Millar \email{colin.millar@jrc.ec.europa.eu}
}

