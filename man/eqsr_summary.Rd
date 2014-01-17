\name{eqsr_summary}
\alias{eqsr_summary}
\title{Summarise predictive distribution of recruitment}
\usage{
  eqsr_summary(fit)
}
\arguments{
  \item{fit}{A list containing fitted MCMC. Returned from
  \code{eqsr_fit}}
}
\value{
  A list containing the following data.frames: \itemize{
  \item data The original ssb-rec estimations \item
  Parameters The estimated alpha (a), beta (b) and cv
  estimates. \item Models The best parameter fit of each
  stock recruitment model and the number (n) and proportion
  (prop) of stochastic recruitment from each model. \item
  srBestFit ssb and rec values for each best fit model,
  facilitates plotting a lines for each model. \item
  recStochastic Contains stochastic recruiment for each
  model. \item recPercentiles The 0.50, 0.05 and 0.95
  percentiles of the stochastic recruitment, grouped along
  10 ssb-groups. }
}
\description{
  The input is a list object that is returned from
  \code{eqsr_fit}. The return object is a list of
  data.frames (see below) that can be used to get e.g. a
  summary plot using the \code{eqsr_plot2} function.

  The function is based on the same algorithm already in
  SRplot, the difference is that here the objects used in
  the plot are returned to the user. The philosophy is to
  separate analysis from plot. Output hopefully facilitates
  user in generating alternative plot to that provided by
  the \code{SRplot} function.
}
\author{
  Einar Hjorleifsson \email{einar.hjorleifsson@gmail.com}
}

