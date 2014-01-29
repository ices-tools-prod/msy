\name{srmsymc_run}
\alias{srmsymc_run}
\title{Run srmsymc}
\usage{
  srmsymc_run(sr, opt_sen = 1, opt_pen = 1, nits = 100,
    path, output = FALSE, do_clean = TRUE, echo = TRUE)
}
\arguments{
  \item{sr}{integer. containg the stock-recruitment model}

  \item{opt_sen}{XXX}

  \item{opt_pen}{XXX}

  \item{nits}{Number of iterations of bootstrapping - if 0,
  does only the deterministic fit}

  \item{path}{characters. Name of the directory for storing
  the output. If missing, the output remains in the working
  directory.}

  \item{output}{boolean. If FALSE (default) no results
  (*.dat) files are returned.}

  \item{do_clean}{XXX}

  \item{echo}{boolean. If TRUE (default) output stuff to
  console. Setting to FALSE is useful when using function
  while knitting.}
}
\description{
  Run the ADM srmsymc program
}
\author{
  Einar Hjorleifsson <einar.hjorleifsson
}

