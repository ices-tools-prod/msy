\name{convertSumSen}
\alias{convertSumSen}
\title{Convert sum and sen files}
\usage{
  convertSumSen(senfilename = NA, indexfilename = NA,
    pfpm = NA, nits = 0, sr = 2, varybiodata = FALSE,
    stockname = "", silent = FALSE, srconstrain = TRUE)
}
\arguments{
  \item{senfilename}{is the name of the senfile, with a
  corresponding sumfile sharing the same name except
  replacing ".sen" with ".sum". Produces an interactive
  window if not specified}

  \item{indexfilename}{is the name of an index file in the
  Lowestoft format pointing to the pf and pm files. If pfpm
  is specified, this is ignored, if neither specified an
  interactive window appears to choose index file}

  \item{pfpm}{is a vector of two values; pm and pf
  (proportion of m and f before spawning}

  \item{nits}{Number of iterations of bootstrapping - if 0,
  does only the deterministic fit}

  \item{sr}{XXX}

  \item{varybiodata}{XXX}

  \item{stockname}{XXX}

  \item{silent}{XXX}

  \item{srconstrain}{XXX}
}
\value{
  something
}
\description{
  Convert sum and sen files
}
\author{
  Tim Earl \email{timothy.earl@cefas.co.uk}
}

