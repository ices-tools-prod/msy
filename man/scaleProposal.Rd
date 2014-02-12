\name{scaleProposal}
\alias{scaleProposal}
\title{Propose a multiplicatively uniform variable.}
\usage{
scaleProposal(A)
}
\arguments{
  \item{A}{the scale of the proposal, 2 is a big jump, 1.1
  is a small jump.}
}
\value{
a draw from a symmetrical (on the multiplicative scale)
random variable
}
\description{
Taken from the proposal used for precisions in the C
library GMRFLib by Havard Rue and described in several
papers and his book on GMRFs with Leonard Held.
}
\author{
Colin Millar \email{colin.millar@jrc.ec.europa.eu}
}

