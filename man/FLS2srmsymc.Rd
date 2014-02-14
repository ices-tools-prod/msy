\name{FLS2srmsymc}
\alias{FLS2srmsymc}
\title{FLR to srmsymc}
\usage{
  FLS2srmsymc(stk, y = 3, name_stock, filename_age,
    opt_sr_model, opt_land = TRUE, opt_sim = TRUE,
    opt_age = TRUE, opt_pen = TRUE)
}
\arguments{
  \item{stk}{FLStock object}

  \item{y}{Years}

  \item{name_stock}{Name of the stock. May later down the
  line be used in table and plot outputs}

  \item{filename_age}{Name of the associated file
  containing biological (weights, maturity and mortalit)
  and fisheries data (selection patter) by age.}

  \item{opt_sr_model}{Recruitment model number (1: Ricker,
  2: Beverton-Holt, 3: Segmented regression).}

  \item{opt_land}{Boolean}

  \item{opt_sim}{0=no simulation, 1=simulation}

  \item{opt_age}{0=error only in recr, 1=error in recr &
  steady-state vectors}

  \item{opt_pen}{0=no SR constraints, 1=apply SR constrain
  (IS THIS ONLY FOR THE SEGMENTED REGRESSION??).}
}
\description{
  XXX
}

