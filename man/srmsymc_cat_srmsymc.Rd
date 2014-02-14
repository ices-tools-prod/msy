\name{srmsymc_cat_srmsymc}
\alias{srmsymc_cat_srmsymc}
\title{Write srmsymc.dat to disk}
\usage{
  srmsymc_cat_srmsymc(name_stock, filename_age, y1, y2, aR,
    aP, opt_sr_model, opt_sim, opt_age, opt_pen, r, ssb,
    year)
}
\arguments{
  \item{name_stock}{Name of the stock. May later down the
  line be used in table and plot outputs}

  \item{filename_age}{Name of the associated file
  containing biological (weights, maturity and mortalit)
  and fisheries data (selection patter) by age.}

  \item{y1}{First year of ssb and recruitment data}

  \item{y2}{Last year of ssb and recruitment data}

  \item{aR}{Recruitment age}

  \item{aP}{Plus group age}

  \item{opt_sr_model}{Recruitment model number (1: Ricker,
  2: Beverton-Holt, 3: Segmented regression).}

  \item{opt_sim}{0=no simulation, 1=simulation}

  \item{opt_age}{0=error only in recr, 1=error in recr &
  steady-state vectors}

  \item{opt_pen}{0=no SR constraints, 1=apply SR constrain
  (IS THIS ONLY FOR THE SEGMENTED REGRESSION??).}

  \item{r}{A vector containing recruitment}

  \item{ssb}{A vector containing spawning stock biomass}

  \item{year}{A vector containinb years (just used in
  comments).}
}
\description{
  The function is used to write srmsymc.dat to the disk.
}
\author{
  Einar Hjorleifsson
}

