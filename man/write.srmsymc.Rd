\name{write.srmsymc}
\alias{write.srmsymc}
\title{write data files to disk for the ADM srmsymc program}
\usage{
write.srmsymc(rby, iba, aR, col_names = c("year", "r", "ssb"),
  opt_sr_model = 1, opt_sim = 1, opt_pen = 1, opt_age = 1, rba,
  filename_age = "age.dat", aP, years)
}
\arguments{
  \item{rby}{A list that contains in its first three rows
  year, recruitment and ssb.}

  \item{iba}{A list - read.sen}

  \item{aR}{\emph{integer}. Age of recruitment. Has to be
  specified. Is used to shift the recruitment data such
  that it aligns with ssb data. Can be set to 0 (zero) if
  ssb-recruitemnt pairs already properly aligned.}

  \item{col_names}{vector of column names for year,
  recruitment and spawning stock biomass (in that order).}

  \item{opt_sr_model}{\emph{integer}. Number (1-3)
  corresponding to recruitment model.}

  \item{opt_sim}{\emph{integer}. If 0 no simulation, if 1
  simulation. #'}

  \item{opt_pen}{\emph{integer}. If 0 then no recruitement
  constraints if 1 then recruitment contraints applied
  (NEED MORE DETAIL)}

  \item{opt_age}{\emph{integer}. If 0 then error only in
  recruitment, if 1 then both error in recruitment and
  steady-state vectors (age based inputs).}

  \item{rba}{data.frame that contains age based input NOT
  YET IMPLEMENTED}

  \item{filename_age}{\emph{character}. Name of the option
  file that contains the age-based-data. If missing, set to
  "age.dat".}

  \item{aP}{\emph{integer}. Plus group age.}

  \item{years}{vector containing years to include in the
  ssb-r}
}
\description{
Some details, including link to the function a) running the
program and creating the data input
}
\author{
Einar Hjorleifsson <einar.hjorleifsson
}

