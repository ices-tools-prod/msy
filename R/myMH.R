

##############################################
# stock recruit formulations
##############################################


#' stock recruitment function
#'
#'
#' @param ab the model parameters
#' @param ssb a vector of ssb
#' @return log recruitment according to model
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
ricker <- function(ab, ssb) {
  log(ab$a) + log(ssb) - ab$b * ssb
}

#' stock recruitment function
#'
#'
#' @param ab the model parameters
#' @param ssb a vector of ssb
#' @return log recruitment according to model
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
segreg <- function(ab, ssb) {
  log(ifelse(ssb >= ab$b, ab$a * ab$b, ab$a * ssb))
}

#' stock recruitment function
#'
#'
#' @param ab the model parameters
#' @param ssb a vector of ssb
#' @return log recruitment according to model
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
bevholt <- function(ab, ssb) {
  log(ab$a * ssb / (1 + ab$b * ssb))
}


#' stock recruitment function
#'
#'
#' @param ab the model parameters
#' @param ssb a vector of ssb
#' @return log recruitment according to model
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
segreg2 <- function(ab, ssb) {
  log(ifelse(ssb >= ab$b, ab$a, ab$a / ab$b * ssb))
}

#' stock recruitment function
#'
#'
#' @param ab the model parameters
#' @param ssb a vector of ssb
#' @return log recruitment according to model
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
bevholt2 <- function(ab, ssb) {
  log(ab$a * ssb / (ab$b + ssb))
}



######
# This function calculates the log likelihood of the stock recruit relationship
######

#' the log-likelihood of the rectuit function
#'
#'
#' @param param the model parameters
#' @param data the rec and ssb data
#' @param model the stock recruit model to use 
#' @param logpar are the parameters on the log scale
#' @return the log-likelihood
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
llik <- function(param, data, model, logpar = FALSE) 
{
  FUN <- match.fun(model)
  if (logpar) {
    pred <- FUN(list(a = exp(param[1]), b = exp(param[2])), data $ ssb)
    sum( dnorm(log(data $ rec), pred, exp(param[3]), log = TRUE) ) #- sum(param) # add on prior so it is uniform
  } else {
    pred <- FUN(list(a = param[1], b = param[2]), data $ ssb)
    sum( dnorm(log(data $ rec), pred, param[3], log = TRUE) )
  }
}



######
# Propose a multiplicatively uniform variable
######
#' Propose a multiplicatively uniform variable.
#'
#' Taken from the proposal used for precisions in the C library GMRFLib by Havard Rue  
#' and described in several papers and his book on GMRFs with Leonard Held.
#'
#' @param A the scale of the proposal, 2 is a big jump, 1.1 is a small jump.
#' @return a draw from a symmetrical (on the multiplicative scale) random variable
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
scaleProposal <- function(A)
{
    if (A <= 1.0) {
        1.0
    } else {
        len = A - 1 / A
        if (runif(1) < len / (len + 2 * log(A))) {
            1 / A + len * runif(1)
        } else {
            A^(2.0 * runif(1) - 1)
        }
    }
}



#' performs a single
#'
#'
#' @param nparam the x location of the bubbles
#' @param param the y location of the bubbles
#' @param ni a vector of bubble sizes
#' @param data a scaling factor 
#' @param llikhood  extra arguments to plot and points
#' @param delta  extra arguments to plot and points
#' @param model  extra arguments to plot and points
#' @return a pointer to the environment in which summaries of the data reside
#' @note \code{extractData} is intended to be used internally
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
updateparam <- function(param, data, llikhood, delta, model)
{
# Keep a record of the current parameter value being updated

  oldparam <- param
  param <- param * replicate(length(param), scaleProposal(delta))  # multiplicative
  newllikhood <- llik(param, data, model)

# MH step:

  if (runif(1) <= exp(newllikhood - llikhood)) { # no need for min(1, ...) 
# Accept the proposed move:
    llikhood <- newllikhood
  } else {
    param <- oldparam
  }

  c(llikhood, param)
}




#' Runs a Metropolis Hasting MCMC on a stock recruit model
#'
#'
#' @param nt The total number of MCMC iterations to run
#' @param nburn The number of initial iterations to remove
#' @param data a list or data.frame with elements \code{ssb} and \code{rec}
#' @param delta a scaling factor for the MH proposal 
#' @param model the recruitment model to use
#' @return a data.frame of nt - nburn rows and columns for the model parameters and the log-likelihood
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
MH <- function(nt, nburn, data, model, delta = 1.3) 
{

# scale data
  sdata <- data
  #sdata $ ssb <- sdata $ ssb / exp(mean(log(sdata $ ssb)))
  #sdata $ rec <- sdata $ rec / exp(mean(log(sdata $ rec)))

# Set initial parameter values - use maximum likelihood estimates:

  opt <- optim(rep(0, 3), llik, data = sdata, model = model, logpar = TRUE, control = list(fnscale = -1)) 
  param <- exp(opt $ par)

# sample is an array in which we put the sample from the posterior distribution.

  sample <- array(0, dim=c(nt, 4))
  sample <- data.frame(sample)
  names(sample) <- c("llik", "a", "b", "cv")

# Calculate log(likelihood) for initial state using a separate function "llik":

  llikhood <- llik(param, sdata, model)

# MCMC updates - MH algorithm:

  for (t in 1:nt) 
  {
    oldparam <- param
    param <- param * replicate(length(param), scaleProposal(delta))  # multiplicative
    newllikhood <- llik(param, sdata, model)

# MH step:

    if (runif(1) <= exp(newllikhood - llikhood)) { # no need for min(1, ...) 
# Accept the proposed move:
      llikhood <- newllikhood
    } else {
      param <- oldparam
    }

# Record the set of parameter values:
    sample[t,] <- c(llikhood, param)
  }

# Calculate the mean and standard deviation of the parameters
# following burn-in:
  subsample <- sample[(nburn+1):nt,]

# rescale parameters and unlog if log transform was used
  # DO IT!

# a crudish approximation to acceptance rate
  cat(model, "acceptance rate:", round(mean(diff(subsample[,1]) != 0), 3), ", try for 0.40\n")

  subsample
}

