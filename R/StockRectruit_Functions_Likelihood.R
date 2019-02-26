##############################################
# stock recruit formulations
##############################################

#' Get starting values for models
#'
#' a quick fix!!
#'
#' @param model XXX
#' @param data XXX
#' @return vector of starting values
#' @export
initial <- function(model, data)
{
  if (model == "Segreg") {
    c(log(stats::median(data$rec/data$ssb, na.rm = TRUE)), b = log(stats::median(data$ssb)), 0)
  } else if (model == "Smooth_hockey") {
    c(log(stats::median(data$rec/data$ssb, na.rm = TRUE)), b = log(stats::median(data$ssb)), 0)
  } else {
    c(0,0,0)
  }
}


#' stock recruitment function
#'
#'
#' @param ab the model parameters
#' @param ssb a vector of ssb
#' @return log recruitment according to model
#' @export
Ricker <- function(ab, ssb) {
  log(ab$a) + log(ssb) - ab$b * ssb
}

#' stock recruitment function
#'
#'
#' @param ab the model parameters
#' @param ssb a vector of ssb
#' @return log recruitment according to model
#' @export
Segreg <- function(ab, ssb) {
  log(ifelse(ssb >= ab$b, ab$a * ab$b, ab$a * ssb))
}

#' stock recruitment function
#'
#'
#' @param ab the model parameters
#' @param ssb a vector of ssb
#' @return log recruitment according to model
#' @export
Bevholt <- function(ab, ssb) {
  log(ab$a * ssb / (1 + ab$b * ssb))
}


#' stock recruitment function
#'
#'
#' @param ab the model parameters
#' @param ssb a vector of ssb
#' @return log recruitment according to model
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
#' @export
bevholt2 <- function(ab, ssb) {
  log(ab$a * ssb / (ab$b + ssb))
}

#' stock recruitment function
#'
#'
#' @param ab the model parameters
#' @param ssb a vector of ssb
#' @param gamma a smoother parameter
#' @return log recruitment according to model
#' @export
Smooth_hockey <- function(ab, ssb, gamma = 0.1) {
  log(ab$a * (ssb + sqrt(ab$b^2 + gamma^2/4) - sqrt((ssb - ab$b)^2 + gamma^2/4)))
}



#segreg3  <- function(ab, ssb, Blim) {
#  log(ifelse(ssb >= Blim, ab$a * Blim, ab$a * ssb))
#}


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
#' @export
llik <- function(param, data, model, logpar = FALSE)
{
  FUN <- match.fun(model)
  if (logpar) {
    pred <- FUN(list(a = exp(param[1]), b = exp(param[2])), data $ ssb)
    sum( stats::dnorm(log(data $ rec), pred, exp(param[3]), log = TRUE) ) #- sum(param) # add on prior so it is uniform
  } else {
    pred <- FUN(list(a = param[1], b = param[2]), data $ ssb)
    sum( stats::dnorm(log(data $ rec), pred, param[3], log = TRUE) )
  }
}
