################################################################################
# Einar 2015-03-13
# Note that most of this stuff that was here originally (prior to wkmsyref2)
# is now in file eqsim.R
# Should consider moving other stuff to more appropriate files (if needed)
# and then delete this file.


##############################################
# utility functions
##############################################



#' A function only available for R 2.15.1
#'
#' for back compatibility
#'
#' @param ... whatever
#' @return result of paste with sep = ""
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
paste0 <- function(...) paste(..., sep = "")


#' @title Progress function
#'
#'
#' @param p Value
#' @return NULL
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
loader <- function(p)
{
  if (p==0) cat("0%                       50%                     100%\n")
 str <- paste0(rep(c("\r[", "=", ">", " ", "]"), c(1, floor(p*50), 1, 50 - floor(p*50), 1)), collapse = "")  
 cat(str); flush.console()
 if (floor(p) == 1) cat("\n")
}





############################################
# more plotting routines
##############################################
#' plots the ordered posterior densities accross MCMC iterations for each model considered in \code{fitModels}
#'
#'
#' @param fit an fitted MCMC returned from \code{fitModels}
#' @return NULL produces a plot
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
LLplot <- function(fit) 
{
  lliks <- sapply(fit $ fits, function(x) sort(x $ llik))

  plot(0, 0, type = "n", 
       main = paste("LL of Bayes model:", fit $ stknam), 
       xlab = "model order", ylab = "log likelihood", 
       ylim = range(lliks), xlim = c(1, nrow(lliks)))
  for (i in 2:ncol(lliks)-1) 
  {
    lines(lliks[,i], lty = i, col = i)
  }
  lines(sort(fit $ fits $ BMA $ llik), lwd = 2)
  legend(x = "bottomright", legend = c(names(fit $ fits)), 
         lty = c(2:ncol(lliks)-1,1), col = c(2:ncol(lliks)-1,1), lwd = c(rep(1, ncol(lliks)-1),2)) 
}


#' plots the fits for each model considered in \code{fitModels}
#'
#'
#' @param fit an fitted MCMC returned from \code{fitModels}
#' @return NULL produces one or several plots
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
LLmodplot <- function(fit)
{
  rec <- fit $ data $ rec[order(fit $ data $ ssb)]
  ssb <- sort(fit $ data $ ssb)

  for (mod in seq(fit $ fits)[2:length(fit $ fits)-1]) 
  { 
    model <- paste(names(fit $ fits)[mod], "fit")
    Rsym <- fit $ preds[[mod]]

    plot(ssb, rec, xlab = 'SSB', ylab = 'Recruits', 
         main = paste(fit $ stknam, model), 
         ylim = c(0, max(rec)), xlim = c(0, max(ssb)), type = "n")
    for (i in sample(1:nrow(Rsym), 1000)) {
      lines(ssb, Rsym[i,], col = paste0(grey(0), "10") )
    }
    for (i in c(0.05, 0.95)) {
      lines(ssb, apply(Rsym, 2, quantile, i), col = 4, lwd = 3)
    }
    lines(ssb, apply(Rsym, 2, quantile, 0.5), col = 7, lwd = 3)
    #TODO lines(ssb, Rsym[which.max(res[[mod]]),], col=1, lwd=2)
    lines(fit $ data $ ssb, fit $ data $ rec, col = "red")
    points(ssb, rec, pch = 19, col = "red")  
  }
}

