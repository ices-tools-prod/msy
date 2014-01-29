#' @title Read key parameter estimates from the srmsymc program
#' 
#' @description Reads the \code{simpar.dat} generate by the \code{srmsymc}.
#' Returns a list containing the stochastic and the deterministic estimates.
#' 
#' @export
#' 
#' @param path character containing the path to the file. If missing it is
#' assumed that the file is in the current working directory.
#' @param file the file name, should not be anything else than \code{simpar.dat}
#' @param srno integer. 1=Ricker, 2=Beverton-Holt, 3=Hockey-stick.
#' @param minSSB numerical. Minimum observed SSB value. Used for the hockey-stick,
#' excluding values below the lowest observed. If set to 0 (default) do not exclude
#' values below the minimum.
#' @param trim boolean. If TRUE (default) come cleaning (NEED TO SPECIFY) is done.
#' @param longformat boolean. If TRUE (default) return a long format.


srmsymc_read_par <- function(path,file="simpar.dat",srno=1,minSSB=0,trim=TRUE,longformat=TRUE) {
  
  if(!missing(path)) file <- paste(path, "/", file, sep="")
  if (!file.exists(file)) stop(paste(file, "not found"))
  
  
  
  d = read.table(file)
  names(d) = c("iter","ap","bp","alpha","beta","sigr","scor",
               "fcrash","fmax","f01","f20","f25","f30","f35","f40","fmsy",
               "msy","bmsy","msypr","bmsypr","msyr",
               "fnfcrash","fnfmax","fnf01","fnf20","fnf15",
               "fnf30","fnf35","fnf40","dydf","penalty","nll","AIC")
  
  d.det <- subset(d, iter == 0)
  d.det$iter <- 0
  
  d <- subset(d, iter != 0)
  d$iter <- c(1:nrow(d))
  # get rid of negative/unresonable parameters
  if(trim) {
    #  
    if (srno==3) {
      i <- d$beta < minSSB | d$alpha <= 0 | d$beta <= 0
    } else {
      i <- d$alpha <= 0 | d$beta <= 0
    }
    cn <- c("ap","bp","alpha","beta","sigr","scor","fcrash","fmsy",
            "msy","bmsy","fnfcrash","dydf","penalty","nll","AIC")
    d[i,cn] <- NA
    
    # Make msy-points as NA if fmsy is greater than 3 
    if (any(!d$fmsy < 3 , na.rm=TRUE)) d[which(!d$fmsy < 3),c("fmsy","bmsy","msy")] <- NA
    # Make Fcrash-points as NA if above 5
    if (any(!d$fcrash<5, na.rm=TRUE)) d[which(!d$fcrash<5),"fcrash"] <- NA
    # Similar for Fmas greater than 3
    if (any(!d$fmax<3 , na.rm=TRUE)) d[which(!d$fmax<3),c("fmax","bmsypr","msypr")] <- NA
    # columns of outputpr[[srtype]] to be output to output[[srtype]]
    indexmap = c("f20","f25","f30","f35","f40","f01") 
    for (i in indexmap) if (any(!d[,i] < 3, na.rm=TRUE)) d[which(!d[,i]<3),i] <- NA
  }
  
  if(longformat) {
    d <- melt(d,id.vars='iter')
    d.det <- melt(d.det,id.vars='iter')
  }
  
  return(list(stochastic=d,
              deterministic=d.det))
}

#' @title Read the yield estimates the srmsymc program
#' 
#' @description Reads the file \code{simpary.dat} generate by the \code{srmsymc}.
#' In the \code{simpary.dat} the first row is the fishing mortality,
#' the second row is the deterministic estimates of yield and the remaining
#' rows represents the catch for each iterations.
#' Returns a list containing two data.frames, one containing the deterministic
#' estimate, the other containing the stochastic simulations.
#' @param path character containing the path to the file. If missing it is
#' assumed that the file is in the current working directory.
#' @param file the file name, should not be anything else than \code{simpary.dat}
#' @param srno integer. 1=Ricker, 2=Beverton-Holt, 3='Hockey-stick.
#' @param trim boolean. If TRUE (default) come cleaning (NEED TO SPECIFY) is done.
#' @param longformat boolean. If TRUE (default) return a long format.
#' @param doPlot boolean. If TRUE (default) returns a ggplot
#' @export

srmsymc_read_yield <- function(path,file="simpary.dat",srno=1,trim=TRUE,longformat=TRUE,doPlot=TRUE) {
  
  if(!missing(path)) file <- paste(path, "/", file, sep="")
  if (!file.exists(file)) stop(paste(file, "not found"))
  
  d = as.matrix(read.table(file))
  
  if(trim) {
    d[-1,] = pmax(d[-1,],0)
    noredlines = prod(dim(d)==0)
  }
  
  colNames <- d[1,]
  d <- as.data.frame(d[-1,])
  names(d) <- colNames
  d.det <- d[1,]
  d.det$iter <- 0
  
  d <- d[-1,]
  d$iter <- c(1:nrow(d))
  
  if(longformat) {
    d <- melt(d,id.vars='iter')
    d$variable <- as.numeric(as.character(d$variable))
    d.det <- melt(d.det,id.vars='iter')
    d.det$variable <- as.numeric(as.character(d.det$variable))
    
    x <- calc.quantiles(d,d.det=d.det)
    
    if(doPlot) {
      gg.plot <- do.ggplot(x,d)
      return(list(deterministic = d.det,
                  stochastic = d,
                  quantiles = x,
                  gg.plot = gg.plot))
    } else {
      return(list(deterministic = d.det,
                  stochastic = d,
                  quantiles = x))
    }
  }
  
  if(!longformat) stop('Not yet implemented for wide format')
  
}


#' @title Read the ssb estimates the srmsymc program
#' 
#' @description Reads the file \code{simpary.dat} generate by the \code{srmsymc}.
#' In the \code{simpary.dat} the first row is the fishing mortality,
#' the second row is the deterministic estimates of ssb and the remaining
#' rows represents the catch for each iterations.
#' Returns a list containing two data.frames, one containing the deterministic
#' estimate, the other containing the stochastic simulations.
#' 
#' @param path character containing the path to the file. If missing it is
#' assumed that the file is in the current working directory.
#' @param file the file name, should not be anything else than \code{simpary.dat}
#' @param srno integer. 1=Ricker, 2=Beverton-Holt, 3='Hockey-stick.
#' @param trim boolean. If TRUE (default) come cleaning (NEED TO SPECIFY) is done.
#' @param longformat boolean. If TRUE (default) return a long format.
#' @param doPlot boolean. If TRUE (default) returns a ggplot
#' @export

srmsymc_read_ssb <- function(path,file="simparssb.dat",srno=1,trim=TRUE,longformat=TRUE,doPlot=TRUE) {
  
  if(!missing(path)) file <- paste(path, "/", file, sep="")
  if (!file.exists(file)) stop(paste(file, "not found"))
  
  d = as.matrix(read.table(file))
  
  if(trim) {
    d[-1,] = pmax(d[-1,],0)
    noredlines = prod(dim(d)==0)
  }
  
  colNames <- d[1,]
  d <- as.data.frame(d[-1,])
  names(d) <- colNames
  d.det <- d[1,]
  d.det$iter <- 0
  
  d <- d[-1,]
  d$iter <- c(1:nrow(d))
  
  if(longformat) {
    d <- melt(d,id.vars='iter')
    d$variable <- as.numeric(as.character(d$variable))
    d.det <- melt(d.det,id.vars='iter')
    d.det$variable <- as.numeric(as.character(d.det$variable))
    
    x <- calc.quantiles(d,d.det=d.det)
    
    if(doPlot) {
      gg.plot <- do.ggplot(x,d)
      return(list(deterministic = d.det,
                  stochastic = d,
                  quantiles = x,
                  gg.plot = gg.plot))
    } else {
      return(list(deterministic = d.det,
                  stochastic = d,
                  quantiles = x))
    }
  }
  
  if(!longformat) stop('Not yet implemented for wide format')
}





#' @title Read all output files from srmsymc program
#' 
#' @description xxx
#' 
#' @export
#' 
#' @param path x
#' @param ... addition arguements

srmsymc_read <- function(path, ...) {
  x <- list()
  x$par   <- srmsymc_read_par(path=path)
  x$yield <- srmsymc_read_yield(path=path)
  x$ssb   <- srmsymc_read_ssb(path=path)
  return(x)
}



#' @title Read the ssb per recruit estimates from the srmsymc program
#' 
#' @description Reads the file \code{simparssbpr.dat} generate by the \code{srmsymc}.
#' In the \code{simparssbpr.dat} the first row is the fishing mortality,
#' the second row is the deterministic estimates of ssb and the remaining
#' rows represents the catch for each iterations.
#' Returns a list containing two data.frames, one containing the deterministic
#' estimate, the other containing the stochastic simulations.
#' 
#' @export
#' 
#' @param file the file name, should not be anything else than \code{simparssbpr.dat}
#' @param trim boolean. If TRUE (default) come cleaning (NEED TO SPECIFY) is done.
#' @param longformat boolean. If TRUE (default) return a long format.
#' @param doPlot boolean. If TRUE (default) returns a ggplot


srmsymc_read_ssbpr <- function(file="simparssbpr.dat",trim=TRUE,longformat=TRUE,doPlot=TRUE) {
  
  if (!file.exists(file)) stop(paste(file, "not found"))
  
  d = as.matrix(read.table(file))
  
  if(trim) {
    d[-1,] = pmax(d[-1,],0)
    noredlines = prod(dim(d)==0)
  }
  
  colNames <- d[1,]
  d <- as.data.frame(d[-1,])
  names(d) <- colNames
  d.det <- d[1,]
  d.det$iter <- 0
  
  d <- d[-1,]
  d$iter <- c(1:nrow(d))
  
  if(longformat) {
    d <- melt(d,id.vars='iter')
    d$variable <- as.numeric(as.character(d$variable))
    d.det <- melt(d.det,id.vars='iter')
    d.det$variable <- as.numeric(as.character(d.det$variable))
    
    x <- calc.quantiles(d,d.det=d.det)
    
    if(doPlot) {
      gg.plot <- do.ggplot(x,d)
      return(list(deterministic = d.det,
                  stochastic = d,
                  quantiles = x,
                  gg.plot = gg.plot))
    } else {
      return(list(deterministic = d.det,
                  stochastic = d,
                  quantiles = x))
    }
  }
  
  if(!longformat) stop('Not yet implemented for wide format')
}

#' @title Read the yield per recruit estimates the srmsymc program
#' 
#' @description Reads the file \code{simparypr.dat} generate by the \code{srmsymc}.
#' In the \code{simparypr.dat} the first row is the fishing mortality,
#' the second row is the deterministic estimates of ssb and the remaining
#' rows represents the catch for each iterations.
#' Returns a list containing two data.frames, one containing the deterministic
#' estimate, the other containing the stochastic simulations.
#' 
#' @param file the file name, should not be anything else than \code{simparypr.dat}
#' @param trim boolean. If TRUE (default) come cleaning (NEED TO SPECIFY) is done.
#' @param longformat boolean. If TRUE (default) return a long format.
#' @param doPlot boolean. If TRUE (default) returns a ggplot
#' @export

srmsymc_read_ypr <- function(file="simparypr.dat",trim=TRUE,longformat=TRUE,doPlot=TRUE) {
  
  if (!file.exists(file)) stop(paste(file, "not found"))
  
  d = as.matrix(read.table(file))
  
  if(trim) {
    d[-1,] = pmax(d[-1,],0)
    noredlines = prod(dim(d)==0)
  }
  
  colNames <- d[1,]
  d <- as.data.frame(d[-1,])
  names(d) <- colNames
  d.det <- d[1,]
  d.det$iter <- 0
  
  d <- d[-1,]
  d$iter <- c(1:nrow(d))
  
  if(longformat) {
    d <- melt(d,id.vars='iter')
    d$variable <- as.numeric(as.character(d$variable))
    d.det <- melt(d.det,id.vars='iter')
    d.det$variable <- as.numeric(as.character(d.det$variable))
    
    x <- calc.quantiles(d,d.det=d.det)
    
    if(doPlot) {
      gg.plot <- do.ggplot(x,d)
      return(list(deterministic = d.det,
                  stochastic = d,
                  quantiles = x,
                  gg.plot = gg.plot))
    } else {
      return(list(deterministic = d.det,
                  stochastic = d,
                  quantiles = x))
    }
  }
  
  if(!longformat) stop('Not yet implemented for wide format')
}




