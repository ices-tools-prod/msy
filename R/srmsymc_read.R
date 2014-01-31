#' @title srmsymc_ricker
#' 
#' @param ssb spawning stock biomass
#' @param alpha alpha
#' @param beta beta
srmsymc_ricker <- function(ssb,alpha,beta) {alpha * ssb*exp(-beta * ssb)}
#' @title srmsymc_bevholt
#' 
#' @param ssb spawning stock biomass
#' @param alpha alpha
#' @param beta beta
srmsymc_bevholt <- function(ssb,alpha,beta) {alpha * ssb /(beta + ssb)}
#' @title srmsymc_segreg
#' 
#' @param ssb spawning stock biomass
#' @param alpha alpha
#' @param beta beta
srmsymc_segreg <- function(ssb,alpha,beta) {alpha*(ssb+sqrt(beta^2+0.001)-sqrt((ssb-beta)^2+0.001))}

#' @title Read the srmsymc setup
#' 
#' @description Reads the setup from the srmsymc.dat file
#' 
#' @export
#' 
#' @param path The path to the directory that stores the srmsymc results

srmsymc_read_setup <- function(path) 
  {
  x <- read.table(paste(path,"srmsymc.dat",sep="/"),header=FALSE,nrows=10)
  res <- list()
  res$name_stock <- x[1,]
  res$filename_age <- x[2,]
  res$y1 <- as.integer(x[3,])
  res$y2 <- as.integer(x[4,])
  res$aR <- as.integer(x[5,])
  res$aP <- as.integer(x[6,])
  res$opt_sr_model <- as.integer(x[7,])
  res$opt_sim <- as.integer(x[8,])
  res$opt_age <- as.integer(x[9,])
  res$opt_pen <- as.integer(x[10,])
  return(res)
}

#' @title Recruitment model
#' 
#' @description  Read recruitment model from the srmsymc.dat file
#' 
#' @param path XXX
srmsymc_readRecModel <- function(path) 
  {
  tmpFile <- paste(path,"/","srmsymc.dat",sep="")
  if(file.exists(tmpFile)) {
    x <- readLines(tmpFile)
    i <- grep("# Ropt:   S-R function type",x)
    srno <- as.integer(substr(x[i],1,1))
  } else {
    stop("srno needs to be specified")
  }
  return(srno)
}


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
#' @param minSSB numerical. Minimum observed SSB value. Used for the hockey-stick,
#' excluding values below the lowest observed. If set to 0 (default) do not exclude
#' values below the minimum.
#' @param trim boolean. If TRUE (default) come cleaning (NEED TO SPECIFY) is done.
#' @param longformat boolean. If TRUE (default) return a long format.

srmsymc_read_par <- function(path,file="simpar.dat",minSSB=0,trim=TRUE,longformat=TRUE) {
  
  if(!missing(path)) file <- paste(path, "/", file, sep="")
  if (!file.exists(file)) stop(paste(file, "not found"))
  
  # Get the recruitment model number used
  srno <- srmsymc_readRecModel(path)
  
  d = read.table(file)
  names(d) = c("iter","ap","bp","alpha","beta","sigr","scor",
               "fcrash","fmax","f01","f20","f25","f30","f35","f40","fmsy",
               "msy","bmsy","msypr","bmsypr","msyr",
               "fnfcrash","fnfmax","fnf01","fnf20","fnf15",
               "fnf30","fnf35","fnf40","dydf","penalty","nll","AIC")
  
  d.det <- subset(d, iter == 0)

  d <- subset(d, iter != 0)
  #d$iter <- c(1:nrow(d))

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
    # Similar for Fmax greater than 3
    if (any(!d$fmax<3 , na.rm=TRUE)) d[which(!d$fmax<3),c("fmax","bmsypr","msypr")] <- NA
    # columns of outputpr[[srtype]] to be output to output[[srtype]]
    columns = c("f20","f25","f30","f35","f40","f01") 
    for (cn in columns) if (any(!d[,cn] < 3, na.rm=TRUE)) d[which(!d[,cn]<3),cn] <- NA
  }
  
  d$srno <- srno
  
  if(longformat) {
    d <- melt(d,id.vars=c('iter','srno'))
    d.det <- melt(d.det,id.vars=c('iter','srno'))
  }
  
  
  return(list(stochastic=d,
              deterministic=d.det,
              srno=srno))
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
#' @param trim boolean. If TRUE (default) come cleaning (NEED TO SPECIFY) is done.
#' @param longformat boolean. If TRUE (default) return a long format.
#' @param doPlot boolean. If TRUE (default) returns a ggplot
#' @export

srmsymc_read_yield <- function(path,file="simpary.dat",trim=TRUE,longformat=TRUE,doPlot=TRUE) {
  
  if(!missing(path)) file <- paste(path, "/", file, sep="")
  if (!file.exists(file)) stop(paste(file, "not found"))
  
  # Get the recruitment number from the srmsymc.dat file
  srno <- srmsymc_readRecModel(path=path)
  
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
  d.det$srno <- srno
  
  d <- d[-1,]
  d$iter <- c(1:nrow(d))
  d$srno <- srno
  
  if(longformat) {
    d <- melt(d,id.vars=c('iter','srno'))
    d$variable <- as.numeric(as.character(d$variable))
    d.det <- melt(d.det,id.vars=c('iter','srno'))
    d.det$variable <- as.numeric(as.character(d.det$variable))

    x <- calc.quantiles(d,d.det=d.det)
    x$srno <- srno
    
    if(doPlot) {
      gg.plot <- do.ggplot(x,d)
      return(list(deterministic = d.det,
                  stochastic = d,
                  quantiles = x,
                  gg.plot = gg.plot,
                  srno=srno))
    } else {
      return(list(deterministic = d.det,
                  stochastic = d,
                  quantiles = x,
                  srno=srno))
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
  
  # Get the recruitment number from the srmsymc.dat file
  srno <- srmsymc_readRecModel(path=path)
  
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
  
  d$srno <- srno
  d.det$srno <- srno
  
  if(longformat) {
    d <- melt(d,id.vars=c('iter','srno'))
    d$variable <- as.numeric(as.character(d$variable))
    d.det <- melt(d.det,id.vars=c('iter','srno'))
    d.det$variable <- as.numeric(as.character(d.det$variable))
    
    x <- calc.quantiles(d,d.det=d.det)
    
    x$srno <- srno
    
    if(doPlot) {
      gg.plot <- do.ggplot(x,d)
      return(list(deterministic = d.det,
                  stochastic = d,
                  quantiles = x,
                  gg.plot = gg.plot,
                  srno=srno))
    } else {
      return(list(deterministic = d.det,
                  stochastic = d,
                  quantiles = x,
                  srno=srno))
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

#' @title Combine the parameters from the three sr-functions
#' 
#' @description XXX
#' 
#' @export
#' 
#' @param srweights a vector of length 3 indicating the weights for each srmodel.
#' If missing then equal weights applied to all models.

srmsymc_combineParameterFiles <- function(srweights=NA) {
  x1 <- srmsymc_read_par("ricker",longformat=FALSE)$stochastic
  x1$srno <- srmsymc_read_par("ricker",longformat=FALSE)$srno
  x2 <- srmsymc_read_par("bevholt",longformat=FALSE)$stochastic
  x2$srno <- srmsymc_read_par("bevholt",longformat=FALSE)$srno
  x3 <- srmsymc_read_par("segreg",longformat=FALSE)$stochastic
  x3$srno <- srmsymc_read_par("segreg",longformat=FALSE)$srno
  if(any(is.na(srweights))) {
    dat <- rbind(x1,x2,x3)
  } else {
    dat <- rbind(x1[1:srweights[1],],
                 x2[1:srweights[2],],
                 x3[1:srweights[3],])
  }
  return(dat)
}

#' @title Estimating weights for each recruitment function
#' 
#' @description XXX
#' 
#' @export
#' 
#' @param srautoweights XXX
#' @param trimming XXX
#' @param nits Number of iterations

srmsymc_srweights <- function(srautoweights=TRUE,
                              trimming=NA,
                              nits=100)
{
  simdata <- list()
  simdata[[1]] <- srmsymc_read_par("ricker",longformat=FALSE)$stochastic
  simdata[[2]] <- srmsymc_read_par("bevholt",longformat=FALSE)$stochastic
  simdata[[3]] <- srmsymc_read_par("segreg",longformat=FALSE)$stochastic
  
  if (srautoweights) {  
    lik <- sapply(simdata, function(x){sort(exp(-x$nll),na.last=FALSE)})   
    if (is.na(trimming)) 
    {
      trim <- 0:(nits/4)
      srwall <- sapply(trim, function(p){1/colMeans(1/lik[(p+1):nits,],na.rm=TRUE)} )*ifelse(is.na(srweights),1,0)
      srwall[is.na(srwall)] <- 0
      srwall <- t(t(srwall)/colSums(srwall))
      rownames(srwall) <- c("Ricker","Beverton-Holt","Smooth hockeystick")
      
      #png(paste0(outputfolder, stockname, "_trim_diag.png"),800,600)
      #plot(100*trim/nits,srwall["Ricker",],ylim=range(0,srwall), lty=1, type='l', 
      #     ylab="Relative weight", xlab="Trimmed percentage", main = paste(stockname,"Trimming the harmonic mean",sep=" - "))
      #lines(100*trim/nits,srwall["Beverton-Holt",],lty=2)
      #lines(100*trim/nits,srwall["Smooth hockeystick",],lty=3)          
      #legend("topright",pch=NA,lty=1:3,legend=srname)
      #dev.off()
      #write.csv(t(rbind(Percentage=100*trim/nits,srwall)),paste0(outputfolder,stockname,"_trim_diag.csv"),row.names=FALSE)
      trimming <- 0
    }
    
    srweightsexact <- sapply(round(trimming*nits), function(p){1/colMeans(1/lik[(p+1):nits,],na.rm=TRUE)} )*ifelse(is.na(srweights),1,0)
    
    srweightsexact[is.na(srweightsexact)] <- 0
    srweightsexact <- srweightsexact/sum(srweightsexact)    
    srweights <- round(srweightsexact*nits)
  } else {
    if (is.na(trimming)) 
    {
      warning("Trimming parameter ignored when applying manual weights")
    }
    srweightsexact <- srweights
    srweights <- round(srweightsexact*nits)
    return(srweights)
  }
}
