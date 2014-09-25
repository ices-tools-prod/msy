#' @title Compile the srmsymc ADMB code
#' @param windose Whether in windows or not
#' @param compile Compile binaries from source. TRUE requires ADMB to be installed, FALSE only works on Windows

#' 
#' @description The \code{plotMSY} approach uses a ADModel Builder code. ...
#' 
#' @export
#' 

srmsymc_compile <- function(windose=TRUE,compile=TRUE) 
{
  if(!(compile | windose)) stop("Must compile in non-Windows operating systems")
  
  if(compile)
  {
    file.copy(paste(path.package("msy"),'tpl/srmsymc.tpl',sep='/'),"srmsymc.tpl")
    file.copy(paste(path.package("msy"),'tpl/srmsymc2.tpl',sep='/'),"srmsymc2.tpl")    
    compile_admb("srmsymc")
    clean_admb("srmsymc")
    compile_admb("srmsymc2")
    clean_admb("srmsymc2")
  } else {
    file.copy(paste(path.package("msy"),'bin/srmsymc.exe',sep='/'),"srmsymc.exe")
    file.copy(paste(path.package("msy"),'bin/srmsymc2.exe',sep='/'),"srmsymc2.exe")
  }
}

#' @title Run srmsymc
#' 
#' @description Run the ADM srmsymc program
#' 
#' @export
#' 
#' @author Einar Hjorleifsson <einar.hjorleifsson@@gmail.com>
#' 
#' @param sr integer. containg the stock-recruitment model
#' @param opt_sen XXX
#' @param opt_pen XXX
#' @param nits Number of iterations of bootstrapping - if 0, does only the deterministic fit
#' @param path characters. Name of the directory for storing the output. If
#' missing, the output remains in the working directory.
#' @param output boolean. If FALSE (default) no results (*.dat) files are returned.
#' @param do_clean XXX
#' @param echo boolean. If TRUE (default) output stuff to console. Setting to FALSE is
#' useful when using function while knitting.
#' @param windose boolean - a quick fix for now

srmsymc_run <- function(sr,opt_sen=1,opt_pen=1,nits=100,path,output=FALSE,
                        do_clean=TRUE,echo=TRUE,windose=TRUE) {
  
  # check if file exists
  if (!file.exists('srmsymc')) srmsymc_compile()
  if (!file.exists('srmsymc.dat')) stop('The srmsymc.dat file must be in the working directory')
  if (!file.exists('age.dat')) stop('The age.dat file must be in the working directory')
  
  # update the .dat file
  # Update the recruitment model number
  x <- readLines("srmsymc.dat")
  i <- grep("# Ropt:   S-R function type",x)
  x[i] <- paste(sr,"  # Ropt:   S-R function type")
  # Update the options
  i <- grep("# senopt", x)
  x[i] <- paste(opt_sen,"# senopt")
  i <- grep("# penopt", x)
  x[i] <- paste(opt_pen,"# penopt")
  write(x,file="srmsymc.dat")
    
  tmpfile = file("sim.dat","w")
  cat('21 1 ', nits, "\n", sep="", file = tmpfile)
  close(tmpfile)
  
  if(windose) {
    programName <- "srmsymc.exe"
  } else {
    programName <- "srmsymc"
  }
  
  cmd <- paste(programName,'-mcmc', (nits+11)*10000, '-mcsave 10000 -nosdmcmc')
  # -nosdmcmc turn off mcmc histogram calcs to make mcsave run faster
  if(!echo) cmd <- paste(cmd,'&> /dev/null')
  system(cmd)
  cmd <- paste(programName," -mceval")
  if(!echo) cmd <- paste(cmd,'&> /dev/null')
  system(cmd)
  
  if(!missing(path)) {
    
    if (!file.exists(path)) {
      cmd <- paste("mkdir",path)
      system(cmd)
    }
    
    cmd <- paste("cp *.dat ",path,"/.",sep="")
    system(cmd)
  }
  
  if(output) {
    x <- list()
    x$par <- srmsymc_read_par()
    x$yield <- srmsymc_read_yield()
    x$ssb <- srmsymc_read_ssb()
  }
  
  # TODO: clean directory
  if(do_clean) {
    clean_admb("srmsymc")
    clean_admb("srmsymc2")
    system("rm biopar.dat sim.dat simpar.dat simparssb.dat simpary.dat srmsymc_pe.dat simparssbpr.dat simparypr.dat")
  }
  
  if(output) return(x)
}

