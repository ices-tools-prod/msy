
#' plotMSY
#'
#'
#' @param    senfilename    is the name of the senfile, with a corresponding sumfile sharing the same name except replacing ".sen" with ".sum". Produces an interactive window if not specified
#' @param    indexfilename  is the name of an index file in the Lowestoft format pointing to the pf and pm files. If pfpm is specified, this is ignored, if neither specified an interactive window appears to choose index file 
#' @param    pfpm           is a vector of two values; pm and pf (proportion of m and f before spawning
#' @param    srweights      is a vector of three values, giving relative weighting of the three stock recruit forms, or a combination of 0 and NA for automatic weighting.
#' @param    trimming       proportion of the simulations to trim before taking harmonic mean. NA causes no trimming, but a diagnostic plot
#' @param    nits           Number of iterations of bootstrapping - if 0, does only the deterministic fit
#' @param    nhair          number of lines to plot on 'hair' plots
#' @param    varybiodata    If TRUE, bootstraps the biological & fleet data (weight, maturity, mortality, selectivity) if FALSE, varies SR relationship only 
#' @param    stockname      Display title for stock used in titles and output path
#' @param    fpa            Value of Fpa to be plotted on output (NA for no value to plot)
#' @param    flim           Value of Flim to be plotted on output (NA for no value to plot)
#' @param bpa Value of Bpa to be plotted on output (NA for no value to plot)
#' @param blim Value of Blim to be plotted on output (NA for no value to plot)
#' @param    outputfolder   Location for output files. Defaults to ".\\output\\[stockname]\\"
#' @param    datfilename    A pre-calculated dat file - if provided, senfilename, indexfilename, varybiodata, srconstrain and pfpm are ignored in preference to values in the dat file.  Data from the sum file will be added to the plots if it can be found
#' @param    silent         Supresses the majority of the output to screen. Default is TRUE
#' @param    onlyYPR        Calculate only the yield per recruit reference points, for stocks where the SRR is unknown or uncertain. Default is FALSE
#' @return mean(5:7) ##Some function
#' @author Tim Earl \email{timothy.earl@@cefas.co.uk}
#' @export

srmsymc_convertsumsen = function(senfilename = NA, 
                   indexfilename = NA, 
                   pfpm = NA, 
                   srweights=c(NA, NA, NA), 
                   trimming = NA, 
                   nits = 100, 
                   nhair = 100, 
                   varybiodata = TRUE, 
                   stockname = "", 
                   fpa=NA, 
                   flim=NA, 
                   bpa=NA, 
                   blim=NA, 
                   outputfolder="", 
                   datfilename=NA, 
                   silent=TRUE, 
                   onlyYPR=FALSE)
{
  
  if (onlyYPR) {
    srname = c("None")
    srsn = c("00")  
    sr = 0
  } else {
    srname = c("Ricker","Beverton-Holt","Smooth hockeystick")
    srsn = c("Ri","BH","HS")
    sr = 1:length(srname)
  }
  
  srconstrain = TRUE  #Should a penalty be applied to keep alpha and beta positive, and hockeystick breakpoint within data.


  #Validate input arguments and ask if none provided
  if (is.na(datfilename))
  {
    #take input from sen and sum files
    #senfilename = tolower(senfilename)
    #indexfilename = tolower(indexfilename)
    if (is.na(senfilename)) stop("sen file needs to be specified") #senfilename = tolower(choose.files("*.sen", "Choose SEN file",multi=FALSE)) 
    if (!file.exists(senfilename)) stop("SEN file not found")
    sumfilename = sub(".sen",".sum",senfilename,fixed=TRUE)
    if (!file.exists(sumfilename)) stop("SUM file not found")  
  
    if (any(is.na(pfpm)))
    {                
      if (is.na(indexfilename)) indexfilename = tolower(choose.files("*.*", "Choose index file",multi=FALSE))                         
      if (!file.exists(indexfilename)) stop("Index file not found")   
    } else {
      if (length(pfpm)!=2) stop("pfpm must be a vector of length 2")
      if (any(pfpm[1]>1,pfpm[1]<0)) stop("pf must be between 0 and 1")
      if (any(pfpm[2]>1,pfpm[2]<0)) stop("pm must be between 0 and 1")      
    }
  
    if (stockname=="") stockname = gsub(".sen", "", basename(senfilename), fixed=TRUE)
    if (outputfolder=="") outputfolder = paste("output/", stockname, "/", sep="")
  } else {
    #datfilename provided
    datfilename = tolower(datfilename)
    if (!file.exists(datfilename)) stop("file not found:", datfilename)
    if (stockname=="") stockname = scan(datfilename,"",comment.char='#',nlines=2,quiet=TRUE)[1]
    if (outputfolder=="") outputfolder = paste(".\\output\\", stockname, "\\", sep="")
    sumfilename = NA 
    if (!is.na(senfilename)) #This can be used for adding points and legends to the yield and SSB plots
    {
       senfilename = tolower(senfilename)
       sumfilename = sub(".sen",".sum",senfilename,fixed=TRUE)
       if (!file.exists(sumfilename)) 
       { 
         sumfilename = NA   
       } else { #convert to space delimited if comma delimited; save as sumf.tmp
          sumf = scan(sumfilename,"",sep="\n",quote=NULL,quiet=silent)
          sumf = gsub(","," ",sumf)
          cat(sumf,file = "sumf.tmp",sep="\n",quote=NULL)
          sumfilename = "sumf.tmp"
       }
    }
  }
  
  if (length(srweights)!=3) stop("srweights must be a vector of length 3")
  if (any(is.na(srweights)))  { 
    if(!all(range(0,srweights,na.rm=TRUE)==c(0,0))) stop("NAs in srweight can only be combined with zeroes")
    srautoweights <- TRUE
  } else {
    srweights <- srweights/sum(srweights)
    if(is.na(sum(srweights))) stop("srweights can't be normalised - possibly infinite values or all zeroes")
    if(sum(srweights)!=1) stop("srweights can't be normalised - possibly infinite values or all zeroes")    
    srautoweights <- FALSE
  }
          
  #Start of plotMSY function proper  
  cat("Stock:", stockname, "\n")
  graphics.off()   #So that graphics output can be sent to files
  dir.create(outputfolder, showWarnings=FALSE, recursive=TRUE)
  outputfilename = paste(outputfolder, stockname, ".txt", sep="")  
  output = list()
  noredlines = simdatadet = simdata = simy = simSSB  = list()     
  
  #Create .dat, run srmsy and read in its output
  #for (srtype in srname)
  #{
    #Create srmsy.dat and out.dat
    srno = sr[srname[sr]==srtype]
    if (is.na(datfilename))
    {
      sumsen = convertSumSen(senfilename, indexfilename, pfpm, nits, srno, varybiodata, stockname,silent=TRUE,srconstrain=srconstrain)
    } else {
      sumsen = convertDat(datfilename, srno) 
    }
} # eof

