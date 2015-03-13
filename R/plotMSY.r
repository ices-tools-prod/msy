
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
#' @param    optLanding     Boolean indicating whether to optimise by landings (default) or total catch
#' @return    Returns a list containing a summary of the data calculated, and produces a folder of files containing plots and data in outputfolder
#' @author Tim Earl \email{timothy.earl@@cefas.co.uk}
#' @export




plotMSY = function(senfilename = NA, indexfilename = NA, pfpm = NA, srweights=c(NA,NA,NA), trimming = NA, nits = 100, nhair = 100, varybiodata = TRUE, stockname = "",
                   fpa=NA, flim=NA, bpa=NA, blim=NA, outputfolder="", datfilename=NA, silent=TRUE, onlyYPR=FALSE, optLanding=TRUE)
{
  ## Subfunctions
##  source("convertSumSen.r") #for convertSumSen() or convertDat()
  
  boxplot2 <- function(x,add=FALSE,...) 
  { 
  #  browser()
    if (sum(!is.na(x))==0) 
    {
      if (add) box()
      return()
    }
    bp = boxplot(x,plot=FALSE,...,na.rm=TRUE)
    bp$stats = as.matrix(quantile(x,c(0.05, 0.25, 0.50, 0.75, 0.95),na.rm=TRUE))
    bp$out = c(x[x<bp$stats[1]],x[x>bp$stats[5]])
    bp$out <- bp$out[!is.na(bp$out)]
    bp$group = rep(1,length(bp$out))
    bxp(bp,add=add,...)
  }
  
  
  
    ## Stock-Recruit functions
  recruitment = function(SSB, alpha, beta, sr)
  {
    if (sr == 1)
    {
      #Ricker
      recruitment = alpha * SSB*exp(-beta * SSB)
    } else if (sr == 2) {
      #Beverton Holt
      recruitment = alpha * SSB /(beta + SSB)
    } else {
      #Smooth Hockeystick
      g = 0.001
      recruitment = alpha*(SSB+sqrt(beta^2+g)-sqrt((SSB-beta)^2+g))
 
    }
    return(recruitment)
  }    
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
    senfilename = tolower(senfilename)
    indexfilename = tolower(indexfilename)
    if (is.na(senfilename)) senfilename = tolower(choose.files("*.sen", "Choose SEN file",multi=FALSE)) 
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
    if (outputfolder=="") outputfolder = paste(".\\output\\", stockname, "\\", sep="")
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
  for (srtype in srname)
  {
    #Create srmsy.dat and out.dat
    srno = sr[srname[sr]==srtype]
    if (is.na(datfilename))
    {
      sumsen = convertSumSen(senfilename, indexfilename, pfpm, nits, srno, varybiodata, stockname,silent=TRUE,srconstrain=srconstrain,optLanding=optLanding)
    } else {
      sumsen = convertDat(datfilename, srno) 
    }
  
    #Run ADMB code in mcmc mode
    cat("Running AD Model Builder code for stock recruitment type ", srtype, '\n')
    command = paste("srmsymc.exe -mcmc", (nits+11)*10000, "-mcsave 10000 -nosdmcmc") ##nits runs +10 for burn in period, +1 for deterministic; 
                                                                                     ##10,000 fairly arbitrary to reduce correlation between successive saved runs
    system(command, show.output.on.console=!silent)
    system("srmsymc.exe -mceval",show.output.on.console=!silent)
    
    #Read in simpar.dat, simpary.dat, simparSSB.dat  and remove unfeasible simulations
    simdata[[srtype]] = read.table(".\\simpar.dat")
    simy[[srtype]] = as.matrix(read.table(".\\simpary.dat"))
    simSSB[[srtype]] = as.matrix(read.table(".\\simparssb.dat")) 
    colnames(simdata[[srtype]]) = c("N","ap","bp","alpha","beta","sigr","scor","fcrash","fmax","f01","f20","f25","f30","f35","f40","fmsy","msy","bmsy","msypr","bmsypr","msyr",
                                    "fnfcrash","fnfmax","fnf01","fnf20","fnf15","fnf30","fnf35","fnf40","dydf","penalty","nll","AIC")
    simdatadet[[srtype]] = subset(simdata[[srtype]], subset=(simdata[[srtype]]$N==0))
    colnames(simdatadet[[srtype]]) = colnames(simdata[[srtype]]) 
    simdata[[srtype]] = subset(simdata[[srtype]], subset=(simdata[[srtype]]$N!=0))
    if (srtype==srname[3]) {
      simdata[[srtype]][simdata[[srtype]]$beta<min(sumsen$SSB) | simdata[[srtype]]$alpha<=0 | simdata[[srtype]]$beta<=0,
                        c("ap","bp","alpha","beta","sigr","scor","fcrash","fmsy","msy","bmsy","fnfcrash","dydf","penalty","nll","AIC")] <- NA
    } else {
      simdata[[srtype]][simdata[[srtype]]$alpha<=0 | simdata[[srtype]]$beta<=0,
                      c("ap","bp","alpha","beta","sigr","scor","fcrash","fmsy","msy","bmsy","fnfcrash","dydf","penalty","nll","AIC")] <- NA
    
    }
               
        
    if (any(!simdata[[srtype]]$fmsy<3,na.rm=TRUE)) simdata[[srtype]][which(!simdata[[srtype]]$fmsy<3),c("fmsy","bmsy","msy")] <- NA
    if (any(!simdata[[srtype]]$fcrash<5,na.rm=TRUE)) simdata[[srtype]][which(!simdata[[srtype]]$fcrash<5),"fcrash"] <- NA 
    if (any(!simdata[[srtype]]$fmax<3,na.rm=TRUE)) simdata[[srtype]][which(!simdata[[srtype]]$fmax<3),c("fmax","bmsypr","msypr")] <- NA     
    indexmap = c("f20","f25","f30","f35","f40","f01") # columns of outputpr[[srtype]] to be output to output[[srtype]]   
    for (i in indexmap) if (any(!simdata[[srtype]][,i]<3,na.rm=TRUE)) simdata[[srtype]][which(!simdata[[srtype]][,i]<3),i] <- NA
    
    simy[[srtype]] = simy[[srtype]][c(1,2,simdata[[srtype]]$N-8),] #N starts at 11, this should correspond to line 3 of simy and simSSB
    simSSB[[srtype]] = simSSB[[srtype]][c(1,2,simdata[[srtype]]$N-8),]    
    simSSB[[srtype]][-1,] = pmax(simSSB[[srtype]][-1,],0) #Set to zero if SSB or y is negative
    simy[[srtype]][-1,] = pmax(simy[[srtype]][-1,],0)
    noredlines[[srtype]] = (prod(dim(simdata[[srtype]]))==0)
    

    #Set up array for ouput
    output[[srtype]] = matrix(NA,nrow =9, ncol=9,dimnames=list(c("Deterministic","Mean","5%ile","25%ile","50%ile","75%ile","95%ile","CV","N"),
                                                               c("Fcrash","Fmsy","Bmsy","MSY","ADMB Alpha","ADMB Beta","Unscaled Alpha","Unscaled Beta","AICc"))) 
    indexmap = c("fcrash","fmsy","bmsy","msy","ap","bp","alpha","beta","AIC") # columns of simdata to be output to output[[srtype]]
  (length(simdata[[srtype]][,1])>0)

    for (i in seq(along=indexmap))
    {
      output[[srtype]]["Deterministic",i] = simdatadet[[srtype]][,indexmap[i]]
      if (length(simdata[[srtype]][,1])>0)
      {
        output[[srtype]]["Mean",i] = mean(simdata[[srtype]][,indexmap[i]],na.rm=TRUE)      
        output[[srtype]][3:7,i] = quantile(simdata[[srtype]][,indexmap[i]],c(0.05, 0.25, 0.50, 0.75, 0.95),na.rm=TRUE)   
        output[[srtype]]["CV",i] = sd(simdata[[srtype]][,indexmap[i]],na.rm=TRUE)/mean(simdata[[srtype]][,indexmap[i]],na.rm=TRUE)  
        output[[srtype]]["N",i] = sum(!is.na(simdata[[srtype]][,indexmap[i]]))
        

      }
  #    output[[srtype]]["N",c("Fmsy","Bmsy","MSY")] = sum(simdata[[srtype]]$fmsy<3)                                                             
  #    output[[srtype]]["N",c("Fcrash")] = sum(simdata[[srtype]]$fcrash<5)   
  #    output[[srtype]]["N",c("ADMB Alpha","ADMB Beta","Unscaled Alpha","Unscaled Beta","AICc")] = length(simdata[[srtype]]$N)   
    }                                                              
  }
  
  
  simSSBpr = as.matrix(read.table(".\\simparssbpr.dat")) 
  simypr = as.matrix(read.table(".\\simparypr.dat"))
  
  outputpr = matrix(NA,nrow =9, ncol=11,dimnames=list(c("Deterministic","Mean","5%ile","25%ile","50%ile","75%ile","95%ile","CV","N"),
                                                     c("F20","F25","F30","F35","F40","F01","Fmax","Bmsypr","MSYpr","Fpa","Flim")))  
# browser()                                                   
  indexmap = c("f20","f25","f30","f35","f40","f01","fmax","bmsypr","msypr") # columns of outputpr[[srtype]] to be output to output[[srtype]]
  for (i in seq(along=indexmap))
  {
    outputpr["Deterministic",i] = simdatadet[[srtype]][,indexmap[i]]
    if (!noredlines[[srtype]])
    {
      outputpr["Mean",i] = mean(simdata[[srtype]][,indexmap[i]],na.rm=TRUE)      
      outputpr[3:7,i] = quantile(simdata[[srtype]][,indexmap[i]],c(0.05, 0.25, 0.50, 0.75, 0.95),na.rm=TRUE)   
      outputpr["CV",i] = sd(simdata[[srtype]][,indexmap[i]],na.rm=TRUE)/mean(simdata[[srtype]][,indexmap[i]],na.rm=TRUE)   
      outputpr["N",i] =  sum(!is.na(simdata[[srtype]][,indexmap[i]]))       
    }
  }      
  outputpr["Deterministic","Fpa"] = fpa
  outputpr["Deterministic","Flim"] = flim 
  
  #Adjust for recruitage so that ssb and recuitment line up
  if (sumsen$recruitage > 0)
  {
    sumsen$Recruits = sumsen$Recruits[-c(1:sumsen$recruitage)]
    sumsen$SSB = sumsen$SSB[seq(along=sumsen$Recruits)]
    sumsen$Year = sumsen$Year[-c(1:sumsen$recruitage)]
  }  
  
  done_fit2 <- FALSE
#  source(".\\noel\\help_example\\fit2.r")
#  b.cm.pred <- fit2(sumsen)
#  done_fit2 <- TRUE  
  
   if (!is.na(sumfilename))
   {
     #Read columns from sumfile
     #sumfilename = "sumf.tmp"
     yearRange =  scan("sumf.tmp",nlines=1,skip=4,quote="\"",quiet=TRUE)
     ncols =  scan("sumf.tmp",nlines=1,skip=1,quote="\"",quiet=TRUE)
       
     if (ncols!=12) stop("Table in sum file should have 12 columns")
     sumData =  t(matrix(scan("sumf.tmp",nlines=yearRange[2]-yearRange[1]+1,skip=27,quote="\"",quiet=TRUE),ncols))
     SSBtext =  gsub("\t"," ",scan("sumf.tmp","",nlines=1,skip=7,quote="\"",sep=",",quiet=TRUE))
     rectext =  gsub("\t"," ",scan("sumf.tmp","",nlines=1,skip=5,quote="\"",sep=",",quiet=TRUE))
     yieldtext =  gsub("\t"," ",scan("sumf.tmp","",nlines=1,skip=13,quote="\"",sep=",",quiet=TRUE))
     if (max(sumData[,10]==0))  #If there's no human consumption, use total
     {  
       sumData[,10] = sumData[,9]
       yieldtext =  gsub("\t"," ",scan("sumf.tmp","",nlines=1,skip=11,quote="\"",sep=",",quiet=TRUE))
     }
     if (max(sumData[,6]==0))  sumData[,6] = sumData[,5]
     fbar = rev(sumData[,10])[1] #Fbar in last year
   } else {
     SSBtext = "SSB"
     rectext = "Recruitment"
     yieldtext = "yield"
     sumData = c()
     fbar = NA
   }
   fbartext = paste("F",rev(sumsen$Year)[1],sep="")
   
   
    xlim = c(0,2)
    if (!onlyYPR)
    {
  
##Plot SRRs
    png(paste(outputfolder, stockname, "_SRR.png",sep=""),height=11.5,width=9,units="in",res=144)
      layout(t(matrix(1:6,2)))
      SSB = max(sumsen$SSB)*(0:105)/100  #x values for plotting

      ##SRRs on first page
      for (srtype in sr)
      {
       
        plot(c(0,sumsen$SSB), c(0,sumsen$Recruits), xlab=SSBtext, ylab=rectext,type='n')
        #if (done_fit2) lines(b.cm.pred$stock.size,b.cm.pred$recruit,col="purple")
        title(paste(sumsen$stock, srname[srtype]))
        points(sumsen$SSB, sumsen$Recruits)  #Data points
        recruits = recruitment(SSB, simdatadet[[srtype]]$alpha[1], simdatadet[[srtype]]$beta[1], srtype)
        if (!noredlines[[srtype]])
        {
          m = matrix(NA,(length(simdata[[srtype]][[1]])),length(SSB))
          for (i in 1:length(simdata[[srtype]][[1]]))
          {
            m[i,] = recruitment(SSB, simdata[[srtype]]$alpha[i], simdata[[srtype]]$beta[i], srtype)
          }
          lines(SSB,apply(m,2,quantile,probs=0.05,na.rm=TRUE),col='red',lty=3)      #Confidence intervals
          lines(SSB,apply(m,2,quantile,probs=0.1,na.rm=TRUE),col='red',lty=2)
          lines(SSB,apply(m,2,quantile,probs=0.5,na.rm=TRUE),col='red')
          lines(SSB,apply(m,2,quantile,probs=0.9,na.rm=TRUE),col='red',lty=2)
          lines(SSB,apply(m,2,quantile,probs=0.95,na.rm=TRUE),col='red',lty=3)
        }
        lines(SSB,recruits, col="blue")   #Deterministic line   
        legend(0,0.7*max(sumsen$Recruits), yjust=0,
               c("5%, 95%", "10%, 90%", "50%", "Deterministic",paste(sum(!is.na(simdata[[srtype]][["alpha"]])),"/",nits,sep="")),
              lty=c(3,2,1,1,1), col=c(rep('red',3), 'blue','white'))
       #----SRR hairplot
          plot(c(0,sumsen$SSB), c(0,sumsen$Recruits), xlab=SSBtext, ylab=rectext,type='n')
          title(paste(sumsen$stock, srname[srtype]))
          recruits = recruitment(SSB, simdatadet[[srtype]]$alpha[1], simdatadet[[srtype]]$beta[1], srtype) 
          srnhair = min(nhair,length(simdata[[srtype]]$fcrash))
          if (srnhair > 1)  for (i in 2+(1:srnhair))  points(SSB,recruitment(SSB, simdata[[srtype]]$alpha[i], simdata[[srtype]]$beta[i], srtype),type='l', col='red')  #Simulations
          lines(SSB,recruits, col="blue")      
      }
    dev.off() #SRR plot
    
     
    #plots of parameter correlation
    png(paste(outputfolder, stockname, "_diagnostics.png",sep=""),height=11.5,width=9,units="in",res=144)     
      layout(t(matrix(c(1:6),2)))
      for (srtype in sr)
      {

        if (srtype==3) 
        {
          plot.new()
        } else {
          plot(c(simdata[[srtype]]$ap,simdatadet[[srtype]]$ap),c(simdata[[srtype]]$bp,simdatadet[[srtype]]$bp),xlab = "Alpha", ylab="Beta")
          points(simdatadet[[srtype]]$ap,simdatadet[[srtype]]$bp,col='blue',pch=19)
          title(paste(stockname,"ADMB parameters",srname[srtype]))         
        }
        plot(c(simdata[[srtype]]$alpha,simdatadet[[srtype]]$alpha),c(simdata[[srtype]]$beta,simdatadet[[srtype]]$beta),xlab = "Alpha", ylab="Beta")
        points(simdatadet[[srtype]]$alpha,simdatadet[[srtype]]$beta,col='blue',pch=19)
        title(paste(stockname,"Unscaled parameters",srname[srtype]))
       }
     dev.off() #Diagnostics plot
     
    
    
    #one plot per SRR

    for (srtype in sr)
    {
      png(paste(outputfolder, stockname, "_Yield_",srsn[srtype],".png",sep=""),height=11.5,width=9,units="in",res=144)
        layout(t(matrix(c(1,1:7),2)), widths=c(5,5), heights=c(1,3,3.5,4))
        par(mai=c(0,0,0.8,0))
        plot.new()
        title(paste(stockname,srname[srtype]))
        srnhair = min(nhair,length(simdata[[srtype]]$fcrash))
        par(mai=c(0.5,1,0,0.5))
      
        ##Boxplots
        for (i in 1:2)
        {
          plot(0,0,axes=FALSE,type='n',xlim=xlim,ylim=c(-2.5,6.5),xlab=" ",ylab="")
          axis(2,(-2):6,c(fbartext,"Fpa","Flim","Fcra","Fmax","F01","Fmsy","F35","F40"),cex.axis=1.3,las=1)
          axis(1)
          box()
          if(!is.na(fbar)) boxplot2(fbar, add=TRUE, at=-2,horizontal=TRUE,axes=FALSE)
          if(!is.na(fpa)) boxplot2(fpa, add=TRUE, at=-1,horizontal=TRUE,axes=FALSE)
          if(!is.na(flim)) boxplot2(flim, add=TRUE, at=0,horizontal=TRUE,axes=FALSE)
          if (!noredlines[[srtype]])
          {
            boxplot2(simdata[[srtype]]$fcrash, add=TRUE, at=1,horizontal=TRUE,axes=FALSE)
            boxplot2(simdata[[srtype]]$fmax,   add=TRUE, at=2,horizontal=TRUE,axes=FALSE)
            boxplot2(simdata[[srtype]]$f01,    add=TRUE, at=3,horizontal=TRUE,axes=FALSE)
            boxplot2(simdata[[srtype]]$fmsy,   add=TRUE, at=4,horizontal=TRUE,axes=FALSE)
            boxplot2(simdata[[srtype]]$f35,    add=TRUE, at=5,horizontal=TRUE,axes=FALSE)
            boxplot2(simdata[[srtype]]$f40,    add=TRUE, at=6,horizontal=TRUE,axes=FALSE)
          }
        }
    
        ##Yield
        plot(t(simy[[srtype]][1,]),1.5*t(simy[[srtype]][2,]),type='n',axes=TRUE, xlab="", ylab = yieldtext, xlim=xlim,ylim = c(0,1.5*max(simy[[srtype]][2,])))
        if (!noredlines[[srtype]])
        {
          #quantiles
          lines(simy[[srtype]][1,],apply(simy[[srtype]][-1,],2,quantile,probs=0.05,na.rm=TRUE),col='red',lty=3)
          lines(simy[[srtype]][1,],apply(simy[[srtype]][-1,],2,quantile,probs=0.1,na.rm=TRUE),col='red',lty=2)
          lines(simy[[srtype]][1,],apply(simy[[srtype]][-1,],2,quantile,probs=0.5,na.rm=TRUE),col='red')
          lines(simy[[srtype]][1,],apply(simy[[srtype]][-1,],2,quantile,probs=0.9,na.rm=TRUE),col='red',lty=2)
          lines(simy[[srtype]][1,],apply(simy[[srtype]][-1,],2,quantile,probs=0.95,na.rm=TRUE),col='red',lty=3)
        }
        points(t(simy[[srtype]][1,]),t(simy[[srtype]][2,]),type='l', col='blue')
        if (!is.na(sumfilename)) points(sumData[,10],sumData[,6],cex=.7,type='b')
        if (!is.na(sumfilename)) text(rev(sumData[,10])[1], rev(sumData[,6])[1],rev(sumData[,1])[1],pos=4 )
        if (!is.na(sumfilename)) points(rev(sumData[,10])[1], rev(sumData[,6])[1],cex=.7,pch=19 )

        plot(t(simy[[srtype]][1,]),1.5*t(simy[[srtype]][2,]),type='n',axes=TRUE, xlab="", ylab = yieldtext, xlim=xlim,ylim = c(0,1.5*max(simy[[srtype]][2,])))
        if (srnhair > 1)  for (i in 2+(1:srnhair))  points(t(simy[[srtype]][1,]),t(simy[[srtype]][i,]),type='l', col='red')
        points(t(simy[[srtype]][1,]),t(simy[[srtype]][2,]),type='l', col='blue')
        
        ##SSB
        par(mai=c(1,1,0,0.5))
        if (!is.na(sumfilename)) plot(t(simSSB[[srtype]][1,]),1.5*t(simSSB[[srtype]][2,]),type='n',axes=TRUE, ylab = SSBtext,xlab="F", xlim=xlim,ylim = c(0,1.5*max(sumData[,3])))  else
                                 plot(t(simSSB[[srtype]][1,]),1.5*t(simSSB[[srtype]][2,]),type='n',axes=TRUE, ylab = SSBtext,xlab="F", xlim=xlim)  
        if (!noredlines[[srtype]])
        {
          #quantiles
          lines(simSSB[[srtype]][1,],apply(simSSB[[srtype]][-1,],2,quantile,probs=0.05),col='red',lty=3)
          lines(simSSB[[srtype]][1,],apply(simSSB[[srtype]][-1,],2,quantile,probs=0.1),col='red',lty=2)
          lines(simSSB[[srtype]][1,],apply(simSSB[[srtype]][-1,],2,quantile,probs=0.5),col='red')
          lines(simSSB[[srtype]][1,],apply(simSSB[[srtype]][-1,],2,quantile,probs=0.9),col='red',lty=2)
          lines(simSSB[[srtype]][1,],apply(simSSB[[srtype]][-1,],2,quantile,probs=0.95),col='red',lty=3)
        }
        points(t(simSSB[[srtype]][1,]),t(simSSB[[srtype]][2,]),type='l',col='blue')
        if (!is.na(sumfilename)) points(sumData[,10],sumData[,3],cex=.7,type='b')   
        if (!is.na(sumfilename)) text(rev(sumData[,10])[1], rev(sumData[,3])[1],rev(sumData[,1])[1],pos=4 )
        if (!is.na(sumfilename)) points(rev(sumData[,10])[1], rev(sumData[,3])[1],cex=.7,type='b',pch=19 )

        if (!is.na(sumfilename)) plot(t(simSSB[[srtype]][1,]),1.5*t(simSSB[[srtype]][2,]),type='n',axes=TRUE, ylab = SSBtext,xlab="F", xlim=xlim,ylim = c(0,1.5*max(sumData[,3])))  else
                                 plot(t(simSSB[[srtype]][1,]),1.5*t(simSSB[[srtype]][2,]),type='n',axes=TRUE, ylab = SSBtext,xlab="F", xlim=xlim)  
        if (srnhair > 1)  for (i in 2+(1:srnhair))  points(t(simSSB[[srtype]][1,]),t(simSSB[[srtype]][i,]),type='l', col='red')
        lines(t(simSSB[[srtype]][1,]),t(simSSB[[srtype]][2,]),col='blue')
      dev.off()
     }
     }
   
      png(paste(outputfolder, stockname, "_YPR.png",sep=""),height=11.5,width=9,units="in",res=144)     
        ##yield per recruit
        srtype = 1 ##Any of them should be the same for fmax, f01, f35, f40
        layout(t(matrix(c(1,1:7),2)), widths=c(5,5), heights=c(1,3,3.5,4))
        par(mai=c(0,0,0.8,0))
        plot.new()
        title(paste(stockname,"- Per recruit statistics"))
  
        ##Boxplots
        srnhair = min(nhair,dim(simypr)[1]-2)
        par(mai=c(0.5,1,0,0.5))
        for (i in 1:2)
        {
          plot(0,0,axes=FALSE,type='n',xlim=xlim,ylim=c(-2.5,4.5),xlab="",ylab="")
          axis(2,(-2):4,c(fbartext,"Fpa","Flim","Fmax","F01","F35","F40"),,cex.axis=1.3,las=1)
          axis(1)
          box()
          if(!is.na(fbar)) boxplot2(fbar, add=TRUE, at=-2,horizontal=TRUE,axes=FALSE)
          if(!is.na(fpa)) boxplot2(fpa, add=TRUE, at=-1,horizontal=TRUE,axes=FALSE)
          if(!is.na(flim)) boxplot2(flim, add=TRUE, at=0,horizontal=TRUE,axes=FALSE)
          if (!noredlines[[1]])
          {
            boxplot2(simdata[[1]]$fmax,   add=TRUE, at=1,horizontal=TRUE,axes=FALSE)
            boxplot2(simdata[[1]]$f01,    add=TRUE, at=2,horizontal=TRUE,axes=FALSE)
            boxplot2(simdata[[1]]$f35,    add=TRUE, at=3,horizontal=TRUE,axes=FALSE)
            boxplot2(simdata[[1]]$f40,    add=TRUE, at=4,horizontal=TRUE,axes=FALSE)
          }
        }
           
        ##YPR
        plot(t(simypr[1,]),1.5*t(simypr[2,]),type='n',xlab="",ylab = "Yield per recruit",xlim=xlim)
        if (prod(dim(simypr)) >0)
        {
          #quantiles
          lines(simypr[1,],apply(simypr[-1,],2,quantile,probs=0.05,na.rm=TRUE),col='red',lty=3)
          lines(simypr[1,],apply(simypr[-1,],2,quantile,probs=0.1,na.rm=TRUE),col='red',lty=2)
          lines(simypr[1,],apply(simypr[-1,],2,quantile,probs=0.5,na.rm=TRUE),col='red')
          lines(simypr[1,],apply(simypr[-1,],2,quantile,probs=0.9,na.rm=TRUE),col='red',lty=2)
          lines(simypr[1,],apply(simypr[-1,],2,quantile,probs=0.95,na.rm=TRUE),col='red',lty=3)
        }
        points(t(simypr[1,]),t(simypr[2,]),type='l', col='blue')   
       
        plot(t(simypr[1,]),1.5*t(simypr[2,]),type='n', xlab="",ylab = "Yield per recruit",axes=TRUE, xlim=xlim)
        if (srnhair > 1)  for (i in 2+(1:srnhair))  points(t(simypr[1,]),t(simypr[i,]),type='l', col='red')
        points(t(simypr[1,]),t(simypr[2,]),type='l', col='blue')
  
        ##SSBPR
        par(mai=c(1,1,0,0.5))
        plot(t(simSSBpr[1,]),1.5*t(simSSBpr[2,]),type='n',axes=TRUE, ylab = "SSB per recruit", xlab="F",xlim=xlim)
        if (prod(dim(simSSBpr))>0)
        {
          #quantiles
          lines(simSSBpr[1,],apply(simSSBpr[-1,],2,quantile,probs=0.05,na.rm=TRUE),col='red',lty=3)
          lines(simSSBpr[1,],apply(simSSBpr[-1,],2,quantile,probs=0.1,na.rm=TRUE),col='red',lty=2)
          lines(simSSBpr[1,],apply(simSSBpr[-1,],2,quantile,probs=0.5,na.rm=TRUE),col='red')
          lines(simSSBpr[1,],apply(simSSBpr[-1,],2,quantile,probs=0.9,na.rm=TRUE),col='red',lty=2)
          lines(simSSBpr[1,],apply(simSSBpr[-1,],2,quantile,probs=0.95,na.rm=TRUE),col='red',lty=3)
        }
        points(t(simSSBpr[1,]),t(simSSBpr[2,]),type='l', col='blue')
         
        plot(t(simSSBpr[1,]),1.5*t(simSSBpr[2,]),type='n', xlab="F",ylab = "SSB per recruit",axes=TRUE, xlim=xlim)
        if (srnhair > 1)  for (i in 2+(1:srnhair))  points(t(simSSBpr[1,]),t(simSSBpr[i,]),type='l', col='red')
        points(t(simSSBpr[1,]),t(simSSBpr[2,]),type='l', col='blue')     
      graphics.off()


           ## Combined distribution
      if(!onlyYPR)           
      {
    #    browser()
        png(paste(outputfolder, stockname, "_Fmsy.png",sep=""),height=11.5,width=9,units="in",res=144) 
        allSimData = rbind(simdata[[1]],simdata[[2]],simdata[[3]])
        f.breaks <- hist(c(allSimData$fmsy,allSimData$fcrash),50, plot=FALSE)$breaks
        allSimDataHist = hist(allSimData$fmsy,f.breaks,plot=FALSE)
        count1 = hist(simdata[[1]]$fmsy,f.breaks,plot=FALSE)$counts
        count2 = hist(simdata[[2]]$fmsy,f.breaks,plot=FALSE)$counts
        count3 = hist(simdata[[3]]$fmsy,f.breaks,plot=FALSE)$counts
        mf = max(allSimDataHist$counts)
        par(mfrow=c(2,1))
        allSimData = rbind(simdata[[1]],simdata[[2]],simdata[[3]])
        tempHist = hist(c(allSimData$fmsy,allSimData$fcrash),50, plot=FALSE)
        my = max(tempHist$counts)
        mx = max(f.breaks)


        hist(rbind(simdata[[1]],simdata[[2]],simdata[[3]])$fmsy, f.breaks, col=grey(0.75), ylim = c(0,1.10*max(tempHist$counts)), xlab="Fmsy", 
        main = paste(stockname,"- combined Fmsy distribution: equally weighted"))
        hist(rbind(simdata[[1]],simdata[[2]])$fmsy, tempHist$breaks, add =TRUE, col=grey(0.5))
        hist(simdata[[1]]$fmsy, tempHist$breaks, add =TRUE, col=grey(0.25))
        qs = cbind(quantile(allSimData$fmsy,c(0.05, 0.25, 0.50, 0.75, 0.95),na.rm=TRUE), 
                   quantile(allSimData$fcrash,c(0.05, 0.25, 0.50, 0.75, 0.95),na.rm=TRUE),
                   quantile(allSimData$msy,c(0.05, 0.25, 0.50, 0.75, 0.95),na.rm=TRUE),
                   quantile(allSimData$bmsy,c(0.05, 0.25, 0.50, 0.75, 0.95),na.rm=TRUE) )
        lines(c(qs[1,1],qs[1,1]),c(0,my),col='red',lty=3)      #Confidence intervals
        lines(c(qs[2,1],qs[2,1]),c(0,my*1.05),col='red',lty=2)
        lines(c(qs[3,1],qs[3,1]),c(0,my),col='red')
        lines(c(qs[4,1],qs[4,1]),c(0,my*1.05),col='red',lty=2)
        lines(c(qs[5,1],qs[5,1]),c(0,my),col='red',lty=3)
        text(qs[1,1],my,"5%",pos=3,col="red")
        text(qs[2,1],my*1.05,"25%",pos=3,col="red")
        text(qs[3,1],my,"50%",pos=3,col="red")
        text(qs[4,1],my*1.05,"75%",pos=3,col="red")
        text(qs[5,1],my,"95%",pos=3,col="red")
        
        legend(0.75*mx,my,c("Ricker","Beverton-Holt","Hockeystick"), fill=grey(0.25*1:3))

        hist(rbind(simdata[[1]],simdata[[2]],simdata[[3]])$fcrash, f.breaks, col=grey(0.75), ylim = c(0,1.10*max(tempHist$counts)), xlab="Fcrash", 
        main = paste(stockname,"- combined Fcrash distribution: equally weighted"))
        hist(rbind(simdata[[1]],simdata[[2]])$fcrash, f.breaks, add =TRUE, col=grey(0.5))
        hist(simdata[[1]]$fcrash, f.breaks, add =TRUE, col=grey(0.25))

#        qs = quantile(allSimData$fcrash,c(0.05, 0.25, 0.50, 0.75, 0.95))
        lines(c(qs[1,2],qs[1,2]),c(0,my),col='red',lty=3)      #Confidence intervals
        lines(c(qs[2,2],qs[2,2]),c(0,my*1.05),col='red',lty=2)
        lines(c(qs[3,2],qs[3,2]),c(0,my),col='red')
        lines(c(qs[4,2],qs[4,2]),c(0,my*1.05),col='red',lty=2)
        lines(c(qs[5,2],qs[5,2]),c(0,my),col='red',lty=3)
        text(qs[1,2],my,"5%",pos=3,col="red")
        text(qs[2,2],my*1.05,"25%",pos=3,col="red")
        text(qs[3,2],my,"50%",pos=3,col="red")
        text(qs[4,2],my*1.05,"75%",pos=3,col="red")
        text(qs[5,2],my,"95%",pos=3,col="red")
     graphics.off()
     
      if (srautoweights)
      {  
        lik <- sapply(simdata, function(x){sort(exp(-x$nll),na.last=FALSE)})   
        if (is.na(trimming)) 
        {
          trim <- 0:(nits/4)
          srwall <- sapply(trim, function(p){1/colMeans(1/lik[(p+1):nits,],na.rm=TRUE)} )*ifelse(is.na(srweights),1,0)
          srwall[is.na(srwall)] <- 0
          srwall <- t(t(srwall)/colSums(srwall)    )
          png(paste0(outputfolder, stockname, "_trim_diag.png"),800,600)
            plot(100*trim/nits,srwall["Ricker",],ylim=range(0,srwall), lty=1, type='l', 
                  ylab="Relative weight", xlab="Trimmed percentage", main = paste(stockname,"Trimming the harmonic mean",sep=" - "))
            lines(100*trim/nits,srwall["Beverton-Holt",],lty=2)
            lines(100*trim/nits,srwall["Smooth hockeystick",],lty=3)          
            legend("topright",pch=NA,lty=1:3,legend=srname)
          dev.off()
          write.csv(t(rbind(Percentage=100*trim/nits,srwall)),paste0(outputfolder,stockname,"_trim_diag.csv"),row.names=FALSE)
          trimming <- 0
          
        }
#        browser()
         srweightsexact <- sapply(round(trimming*nits), function(p){1/colMeans(1/lik[(p+1):nits,],na.rm=TRUE)} )*ifelse(is.na(srweights),1,0)
#        srweightsexact <- sapply(simdata,function(x){1/mean(1/exp(-x$nll),na.rm=TRUE)})*ifelse(is.na(srweights),1,0)
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
        
      }
     
      # browser()
      png(paste0(outputfolder, stockname, "_Fmsy2.png"),height=11.5,width=9,units="in",res=144) 
        #par(mfrow=c(2,1))
        allSimData = rbind(simdata[[1]][seq(length=srweights[1]),],
                           simdata[[2]][seq(length=srweights[2]),],
                           simdata[[3]][seq(length=srweights[3]),])
                           
        allSimData <- allSimData[!is.na(rowSums(allSimData[,!(colnames(allSimData) %in% c("fmax","msypr","bmsypr"))])),]
        allSimSSB = rbind(simSSB[[1]][seq(length=srweights[1]+1),],
                          simSSB[[2]][seq(length=srweights[2]),],
                          simSSB[[3]][seq(length=srweights[3]),])
        allSimSSB <- allSimSSB[!is.na(rowSums(allSimSSB)),]        
        allSimDataHist = hist(allSimData$fmsy,f.breaks,plot=FALSE)
        qs = cbind(qs,quantile(allSimData$fmsy,c(0.05, 0.25, 0.50, 0.75, 0.95)), 
                   quantile(allSimData$fcrash,c(0.05, 0.25, 0.50, 0.75, 0.95),na.rm=TRUE),
                   quantile(allSimData$msy,c(0.05, 0.25, 0.50, 0.75, 0.95),na.rm=TRUE),
                   quantile(allSimData$bmsy,c(0.05, 0.25, 0.50, 0.75, 0.95),na.rm=TRUE) )

        tempHist = hist(c(allSimData$fmsy,allSimData$fcrash),f.breaks, plot=FALSE)
        my = max(tempHist$counts)
        mx = max(tempHist$breaks)
        par(mfrow=c(2,1))
        hist(rbind(simdata[[1]][seq(length=srweights[1]),],simdata[[2]][seq(length=srweights[2]),],simdata[[3]][seq(length=srweights[3]),])$fmsy, f.breaks, col=grey(0.75),
             ylim = c(0,1.10*max(tempHist$counts)), xlab="Fmsy", 
             main = paste(stockname,"- combined Fmsy distribution:",ifelse(srautoweights,"automatically weighted","manually weighted")))
        hist(rbind(simdata[[1]][seq(length=srweights[1]),],simdata[[2]][seq(length=srweights[2]),])$fmsy, f.breaks, add =TRUE, col=grey(0.5))
        hist(simdata[[1]][seq(length=srweights[1]),]$fmsy, f.breaks, add =TRUE, col=grey(0.25))
         

        lines(c(qs[1,5],qs[1,5]),c(0,my),col='red',lty=3)      #Confidence intervals
        lines(c(qs[2,5],qs[2,5]),c(0,my*1.05),col='red',lty=2)
        lines(c(qs[3,5],qs[3,5]),c(0,my),col='red')
        lines(c(qs[4,5],qs[4,5]),c(0,my*1.05),col='red',lty=2)
        lines(c(qs[5,5],qs[5,5]),c(0,my),col='red',lty=3)
        text(qs[1,5],my,"5%",pos=3,col="red")
        text(qs[2,5],my*1.05,"25%",pos=3,col="red")
        text(qs[3,5],my,"50%",pos=3,col="red")
        text(qs[4,5],my*1.05,"75%",pos=3,col="red")
        text(qs[5,5],my,"95%",pos=3,col="red")
        
        
        legtext <- paste(c("Ricker","Beverton-Holt","Hockeystick")," (",round(100*srweightsexact),"%)",sep="")
        legend(0.75*mx,my,legtext, fill=grey(0.25*1:3))

        hist(rbind(simdata[[1]][seq(length=srweights[1]),],simdata[[2]][seq(length=srweights[2]),],simdata[[3]][seq(length=srweights[3]),])$fcrash, f.breaks, col=grey(0.75), ylim = c(0,1.10*max(tempHist$counts)), xlab="Fcrash",
                     main = paste(stockname,"- combined Fcrash distribution:",ifelse(srautoweights,"automatically weighted","manually weighted")))
        hist(rbind(simdata[[1]][seq(length=srweights[1]),],simdata[[2]][seq(length=srweights[2]),])$fcrash, f.breaks, add =TRUE, col=grey(0.5))
        hist(simdata[[1]][seq(length=srweights[1]),]$fcrash, f.breaks, add =TRUE, col=grey(0.25))


#        qs = quantile(allSimData$fcrash,c(0.05, 0.25, 0.50, 0.75, 0.95))
        lines(c(qs[1,6],qs[1,6]),c(0,my),col='red',lty=3)      #Confidence intervals
        lines(c(qs[2,6],qs[2,6]),c(0,my*1.05),col='red',lty=2)
        lines(c(qs[3,6],qs[3,6]),c(0,my),col='red')
        lines(c(qs[4,6],qs[4,6]),c(0,my*1.05),col='red',lty=2)
        lines(c(qs[5,6],qs[5,6]),c(0,my),col='red',lty=3)
        text(qs[1,6],my,"5%",pos=3,col="red")
        text(qs[2,6],my*1.05,"25%",pos=3,col="red")
        text(qs[3,6],my,"50%",pos=3,col="red")
        text(qs[4,6],my*1.05,"75%",pos=3,col="red")
        text(qs[5,6],my,"95%",pos=3,col="red")
     graphics.off()    
     
    #  browser()
      if(!is.na(blim))
      {
        png(paste(outputfolder, stockname, "_Fmsy3.png",sep=""),height=11.5,width=9,units="in",res=144) 
          plot(NA,xlim=range(allSimSSB[1,]),ylim=c(0,1),xlab="F",ylab="Probability SSB < Blim",
               main = paste(stockname,"- Probability SSB < Blim"))
          probs <- apply(allSimSSB[-1,],2,function(x){sum(x<blim)/length(x)})
          lines(allSimSSB[1,], probs,lty=1)         
          abline(h=0.05,lty=2)
          est.fpa <- approx(probs, allSimSSB[1,], 0.05,ties="ordered")$y  
          abline(v=est.fpa,lty=2) 
          text(est.fpa,1,round(est.fpa,2))
          text(max(allSimSSB[1,]),0.05,"0.05")
        graphics.off()      
      }
     
      }

     ##Text output
     cat("Stock name\n", stockname, '\n', file=outputfilename,sep="")
     cat("Sen filename\n", senfilename, '\n', file=outputfilename,append=TRUE,sep="")
     if (any(is.na(pfpm)))  cat("index filename\n", indexfilename, '\n', file=outputfilename,append=TRUE,sep="") else
                            cat("pf, pm\n", pfpm[1], '\t', pfpm[2],  '\n', file=outputfilename,append=TRUE,sep="") 
     cat("Number of iterations\n", nits, '\n', file=outputfilename,append=TRUE,sep="")
     cat("Simulate variation in Biological parameters\n", varybiodata, '\n', file=outputfilename,append=TRUE,sep="")  
     cat("SR relationship constrained\n", srconstrain, '\n', file=outputfilename,append=TRUE,sep="")     
     if (!onlyYPR)       
       for (srtype in sr)
       {
              cat('\n',srname[srtype], '\n', file=outputfilename,append=TRUE)
              cat(sum(!is.na(simdata[[srtype]]$alpha)),"/",nits,' Iterations resulted in feasible parameter estimates\n',sep="", file=outputfilename,append=TRUE)
              cat('',colnames(output[[srtype]]), '\n', file=outputfilename,append=TRUE,sep='\t')
              write.table(output[[srtype]],file=outputfilename,append=TRUE,quote=FALSE,col.names=FALSE,sep='\t',na="")
       }     

     cat('\n',"Per recruit", '\n', file=outputfilename,append=TRUE)
     cat('',colnames(outputpr), '\n', file=outputfilename,append=TRUE,sep='\t')
     write.table(outputpr,file=outputfilename,append=TRUE,quote=FALSE,col.names=FALSE,sep='\t',na="")
     if(!onlyYPR)     
     {
       cat('\n',"Combining all SRRs\n", file=outputfilename,append=TRUE)
       if (srautoweights) {
         cat('Automatically specified weights\n',names(srweights),'\n',srweightsexact,'\n', sep='\t', file=outputfilename,append=TRUE)
       } else {
         cat('Manually specified weights\n',names(srweights),'\n',srweightsexact,'\n', sep='\t', file=outputfilename,append=TRUE)       
       }
       cat('\n','Percentage\tFmsy\tFcrash\tMSY\tBmsy\tFmsy_w\tFcrash_w\tMSY_w\tBmsy_w\n', file=outputfilename,append=TRUE)       
       write.table(qs,file=outputfilename,append=TRUE,quote=FALSE,col.names=FALSE,sep='\t',na="")     
     }
    
     #Tidy up files
     #unlink(c("senf.tmp","sumf.tmp"))
     
     #return data to R as invisible list
     output$perRecruit = outputpr
     names(output) = srname
#     output$simData = simdata
     output$simData = list() 
     for (sr in srname) output$simData[[sr]] = rbind(simdatadet[[sr]],simdata[[sr]])
     output$simData[["combined"]] = rbind(NA,allSimData)
     output$simSSB = simSSB
     output$simSSBpr = simSSBpr
     output$simY = simy
     output$simYpr = simypr     
     invisible(output)
}

#' @title XXX
#' 
#' @description XXX
#' 
#' @export
#' 
#' @param sep XXX

paste0 <- function(... , sep='') paste(... , sep=sep)