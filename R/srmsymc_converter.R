#' @title FLR to srmsymc
#' 
#' @description XXX
#' 
#' @export
#' 
#' @param stk FLStock object
#' @param y Years
#' @param name_stock Name of the stock. May later down the line be used in 
#' table and plot outputs
#' @param filename_age Name of the associated file containing biological (weights,
#' maturity and mortalit) and fisheries data (selection patter) by age.
#' @param opt_sr_model Recruitment model number (1: Ricker, 2: Beverton-Holt,
#' 3: Segmented regression).
#' @param opt_land Boolean
#' @param opt_sim 0=no simulation, 1=simulation
#' @param opt_age 0=error only in recr, 1=error in recr & steady-state vectors
#' @param opt_pen 0=no SR constraints, 1=apply SR constrain (IS THIS ONLY FOR
#' THE SEGMENTED REGRESSION??).
FLS2srmsymc <- function(stk,y=3,name_stock,filename_age,opt_sr_model,
                        opt_land=TRUE,opt_sim=TRUE,opt_age=TRUE,opt_pen=TRUE) {
  x <- fls2list(stk,y=y,optLand=opt_land)
  y <- srmsymc_cat_srmsymc("codNS","age.dat",min(x$rby$year),max(x$rby$year),
                           aR=min(as.numeric((x$dims$age))),
                           aP=max(as.numeric((x$dims$age))),
                           1,1,1,1,
                           r=x$rby$rec,
                           ssb=x$rby$ssb,
                           year=x$rby$year)
  z <- srmsymc_cat_age(filename_age,n_fleets=2,pf=0,pm=0,
                       sel=cbind(x$sH[,1],x$sD[,1]),
                       sel_cv=cbind(x$sH[,2],x$sD[,2]),
                       w=cbind(x$wH[,1],x$wD[,1]),
                       w_cv=cbind(x$wH[,2],x$wD[,2]),
                       bio=cbind(x$M[,1],x$MT[,1],x$wS[,1]),
                       bio_cv=cbind(x$M[,2],x$MT[,2],x$wS[,1]))
}

#' @title XXX
#' 
#' @description From Tim
#' 
#' @export
#' 
#' @param stk FLStock object
#' @param y year
#' @param optLand Boolean, if TRUE (default) then ...

fls2list <- function(stk, y, optLand=TRUE)
{                                   
  
  d.flag <- (sum(discards(stk)) > 0) 
  ret <- vector("list",9)
  names(ret) <- c("rby","sH","sD","wH","wD","wS","M","MT","dims")
  
  ret$rby <- data.frame(year=as.numeric(dimnames(stk@stock.n)$year),
                    rec= c(stock.n(stk)[1,]),
                    ssb=c(ssb(stk)),
                    fbar=c(fbar(stk)),
                    yield=c(catch(stk)))
  
  years <- rev(rev(dimnames(stk@stock.n)$year)[1:y])
  ages <- dimnames(stk@stock.n)$age
  my_sum <- function(x) {r=cbind(mean=apply(x[,years,],1,mean),CV=apply(x[,years,],1,sd)/apply(x[,years,],1,mean));r[is.na(r)]<- 0; r}
  
  
  if (d.flag & optLand) {
    sC <- harvest(stk)[,years]
    sH <- harvest(stk)[,years] * (landings.n(stk)/catch.n(stk))[,years]
    sD <- harvest(stk)[,years] * (discards.n(stk)/catch.n(stk))[,years]    
    for (yy in years) sH[,yy] <- sH[,yy]/mean(sC[range(stk)["minfbar"]:range(stk)["maxfbar"],yy])
    for (yy in years) sD[,yy] <- sD[,yy]/mean(sC[range(stk)["minfbar"]:range(stk)["maxfbar"],yy])    
    ret$sH <- my_sum(sH)  
    ret$sD <- my_sum(sD)      
  } else {
    sH <- harvest(stk)[,years]
    for (yy in years) sH[,yy] <- sH[,yy]/mean(sH[range(stk)["minfbar"]:range(stk)["maxfbar"],yy])
    ret$sH <- my_sum(sH)
    ret$sD <- ret$sH*0  
  }
  
  ret$wH <- my_sum(landings.wt(stk))
  ret$wD <- my_sum(discards.wt(stk))
  ret$wD[is.na(ret$wD)] <- 0
  ret$wS <- my_sum(stock.wt(stk))  
  ret$M <- my_sum(m(stk)) 
  ret$MT <- my_sum(mat(stk)) 
  ret$MT[is.na(ret$MT)] <- 0 
  ret$dims <- list(age=ages,year=years,fbarage=range(stk)["minfbar"]:range(stk)["maxfbar"])
  
  #Set CVs to a more realistic value
  #i <- (ret$MT[,"mean"] >=0.05 & ret$MT[,"mean"] <=0.95)
  #ret$MT[i,"CV"] <- max(ret$MT[i,"CV"],0.1)
  
  #ret$M[,"CV"] <- max(ret$M[,"CV"],0.1)
  return (ret)
}



#' @title Write srmsymc.dat to disk
#' 
#' @description The function is used to write srmsymc.dat to the disk.
#' 
#' @author Einar Hjorleifsson
#' 
#' @export
#' 
#' @param name_stock Name of the stock. May later down the line be used in 
#' table and plot outputs
#' @param filename_age Name of the associated file containing biological (weights,
#' maturity and mortalit) and fisheries data (selection patter) by age.
#' @param y1 First year of ssb and recruitment data
#' @param y2 Last year of ssb and recruitment data
#' @param aR Recruitment age
#' @param aP Plus group age
#' @param opt_sr_model Recruitment model number (1: Ricker, 2: Beverton-Holt,
#' 3: Segmented regression).
#' @param opt_sim 0=no simulation, 1=simulation
#' @param opt_age 0=error only in recr, 1=error in recr & steady-state vectors
#' @param opt_pen 0=no SR constraints, 1=apply SR constrain (IS THIS ONLY FOR
#' THE SEGMENTED REGRESSION??).
#' @param r A vector containing recruitment
#' @param ssb A vector containing spawning stock biomass
#' @param year A vector containinb years (just used in comments).

srmsymc_cat_srmsymc <- function(name_stock, filename_age, y1, y2, aR, aP,
                                opt_sr_model, opt_sim, opt_age, opt_pen,
                                r, ssb, year) {
  tmpfile <- file('srmsymc.dat',open='w')
  cat('# Header: Some nice description\n',file=tmpfile,append=TRUE)
  cat(name_stock, '# stkname: Name of the stock\n',file=tmpfile,append=TRUE)
  cat(filename_age,'# filname: Name of the option file (2nd file\n',file=tmpfile,append=TRUE)
  cat(y1,    ' # ybeg:   First year                                                   (yearRange[1])\n',file=tmpfile,append=TRUE)
  cat(y2,    ' # yend:   Last year                                                    (yearRange[2])\n',file=tmpfile,append=TRUE)
  cat(aR,    ' # r:      Recruitment age                                              (senhead[1])\n',file=tmpfile,append=TRUE)
  cat(aP,    ' # A:      Plus group age                                              (senhead[2])\n',file=tmpfile,append=TRUE)
  cat(opt_sr_model,  ' # Ropt:   S-R function type                                            (sr)\n',file=tmpfile,append=TRUE)
  cat(opt_sim,' # simopt: 0=no simulation, 1=simulation                                (ifelse(nits==0,0,1)) \n',file=tmpfile,append=TRUE)
  cat(opt_age,' # senopt: 0=error only in recr, 1=error in recr & steady-state vectors (ifelse(varybiodata,1,0))\n',file=tmpfile,append=TRUE)
  cat(opt_pen,' # penopt: 0=no SR constraints, 1=apply SR constraints                  (ifelse(srconstrain,1,0))\n',file=tmpfile,append=TRUE)
  cat('# r ssb\n', file=tmpfile, append=TRUE)
  cat(paste(r,ssb,'#',year), file = tmpfile,append = TRUE,sep="\n")
  close(tmpfile)
}


#' @title Write age.dat to disk
#' 
#' @description XXX
#' 
#' @author Einar Hjorleifsson
#' 
#' @export
#' 
#' @param filename_age XXX
#' @param n_fleets XXX
#' @param pf XXX
#' @param pm XXX
#' @param sel XXX
#' @param sel_cv XXX
#' @param w XXX
#' @param w_cv XXX
#' @param bio XXX
#' @param bio_cv XXX
#' 
srmsymc_cat_age <- function(filename_age,n_fleets,pf,pm,sel,sel_cv,w,w_cv,bio,bio_cv) {
  
  n <- 4 # number of digits to write
  
  tmpfile <- file(filename_age,open='w')
  cat('#Header: Some nice description\n',file=tmpfile,append=TRUE)
  cat(n_fleets,  '# fno: Number of fleets (nstocks)\n',file=tmpfile,append=TRUE)
  cat(1,      '# sno: Fleets for yield per recruit stats (always 1)\n',file=tmpfile,append=TRUE)
  cat(pf,   '# f: proportional fishing mortality before spawning time (pf)\n',file=tmpfile,append=TRUE)
  cat(pm,   '# m: proportional natural mortality before spawning time (pm)\n',file=tmpfile,append=TRUE)
  cat('# Selection pattern\n',file=tmpfile,append=TRUE)
  for(i in 1:nrow(sel))   cat(format(round(sel[i,],n),n),'\n',file=tmpfile,append=TRUE)
  cat('# cv Selection pattern\n',file=tmpfile,append=TRUE)
  for(i in 1:nrow(sel_cv)) cat(format(round(sel_cv[i,],n)),'\n',file=tmpfile,append=TRUE)
  cat('# Weight at age\n',file=tmpfile,append=TRUE)
  for(i in 1:nrow(w))   cat(format(round(w[i,],n),n),'\n',file=tmpfile,append=TRUE)
  cat('# cv Weight at age\n',file=tmpfile,append=TRUE)
  for(i in 1:nrow(w_cv)) cat(format(round(w_cv[i,],n),n),'\n',file=tmpfile,append=TRUE)
  cat('# Biological data\n',file=tmpfile,append=TRUE)
  cat('#    M,     mat,   wSSB\n',file=tmpfile,append=TRUE)
  for(i in 1:nrow(bio))   cat(format(round(bio[i,],n),n),'\n',file=tmpfile,append=TRUE)
  cat('# cv Biological data\n',file=tmpfile,append=TRUE)
  cat('#  cvM,   cvmat, cvwSSB\n',file=tmpfile,append=TRUE)
  for(i in 1:nrow(bio_cv)) cat(format(round(bio_cv[i,],n),n),'\n',file=tmpfile,append=TRUE)
  close(tmpfile)
}




#' @title write data files to disk for the ADM srmsymc program
#' 
#' @description Some details, including link to the function a) running the
#' program and creating the data input
#' 
#' @author Einar Hjorleifsson <einar.hjorleifsson@gmail.com>
#' 
#' @export
#' 
#' @param rby A list that contains in its first three rows year, recruitment
#' and ssb.
#' @param iba A list - read.sen
#' @param aR \emph{integer}. Age of recruitment. Has to be specified. Is used to
#' shift the recruitment data such that it aligns with ssb data. Can be set to
#' 0 (zero) if ssb-recruitemnt pairs already properly aligned.
#' @param col_names vector of column names for year, recruitment and spawning stock
#' biomass (in that order).
#' @param opt_sr_model \emph{integer}. Number (1-3) corresponding to recruitment model.
#' @param opt_sim \emph{integer}. If 0 no simulation, if 1 simulation.
#' #' @param opt_pen \emph{integer}. If 0 then no recruitement constraints if 1
#' then recruitment contraints applied (NEED MORE DETAIL)
#' @param opt_age \emph{integer}. If 0 then error only in recruitment, if 1 then
#' both error in recruitment and steady-state vectors (age based inputs).
#' @param rba data.frame that contains age based input NOT YET IMPLEMENTED
#' @param filename_age \emph{character}. Name of the option file that contains the
#' age-based-data. If missing, set to "age.dat".
#' @param aP \emph{integer}. Plus group age.
#' @param years vector containing years to include in the ssb-r

write.srmsymc <- function(rby,iba,aR,col_names=c("year","r","ssb"),opt_sr_model=1,
                          opt_sim=1,opt_pen=1,opt_age=1,rba,
                          filename_age="age.dat",aP,years) {
  
  
  
  value <- NULL
  
  ## rby tests
  if(missing(rby)) stop("summary (year) data must be specified")    
  class2 <- rby$info$class2
  if(is.null(class2)) stop("Not implemented yet for provided file type")
  if(class2 != "rby") stop("Not implemented yet for provided file type")
  
  x <- rby$data[,col_names]
  y <- rby$info
  ## rba test
  #if(missing(rba)) stop("prediction (age) data must be specified")
  #class2 <- attributes(rba)$class2
  #if(is.null(class2)) stop("Not implemented yet for provided file type")
  #if(class2 != "rby") stop("Not implemented yet for provided file type")
  
  if(missing(filename_age)) filename_age  <- "age.dat"
  
  ### A. Setting up the srmsymc.dat
  
  ## setting stock name
  name_stock  <- y$name_stock
  if(is.na(name_stock)) name_stock <- "nn"
  name_stock <- y$name_stock <- str_replace_all(name_stock," ","")
  
  aR <- y$time[3]
  if(is.na(aR)) stop("Recruitment age (aR) must be specified")
  
  ## Align recruitment
  x <- align_ssb_r(x,col.names=names(x),aR)
  x <- x[!is.na(x$r),]
  
  if(missing(years)) {
    y1 <- y$time[1]
    y2 <- y$time[2]
    years <- c(y1:y2)
  } else {
    y$time[1] <- y1 <- min(years)
    y$time[2] <- y2 <- max(years)
  }
  
  x <- x[x$year %in% years,]
  
  if(missing(aP)) stop("Plus group age (aP) must be specified")
  y$time[6] <- aP
  
  # return file
  x <- data.frame(r=x$r,ssb=x$ssb,year=x$year)
  y$filename_age = filename_age
  y$opt_sr_model = opt_sr_model
  y$opt_sim = opt_sim
  y$opt_age = opt_age
  y$opt_pen = opt_pen
  rby <- list(data=x,info=y)
  
  srmsymc_cat_srmsymc(name_stock, filename_age, y1, y2, aR, aP,
                  opt_sr_model, opt_sim, opt_age, opt_pen,
                  r=x$r, ssb=x$ssb, year=x$year)
  
  ### B. Setting up the age.dat
  
  if(iba$info$creator == "created from function fishvise:::read.sen") {
    x <- iba$data
    y <- iba$info
    nFleets <- sum(y$fleets[2:4])
    fleet_Fs <- y$mort[,c(2:4)]
    fleet_Fs_names <- colnames(fleet_Fs)
    weight_names <- str_replace(fleet_Fs_names,"F","W")
    ages <- c(y$time[4]:y$time[5])
    sel <- cvsel <- wgt <- cvwgt <- matrix(NA,nrow=length(ages),ncol=nFleets)
    for (i in 1:nFleets) {
      x1 <- x[x$id %in% fleet_Fs_names[i],]
      x1$value <- x1$value/mean(x1$value[x1$age %in% c(fleet_Fs[1,i]:fleet_Fs[2,i])])
      sel[,i] <- x1[,'value']
      cvsel[,i] <- x1[,'cv']
      x1 <- x[x$id %in% weight_names[i],]
      wgt[,i] <- x1[,"value"]
      cvwgt[,i] <- x1[,"value"]
    }
    bio <- cvbio <- matrix(NA,nrow=length(ages),ncol=3)
    bio[,1]   <- x[x$id %in% 'M','value']
    bio[,2]   <- x[x$id %in% 'xN','value']
    bio[,3]   <- x[x$id %in% 'wN','value']
    cvbio[,1] <- x[x$id %in% 'M','cv']
    cvbio[,2] <- x[x$id %in% 'xN','cv']
    cvbio[,3] <- x[x$id %in% 'wN','cv']
    cat_age.dat(filename_age,nFleets,y$pF,y$pM,sel,cvsel,wgt,cvwgt,bio,cvbio)
    iba <- list(sel=sel,sel_cv=cvsel,w=wgt,w_cv=cvwgt,bio=bio,bio_cv=cvbio,info=y)
  }
  
  if(iba$info$creator == "fishvise::read.rby") {
    
    x <- iba$data[,c("year","age","tF","wC","wX","xN")]
    
    tF <- ddply(x[x$age %in% 5:10,],"year",summarise,refF=mean(tF))
    x <- join(x,tF)
    x$sF <- x$tF/x$refF
    
    d <- melt(x[,c("year","age","sF","wC","wX","xN")],id.vars=c("year","age"))
    d <- ddply(d,c("variable","age"),summarise,ave=mean(value,na.rm=T),cv=sqrt(var(value,na.rm=T))/ave)
    
    nFleets <- 1
    fleet_Fs <- matrix(c(5,10),ncol=1,nrow=2)
    colnames(fleet_Fs) <- "sF"
    fleet_Fs_names <- colnames(fleet_Fs)
    weight_names <- "wC"
    
    ages <- c(min(x$age):max(x$age))
    cat("Line 621")
    sel <- cvsel <- wgt <- cvwgt <- matrix(NA,nrow=length(ages),ncol=nFleets)
    x <- d
    names(x) <- c("id","age","value","cv")
    for (i in 1:nFleets) {
      x1 <- x[x$id %in% fleet_Fs_names[i],]
      x1$value <- x1$value/mean(x1$value[x1$age %in% c(fleet_Fs[1,i]:fleet_Fs[2,i])])
      sel[,i] <- x1[,'value']
      cvsel[,i] <- x1[,'cv']
      x1 <- x[x$id %in% weight_names[i],]
      wgt[,i] <- x1[,"value"]
      cvwgt[,i] <- x1[,"value"]
    }
    bio <- cvbio <- matrix(NA,nrow=length(ages),ncol=3)
    bio[,1]   <- 0.2
    bio[,2]   <- x[x$id %in% 'xN','value']
    bio[,3]   <- x[x$id %in% 'wX','value']
    cvbio[,1] <- 0.0
    cvbio[,2] <- x[x$id %in% 'xN','cv']
    cvbio[,3] <- x[x$id %in% 'wX','cv']
    
    cat_age.dat("age.dat",1,0,0,sel,cvsel,wgt,cvwgt,bio,cvbio) 
    iba <- list(sel=sel,sel_cv=cvsel,w=wgt,w_cv=cvwgt,bio=bio,bio_cv=cvbio,info=y)
  }
  
  
  
  return(list(rby=rby,iba=iba))
}




#' @title Write age.dat to disk
#' 
#' @description XXX
#' 
#' @author Einar Hjorleifsson
#' 
#' @export
#' 
#' @param filename_age XXX
#' @param n_fleets XXX
#' @param pf XXX
#' @param pm XXX
#' @param sel XXX
#' @param sel_cv XXX
#' @param w XXX
#' @param w_cv XXX
#' @param bio XXX
#' @param bio_cv XXX
#' 
cat_age.dat <- function(filename_age,n_fleets,pf,pm,sel,sel_cv,w,w_cv,bio,bio_cv) {
  
  n <- 4 # number of digits to write
  
  tmpfile <- file(filename_age,open='w')
  cat('#Header: Some nice description\n',file=tmpfile,append=TRUE)
  cat(n_fleets,  '# fno: Number of fleets (nstocks)\n',file=tmpfile,append=TRUE)
  cat(1,      '# sno: Fleets for yield per recruit stats (always 1)\n',file=tmpfile,append=TRUE)
  cat(pf,   '# f: proportional fishing mortality before spawning time (pf)\n',file=tmpfile,append=TRUE)
  cat(pm,   '# m: proportional natural mortality before spawning time (pm)\n',file=tmpfile,append=TRUE)
  cat('# Selection pattern\n',file=tmpfile,append=TRUE)
  for(i in 1:nrow(sel))   cat(format(round(sel[i,],n),n),'\n',file=tmpfile,append=TRUE)
  cat('# cv Selection pattern\n',file=tmpfile,append=TRUE)
  for(i in 1:nrow(sel_cv)) cat(format(round(sel_cv[i,],n)),'\n',file=tmpfile,append=TRUE)
  cat('# Weight at age\n',file=tmpfile,append=TRUE)
  for(i in 1:nrow(w))   cat(format(round(w[i,],n),n),'\n',file=tmpfile,append=TRUE)
  cat('# cv Weight at age\n',file=tmpfile,append=TRUE)
  for(i in 1:nrow(w_cv)) cat(format(round(w_cv[i,],n),n),'\n',file=tmpfile,append=TRUE)
  cat('# Biological data\n',file=tmpfile,append=TRUE)
  cat('#    M,     mat,   wSSB\n',file=tmpfile,append=TRUE)
  for(i in 1:nrow(bio))   cat(format(round(bio[i,],n),n),'\n',file=tmpfile,append=TRUE)
  cat('# cv Biological data\n',file=tmpfile,append=TRUE)
  cat('#  cvM,   cvmat, cvwSSB\n',file=tmpfile,append=TRUE)
  for(i in 1:nrow(bio_cv)) cat(format(round(bio_cv[i,],n),n),'\n',file=tmpfile,append=TRUE)
  close(tmpfile)
}







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

  
  
  
  
  
  
