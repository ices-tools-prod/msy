# convertSumSen.r: R code to convert sen and sum files to out.dat and srmsy.dat used for
#                  input to srmsy.exe
# developed by Tim Earl at Cefas, UK, April 2010,
# in collaboration with Chris Darby and José De Oliveira, also at Cefas.
#
# Requires the following file in the same directory:
# - readLow.r if pf anf pm are to be extracted from Lowestoft format files
#   
# plotMSY will often be used as a wrapper to this function
#


#call convertfiles(senfilename=NA,indexfilename=NA,pfpm=NA) for interactive version
#     senfilename    is the name of the senfile, with a corresponding sumfile sharing the 
#                    same name except replacing ".sen" with ".sum".
#                    Produces an interactive window if not specified
#     pfpm           is a vector of two values; pm and pf (proportion of m and f before spawning
#     indexfilename  is the name of an index file in the Lowestoft format pointing to the pf and pm files.
#                    If pfpm is specified, this is ignored, if neither specified an interactive window appears to choose index file 
#     nits           Number of iterations of bootstrapping - if 0, does only the deterministic fit
#     sr             Stock recruit relationship, 1 = Ricker, 2= BevHolt, 3=Hockeystick
#     varybiodata    If TRUE, bootstraps the biological data (weight, maturity, mortality) 
#                    if FALSE, varies SR relationship only 
#     stockname      Display title for stock
#     silent         If TRUE output to screen then routine output to screen is supressed
#     srconstrain    If TRUE, SR parameters must be greater than zero (and breakpoint within data for hockeystick), otherwise parameters may take any value

##Version 2 - reads comma separated sen and sum files

###Only tested on windows - uses 'choose.files' which may not work on other operating systems





convertSumSen <- function(senfilename=NA,indexfilename=NA,pfpm=NA, nits=0, sr=2, varybiodata=FALSE,stockname="",silent=FALSE,srconstrain=TRUE)
{

  #filenames
  senfilename = tolower(senfilename)
  indexfilename = tolower(indexfilename)
  if (is.na(senfilename)) senfilename = tolower(choose.files("*.sen", "Choose SEN file",multi=FALSE)) 
  if (!file.exists(senfilename)) stop("SEN file not found")
  sumfilename = sub(".sen",".sum",senfilename,fixed=TRUE)
  if (!file.exists(sumfilename)) stop("SUM file not found")  
 # path = gsub("/","\\\\",senfilename) #converts "c:/" to "c:\\"  
 # path = paste(rev(rev(strsplit(senfilename, "\\\\")[[1]])[-1]),collapse="\\")

  opfilename = ".\\out.dat"
  srfilename = ".\\srmsymc.dat"
  
  # create temp file with commas replaced by spaces
  senf = scan(senfilename,"",sep="\n",quote=NULL,quiet=silent)
  senf = gsub(","," ",senf) 
  cat(senf,file = "senf.tmp",sep="\n",quote=NULL)
  senfilename = "senf.tmp"
  sumf = scan(sumfilename,"",sep="\n",quote=NULL,quiet=silent)
  sumf = gsub(","," ",sumf)
  cat(sumf,file = "sumf.tmp",sep="\n",quote=NULL)
  sumfilename = "sumf.tmp"
  
  
                  
  if (any(is.na(pfpm)))
  {                
    if (is.na(indexfilename)) indexfilename = tolower(choose.files("*.*", "Choose index file",multi=FALSE))                         
    if (!file.exists(indexfilename)) stop("Index file not found")   
    source("ReadLow.r")     
    stock = read.Lowestoft(indexfilename, 7:8,expand=FALSE,silent=silent)
  } else {
    if (length(pfpm)!=2) stop("pfpm must be a vector of length 2")
    stock = list(pF = pfpm[1],pM = pfpm[2])
  }

  senhead = scan(senfilename,nlines=2,skip=1,quiet=TRUE)
  ages = senhead[1]:senhead[2]
  nstocks = sum(senhead[5:7])
  
    ##Set up matrix to contain parameters
    #ages = 1:7
    Vars = c("SH","SD","SI","WH","WD","WI","M","MT","Wx","WS")       
    Params = t(matrix(paste("\'",outer(Vars,ages,paste,sep=""),"\'",sep=""),length(Vars),length(ages)))
    colnames(Params) = Vars
    rownames(Params) = ages  
    Params[,"Wx"] = 0
    if (senhead[5] == 0) exit("No Human consumption data")
    if (senhead[6] == 0) Params[,c("SD","WD")] = 0      #Discards
    if (senhead[7] == 0) Params[,c("SI","WI")] = 0      #Industrial
    cv = Params
    
                                                                     
    #Read in data from SEN file
    headlines = 3
    reading_parameters = TRUE
    while(reading_parameters)
    {
      line = scan(senfilename,what="",nlines=1,skip=headlines,quote="\"",quiet=TRUE)     
      
      if (length(line)< 3) 
      {
        reading_parameters = FALSE 
      } else {
      if(substr(line[1],1,1) == "\'")
      {
        if (line[1] == "\'")
        {
          line[1] = paste(line[1],line[2],sep="")   #to get rid of any space in the quotes
          line[2] = line[3]
          line[3] = line[4]          
        }
        cv[Params==toupper(line[1])]=line[3]
        Params[Params==toupper(line[1])]=line[2]
        headlines = headlines + 1
      } else {
        line[1] = paste('\'',line[1],'\'',sep="")
        cv[Params==toupper(line[1])]=line[3]
        Params[Params==toupper(line[1])]=line[2]
        headlines = headlines + 1
        
        
        }
      }
    
    }  

    #Read in data from SUM file
    yearRange =  scan(sumfilename,nlines=1,skip=4,quote="\"",quiet=TRUE)
    ncols =  scan(sumfilename,nlines=1,skip=1,quote="\"",quiet=TRUE)
    sumData =  t(matrix(scan(sumfilename,nlines=yearRange[2]-yearRange[1]+1,skip=27,quote="\"",quiet=TRUE),ncols))
    Params = matrix(as.numeric(Params),length(ages),length(Vars),dimnames=dimnames(Params))
    cv = matrix(as.numeric(cv),length(ages),length(Vars),dimnames=dimnames(Params))
    #normalise selectivity
    agerange = scan(sumfilename,nlines=1,skip=20,quiet=TRUE)     
    
    normfactor = sum(mean(Params[ages %in% (agerange[1]:agerange[2]),"SD"] ),mean(Params[ages %in% (agerange[1]:agerange[2]),"SH"] ),mean(Params[ages %in% (agerange[1]:agerange[2]),"SI"] ))

    if(mean(Params[ages %in% (agerange[1]:agerange[2]),"SD"])>0) Params[,"SD"] = Params[,"SD"]/ normfactor
    if(mean(Params[ages %in% (agerange[1]:agerange[2]),"SH"])>0) Params[,"SH"] = Params[,"SH"]/ normfactor
    if(mean(Params[ages %in% (agerange[1]:agerange[2]),"SI"])>0) Params[,"SI"] = Params[,"SI"]/ normfactor
    
    Params[,"Wx"] = (Params[,"WD"]*Params[,"SD"]+Params[,"WH"]*Params[,"SH"]+Params[,"WI"]*Params[,"SI"])/(Params[,"SH"]+Params[,"SD"]+Params[,"SI"])
    Params[is.na(Params[,"Wx"]),"Wx"] =0
    cv[,"Wx"] = (cv[,"WD"]*Params[,"SD"]+cv[,"WH"]*Params[,"SH"]+cv[,"WI"]*Params[,"SI"])/(Params[,"SH"]+Params[,"SD"]+Params[,"SI"])

    
    
    #Read in data from pf and pm file
    pf = stock$pF[1]
    pm = stock$pM[1]
    if (all((prod(dim(stock$pF))>1),sd(as.vector(stock$pF))!=0)) warning("pf not a single value")
    if (all((prod(dim(stock$pM))>1),sd(as.vector(stock$pM))!=0)) warning("pm not a single value")    

##  output  to \\out.dat
    s=t(matrix(rep(senhead[5:7],times=length(ages)),3))
    opfile = file(opfilename,"w")
      cat('#fno, sno, f, m // fno=nr fleets; sno=fleet for ypr stats; f=F before spwn; m=M before spwn\n', file = opfile)
      cat(nstocks,'1', pf, pm, '\n', file = opfile,append = TRUE,sep=" ")        ##Always maximises stock 1 i.e. Human consumption
      cat('#sdat // col1=sel fleet 1 (landings); col2=sel fleet 2\n', file = opfile,append = TRUE,sep="")
      cat(paste(ifelse(s[,1]==1,Params[,"SH"],""),ifelse(s[,2]==1,Params[,"SD"],""),ifelse(s[,3]==1,Params[,"SI"],"")), file = opfile,append = TRUE,sep="\n")
      cat('#sdat cv // col1=sel fleet 1 (landings); col2=sel fleet 2\n', file = opfile,append = TRUE,sep="")
      cat(paste(ifelse(s[,1]==1,cv[,"SH"],""),ifelse(s[,2]==1,cv[,"SD"],""),ifelse(s[,3]==1,cv[,"SI"],"")), file = opfile,append = TRUE,sep="\n")
      cat('#wdat // col1=wght fleet 1 (landings); col2=wght fleet 2\n', file = opfile,append = TRUE,sep="")
      cat(paste(ifelse(s[,1]==1,Params[,"WH"],""),ifelse(s[,2]==1,Params[,"WD"],""),ifelse(s[,3]==1,Params[,"WI"],"")), file = opfile,append = TRUE,sep="\n")
      cat('#wdat cv // col1=wght fleet 1 (landings); col2=wght fleet 2\n', file = opfile,append = TRUE,sep="")
      cat(paste(ifelse(s[,1]==1,cv[,"WH"],""),ifelse(s[,2]==1,cv[,"WD"],""),ifelse(s[,3]==1,cv[,"WI"],"")), file = opfile,append = TRUE,sep="\n")
      cat('#biodat // col1=natural mortality; col2=maturity; col3=wght in stock\n', file = opfile,append = TRUE,sep="")
      cat(paste(Params[,"M"],Params[,"MT"],Params[,"WS"]), file = opfile,append = TRUE,sep="\n")
      cat('#biodat cv // col1=natural mortality; col2=maturity; col3=wght in stock\n', file = opfile,append = TRUE,sep="")
      cat(paste(cv[,"M"],cv[,"MT"],cv[,"WH"]), file = opfile,append = TRUE,sep="\n")
    close(opfile)
    

 
## output to \\srmsy.dat
    opfile = file(srfilename,"w")
      cat('#stkname, filname // stkname=stock dealing with; filname=name of 2nd file\n', file = opfile)
      cat(sub(".dat","", rev(strsplit(opfilename,"\\\\")[[1]])[1]),rev(strsplit(opfilename,"\\\\")[[1]])[1],'\n', file = opfile,append = TRUE,sep=" ") 
      cat('#ybeg, yend, r, A, Ropt, simopt, senopt,penopt // ybeg=1st yr; yend=last yr; r=recr age; A=plusgroup; Ropt=S-R function type, simopt (0=no sim, 1=do sim); senopt (0=error only in recr, 1=error in recr & steady-state vectors);  penopt (0=no SR constraints, 1=apply SR constraints)\n', file = opfile)
      cat(yearRange, senhead[1], senhead[2], sr, ifelse(nits==0,0,1), ifelse(varybiodata,1,0), ifelse(srconstrain,1,0),'\n', file = opfile,append = TRUE,sep=" ")  
      cat('#R, Bssb // R=recr; Bssb=SSB\n', file = opfile,append = TRUE,sep="")
      cat(paste(sumData[,2],sumData[,3]), file = opfile,append = TRUE,sep="\n")
    close(opfile)
    
## output to \\sim.dat
    opfile = file(".\\sim.dat","w")
      cat('21 1 ', nits, "\n", sep="", file = opfile)
    close(opfile)    

  if (!silent) cat("Finished writing files\n")
  if (!silent) cat (opfilename,"\n")
  if (!silent) cat (srfilename,"\n")  
  
#  unlink("senf.tmp")
#  unlink("sumf.tmp")  
  ##Return variables for debugging
  return(list(Year = sumData[,1],Recruits = sumData[,2],SSB=sumData[,3],recruitage=senhead[1],stock=stockname))
}





#tmp = convertSumSen("cod_vmcm_late.sen",pfpm=c(0,0))
#tmp = convertfiles()
                                           
 
                                           
convertDat = function(datfilename, srtype)                                           
{
##
#  datfilename = ".\\data\\cod\\srmsy.dat"
#  srtype = 1
##
  datfilename = tolower(datfilename)
  if (!file.exists(datfilename)) stop("file not found:", datfilename)
  dathead1 = scan(datfilename,"",comment.char='#',nlines=2,quiet=TRUE)
  dathead2 = scan(datfilename,0,comment.char='#',nlines=2,skip=2,quiet=TRUE)  
  dathead2[5] = srtype
  datbody = scan(datfilename,list(0,0),comment.char='#',skip=5,quiet=TRUE) 
  outfilename = paste(dirname(datfilename),dathead1[[2]],sep="/")
  if (!file.exists(outfilename)) stop("file not found:", outfilename)
  outtext = scan(outfilename,"",quiet=TRUE,sep='\n')
  dattext = scan(datfilename,"",quiet=TRUE,sep='\n')  
  
#  write to out.dat and srmsymc.dat
  cat(outtext, sep='\n', file=".\\out.dat")
  cat(dattext[1], sep='\n', file=".\\srmsymc.dat")
  cat(dathead1[1],'out.dat\n', sep=' ', file=".\\srmsymc.dat", append=TRUE)    
  cat(dattext[3], sep='\n', file=".\\srmsymc.dat", append=TRUE)  
  cat(dathead2,'\n', sep=' ', file=".\\srmsymc.dat", append=TRUE)    
  cat(dattext[-(1:4)], sep='\n', file=".\\srmsymc.dat", append=TRUE)    

  return(list(Year = dathead2[1]:dathead2[2],Recruits = datbody[[1]],SSB=datbody[[2]],recruitage=dathead2[3],stock=dathead1[1]))
}


                        
                        